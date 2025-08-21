from __future__ import annotations

import os
from collections.abc import Callable
from pathlib import Path

from flask import Response, jsonify, make_response, request

from rnaprobes.util import validate_arg
from .usage_tracker import add_run_to_db
from .Program import Program
from .DelayedProgram import DelayedProgram
from ..rnaprobes.util import optional_argument, safe_remove_tree
from ..rnaprobes.TFOFinder import tfofinder
from ..rnaprobes.PinMol import pinmol
from ..rnaprobes.smFISH import smFISH

from werkzeug.datastructures import FileStorage
from werkzeug.utils import secure_filename
from functools import partial
import time, io, zipfile

root = Path(os.getcwd())
output_dir = root / "user-files"
pinmol_output_dir, sm_fish_output_dir = output_dir / "pinmol", output_dir / "smfish"

def set_root(root_path: Path, output_file_path: str | Path = "user-files"):
    global root, output_dir, pinmol_output_dir, sm_fish_output_dir
    root = root_path
    output_dir = (root / output_file_path).resolve()
    assert root_path.is_absolute(), "Given root must be an absolute path"
    assert root in output_dir.parents, "output_file_path must be a child of root"
    pinmol_output_dir, sm_fish_output_dir = output_dir / "pinmol", output_dir / "smfish"

def close_file(func: Callable, *, file_arg_name: str = "filein", **kwargs):
    """
    Execute a function and close the file afterward.
    :param func: the function to execute
    :param file_arg_name: the name of the argument representing the file to be closed
    :param kwargs: arguments including one corresponding to file_arg_name
    :return:
    """
    with kwargs[file_arg_name]:
        return func(**kwargs)

def get_zip_bytes(filename: str, **kwargs) -> bytes:
    file_obj = io.BytesIO()
    with zipfile.ZipFile(file_obj, 'w', compression=zipfile.ZIP_DEFLATED) as zip_file:
        for name, value in kwargs.items():
            zip_info = zipfile.ZipInfo(name)
            zip_info.date_time = time.localtime(time.time())[:6]
            zip_file.writestr(f"{filename}-{name.replace('_', '-')}{value[0]}", value[1])
            # print(f"{filename}-{name.replace("_", "-")}{value[0]}")
    file_obj.seek(0)
    return file_obj.getvalue()

def save_to_file(file_storage: FileStorage, path: Path, parent_exist_ok = True):
    path.parent.mkdir(parents=True, exist_ok=parent_exist_ok) #not needed in this case
    file_storage.save(path)


def smFISH_get_args(req, output_dir: Path) -> dict:
    file_path = output_dir / secure_filename(req.files.get("ct-file").filename)
    save_to_file(req.files.get("ct-file"), file_path)
    return dict(file_path = file_path,
        output_dir = output_dir,
        arguments = smFISH.parse_arguments("-d " + ("-i" if req.form.get("smFISH-intermolecular") else "-ni"), from_command_line=False)
         )

program_dict = { #get args, validate args, return value
    'tfofinder': Program("TFOFinder",
                         lambda req, _: dict(filein = req.files.get("ct-file").stream,
                                             probe_lengths = req.form.get("tfofinder-probe-length"),  #validation is in validate_arguments
                                             filename = secure_filename(req.files.get("ct-file").filename),
                                             arguments = tfofinder.parse_arguments("", from_command_line=False)),
                         tfofinder.validate_arguments,
                         partial(close_file, tfofinder.calculate_result), root_dir=output_dir, folder_not_needed=True),
    'pinmol': Program("PinMol", lambda req, output_dir: dict(filein= req.files.get("ct-file").stream,
                                   probe_length = req.form.get("pinmol-probe-length", type=int),
                                   output_dir = output_dir,
                                   filename = secure_filename(req.files.get("ct-file").filename),
                                    arguments = pinmol.parse_arguments(f"-nb -w"
                                                                        f"{optional_argument(req, 'pinmol-start-base', '-s', default_value=1)}"
                                                                        f"{optional_argument(req, 'pinmol-end-base', '-e', default_value=-1)}",
                                                                       from_command_line=False)
), pinmol.validate_arguments, partial(close_file, pinmol.calculate_result), output_dir=pinmol_output_dir, root_dir=pinmol_output_dir),
    'smfish': DelayedProgram("smFISH", smFISH_get_args, smFISH.validate_arguments, smFISH.calculate_result, output_dir=sm_fish_output_dir, root_dir=sm_fish_output_dir)
    .set_extra_notification_string_callback(lambda args: smFISH.get_size_warning(args["nucleotide_length"]))
}

def get_program_object(prog_name: str) -> Program:
    return program_dict[prog_name.lower()]

def run_program(prog_name: str, user_id: str, error_message_validation: str ="Something went wrong", error_message_program: str="Something went wrong"):
    # import time
    # prev = time.time_ns()
    program = get_program_object(prog_name)
    result =  program.run(request, error_message_validation, error_message_program)
    response = get_program_response(result, program.name)
    if response.status_code == 200: log_program_success(program.name, user_id) #OK for DelayedProgram since it still verifies the arguments
    # print(((time.time_ns() - prev) // 1_000) / 1_000)
    return response

def get_program_response(result: tuple[str, int] | Response | dict, program: str) -> Response:
    if type(result) == Response: return result
    if type(result) == dict: return jsonify(**result)
    return make_response(*result) #tuple

def log_program_success(program: str, user_id):
    add_run_to_db(user_id, program)