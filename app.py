#python version >=3.9
from __future__ import annotations

import os
import uuid
from collections.abc import Callable

from flask import Flask, render_template, request, redirect, url_for, jsonify, Response
from sys import argv

from werkzeug.datastructures import FileStorage

from DelayedProgram import DelayedProgram
from Program import Program

from src.TFOFinder import tfofinder
from src.PinMol import pinmol
from src.smFISH import smFISH
from src.util import optional_argument, get_folder_as_zip, ValidationError, safe_remove_tree
from werkzeug.utils import secure_filename
from functools import partial
import time, io, zipfile, base64
from pathlib import Path

app = Flask(__name__)
app.config['TEMPLATES_AUTO_RELOAD'] = True\

root = Path(__file__).parent
output_dir = root / "user-files"
pinmol_output_dir, sm_fish_output_dir = output_dir / "pinmol", output_dir / "smfish"


@app.errorhandler(500)
def application_error(error):
    return str(error), 500

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

def delete_folder_tree(dir_arg_name: str = "output_dir", root_dir: Path = root, **kwargs):
    """
    Delete the folder tree located in the kwargs. For use in as Program object's run_finally function
    :param dir_arg_name: the name of the argument representing the folder to be deleted
    :param root_dir: the root of the directory. Will fail if the directory is not a child of the root
    :param kwargs: arguments including one corresponding to dir_arg_name
    :return:
    """
    safe_remove_tree(kwargs[dir_arg_name], root_dir)


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

def get_result_temp(program: str, result: str, **kwargs) -> tuple[bytes, str]:
    zip_file = get_zip_bytes("tempfile", **{program: (".txt", result)})
    filename = f"TempResult-{program}.zip"
    return zip_file, filename

def save_to_file(file_storage: FileStorage, path: Path, parent_exist_ok = True):
    path.parent.mkdir(parents=True, exist_ok=parent_exist_ok) #not needed in this case
    file_storage.save(path)


def smFISH_get_args(req, output_dir: Path) -> dict:
    file_path = output_dir / secure_filename(req.files.get("ct-file").filename)
    save_to_file(req.files.get("ct-file"), file_path)
    return dict(file_path = file_path,
        output_dir = output_dir,
        arguments = smFISH.parse_arguments("-d "+ ("-i" if req.form.get("smFISH-intermolecular") else "-ni"), from_command_line=False)
         )

program_dict = { #get args, validate args, return value
    'tfofinder': Program("TFOFinder",
                        lambda req, _: dict(filein = req.files.get("ct-file").stream,
                                      probe_lengths = req.form.get("tfofinder-probe-length"), #validation is in validate_arguments
                                      filename = secure_filename(req.files.get("ct-file").filename),
                                         arguments = tfofinder.parse_arguments("", from_command_line=False)),
                tfofinder.validate_arguments,
                  partial(close_file, tfofinder.calculate_result), root_dir=root, folder_not_needed=True),
    'pinmol': Program("PinMol", lambda req, output_dir: dict(filein= req.files.get("ct-file").stream,
                                   probe_length = req.form.get("pinmol-probe-length", type=int),
                                   output_dir = output_dir,
                                   filename = secure_filename(req.files.get("ct-file").filename),
                                    arguments = pinmol.parse_arguments(f"-nb -w"
                                                                        f"{optional_argument(req, 'pinmol-start-base', '-s', default_value=1)}"
                                                                        f"{optional_argument(req, 'pinmol-end-base', '-e', default_value=-1)}",
                                                                        from_command_line=False)
), pinmol.validate_arguments, partial(close_file, pinmol.calculate_result), output_dir=pinmol_output_dir, root_dir=root),
    'smfish': DelayedProgram("smFISH", smFISH_get_args, smFISH.validate_arguments, smFISH.calculate_result, output_dir=sm_fish_output_dir, root_dir=root)
}

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/form', methods=['GET', 'POST'])
def form():
    if request.method == 'POST':
        programs = request.form.getlist('programs')
        if not programs:
            # Redirect back to home if no programs were selected. Needed if you navigate directly to it without using index.html first.
            return redirect(url_for('index'))
        return render_template('request.html', programs=programs, exports=get_exports())
    # If accessed via GET, redirect to home
    return redirect(url_for('index'))

@app.route('/send-request', methods=['POST'])
def send_request():
    program = request.args.get("program")
    return run_program(program, "The given arguments are invalid",
                         "Something went wrong when calculating your result")

@app.route('/query-result', methods=['GET'])
def query_result():
    program_name = secure_filename(request.args.get("program")) #prevents injection attacks
    id = uuid.UUID(request.args.get("id")) #prevents injection attacks
    program = program_dict[program_name.lower()]
    output_dir = program.output_dir / str(id)
    if isinstance(program, DelayedProgram):
        return program.get_current_result(output_dir, id)
    else:
        return f"Can only query a DelayedProgram. {program_name} is not a DelayedProgram", 400

def get_exports():
    return {"TFOFinder": tfofinder.exported_values, "PinMol": pinmol.exported_values, "smFISH": smFISH.exported_values}

def run_program(prog_name, error_message_validation="Something went wrong", error_message_program="Something went wrong"):
    # import time
    # prev = time.time_ns()
    program = program_dict[prog_name.lower()]
    result =  program.run(request, error_message_validation, error_message_program)
    # print(((time.time_ns() - prev) // 1_000) / 1_000)
    return result

if __name__ == '__main__':
    app.run(debug= (True if len(argv) > 1 else False)) #use debug if any commmand line arguments are inputted
