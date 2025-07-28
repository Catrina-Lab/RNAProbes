#python version >=3.9
from __future__ import annotations
import uuid
from collections.abc import Callable

from flask import Flask, render_template, request, redirect, url_for, jsonify, Response
from sys import argv

from src.RNASuiteUtil import ProgramObject
from src.TFOFinder import tfofinder
from src.PinMol import pinmol
from src.util import optional_argument, get_folder_as_zip, ValidationError, safe_remove_tree
from werkzeug.utils import secure_filename
from functools import partial
import time, io, zipfile, base64
from pathlib import Path

app = Flask(__name__)
app.config['TEMPLATES_AUTO_RELOAD'] = True

root = Path(__file__).parent
output_dir = root / "user-files"
pinmol_output_dir, sm_fish_output_dir = output_dir / "pinmol", output_dir / "smfish"


def send_error_response(error: Exception, **kwargs):
    if isinstance(error, ValidationError):
        return str(error), 400
    else:
        return str(error), 500


class Program:
    def _get_args(self, request, job_id, error_message: str="Something went wrong"):
        try:
            return self._get_args_raw(request, job_id)
        except Exception as e:
            raise ValidationError(f"{error_message}: {str(e)}") from e

    def _validate_args(self, kwargs: dict, error_message: str="Something went wrong") -> dict:
        try:
            return self._validate_args_raw(**kwargs) or dict()
        except Exception as e:
            raise ValidationError(f"{error_message}: {str(e)}") from e

    def _run_program(self, kwargs: dict, error_message: str="Something went wrong", validate_err_msg: str=None):
        try:
            return self._run_program_raw(**kwargs)
        except ValidationError as e:
            raise ValidationError(f"{validate_err_msg or error_message}: {str(e)}") from e
        except Exception as e:
            raise Exception(f"{error_message}: {str(e)}") from e

    def run(self, kwargs: dict, validate_err_msg: str, runtime_err_msg: str) -> Response | tuple[str, int]:
        try:
            kwargs = self._get_args(request, uuid.uuid4())
            modified_kwargs = self._validate_args(kwargs, validate_err_msg)
            result = self._run_program(kwargs | modified_kwargs, runtime_err_msg, validate_err_msg=validate_err_msg)
            return self._get_response(result, **kwargs)
        except Exception as e:
            return send_error_response(e, **kwargs)
        finally:
            if self._run_finally is not None: self._run_finally(**kwargs)

    def __init__ (self, name: str, get_args: Callable, validate_args: Callable, run_program: Callable, get_response: Callable = None, run_finally: Callable = None):
        """
        Create a Program object
        :param get_args: Get all the args needed for validating and running the program. If some are not needed in validate or run_program,
        add **ignore to the function signature
        :param validate_args: the function ran to make sure that all the arguments are valid
        :param run_program: run the program itself, returning a result
        :param get_response: create a valid flask response from the program's response
        :param run_finally: an optional argument to run any needed cleanup after running the program (e.g. file deletion)
        """
        self.name = name
        self._get_args_raw = get_args
        self._validate_args_raw = validate_args
        self._run_program_raw = run_program
        self._run_finally = run_finally
        self._get_response = get_response or self._get_response_default #temporary, while there still are temp methods

    def _get_response_default(self, result: ProgramObject, **kwargs) -> Response:
        zip_bytes, zipfile_name = result.to_zip(f"{self.name}ResultsFor-[fname].zip")
        zip_file_encoded = base64.b64encode(zip_bytes).decode("utf-8")
        return jsonify(zip=zip_file_encoded, html=render_template(
            'request-completed.html', program=self.name, filename=zipfile_name))

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

def get_result_temp(program: str, result: str):
    zip_file = get_zip_bytes("tempfile", **{program: (".txt", result)})
    zip_file_encoded = base64.b64encode(zip_file).decode("utf-8")
    filename = f"TempResult-{program}.zip"
    return jsonify(zip=zip_file_encoded, html=render_template(
        'request-completed.html', result=result, program = program, filename=filename))

program_dict = { #get args, validate args, return value
    'tfofinder': Program("TFOFinder",
                        lambda req, _: {"filein": req.files.get("ct-file").stream,
                                      "probe_lengths": req.form.get("tfofinder-probe-length"), #validation is in validate_arguments
                                      "filename": secure_filename(req.files.get("ct-file").filename),
                                         "arguments": tfofinder.parse_arguments("", from_command_line=False)},
                tfofinder.validate_arguments,
                  partial(close_file, tfofinder.calculate_result)), #temp functions
    'pinmol': Program("PinMol", lambda req, id: {"filein": req.files.get("ct-file").stream,
                                   "probe_length": req.form.get("pinmol-probe-length", type=int),
                                   "output_dir": pinmol_output_dir / str(id),
                                   "filename": secure_filename(req.files.get("ct-file").filename),
                                   "probe_count_max": optional_argument(req, "pinmol-probe-count-max", default_value=50, type=int),
                                    "arguments": pinmol.parse_arguments(f"-nb -w"
                                                                        f"{optional_argument(req, 'pinmol-start-base', '-s', default_value=1)}"
                                                                        f"{optional_argument(req, 'pinmol-end-base', '-e', default_value=-1)}",
                                                                        from_command_line=False)
    }, pinmol.validate_arguments, partial(close_file, pinmol.calculate_result), run_finally=delete_folder_tree),
    'smfish': Program("smFISH", lambda x: {"x": "test"}, lambda x: True, lambda x: x, partial(get_result_temp, "smFISH"))
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

def get_exports():
    return {"TFOFinder": tfofinder.exported_values, "PinMol": pinmol.exported_values}

def run_program(prog_name, error_message_validation="Something went wrong", error_message_program="Something went wrong"):
    # import time
    # prev = time.time_ns()
    program = program_dict[prog_name.lower()]
    result =  program.run(request, error_message_validation, error_message_program)
    # print(((time.time_ns() - prev) // 1_000) / 1_000)
    return result

if __name__ == '__main__':
    app.run(debug= (True if len(argv) > 1 else False)) #use debug if any commmand line arguments are inputted
