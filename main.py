from flask import Flask, render_template, request, redirect, url_for, jsonify
from sys import argv
from src.TFOFinder import tfofinder
from werkzeug.utils import secure_filename
from functools import partial
import time, io, zipfile, base64
from pathlib import Path

app = Flask(__name__)

class Program:
    def _validate_args(self, kwargs, error_message="Something went wrong"):
        try:
            return self._validate_args_raw(**kwargs)
        except TypeError as e:
            return TypeError(f"{error_message}: {str(e)}")

    def _run_program(self, kwargs, error_message="Something went wrong"):
        try:
            return self._run_program_raw(**kwargs)
        except Exception as e:
            return Exception(f"{error_message}: {str(e)}")

    def run(self, kwargs, validate_err_msg, runtime_err_msg):
        self._validate_args(kwargs, validate_err_msg)
        result =  self._run_program(kwargs, runtime_err_msg)
        return self._get_response(result, **kwargs)

    def __init__ (self, get_args, validate_args, run_program, get_response):
        self.get_args = get_args
        self._validate_args_raw = validate_args
        self._run_program_raw = run_program
        self._get_response = get_response

def close_file(func, *, file_index = "filein", **kwargs):
    """
    Execute a function and close the file afterward.
    :param func: the function to execute
    :param args: the arguments to pass to that function.
    :param file_index: the index of the argument representing the file to be closed
    :return:
    """
    with kwargs[file_index]:
        return func(**kwargs)

def get_zip(**kwargs):
    file_obj = io.BytesIO()
    with zipfile.ZipFile(file_obj, 'w', compression=zipfile.ZIP_DEFLATED) as zip_file:
        for name, value in kwargs.items():
            zip_info = zipfile.ZipInfo(name)
            zip_info.date_time = time.localtime(time.time())[:6]
            zip_file.writestr(name + value[0], value[1])
    file_obj.seek(0)
    return file_obj

def get_tfofinder_result(result, **kwargs):
    zip_file = get_zip(TFOProbes = (".txt", result))
    zip_file_encoded = base64.b64encode(zip_file.getvalue()).decode("utf-8")
    filename = Path(kwargs["filename"]).stem
    filename = f"TFOFinderResultsFor-{filename}.zip"
    return jsonify(zip=zip_file_encoded, html=render_template(
        'request-completed.html', result=result, program = "TFOFinder", filename=filename))

def get_result_temp(program, result, **kwargs):
    zip_file = get_zip(**{program: (".txt", result)})
    zip_file_encoded = base64.b64encode(zip_file.getvalue()).decode("utf-8")
    filename = f"TempResult-{program}.zip"
    return jsonify(zip=zip_file_encoded, html=render_template(
        'request-completed.html', result=result, program = program, filename=filename))

program_dict = { #get args, validate args, return value
    'tfofinder': Program(lambda req: {"filein": req.files.get("ct-file").stream,
                                      "probe_length": req.form.get("tfofinder-probe-length", type=int),
                                      "filename": secure_filename(req.files.get("ct-file").filename)},
                tfofinder.validate_arguments,
                  partial(close_file, tfofinder.calculate_result),
                         get_tfofinder_result), #temp functions
    'smfish': Program(lambda x: {"x": "test"}, lambda x: True, lambda x: x, partial(get_result_temp, "smFISH")),
    'pinmol': Program(lambda x: {"x": "test"}, lambda x: True, lambda x: x, partial(get_result_temp, "PinMol"))
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
        return render_template('request.html', programs=programs)
    # If accessed via GET, redirect to home
    return redirect(url_for('index'))

@app.route('/send-request', methods=['POST'])
def send_request():
    program = request.args.get("program")
    return run_program(program, "The given arguments are invalid",
                         "Something went wrong when calculating your result")


def run_program(prog_name, error_message_validation="Something went wrong", error_message_program="Something went wrong"):
    program = program_dict[prog_name.lower()]
    kwargs = program.get_args(request)
    return program.run(kwargs, error_message_validation, error_message_program)

if __name__ == '__main__':
    app.run(debug= (True if len(argv) > 1 else False)) #use debug if any commmand line arguments are inputted
