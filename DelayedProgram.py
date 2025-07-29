from __future__ import annotations

import threading
from pathlib import Path
from uuid import uuid4

from flask import Response, jsonify, render_template, request
from Program import Program
import json

from src.RNASuiteUtil import ProgramObject
from src.run_program import programs
from src.util import ValidationError, is_empty

result_dir_name = "program-result"
error_file = "error.json"

def write_error_response(error: BaseException, file: Path, **kwargs):
    to_serialize = dict(message = str(error))
    to_serialize["code"] = 400 if isinstance(error, ValidationError) else 500
    msg = json.dumps(to_serialize)
    with open(file, "w") as f:
        f.write(msg)

def send_error_from_file(path: Path):
    with open(path, "r") as file:
        data = json.load(file)
        return data['message'], data['code']

class DelayedProgram(Program):
    #result = self._run_program(kwargs | modified_kwargs, runtime_err_msg, validate_err_msg=validate_err_msg)
            #return self._get_response(result, **kwargs)
    def _run_program(self, kwargs: dict, job_id: uuid4, error_message: str="Something went wrong", validate_err_msg: str=None):
        threading.Thread(target=self._run_in_background, args=(kwargs, job_id, error_message, validate_err_msg)).run()
        Path(self.output_dir / str(job_id) / result_dir_name).mkdir(parents=True, exist_ok=True)
        return job_id

    def _run_in_background(self, kwargs, job_id: uuid4, error_message: str="Something went wrong", validate_err_msg: str=None):
        output_dir = self.output_dir / str(job_id) / result_dir_name
        output_dir.mkdir(exist_ok=False, parents=True) #if it exists, we have a collision
        try:
            result = super()._run_program(kwargs, error_message, validate_err_msg)
            bytes, name = self._result_bytes(result)
            with open(output_dir / name, "wb") as zip_file:
                zip_file.write(bytes)
        except BaseException as e:
            write_error_response(e, output_dir / error_file, **kwargs)
        finally:
            if self._run_finally is not None: self._run_finally(**kwargs)

    def _get_response(self, id: uuid4, **kwargs) -> Response:
        return jsonify(id=str(id), status="running", html=render_template(
            'request-received.html', program=self.name, delayed=True))

    def get_current_result(self, output_dir: Path, id: uuid4):
        if not output_dir.exists():
            return 400, "Your output has been deleted from the system (or was never ran). Please run it again"
        if is_empty(output_dir):
            return json.dumps(dict(status="running", id=id))
        elif (output_dir / error_file).exists():
            return send_error_from_file(output_dir / error_file)
        else:
            path = next(output_dir.iterdir())
            with open(path, "rb") as file:
                return self._get_response_from_zip(file.read(), path.name)
