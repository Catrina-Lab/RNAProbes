from __future__ import annotations

import threading
from collections.abc import Callable
from pathlib import Path
from uuid import UUID

from flask import render_template
from .Program import Program
import json

from ..rnaprobes.util import is_empty
from ..rnaprobes.util import ValidationError

result_dir_name = "program-result"
error_file = "error.json"

def write_error_response(error: BaseException, file: Path, **kwargs):
    to_serialize = dict(message = str(error) if isinstance(error, Exception) else "")
    to_serialize["code"] = 400 if isinstance(error, ValidationError) else 500
    msg = json.dumps(to_serialize)
    with open(file, "w") as f:
        f.write(msg)

def send_error_from_file(path: Path):
    with open(path, "r") as file:
        data = json.load(file)
        return data['message'], data['code'] if data['message'] != "" else data['code']

class DelayedProgram(Program):
    get_extra_notification = lambda x: ""
    #result = self._run_program(kwargs | modified_kwargs, runtime_err_msg, validate_err_msg=validate_err_msg)
            #return self._get_response(result, **kwargs)
    def _remove_directory(self, output_dir):
        pass

    def _run_program(self, kwargs: dict, job_id: UUID, error_message: str="Something went wrong", validate_err_msg: str=None):
        threading.Thread(target=self._run_in_background, args=(kwargs, job_id, error_message, validate_err_msg)).start()
        # Path(self.output_dir / str(job_id) / result_dir_name).mkdir(parents=True, exist_ok=True)
        return job_id

    def _run_in_background(self, kwargs, job_id: UUID, error_message: str="Something went wrong", validate_err_msg: str=None):
        output_dir = self.output_dir / str(job_id) / result_dir_name
        output_dir.mkdir(exist_ok=False, parents=True) #if it exists, we have a collision
        try:
            result = super()._run_program(kwargs, job_id, error_message, validate_err_msg)
            bytes, name = self._result_bytes(result)
            with open(output_dir / name, "wb") as zip_file:
                zip_file.write(bytes)
            result.cleanup()
        except BaseException as e:
            write_error_response(e, output_dir / error_file, **kwargs)

    def _get_response(self, _, job_id: UUID, **kwargs) -> dict:
        return dict(id=str(job_id), status="running", html=render_template(
            'request-results/request-received.html', program=self.name, delayed=True, extra_notification = self.get_extra_notification(kwargs)))

    def get_current_result(self, output_dir: Path, id: UUID):
        if not output_dir.exists():
            return 400, "Your output has been deleted from the system (or was never ran). Please run it again"
        if not (output_dir / result_dir_name).exists() or is_empty(output_dir / result_dir_name):
            return json.dumps(dict(status="running", id=str(id)))

        return self.send_final_result(output_dir, id)

    def send_final_result(self, output_dir: Path, id: UUID):
        result_dir = output_dir / result_dir_name
        try:
            if (result_dir / error_file).exists():
                return send_error_from_file(result_dir / error_file)
            else:
                path = next(result_dir.iterdir())
                with open(path, "rb") as file:
                    return self._get_response_from_zip(id, file.read(), path.name)
        finally:
            super()._remove_directory(output_dir)

    def set_extra_notification_string_callback(self, func: Callable[dict, str]):
        self.get_extra_notification = func
        return self