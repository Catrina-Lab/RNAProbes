from __future__ import annotations
import base64
import json
import threading
import uuid
from collections.abc import Callable
from copyreg import constructor
from keyword import kwlist
from pathlib import Path
from uuid import UUID

from flask import Response, jsonify, render_template, abort

from ..rnaprobes.RNAProbesUtil import ProgramObject
from ..rnaprobes.util import safe_remove_tree, is_empty
from ..rnaprobes.util import ValidationError
import traceback

IS_DELAYED = "_delayed_"

def send_error_response(error: BaseException, **kwargs):
    print(traceback.print_exc())
    if isinstance(error, ValidationError):
        return str(error), 400
    elif isinstance(error, Exception):
        return str(error), 500
    else:
        abort(500)

class Program:
    def _get_args(self, request, output_dir: Path, error_message: str="Something went wrong") -> tuple[dict, bool]:
        try:
            kwargs : dict = self._get_args_raw(request, output_dir)
            is_delayed = IS_DELAYED in kwargs
            kwargs.pop(IS_DELAYED, None)
            return kwargs, is_delayed
        except Exception as e:
            raise ValidationError(f"{error_message}: {str(e)}") from e

    def _validate_args(self, kwargs: dict, error_message: str="Something went wrong") -> dict:
        try:
            return kwargs | self._validate_args_raw(**kwargs) or dict()
        except Exception as e:
            raise ValidationError(f"{error_message}: {str(e)}") from e

    def _run_program(self, kwargs: dict, job_id: UUID, error_message: str= "Something went wrong", validate_err_msg: str=None, is_delayed: bool = False):
        try:
            return self._run_program_raw(**kwargs)
        except ValidationError as e:
            raise ValidationError(f"{validate_err_msg or error_message}: {str(e)}") from e
        except Exception as e:
            raise Exception(f"{error_message}: {str(e)}") from e


    def _get_response(self, result: ProgramObject, job_id, **kwargs) -> dict:
        zip_bytes, zipfile_name = self.get_response_bytes(result, **kwargs)
        return self._get_response_from_zip(job_id, zip_bytes, zipfile_name)

    def remove_directory(self, output_dir):
        safe_remove_tree(output_dir, self.root_dir)

    def _get_program_object(self, is_delayed, kwargs: dict, job_id: UUID, output_dir: Path, validate_err_msg: str, runtime_err_msg: str):
        constructor = DelayedRunnableProgram if is_delayed else RunnableProgram
        return constructor(self, kwargs, job_id, output_dir, validate_err_msg, runtime_err_msg)

    def run(self, request, validate_err_msg: str, runtime_err_msg: str) -> dict | tuple[str, int]:
        kwargs = dict()
        output_dir, runnable = None, None
        try:
            job_id, output_dir = self.set_id()
            kwargs, is_delayed = self._get_args(request, output_dir)
            kwargs = self._validate_args(kwargs, validate_err_msg) #join the result with kwargs
            runnable = self._get_program_object(is_delayed, kwargs, job_id, output_dir, validate_err_msg, runtime_err_msg)
            return runnable.run(self._run_program, self._get_response)
        except BaseException as e:
            return send_error_response(e, **kwargs)
        finally:
            if not (output_dir is None or isinstance(runnable, DelayedRunnableProgram)): self.remove_directory(output_dir)

    def __init__ (self, name: str, get_args: Callable, validate_args: Callable, run_program: Callable, root_dir, output_dir : Path = None, folder_not_needed = False):
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
        self.root_dir = root_dir
        self.output_dir = output_dir
        self.folder_not_needed = folder_not_needed
        self.get_extra_notification = lambda x: ""

    def set_id(self) -> tuple[UUID, Path]:
        job_id = uuid.uuid4()
        output_dir = None
        if not self.folder_not_needed:
            output_dir = self.output_dir / str(job_id)
            output_dir.mkdir(parents=True)
        return job_id, output_dir

    def get_response_bytes(self, result: ProgramObject, **kwargs) -> tuple[bytes, str]:
        return result.to_zip(f"{self.name}ResultsFor-[fname].zip")

    def _get_response_from_zip(self, job_id: UUID, zip_bytes: bytes, zipfile_name: str):
        zip_file_encoded = base64.b64encode(zip_bytes).decode("utf-8")
        return dict(zip=zip_file_encoded, status="complete", html=render_template(
            'request-results/request-completed.html', program=self.name,
            filename=zipfile_name, id=str(job_id)))
    def set_extra_notification_string_callback(self, func: Callable[dict, str]):
        self.get_extra_notification = func
        return self

    def get_delayed_response(self, output_dir: Path, job_id: UUID):
        return DelayedRunnableProgram.get_current_result(self, self._get_response_from_zip, output_dir, job_id)

class RunnableProgram:
    def __init__(self, program: Program, kwargs: dict, job_id: UUID, output_dir: Path, validate_err_msg: str, runtime_err_msg: str):
        self.program = program
        self.kwargs = kwargs
        self.job_id = job_id
        self.output_dir = output_dir
        self.runtime_err_msg = runtime_err_msg
        self.validate_err_msg=validate_err_msg

    def run(self, run_program, get_response):
        result = run_program(self.kwargs, self.job_id, error_message=self.runtime_err_msg, validate_err_msg=self.validate_err_msg)
        return get_response(result, self.job_id, **self.kwargs)


result_dir_name = "program-result"
error_file_name = "error.json"

class DelayedRunnableProgram(RunnableProgram):
    def run(self, run_program, get_response):
        threading.Thread(target=self._run_in_background, args=(run_program, get_response)).start()

        return self._get_running_response()

    def _get_running_response(self):
        return dict(id=str(self.job_id), status="running", html=render_template(
            'request-results/request-received.html', program=self.program.name, delayed=True, extra_notification = self.program.get_extra_notification(self.kwargs))), 202

    def _run_in_background(self, run_program, get_response):
        output_dir = self.output_dir / result_dir_name
        output_dir.mkdir(exist_ok=False, parents=True) #if it exists, we have a collision
        try:
            result = run_program(self.kwargs, self.job_id, error_message=self.runtime_err_msg, validate_err_msg=self.validate_err_msg)
            bytes, name = self.program.get_response_bytes(result)
            with open(output_dir / name, "wb") as zip_file:
                zip_file.write(bytes)
            result.cleanup()
        except BaseException as e:
            write_error_response(e, output_dir / error_file_name, **self.kwargs)

    @staticmethod
    def get_current_result(program: Program, get_response_from_zip, output_dir: Path, id: UUID):
        if not output_dir.exists():
            return 400, "Your output has been deleted from the system (or was never ran). Please run it again"
        if not (output_dir / result_dir_name).exists() or is_empty(output_dir / result_dir_name):
            return json.dumps(dict(status="running", id=str(id))), 202

        return DelayedRunnableProgram.send_final_result(program, get_response_from_zip, output_dir, id)

    @staticmethod
    def send_final_result(program: Program, get_response_from_zip, output_dir: Path, id: UUID):
        result_dir = output_dir / result_dir_name
        try:
            if (result_dir / error_file_name).exists():
                return send_error_from_file(result_dir / error_file_name)
            else:
                path = next(result_dir.iterdir())
                with open(path, "rb") as file:
                    return get_response_from_zip(id, file.read(), path.name)
        finally:
            program.remove_directory(output_dir)

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