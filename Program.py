from __future__ import annotations
import base64
import uuid
from collections.abc import Callable
from pathlib import Path
from uuid import uuid4

from flask import Response, jsonify, render_template, request

from src.RNASuiteUtil import ProgramObject
from src.util import ValidationError
import traceback

def send_error_response(error: BaseException, **kwargs):
    print(traceback.print_exc())
    if isinstance(error, ValidationError):
        return str(error), 400
    else:
        return str(error), 500

class Program:
    def _get_args(self, request, output_dir: Path, error_message: str="Something went wrong"):
        try:
            return self._get_args_raw(request, output_dir)
        except Exception as e:
            raise ValidationError(f"{error_message}: {str(e)}") from e

    def _validate_args(self, kwargs: dict, error_message: str="Something went wrong") -> dict:
        try:
            return self._validate_args_raw(**kwargs) or dict()
        except Exception as e:
            raise ValidationError(f"{error_message}: {str(e)}") from e

    def _run_program(self, kwargs: dict, job_id: uuid4, error_message: str= "Something went wrong", validate_err_msg: str=None):
        try:
            return self._run_program_raw(**kwargs)
        except ValidationError as e:
            raise ValidationError(f"{validate_err_msg or error_message}: {str(e)}") from e
        except Exception as e:
            raise Exception(f"{error_message}: {str(e)}") from e

    def _get_response(self, result: ProgramObject, **kwargs) -> Response:
        zip_bytes, zipfile_name = self._result_bytes(result, **kwargs)
        return self._get_response_from_zip(zip_bytes, zipfile_name)

    def run(self, kwargs: dict, validate_err_msg: str, runtime_err_msg: str) -> Response | tuple[str, int]:
        try:
            job_id = uuid.uuid4()
            kwargs = self._get_args(request, self.output_dir / str(job_id))
            modified_kwargs = self._validate_args(kwargs, validate_err_msg)
            result = self._run_program(kwargs | modified_kwargs, job_id, runtime_err_msg, validate_err_msg=validate_err_msg)
            return self._get_response(result, **kwargs)
        except Exception as e:
            return send_error_response(e, **kwargs)
        finally:
            if self._run_finally is not None: self._run_finally(**kwargs)

    def __init__ (self, name: str, get_args: Callable, validate_args: Callable, run_program: Callable, get_response_bytes: Callable = None, run_finally: Callable = None, output_dir : Path = None):
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
        self._result_bytes = get_response_bytes or self._get_response_bytes #temporary, while there still are temp methods
        self.output_dir = output_dir

    def _get_response_bytes(self, result: ProgramObject, **kwargs) -> tuple[bytes, str]:
        return result.to_zip(f"{self.name}ResultsFor-[fname].zip")

    def _get_response_from_zip(self, zip_bytes: bytes, zipfile_name: str):
        zip_file_encoded = base64.b64encode(zip_bytes).decode("utf-8")
        return jsonify(zip=zip_file_encoded, status="complete", html=render_template(
            'request-completed.html', program=self.name, filename=zipfile_name))