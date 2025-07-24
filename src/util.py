import io
import zipfile
from collections import namedtuple
from collections.abc import Callable
from pathlib import Path
import argparse
import os

class ValidationError(Exception):
    pass

def input_int(msg : str, int_predicate = lambda n: True, fail_message : str = None, initial_value = None, retry_if_fail=False) -> int:
    fail_message = fail_message or "Please input an integer"
    while True:
        try:
            x = initial_value or int(input(msg))
            initial_value = None
            if not int_predicate(x):
                raise ValueError()
            return x
        except (ValueError, EOFError) as e:
            if retry_if_fail: print(fail_message)
            else: raise ValueError(fail_message) from e

def input_int_in_range(min: int = None, max: int = None, msg : str = None, fail_message : str = None,
                       initial_value = None, extra_predicate=lambda x: False, retry_if_fail=False) -> int:
    fail_message = fail_message or f"Please input an integer between {min} (inclusive) and {max} (exclusive)."
    msg = msg or f"Input an integer between {min} (inclusive) and {max} (exclusive): "
    return input_int(msg, initial_value=initial_value, int_predicate = lambda x: ((min is None or x >= min) and (max is None or x < max)) or extra_predicate(x),
                     fail_message = fail_message, retry_if_fail=retry_if_fail)

def remove_files(*files):
    for file in files:
        remove_if_exists(file)

def remove_if_exists(path: Path):
    if path.exists():
        os.remove(path)

def validate_arg(boolean: bool, msg):
    if not boolean: raise ValidationError(msg)

def validate_range_arg(input: int, min = None, max = None, name = "input", overwrite_msg = None, extra_predicate = lambda x: False):
    """
    Validate an argument that conforms to a range
    :param input:
    :param min: the max value of the argument, exclusive
    :param max:
    :param msg:
    :return:
    """
    validate_arg(type(input) is int, overwrite_msg or f"The given {name} is not an integer! It must be an integer between {min or "negative Infinity"} inclusive and {max or "Infinity"} exclusive")
    min_OK = (min is None or input >= min)
    max_OK = (max is None or input < max)
    validate_arg((max_OK and min_OK) or extra_predicate(input), overwrite_msg or f"The given {name} ({input}) is {"too small" if max_OK else "too large"}! It must be an integer between {min or "negative Infinity"} inclusive and {max or "Infinity"} exclusive")

def optional_argument(request, arg_name: str, cmd_line_name: str = None, default_value = None, type: Callable =None):
    if cmd_line_name is not None:
        if request.form[arg_name] != "":
            return f" {cmd_line_name} {request.form[arg_name]}"
        else:
            return f" {cmd_line_name} {default_value}" if default_value is not None else ""
    elif default_value is not None:
        return default_value if request.form[arg_name] == "" else type(request.form[arg_name])
    else:
        raise ValueError("Must include either cmd_line_name or default_value (or both)")

def get_or_none(obj, attr):
    return getattr(obj, attr, None)
def range_type(string, min=0, max=1000):
       value = int(string)
       if min <= value <= max:
           return value
       else:
           raise argparse.ArgumentTypeError(
               f'value not in range {min}-{max}. Please either keep it in range or leave it out.')

def path_string(string, suffix=".ct"):
    path = Path(string).resolve()
    if path.is_file() and path.suffix in suffix:  # in returns true if string ==
        return str(path)
    else:
        raise argparse.ArgumentTypeError(f'Invalid file given. File must be an existing {suffix} file')

def path_arg(string, suffix=".ct"):
    path = Path(string).resolve()
    if path.is_file() and path.suffix in suffix:  # in returns true if string ==
        return path
    else:
        raise argparse.ArgumentTypeError(f'Invalid file given. File must be an existing {suffix} file')

def get_folder_as_zip(folder_path: Path) -> bytes:
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                abs_path = os.path.join(root, file)
                rel_path = os.path.relpath(abs_path, start=folder_path)
                zipf.write(abs_path, arcname=rel_path)
    zip_buffer.seek(0)
    return zip_buffer.getvalue()

ParsedFile = namedtuple("ParsedFile", ["parent", "stem", "suffix"])
def parse_file_input(filein, output_dir = None) -> ParsedFile:
    filein_path = Path(filein).resolve()  # Convert the filein to a Path object and resolve its full path
    mb_userpath = output_dir or filein_path.parent  # Use the parent directory of the input file to save all files
    fname = filein_path.stem
    return ParsedFile(mb_userpath, fname, filein_path.suffix)


if __name__ == "__main__":
    print("Debug")