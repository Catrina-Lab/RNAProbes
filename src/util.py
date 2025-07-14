from pathlib import Path
import argparse

def input_int(msg : str, int_predicate = lambda n: True, fail_message : str = None) -> int:
    fail_message = fail_message or "Please input an integer"
    while True:
        try:
            x = int(input(msg))
            if not int_predicate(x):
                raise ValueError()
            return x
        except (ValueError, EOFError):
            print(fail_message)

def input_int_in_range(min: int = None, max: int = None, msg : str = None, fail_message : str = None) -> int:
    fail_message = fail_message or f"Please input an integer between {min} (inclusive) and {max} (exclusive)."
    msg = msg or f"Input an integer between {min} (inclusive) and {max} (exclusive): "
    return input_int(msg, int_predicate = lambda x: (min is None or x >= min) and (max is None or x < max), fail_message = fail_message)

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