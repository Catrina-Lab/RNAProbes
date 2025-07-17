from pathlib import Path
from collections.abc import Callable
import argparse
import os

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

def remove_if_exists(path: Path):
    if path.exists():
        os.remove(path)


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

class ArgparseClassBuilder:
    def __init__(self, parse, matches=lambda x: True, to_return=lambda x: x):
        self.to_return = to_return
        self.parse_msg = None
        self.msg = ""
        self.parse = parse
        self.matches = matches

    def set_default_match_fail_msg(self, msg):
        self.msg = msg
        return self
    
    def set_default_parse_fail_msg(self, msg):
        self.parse_msg = msg
        return self

    def _parse_msg(self, msg, *args, **kwargs):
        if type(msg) == str:
            return msg
        elif type(msg) == Callable: #fails for now
            return msg(*args, **kwargs)
        else:
            raise TypeError("Incompatible message type. Only allowed strings or functions, given " + str(type(msg)))


    def __call__(self, *args, err_msg=None, parse_msg=None, **kwargs):
        def instance_class(*cmdline_args):
            try:
                to_check = self.parse(*cmdline_args)
                if self.matches(to_check, *args, **kwargs):  # in returns true if string ==
                    return self.to_return(to_check)
                else:
                    raise argparse.ArgumentTypeError(self._parse_msg(err_msg or self.msg, *args, **kwargs))
            except argparse.ArgumentTypeError:
                raise
            except Exception as e:
                raise argparse.ArgumentTypeError(self._parse_msg(parse_msg or self.parse_msg, *args, **kwargs), e)
        return instance_class

range_class = (ArgparseClassBuilder(int, lambda n, min=0, max=100: min <= n <= max)
               .set_default_match_fail_msg(lambda min=0, max=100:f"Must be between {min} and {max}")
               .set_default_parse_fail_msg(lambda min=0, max=100:f"Must be a valid integer between {min} and {max}"))
fiveToTen = range_class(min=5, max = 10)

if __name__ == "__main__":
    import sys
    fiveToTen(6)