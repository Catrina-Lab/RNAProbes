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
