from __future__ import annotations
import sys

from src.PinMol import pinmol
from src.RNASuiteUtil import run_command_line
from src.TFOFinder import tfofinder
from src.smFISH import smFISH
from src.util import input_value

programs = {
    "tfofinder": tfofinder.run,
    "pinmol": pinmol.run,
    "smfish": smFISH.run
}
if __name__ == "__main__":
    program = input_value("Input a program (either tfofinder, pinmol, or smfish): ", str.lower,
                            lambda program: program in programs.keys(), retry_if_fail=len(sys.argv) >= 2,
                          initial_value=sys.argv[1].lower() if len(sys.argv) >= 2 else None)
    run = programs[program]
    run_command_line(run, sys.argv[2:])