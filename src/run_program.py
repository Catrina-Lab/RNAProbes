import sys

from src.PinMol import pinmol
from src.TFOFinder import tfofinder
from src.smFISH import smFISH

programs = {
    "tfofinder": tfofinder.run,
    "pinmol": pinmol.run,
    "smfish": smFISH.run
}
if __name__ == "__main__":
    assert len(sys.argv) >= 2, "You must choose a program: either tfofinder, smFISH, or PinMol (case insensitive)"
    run = programs[sys.argv[1].lower()]
    run(sys.argv[2:])