
---

# RNAProbes

**RNAProbes** is set of powerful tools for RNA detection.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)  
  - [Using pip](#using-pip)  
  - [Using Poetry](#Using-Poetry-recommended-for-linux)
  - [Mac Additions](#Mac-Additions)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

---

## Overview
RNAProbes consists of 3 programs: smFISH, TFOFinder, and PinMol.
 - TFOFinder: Designs triplex-forming oligonucleotides (TFO) probes.
 - smFISH: Designs molecular beacons for live cell imaging of mRNA targets.
 - PinMol: For designing probes for single-molecule fluorescence in situ hybridization (smFISH), taking into consideration target structure.

The repository contains two ways to run the programs:
1. A webapp located in app.py
2. CLI tools located in src/ (referenced by app.py)
   1. You can also find these in a read-only repository [here](https://github.com/Catrina-Lab/RNAProbesSource), without the webapp

## Installation

### Clone the repository
```commandline
git clone https://github.com/Catrina-Lab/RNAProbes.git\
cd RNAProbes
```

### Using pip

#### Windows, macOS, or linux with venv with `pip`:
```commandline
pip install -r requirements.txt
````

### Using Poetry (recommended for linux)

Ensure you have [Poetry](https://python-poetry.org/docs/#installation) installed (make sure to follow the instructions there to install poetry >= 2.0).

#### Clone the repo and install dependencies:

```commandline
poetry install
```

---

### Mac Additions
We use [RNAStructure text interface](http://rna.urmc.rochester.edu/RNAstructure.html) (slightly improved) by Dave Mathews at the University of Rochester.
This is included for Windows and Linux (64 bit), but must be installed for Mac (for now).

After installing, add the following lines to the .commandline_profile or .profile file:

```commandline
export PATH="[RNAstructure_directory]/RNAstructure/exe:$PATH"
export DATAPATH="[RNAstructure_path]/RNAstructure/data_tables/"
```

## Usage

Once installed, you can run the main program:

```commandline
python -m flask run
```

Or, if using Poetry:

```commandline
poetry run python -m flask run
```

---
Run the CLI tools like so:

```commandline
cd src
python -m rnaprobes
```
Or, if using poetry:
```commandline
cd src
poetry run python -m rnaprobes
```

## Contributing

We welcome any contributions. To contribute:
1. Fork the [RNAProbes repository](https://github.com/Catrina-Lab/RNAProbes)
2. Clone your fork and make desired changes.
3. Submit a pull request against the "main" branch. 
   1. If it is accepted and modifies the RNAProbes source code, it will automatically be changed in that repository also.

Pull requests should include:

* Clear description of changes
* Updated tests (if applicable)
* Documentation updates (as needed)
---

## License

GNU GPL v2

---
