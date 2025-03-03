import json
import sys
import pathlib

# define path to stuff
script_dirpath = pathlib.Path(__file__).parent
dat_filepath = script_dirpath.parent.parent / "chem" / "AD4_parameters.dat"
with open(dat_filepath, "r") as f:
    lines = f.readlines()
    data = {}
    for line in lines:
        if "FE_coeff" in line and not "#" in line:
            parts = line.split()
            data[parts[0]] = parts[1]

    print(json.dumps(data, indent=4))
