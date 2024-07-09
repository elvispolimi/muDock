import json
import sys
import pathlib

# define path to stuff
script_dirpath = pathlib.Path(__file__).parent
dat_filepath = script_dirpath.parent.parent / "chem" / "AD4_parameters.dat"
with open(dat_filepath, "r") as f:
    lines = f.readlines()
    data = []
    for line in lines:
        if line.startswith("atom_par"):
            parts = line.split()
            data.append(
                {
                    "value": parts[1],
                    "name": parts[1],
                    "Rii": parts[2],
                    "epsii": parts[3],
                    "vol": parts[4],
                    "solpar": parts[5],
                    "Rij_hb": parts[6],
                    "epsij_hb": parts[7],
                    "hbond": parts[8],
                }
            )

    print(json.dumps(data, indent=4))
