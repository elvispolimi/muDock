import json
import pathlib
import re


LINE_RE = re.compile(
    r"^\s*(?P<index>\d+)\s+"
    r"(?P<atom_type>\S+)\s+"
    r"(?P<atomic_weight>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<vdw_radius>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<vdw_potential>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<par_charge>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<hbond>\S+)"
    r"(?:\s+.*)?$"
)


def to_enum_name(atom_type: str) -> str:
    return atom_type.replace(".", "").replace("+", "plus")

# define path to stuff
script_dirpath = pathlib.Path(__file__).parent
dat_filepath = script_dirpath.parent.parent / "chem" / "ATOM_DEF_XTOOL.dat"
with open(dat_filepath, "r") as f:
    lines = f.readlines()
    data = []
    for line in lines:
        match = LINE_RE.match(line)
        if not match:
            continue

        atom_type = match.group("atom_type")
        data.append(
            {
                "value": to_enum_name(atom_type),
                "name": to_enum_name(atom_type),
                "atomic_weight": match.group("atomic_weight"),
                "vdw_radius": match.group("vdw_radius"),
                "vdw_potential": match.group("vdw_potential"),
                "par_charge": match.group("par_charge"),
                "hbond": match.group("hbond"),
            }
        )

    print(json.dumps(data, indent=4))
