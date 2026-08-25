import json
import pathlib
import re


LINE_RE = re.compile(
    r"^\s*(?P<index>\d+)\s+(?P<atom_type>\S+)\s+(?P<hbond>\S+)\s+(?P<hydrophobic_scale>[+-]?\d+(?:\.\d+)?)\s*$"
)


def to_enum_name(atom_type: str) -> str:
    name = atom_type
    name = re.sub(r"=", "_eq_", name)
    name = re.sub(r">", "_gt_", name)
    name = re.sub(r"<", "_lt_", name)
    name = re.sub(r"\.", "_", name)
    name = re.sub(r"[()]+", "_", name)
    name = re.sub(r"[^A-Za-z0-9_]+", "_", name)
    name = re.sub(r"_+", "_", name)
    return name.strip("_")

# define path to stuff
script_dirpath = pathlib.Path(__file__).parent
dat_filepath = script_dirpath.parent.parent / "chem" / "ATOM_DEF_XLOGP.dat"
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
                "hbond": match.group("hbond"),
                "hydrophobic_scale": match.group("hydrophobic_scale"),
            }
        )

    print(json.dumps(data, indent=4))
