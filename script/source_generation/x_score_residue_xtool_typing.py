import json
import pathlib
import re


RESI_RE = re.compile(
    r"^\s*RESI\s+(?P<name>\S+)\s+(?P<num_atoms>\S+)\s+(?P<num_bonds>\S+)\s+(?P<total_charge>\S+)\s+!\s*(?P<description>.*?)\s*$"
)
ATOM_RE = re.compile(
    r"^\s*ATOM\s+(?P<atom_name>\S+)\s+"
    r"(?P<basic_atom_type>\S+)\s+"
    r"(?P<x_tool_atom_type>\S+)\s+"
    r"(?P<vdw_radius>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<vdw_potential>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<partial_charge>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<hbond>\S+)\s+"
    r"(?P<hydrophobic_scale>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<sas_parameter>[+-]?\d+(?:\.\d+)?)\s+"
    r"(?P<ring_indicator>\d+)\s+"
    r"(?P<pmf_atom_type>\S+)\s*$"
)
BOND_RE = re.compile(
    r"^\s*BOND\s+(?P<atom_1>\S+)\s+(?P<atom_2>\S+)\s+(?P<bond_type>\S+)\s*$"
)


def to_enum_name(text: str) -> str:
    name = text
    name = re.sub(r"\+", "plus", name)
    name = re.sub(r"-", "minus", name)
    name = re.sub(r"=", "eq", name)
    name = re.sub(r">", "gt", name)
    name = re.sub(r"<", "lt", name)
    name = re.sub(r"\.", "", name)
    name = re.sub(r"[()]+", "_", name)
    name = re.sub(r"[^A-Za-z0-9_]+", "_", name)
    name = re.sub(r"_+", "_", name)
    return name.strip("_")


def parse_residue_definition(file_path: pathlib.Path) -> dict:
    residues = []
    current_residue = None

    with open(file_path, "r") as input_file:
        for line in input_file:
            residue_match = RESI_RE.match(line)
            if residue_match:
                if current_residue is not None:
                    residues.append(current_residue)

                residue_name = residue_match.group("name")
                current_residue = {
                    "value": to_enum_name(residue_name),
                    "name": to_enum_name(residue_name),
                    "num_atoms": residue_match.group("num_atoms"),
                    "num_bonds": residue_match.group("num_bonds"),
                    "total_charge": residue_match.group("total_charge"),
                    "description": residue_match.group("description").strip(),
                    "atoms": [],
                    "bonds": [],
                }
                continue

            if current_residue is None:
                continue

            atom_match = ATOM_RE.match(line)
            if atom_match:
                atom_name = atom_match.group("atom_name")
                current_residue["atoms"].append(
                    {
                        "value": to_enum_name(atom_name),
                        "name": to_enum_name(atom_name),
                        "basic_atom_type": to_enum_name(atom_match.group("basic_atom_type")),
                        "x_tool_atom_type": to_enum_name(atom_match.group("x_tool_atom_type")),
                        "vdw_radius": atom_match.group("vdw_radius"),
                        "vdw_potential": atom_match.group("vdw_potential"),
                        "par_charge": atom_match.group("partial_charge"),
                        "hbond": to_enum_name(atom_match.group("hbond")),
                        "hydrophobic_scale": atom_match.group("hydrophobic_scale"),
                        "sas_parameter": atom_match.group("sas_parameter"),
                        "ring_indicator": atom_match.group("ring_indicator"),
                        "pmf_atom_type": to_enum_name(atom_match.group("pmf_atom_type")),
                    }
                )
                continue

            bond_match = BOND_RE.match(line)
            if bond_match:
                current_residue["bonds"].append(
                    {
                        "atom_1": bond_match.group("atom_1"),
                        "atom_2": bond_match.group("atom_2"),
                        "bond_type": bond_match.group("bond_type"),
                    }
                )

    if current_residue is not None:
        residues.append(current_residue)

    return {"residue_type": residues}


def main() -> None:
    script_dirpath = pathlib.Path(__file__).parent
    dat_filepath = script_dirpath.parent.parent / "chem" / "RESIDUE_DEF_XTOOL.dat"
    json_filepath = script_dirpath.parent.parent / "chem" / "x_score_residue_xtool_types.json"

    data = parse_residue_definition(dat_filepath)

    with open(json_filepath, "w") as output_file:
        json.dump(data, output_file, indent=4)
        output_file.write("\n")


if __name__ == "__main__":
    main()