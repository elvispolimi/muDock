import json
import pathlib
import sys
import jinja2

# the generated cpp has to be heavily modified
# - HET Bonds (No candidate for construction) -> empty
# - Xtool_residue --> residue (renaming)
# - moved enum map from hpp to cpp, and maybe we should move that map to a dedicated residue.cpp

def main() -> None:
    script_dirpath = pathlib.Path(__file__).parent
    root_dirpath = script_dirpath.parent.parent
    
    json_filepath = root_dirpath / "chem" / "x_score_residue_xtool_types.json"
    template_hpp_filepath = script_dirpath / "templates" / "x_score_residue_xtool_types.template.hpp"
    template_cpp_filepath = script_dirpath / "templates" / "x_score_residue_xtool_types.template.cpp"
    
    output_hpp_filepath = root_dirpath / "mudock" / "include" / "mudock" / "chem" / "x_score_residue_xtool_types.hpp"
    output_cpp_filepath = root_dirpath / "mudock" / "src" / "chem" / "x_score_residue_xtool_types.cpp"

    if not json_filepath.exists():
        print(f"Error: {json_filepath} not found.")
        sys.exit(1)

    with open(json_filepath, "r") as f:
        data = json.load(f)

    # Prepare data for templates
    # data["residue_type"] is already a list of residues with atoms and bonds
    # Map bond_type 'am' and 'ar' to integer constants
    for residue in data["residue_type"]:
        for bond in residue["bonds"]:
            if bond["bond_type"] == "am":
                bond["bond_type"] = "3"
            elif bond["bond_type"] == "ar":
                bond["bond_type"] = "4"

    context = {
        "residue_type": data["residue_type"],
        "num_residues": len(data["residue_type"])
    }

    # Setup Jinja2 environment
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(script_dirpath / "templates")),
        variable_start_string="{@",
        variable_end_string="@}",
    )

    # Render Header
    template_hpp = env.get_template("x_score_residue_xtool_types.template.hpp")
    with open(output_hpp_filepath, "w") as f:
        f.write(template_hpp.render(context))
    
    # Render Source
    template_cpp = env.get_template("x_score_residue_xtool_types.template.cpp")
    with open(output_cpp_filepath, "w") as f:
        f.write(template_cpp.render(context))

    print(f"Generated {output_hpp_filepath}")
    print(f"Generated {output_cpp_filepath}")

if __name__ == "__main__":
    main()
