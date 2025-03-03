import pathlib
import jinja2

# define path to stuff
script_dirpath = pathlib.Path(__file__).parent
dat_filepath = script_dirpath.parent.parent / "chem" / "AD4_parameters.dat"
template_dirpath = script_dirpath / "templates"
header_source_filepath = (
    script_dirpath.parent.parent
    / "mudock"
    / "include"
    / "mudock"
    / "chem"
    / "autodock_parameters.hpp"
)

# initialize the jinja environement
env = jinja2.Environment(
    loader=jinja2.FileSystemLoader([str(template_dirpath)]),
    variable_start_string="{@",
    variable_end_string="@}",
)

# load the knowledge
with open(dat_filepath, "r") as f:
    lines = f.readlines()
    data = {}
    for line in lines:
        if "FE_coeff" in line and not "#" in line:
            parts = line.split()
            data[parts[0]] = parts[1]

# render the elements sources
template = env.get_template("autodock_parameters.template.hpp")
with open(header_source_filepath, "w") as outfile:
    outfile.write(template.render(data))