import json
import pathlib
import jinja2

# define path to stuff
script_dirpath = pathlib.Path(__file__).parent
json_filepath = script_dirpath.parent.parent / "chem" / "autodock_types.json"
template_dirpath = script_dirpath / "templates"
header_source_filepath = (
    script_dirpath.parent.parent
    / "mudock"
    / "include"
    / "mudock"
    / "chem"
    / "autodock_types.hpp"
)
source_source_filepath = (
    script_dirpath.parent.parent / "mudock" / "src" / "chem" / "autodock_types.cpp"
)


# initialize the jinja environement
env = jinja2.Environment(
    loader=jinja2.FileSystemLoader([str(template_dirpath)]),
    variable_start_string="{@",
    variable_end_string="@}",
)

# load the knowledge
with open(json_filepath, "r") as json_file:
    data = json.load(json_file)

# add the index for random access to the dictionary
for index, element in enumerate(data["ff_type"]):
    element["index"] = index

# render the elements sources
context = {"data": data["ff_type"], "num_elements": len(data["ff_type"])}
template = env.get_template("autodock_types.template.hpp")
with open(header_source_filepath, "w") as outfile:
    outfile.write(template.render(context))
template = env.get_template("autodock_types.template.cpp")
with open(source_source_filepath, "w") as outfile:
    outfile.write(template.render(context))
