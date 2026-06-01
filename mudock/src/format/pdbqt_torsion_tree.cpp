#include <mudock/format/pdbqt_torsion_tree.hpp>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>

namespace mudock {

std::uint32_t atom_id_from_pdbqt_line(const std::string& line) {
  return static_cast<std::uint32_t>(std::stoul(line.substr(6, 5)));
}

std::size_t register_atom_id(std::unordered_map<std::uint32_t, std::size_t>& atom_id_to_index,
                             const std::uint32_t atom_id) {
  const auto index = atom_id_to_index.size();
  const auto [_, inserted] = atom_id_to_index.emplace(atom_id, index);
  if (!inserted) {
    throw std::runtime_error("Duplicate PDBQT atom id");
  }
  return index;
}

// Recursively parse branches and build the representation of the torsion tree.
pdbqt_torsion_branch parse_branch_tree(const std::vector<std::string>& lines,
                                       std::size_t& pos,
                                       const std::uint32_t from,
                                       const std::uint32_t to,
                                       std::vector<pdbqt_torsion_branch_node>& storage,
                                       std::vector<pdbqt_rotor>& rotors,
                                       std::unordered_map<std::uint32_t, std::size_t>& atom_id_to_index) {

  storage.emplace_back(); //Allocate space for the new node
  auto& node = storage.back();
  rotors.push_back(pdbqt_rotor{from, to});

  while (pos < lines.size()) {
    const auto& line = lines[pos];
    pos++;
    if (line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0) { //(line.starts_with("ATOM") || line.starts_with("HETATM"))
      const auto atom_id = atom_id_from_pdbqt_line(line);
      register_atom_id(atom_id_to_index, atom_id);
      node.atom_indices.push_back(atom_id);
    } else if (line.rfind("BRANCH", 0) == 0) {
      //We have found a new recursive branch , we need to parse:
      // Branch from(parent) to(child)
      std::istringstream stream{line};
      std::string token;
      std::uint32_t child_from;
      std::uint32_t child_to;
      stream >> token >> child_from >> child_to;
      auto child = parse_branch_tree(lines, pos, child_from, child_to, storage, rotors, atom_id_to_index);
      node.atom_indices.push_back(child.to_atom_index);
      node.children.push_back(child);
    } else if (line.rfind("ENDBRANCH", 0) == 0) {
      //We have found the closing keyword of the current branch we can return to the parent
      break;
    }
  }
  // The child-side axis atom is kept in the branch metadata, not as a normal child-node atom.
  node.atom_indices.erase(std::remove(node.atom_indices.begin(), node.atom_indices.end(), to), node.atom_indices.end());
  return pdbqt_torsion_branch{from, to, &node}; //return the branch with the pointer to the node that contains the children and the atoms in the branch
}

void convert_ids_to_indices(pdbqt_torsion_branch_node& node,
                            const std::unordered_map<std::uint32_t, std::size_t>& atom_id_to_index) {
  for (auto& atom_index: node.atom_indices) {
    atom_index = atom_id_to_index.at(static_cast<std::uint32_t>(atom_index));
  }
  for (auto& child: node.children) {
    child.from_atom_index = atom_id_to_index.at(static_cast<std::uint32_t>(child.from_atom_index));
    child.to_atom_index   = atom_id_to_index.at(static_cast<std::uint32_t>(child.to_atom_index));
    convert_ids_to_indices(*child.child, atom_id_to_index);
  }
}

void convert_ids_to_indices(pdbqt_torsion_tree& tree,
                            const std::unordered_map<std::uint32_t, std::size_t>& atom_id_to_index) {
  convert_ids_to_indices(tree.root, atom_id_to_index);
  for (auto& rotor: tree.rotors) {
    rotor.from_atom_index = atom_id_to_index.at(static_cast<std::uint32_t>(rotor.from_atom_index));
    rotor.to_atom_index   = atom_id_to_index.at(static_cast<std::uint32_t>(rotor.to_atom_index));
  }
}

//Same parsing as in smina
pdbqt_torsion_tree parse_pdbqt_torsion_tree(const std::vector<std::string>& lines) {
  pdbqt_torsion_tree tree;
  tree.storage.reserve(lines.size());
  tree.rotors.reserve(lines.size());
  std::unordered_map<std::uint32_t, std::size_t> atom_id_to_index;

  std::size_t pos = 0;
  while (pos < lines.size() && lines[pos].rfind("ROOT", 0) != 0) {
    ++pos;
  }
  if (pos == lines.size()) {
    throw std::runtime_error("Missing ROOT in PDBQT file");
  }
  ++pos;
  //Here we parse the root block
  while (pos < lines.size()) {
    const auto& current = lines[pos];
    pos++;
    if (current.rfind("ATOM", 0) == 0 || current.rfind("HETATM", 0) == 0) {
      const auto atom_id = atom_id_from_pdbqt_line(current);
      register_atom_id(atom_id_to_index, atom_id);
      tree.root.atom_indices.push_back(atom_id);
    } else if (current.rfind("ENDROOT", 0) == 0) {
      break;
    }
  }

  //Now we need to parse the branch section
  while (pos < lines.size()) {
    const auto& current = lines[pos];
    pos++;
    if (current.rfind("BRANCH", 0) != 0) {
      continue;
    }
    std::istringstream stream{current};
    std::string token;
    std::uint32_t from;
    std::uint32_t to;
    stream >> token >> from >> to;
    auto child = parse_branch_tree(lines, pos, from, to, tree.storage, tree.rotors, atom_id_to_index);
    //The invariant here is that we maintain the "to" atom in the parent,
    //from and to are maintained also in the struct as they are needed to apply the mobility semantics
    tree.root.atom_indices.push_back(child.to_atom_index);
    tree.root.children.push_back(child);
  }

  // During parsing these fields temporarily store PDBQT atom ids.
  // Before returning, they are converted in-place to molecule atom indices.
  convert_ids_to_indices(tree, atom_id_to_index);
  return tree;
}
//Overloaded so now we are able to pass directly the string as description instead of the file
pdbqt_torsion_tree parse_pdbqt_torsion_tree(std::string_view pdbqt_description) {
  std::vector<std::string> lines;
  std::string line;
  std::istringstream input{std::string{pdbqt_description}};
  while (std::getline(input, line)) {
    lines.push_back(line);
  }
  return parse_pdbqt_torsion_tree(lines);
}

pdbqt_torsion_tree parse_pdbqt_torsion_tree(const std::filesystem::path& pdbqt_path) {
  std::ifstream input{pdbqt_path};
  if (!input.good()) {
    throw std::runtime_error("Cannot open PDBQT file");
  }

  std::vector<std::string> lines;
  std::string line;
  while (std::getline(input, line)) {
    lines.push_back(line);
  }
  return parse_pdbqt_torsion_tree(lines);
}

} // namespace mudock
