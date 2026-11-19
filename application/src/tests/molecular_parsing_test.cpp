#include <cassert>
#include <cmath>
#include <iomanip>
#include <mudock/format/pdbqt.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/splitter.hpp>
#include <future>
#include <sstream>
#include <stdexcept>

namespace {
  std::string atom_record(const int serial,
                          const std::string_view name,
                          const double x,
                          const double y,
                          const double z,
                          const double charge,
                          const std::string_view type) {
    std::ostringstream line;
    line << std::left << std::setw(6) << "ATOM" << std::right << std::setw(5) << serial << ' '
         << std::left << std::setw(4) << name << ' ' << std::setw(3) << "LIG" << " A" << std::right
         << std::setw(4) << 1 << "    " << std::fixed << std::setprecision(3) << std::setw(8) << x
         << std::setw(8) << y << std::setw(8) << z << std::setw(6) << 1.0 << std::setw(6) << 0.0
         << "    " << std::setw(6) << std::setprecision(3) << charge << ' ' << std::setw(2) << type;
    return line.str();
  }
} // namespace

int main() {
  const std::string pdbqt = "REMARK  Name = ethane\n" + atom_record(1, "C1", 0.0, 0.0, 0.0, -0.1, "C") +
                            "\n" + atom_record(2, "C2", 1.5, 0.0, 0.0, 0.1, "C") +
                            "\nROOT\nENDROOT\nBRANCH 1 2\nENDBRANCH 1 2\nTORSDOF 1\n";
  auto molecule = mudock::parser<mudock::supported_format::PDBQT, mudock::static_molecule>(
      pdbqt);
  assert(molecule.num_atoms() == 2);
  assert(molecule.num_bonds() == 1);
  assert(molecule.num_rotamers() == 1);
  assert(molecule.torsdof() == 1);
  assert(molecule.elements(0) == mudock::element::C);
  assert(molecule.autodock_type(0) == mudock::autodock_ff::C);
  assert(molecule.atom_name(1) == "C2");
  assert(molecule.properties.get(mudock::property_type::NAME) == "ethane");
  assert(std::abs(molecule.x(1) - 1.5F) < 1e-6F);
  assert(std::abs(molecule.charge(0) + 0.1F) < 1e-6F);
  assert(molecule.bonds(0).can_rotate);
  mudock::autodock_ligand prepared_ligand{molecule};
  assert(molecule.autodock_type(0) == mudock::autodock_ff::C);
  assert(std::abs(molecule.charge(0) + 0.1F) < 1e-6F);
  assert(prepared_ligand.Rii(0) == mudock::get_description(mudock::autodock_ff::C).Rii);

  bool rejected_unprepared = false;
  try {
    mudock::static_molecule unprepared;
    (void)mudock::autodock_ligand{unprepared};
  } catch (const std::invalid_argument&) {
    rejected_unprepared = true;
  }
  assert(rejected_unprepared);

  bool rejected_empty = false;
  try {
    (void)mudock::parser<mudock::supported_format::PDBQT, mudock::static_molecule>(std::string_view{});
  } catch (const std::runtime_error&) {
    rejected_empty = true;
  }
  assert(rejected_empty);

  const std::string mol2_record = "@<TRIPOS>MOLECULE\nethane\n2 1 0 0 0\nSMALL\nUSER_CHARGES\n"
                                  "@<TRIPOS>ATOM\n1 C1 0.0 0.0 0.0 C.3 1 LIG -0.25\n"
                                  "2 C2 1.5 0.0 0.0 C.3 1 LIG 0.25\n"
                                  "@<TRIPOS>BOND\n1 1 2 1\n";
  const auto mol2_molecule =
      mudock::parser<mudock::supported_format::MOL2, mudock::static_molecule>(mol2_record);
  assert(mol2_molecule.num_atoms() == 2 && mol2_molecule.num_bonds() == 1);
  assert(mol2_molecule.atom_name(0) == "C1");
  assert(mol2_molecule.sybyl_type(0) == mudock::sybyl_atom_type::C_3);
  assert(std::abs(mol2_molecule.charge(1) - 0.25F) < 1e-6F);
  assert(std::abs(mol2_molecule.x(1) - 1.5F) < 1e-6F);

  auto pdb_line = atom_record(1, "C1", 2.0, 3.0, 4.0, 0.0, "C");
  pdb_line.replace(76, 2, " C");
  const auto pdb_molecule =
      mudock::parser<mudock::supported_format::PDB, mudock::static_molecule>(pdb_line + "\nEND\n");
  assert(pdb_molecule.num_atoms() == 1 && pdb_molecule.num_bonds() == 0);
  assert(pdb_molecule.elements(0) == mudock::element::C);
  assert(std::abs(pdb_molecule.z(0) - 4.0F) < 1e-6F);

  bool rejected_malformed = false;
  try {
    (void)mudock::parser<mudock::supported_format::PDBQT, mudock::static_molecule>("ATOM      1  BAD");
  } catch (const std::runtime_error&) {
    rejected_malformed = true;
  }
  assert(rejected_malformed);

  const std::string sdf_record = "sample\nmuDock\n\n  1  0  0  0  0  0            999 V2000\n"
                                "    1.0000    2.0000    3.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                                "M  END\n$$$$\n";
  mudock::splitter<mudock::sdf> sdf_splitter;
  auto records = sdf_splitter(sdf_record + sdf_record);
  auto trailing = sdf_splitter.flush();
  if (!trailing.empty()) records.emplace_back(std::move(trailing));
  assert(records.size() == 2);
  const auto sequential = mudock::parser<mudock::supported_format::SDF, mudock::static_molecule>(records[0]);
  auto parallel = std::async(std::launch::async, [&] {
    return mudock::parser<mudock::supported_format::SDF, mudock::static_molecule>(records[1]);
  });
  const auto parallel_result = parallel.get();
  assert(sequential.num_atoms() == 1 && sequential.num_bonds() == 0);
  assert(sequential.elements(0) == mudock::element::C);
  assert(std::abs(sequential.x(0) - 1.0F) < 1e-6F);
  assert(parallel_result.num_atoms() == sequential.num_atoms());
  assert(parallel_result.elements(0) == sequential.elements(0));
  assert(std::abs(parallel_result.z(0) - sequential.z(0)) < 1e-6F);

  const mudock::pdbqt splitter;
  assert(splitter.next_molecule_start_index(pdbqt + "REMARK  Name = next\n") != std::string_view::npos);
  assert(splitter.next_molecule_start_index("MODEL 1\n" + pdbqt + "ENDMDL\nMODEL 2\n" + pdbqt) !=
         std::string_view::npos);
  const mudock::pdb pdb_splitter;
  assert(pdb_splitter.next_molecule_start_index("MODEL 1\n" + pdb_line + "\nENDMDL\nMODEL 2\n" + pdb_line) !=
         std::string_view::npos);
}
