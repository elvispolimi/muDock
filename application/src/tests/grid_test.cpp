#include <boost/program_options.hpp>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/format/autogrid.hpp>
#include <stdexcept>
#include <string>

int main(int argc, char* argv[]) {
  namespace po = boost::program_options;
  std::filesystem::path dpf_path{"rec.dpf"};
  po::options_description options("Available options");
  options.add_options()("help", "print this help message");
  options.add_options()("dpf", po::value(&dpf_path)->default_value(dpf_path), "Path to the AutoDock DPF");
  po::variables_map values;
  po::store(po::command_line_parser(argc, argv).options(options).run(), values);
  po::notify(values);

  const auto grid = mudock::load_autogrid_maps(dpf_path);
  if (grid.get_size_x() != 97 || grid.index.size_y() != 91 || grid.index.size_z() != 93)
    throw std::runtime_error("Unexpected 1fkb AutoGrid dimensions");
  const auto center = grid.get_center();
  if (std::abs(center.x() - mudock::fp_type{3}) > mudock::fp_type{1e-6f} ||
      std::abs(center.y() - mudock::fp_type{7}) > mudock::fp_type{1e-6f} ||
      std::abs(center.z() - mudock::fp_type{12.5}) > mudock::fp_type{1e-6f})
    throw std::runtime_error("Unexpected 1fkb AutoGrid center");
  if (!grid.has_map(mudock::autodock_grid_type::A) || !grid.has_map(mudock::autodock_grid_type::ELEC) ||
      !grid.has_map(mudock::autodock_grid_type::DESOLV))
    throw std::runtime_error("AutoGrid loader did not register expected channels");
  const auto affinity = grid.get_atom_map(mudock::autodock_grid_type::A);
  if (std::abs(affinity.get(0, 0, 0)) > mudock::fp_type{1e-6f})
    throw std::runtime_error("Unexpected first value in 1fkb A affinity map");

  const auto test_dir = std::filesystem::temp_directory_path() /
                        ("mudock_invalid_grid_test_" +
                         std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  std::filesystem::create_directory(test_dir);
  const auto invalid_dpf = test_dir / "invalid.dpf";
  {
    std::ofstream output(invalid_dpf);
    output << "autodock_parameter_version 4.2\n";
  }
  bool rejected_missing_fld = false;
  try {
    (void) mudock::load_autogrid_maps(invalid_dpf);
  } catch (const std::runtime_error&) {
    rejected_missing_fld = true;
  }
  if (!rejected_missing_fld) throw std::runtime_error("AutoGrid loader accepted a DPF without fld");

  {
    std::ofstream output(invalid_dpf);
    output << "fld missing_maps.fld\n";
  }
  bool rejected_missing_maps = false;
  try {
    (void) mudock::load_autogrid_maps(invalid_dpf);
  } catch (const std::runtime_error&) {
    rejected_missing_maps = true;
  }
  if (!rejected_missing_maps) throw std::runtime_error("AutoGrid loader accepted a missing FLD file");

  {
    std::ofstream dpf(invalid_dpf);
    dpf << "fld truncated.fld\n";
    std::ofstream fld(test_dir / "truncated.fld");
    fld << "dim1=2\ndim2=2\ndim3=2\n"
           "coord 1 file=coords.xyz filetype=ascii offset=0\n"
           "label=A-affinity\nlabel=Electrostatics\nlabel=Desolvation\n"
           "variable 1 file=A.map filetype=ascii skip=6\n"
           "variable 2 file=e.map filetype=ascii skip=6\n"
           "variable 3 file=d.map filetype=ascii skip=6\n";
    std::ofstream xyz(test_dir / "coords.xyz");
    xyz << "0 1\n0 1\n0 1\n";
    std::ofstream map(test_dir / "A.map");
    for (int header = 0; header < 6; ++header) map << "header\n";
    map << "0 0 0 0 0 0 0\n";
  }
  bool rejected_truncated_map = false;
  try {
    (void) mudock::load_autogrid_maps(invalid_dpf);
  } catch (const std::runtime_error&) {
    rejected_truncated_map = true;
  }
  std::filesystem::remove_all(test_dir);
  if (!rejected_truncated_map) throw std::runtime_error("AutoGrid loader accepted a truncated map");
  return EXIT_SUCCESS;
}
