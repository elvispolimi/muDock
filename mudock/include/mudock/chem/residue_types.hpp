#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <string_view>

namespace mudock {

  // Compact classifier for residue/substructure names.
  // This is NOT intended to replace the original MOL2 subst_name token.
  // Keep the original subst_name in molecule.atom_residue_name/residue_name for round-tripping.
  enum class residue_type : int {
    UNKNOWN = 0,

    // Standard amino acids
    ALA,
    ARG,
    ASN,
    ASP,
    CYS,
    GLN,
    GLU,
    GLY,
    HIS,
    HSD,
    HSE,
    HSP,
    HID,
    HIE,
    HIP,
    ILE,
    LEU,
    LYS,
    MET,
    PHE,
    PRO,
    SER,
    THR,
    TRP,
    TYR,
    VAL,

    // Common protonation/terminal variants
    ASH,
    CYM,
    CYX,
    GLH,
    LYN,
    MSE,
    NME,
    ACE,

    // Nucleic-acid bases/residues
    A,
    C,
    G,
    T,
    U,
    DA,
    DC,
    DG,
    DT,
    DU,
    RA,
    RC,
    RG,
    RT,
    RU,
    ADE,
    CYT,
    GUA,
    THY,
    URA,

    // Solvent / ions / generic ligand-like names
    HOH,
    WAT,
    SOL,
    TIP3,
    SPC,
    NA,
    K,
    CL,
    CA,
    MG,
    MN,
    FE,
    CO,
    NI,
    CU,
    ZN,
    CD,
    HG,

    LIG,
    UNL,
    UNK,
    MOL,
    DRG,
    INH,

    // Common cofactors / biochemical groups
    ATP,
    ADP,
    AMP,
    GTP,
    GDP,
    GMP,
    NAD,
    NAP,
    NDP,
    FAD,
    FMN,
    HEM,
    HEC,
    HEB,
    NAG,
    NDG,
    BMA,
    MAN,
    GAL,
    GLC,
    SO4,
    PO4,
    ACT,
    EDO,
    GOL,
    DMS,
    DMSO
  };

  struct residue_type_description {
    residue_type value;
    std::string_view name;
    bool is_standard_amino_acid;
    bool is_nucleic_acid;
    bool is_water;
    bool is_ion;
    bool is_ligand_like;
  };

  constexpr std::size_t num_residue_types() { return static_cast<std::size_t>(residue_type::DMSO) + 1; }

  extern const std::array<residue_type_description, num_residue_types()> RESIDUE_TYPE_DICTIONARY;

  inline const residue_type_description& get_description(const residue_type type) {
    assert(RESIDUE_TYPE_DICTIONARY[static_cast<std::size_t>(type)].value == type);
    return RESIDUE_TYPE_DICTIONARY[static_cast<std::size_t>(type)];
  }

  // Returns residue_type::UNKNOWN for arbitrary ligand/substructure names such as C20.
  residue_type parse_residue_type(std::string_view token);

  // Throws if token is not present in RESIDUE_TYPE_DICTIONARY.
  residue_type parse_residue_type_or_throw(std::string_view token);

  inline std::string_view to_string(const residue_type type) { return get_description(type).name; }

  inline bool is_standard_amino_acid(const residue_type type) {
    return get_description(type).is_standard_amino_acid;
  }

  inline bool is_nucleic_acid(const residue_type type) { return get_description(type).is_nucleic_acid; }

  inline bool is_water(const residue_type type) { return get_description(type).is_water; }

  inline bool is_ion(const residue_type type) { return get_description(type).is_ion; }

  inline bool is_ligand_like(const residue_type type) { return get_description(type).is_ligand_like; }

} // namespace mudock
