
#include <algorithm>
#include <mudock/chem/drugscore/drugscore_types.hpp>
#include <stdexcept>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/drugscore_types.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  ds_type parse_ds_type_symbol(const std::string_view symbol) {
    const auto element_it = std::find_if(std::begin(DS_TYPE_DICTIONARY),
                                         std::end(DS_TYPE_DICTIONARY),
                                         [&symbol](const auto& e) { return e.symbol == symbol; });
    if (element_it != std::end(DS_TYPE_DICTIONARY))
      return element_it->value;
    else
      throw std::runtime_error("Parsing drugscore type symbol error: Missing element: " +
                               std::string(symbol));
  }

  const std::array<ds_type_description, num_ds_types()> DS_TYPE_DICTIONARY = {{
      {ds_type::Br, "Br"},       {ds_type::C_2, "C.2"},     {ds_type::C_3, "C.3"},   {ds_type::C_ar, "C.ar"},
      {ds_type::C_cat, "C.cat"}, {ds_type::N_2, "N.2"},     {ds_type::N_4, "N.4"},   {ds_type::N_am, "N.am"},
      {ds_type::N_ar, "N.ar"},   {ds_type::N_pl3, "N.pl3"}, {ds_type::O_2, "O.2"},   {ds_type::O_3, "O.3"},
      {ds_type::O_co2, "O.co2"}, {ds_type::S_3, "S.3"},     {ds_type::C_1, "C.1"},   {ds_type::N_3, "N.3"},
      {ds_type::Cl, "Cl"},       {ds_type::F, "F"},         {ds_type::I, "I"},       {ds_type::N_1, "N.1"},
      {ds_type::P_3, "P.3"},     {ds_type::S_2, "S.2"},     {ds_type::S_o2, "S.o2"}, {ds_type::S_o, "S.o"},
      {ds_type::METAL, "METAL"}, {ds_type::ND, "ND"},
  }};

} // namespace mudock
