#include <mudock/tbb_implementation/parser_filter.hpp>
#include <mudock/format/adt_mol2.hpp>
#include <mudock/format/reader.hpp>
#include <string_view>

namespace mudock {

    static constexpr std::string_view molecule_token = adt_mol2_tokens::MOLECULE_TOKEN;
    
    parser_filter::mol_vec parser_filter::operator()(std::string_view sv) const {
        mol_vec result;
    
        while (!sv.empty()) {
            // Assume the format always starts with the molecule token
            size_t next = sv.find(molecule_token, molecule_token.size());

            std::string_view mol =
                (next == std::string_view::npos) ? sv : sv.substr(0, next);

            result.push_back(std::make_unique<static_molecule>(
                mudock::parser<supported_format::ADTMOL2, static_molecule>(mol)));

            if (next == std::string_view::npos) break;
            sv.remove_prefix(next);
        }

        return result;
    }
    
} // namespace mudock    
