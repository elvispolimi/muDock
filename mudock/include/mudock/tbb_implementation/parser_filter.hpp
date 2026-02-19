#pragma once

#include <memory>
#include <string_view>
#include <vector>

#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>

namespace mudock {
    
    class parser_filter {
        public:
            using molVec = std::vector<std::unique_ptr<static_molecule>>;

            molVec operator()(std::string_view sv) const;
    };

} // namespace mudock
