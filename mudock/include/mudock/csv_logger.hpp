#pragma once

#include <filesystem>
#include <vector>
#include <fstream>
#include <utility>

namespace mudock {

    class csv_logger {
    public:
        csv_logger(const std::string& path, const std::vector<std::string>& headers) {
            const bool exists =
                std::filesystem::exists(path) &&
                std::filesystem::file_size(path) > 0;

            file_.open(path, std::ios::app);

            if (!exists) {
                for (std::size_t i = 0; i < headers.size(); ++i) {
                    if (i) file_ << ",";
                    file_ << headers[i];
                }
                file_ << "\n";
            }
        }

        template<typename... Args>
        void log(Args&&... args) {
            bool first = true;

            ((file_ << (std::exchange(first, false) ? "" : ",")
                    << args), ...);

            file_ << "\n";
        }

    private:
        std::ofstream file_;
    };
}