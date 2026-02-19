#include <mudock/tbb_implementation/stream_filter.hpp>
#include <mudock/format/adt_mol2.hpp>
#include <string_view>

namespace mudock {
    
    stream_filter::stream_filter(std::istream& in): stream_(in) {}

    std::string stream_filter::operator()(oneapi::tbb::flow_control& fc, 
                                          std::size_t max_bytes) const {
        std::string buf; 
        buf.resize(max_bytes);

        stream_.read(buf.data(), buf.size());
        if (stream_.gcount() <= 0) {
            fc.stop();
            return {};
        }

        buf.resize(static_cast<size_t>(stream_.gcount()));

        if (!stream_.eof()) {
            // Supported format: adtmol2 (thread-safe)
            // TODO make it more generic to support other formats as well
            size_t cut = buf.rfind(adt_mol2_tokens::MOLECULE_TOKEN);
            if (cut != std::string::npos && cut != 0) {
                std::streamoff unread = static_cast<std::streamoff>(buf.size() - cut);
                stream_.clear();
                stream_.seekg(-unread, std::ios::cur);
                buf.resize(cut);
            }
        }

        return buf;
    }
    
} // namespace mudock    
