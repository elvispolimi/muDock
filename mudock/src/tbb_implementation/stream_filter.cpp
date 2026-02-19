#include <mudock/tbb_implementation/stream_filter.hpp>
#include <mudock/format/adt_mol2.hpp>

namespace mudock {
    
    stream_filter::stream_filter(std::istream& in, std::size_t end)
    : stream_(in)
    {
        if (end != std::numeric_limits<std::size_t>::max()) {
            end_ = end;
            return;
        }

        // default case (no MPI involved)
        std::streampos cur = stream_.tellg();
        stream_.seekg(0, std::ios::end);
        end_ = static_cast<std::size_t>(stream_.tellg());
        stream_.seekg(cur);
    }

    std::string stream_filter::operator()(oneapi::tbb::flow_control& fc,
                                        std::size_t max_bytes) const {
        // Check empty/invalid file
        const std::streampos pos = stream_.tellg();
        if (pos == std::streampos(-1)) { 
            fc.stop();
            return {};
        }

        // Check if we've reached the end of the assigned range
        const std::size_t cur = static_cast<std::size_t>(pos);
        if (cur >= end_) {               
            fc.stop();
            return {};
        }

        // Ensure we don't read past the end of the assigned range
        const std::size_t remaining = end_ - cur;
        if (max_bytes > remaining) max_bytes = remaining;  

        std::string buf;
        buf.resize(max_bytes);

        stream_.read(buf.data(), buf.size());
        if (stream_.gcount() <= 0) {
            fc.stop();
            return {};
        }

        buf.resize(static_cast<size_t>(stream_.gcount()));

        if (!stream_.eof()) {
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