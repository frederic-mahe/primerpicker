// libFuzzer harness: the command line (src/cli/cli.cpp). The input is a
// list of arguments separated by newlines (the program name is added
// first). An argument containing '/' is replaced by a recorded molecule
// file, one containing '\\' by a recorded settings file (pd2.cfg), and
// the harness runs in an empty directory: commands read those files or
// find no file, never other files of the machine. Output is discarded.

#include "cli/cli.hpp"
#include "fuzz_input.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <span>
#include <streambuf>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

constexpr auto max_arguments = std::size_t{8};
constexpr auto program_name = std::string_view{"primerpicker"};
constexpr auto molecule_file = std::string_view{PRIMERPICKER_LEGACY_DIR "/pairs/TEST600.GB"};
constexpr auto settings_file = std::string_view{PRIMERPICKER_LEGACY_DIR "/config/PD2_RELAXED.CFG"};

// A stream buffer that discards what is written
class NullBuffer : public std::streambuf {
protected:
    auto overflow(int_type const character) -> int_type override { return traits_type::not_eof(character); }
};

auto arguments(std::string_view text) -> std::vector<std::string> {
    auto words = std::vector<std::string>{std::string{program_name}};
    while (words.size() <= max_arguments) {
        auto const end = text.find('\n');
        auto word = std::string{text.substr(0, end)};
        if (word.find('/') != std::string::npos) {
            word = molecule_file;
        } else if (word.find('\\') != std::string::npos) {
            word = settings_file;
        }
        words.push_back(std::move(word));
        if (end == std::string_view::npos) {
            break;
        }
        text.remove_prefix(end + 1);
    }
    return words;
}

}  // namespace

extern "C" auto LLVMFuzzerInitialize(int* /*argc*/, char*** /*argv*/) -> int {
    static auto discarded = NullBuffer{};
    std::cout.rdbuf(&discarded);
    std::cerr.rdbuf(&discarded);
    auto const empty = std::filesystem::temp_directory_path() / "primerpicker_fuzz_cli";
    std::filesystem::create_directories(empty);
    std::filesystem::current_path(empty);
    return 0;
}

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const* const data, std::size_t const size) -> int {
    auto const words = arguments(primerpicker::fuzz::text(primerpicker::fuzz::bytes(data, size)));
    auto pointers = std::vector<char const*>{};
    pointers.reserve(words.size());
    std::ranges::transform(words, std::back_inserter(pointers), [](std::string const& word) { return word.c_str(); });
    auto const status = primerpicker::cli::run(std::span{pointers});
    primerpicker::fuzz::check(status == EXIT_SUCCESS or status == EXIT_FAILURE);
    return 0;
}
