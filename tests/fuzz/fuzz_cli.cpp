// libFuzzer harness: the command line (src/cli/cli.cpp). The input is a
// list of arguments separated by newlines (the program name is added
// first). An argument containing '/' is replaced by a recorded molecule
// file (the number of '/' chooses which), one containing '\\' by a
// recorded settings file (pd2.cfg), and the argument after "--report"
// by a file name in the harness's directory, which starts empty:
// commands read those files or find no file, never other files of the
// machine, and write reports only there (never over a recorded file).
// Output is discarded.
// The commands get copies of the recorded files, made in the harness's
// directory at start-up, and each argument ends at its first NUL byte,
// as the command line sees it: an argument "--report<NUL>x" made the
// next argument a recorded file's path, and a report was written over
// that file (found 2026-10-20).

#include "cli/cli.hpp"
#include "fuzz_input.hpp"

#include <algorithm>
#include <array>
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
// TEST600 (one '/'), then molecules with homology sites, a short one
// and a circular one
constexpr auto molecule_files = std::array{
    std::string_view{PRIMERPICKER_LEGACY_DIR "/pairs/TEST600.GB"},
    std::string_view{PRIMERPICKER_LEGACY_DIR "/homology/TREP2.GB"},
    std::string_view{PRIMERPICKER_LEGACY_DIR "/homology/map/TMAP660F.GB"},
    std::string_view{PRIMERPICKER_LEGACY_DIR "/homology/map/TMAP100.GB"},
    std::string_view{PRIMERPICKER_LEGACY_DIR "/enzymes/TENZC.GB"}};
constexpr auto settings_file = std::string_view{PRIMERPICKER_LEGACY_DIR "/config/PD2_RELAXED.CFG"};
constexpr auto report_option = std::string_view{"--report"};
constexpr auto report_file = std::string_view{"report.prn"};  // in the harness's directory

// A stream buffer that discards what is written
class NullBuffer : public std::streambuf {
protected:
    auto overflow(int_type const character) -> int_type override { return traits_type::not_eof(character); }
};

// The name of the copy of a recorded file, in the harness's directory
auto copy_name(std::string_view const recorded) -> std::string {
    return std::filesystem::path{recorded}.filename().string();
}

auto arguments(std::string_view text) -> std::vector<std::string> {
    auto words = std::vector<std::string>{std::string{program_name}};
    while (words.size() <= max_arguments) {
        auto const end = text.find('\n');
        auto const line = text.substr(0, end);
        auto word = std::string{line.substr(0, line.find('\0'))};
        auto const slashes = static_cast<std::size_t>(std::ranges::count(word, '/'));
        if (words.back() == report_option) {
            word = report_file;
        } else if (slashes > 0) {
            word = copy_name(molecule_files.at((slashes - 1) % molecule_files.size()));
        } else if (word.find('\\') != std::string::npos) {
            word = copy_name(settings_file);
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
    for (auto const recorded : molecule_files) {
        std::filesystem::copy_file(recorded, copy_name(recorded), std::filesystem::copy_options::overwrite_existing);
    }
    std::filesystem::copy_file(settings_file, copy_name(settings_file),
                               std::filesystem::copy_options::overwrite_existing);
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
