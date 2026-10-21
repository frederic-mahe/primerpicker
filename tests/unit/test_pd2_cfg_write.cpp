#include "core/io/pd2_cfg.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using primerpicker::io::LegacySystemDefaults;
using primerpicker::io::parse_pd2_cfg;
using primerpicker::io::write_pd2_cfg;

auto read_bytes(std::filesystem::path const& path) -> std::vector<std::byte> {
    auto file = std::ifstream{path, std::ios::binary};
    if (not file) {
        throw std::runtime_error{"cannot read " + path.string()};
    }
    auto text = std::ostringstream{};
    text << file.rdbuf();
    auto const contents = std::move(text).str();
    auto bytes = std::vector<std::byte>{};
    std::ranges::transform(contents, std::back_inserter(bytes),
                           [](char const character) { return static_cast<std::byte>(character); });
    return bytes;
}

// every settings file of tests/legacy: written by PRIMER.exe (System
// Defaults saved) or written as inputs of the probes
auto settings_files() -> std::vector<std::filesystem::path> {
    auto files = std::vector<std::filesystem::path>{};
    for (auto const& entry : std::filesystem::recursive_directory_iterator{PRIMERPICKER_LEGACY_DIR}) {
        auto extension = entry.path().extension().string();
        std::ranges::transform(extension, extension.begin(), [](char const letter) {
            return letter >= 'a' and letter <= 'z' ? static_cast<char>(letter - 'a' + 'A') : letter;
        });
        if (entry.is_regular_file() and extension == ".CFG") {
            files.push_back(entry.path());
        }
    }
    return files;
}

}  // namespace

// Reference: every pd2.cfg file of tests/legacy (tests/legacy/README.md,
// "Settings file pd2.cfg"; PRIMER.exe keeps bytes 0-3 and 47-54 as read):
// writing what was read gives the same bytes
TEST_CASE("write_pd2_cfg: every recorded settings file, byte for byte", "[io][legacy]") {
    auto const files = settings_files();
    CHECK(files.size() >= 40);
    for (auto const& file : files) {
        INFO(file.string());
        auto const bytes = read_bytes(file);
        auto const parsed = parse_pd2_cfg(bytes);
        REQUIRE(parsed.defaults.has_value());
        CHECK(write_pd2_cfg(parsed.defaults.value_or(LegacySystemDefaults{})) == bytes);
    }
}

// Reference: config/PD2_LENGTH25.CFG, saved by PRIMER.exe with the
// built-in defaults but the primer length (tests/legacy/README.md,
// "Settings file pd2.cfg")
TEST_CASE("write_pd2_cfg: the built-in defaults with a primer length of 25", "[io][legacy]") {
    auto defaults = LegacySystemDefaults{};
    constexpr auto length = std::size_t{25};
    defaults.primer_length = length;
    CHECK(write_pd2_cfg(defaults) ==
          read_bytes(std::filesystem::path{PRIMERPICKER_LEGACY_DIR} / "config" / "PD2_LENGTH25.CFG"));
}
