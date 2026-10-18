// primerpicker: entry point of the command-line front-end (cli.cpp)

#include "cli/cli.hpp"

#include <cstddef>
#include <span>

auto main(int const argc, char const* const* const argv) -> int {
    return primerpicker::cli::run(std::span{argv, static_cast<std::size_t>(argc)});
}
