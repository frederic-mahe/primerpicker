#pragma once

#include <span>

namespace primerpicker::cli {

// Runs the command line: arguments as main() receives them (program
// name first); results go to std::cout, errors to std::cerr; returns
// the exit status. Exceptions (std::bad_alloc, std::format_error) are
// reported as errors, not thrown.
[[nodiscard]] auto run(std::span<char const* const> arguments) -> int;

}  // namespace primerpicker::cli
