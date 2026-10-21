#pragma once

#include <string_view>

namespace primerpicker::glossary {

// A built-in definition (data/glossary/glossary.md): its key (the name
// `primerpicker explain` takes), other accepted names separated by
// spaces, its title, and its text (Markdown, UTF-8)
struct Entry {
    std::string_view key;
    std::string_view aliases;
    std::string_view title;
    std::string_view text;
};

}  // namespace primerpicker::glossary
