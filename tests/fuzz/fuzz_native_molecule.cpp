// libFuzzer harness: Primer Designer's own molecule format
// (core/io/legacy_native_molecule.hpp). A file that is read is written
// back; what is written must be read again, and written again
// identically (the writer's output is a fixed point); it must also load
// through parse_molecule_file as the native format.

#include "core/io/legacy_native_molecule.hpp"
#include "core/io/molecule_file.hpp"
#include "fuzz_input.hpp"

#include <cstddef>
#include <cstdint>
#include <string>

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const* const data, std::size_t const size) -> int {
    auto const text = primerpicker::fuzz::text(primerpicker::fuzz::bytes(data, size));
    auto const result = primerpicker::io::parse_native_molecule(text);
    using primerpicker::fuzz::check;
    check(result.molecule.has_value() == (result.error == primerpicker::io::NativeError::none));
    if (not result.molecule.has_value()) {
        return 0;
    }
    auto const written = primerpicker::io::write_native_molecule(*result.molecule);
    check(written.ends_with("\r\n"));
    auto const reread = primerpicker::io::parse_native_molecule(written);
    check(reread.molecule.has_value());
    if (not reread.molecule.has_value()) {
        return 0;
    }
    check(primerpicker::io::write_native_molecule(*reread.molecule) == written);
    check(reread.molecule->genes.size() == result.molecule->genes.size());
    check(reread.molecule->markers.size() == result.molecule->markers.size());
    check(reread.molecule->linear_ends.has_value() == result.molecule->linear_ends.has_value());
    check(reread.molecule->sequence == result.molecule->sequence);
    // the general reader takes it as the native format (unless its name
    // line looks like GenBank's or FASTA's first line): a molecule, or
    // one of the native format's errors
    auto const& name = result.molecule->name;
    if (name.starts_with("LOCUS") or name.starts_with('>')) {
        return 0;
    }
    auto const loaded = primerpicker::io::parse_molecule_file({.text = written, .name = "FUZZ"});
    using primerpicker::io::MoleculeError;
    check(loaded.molecule.has_value() or loaded.error == MoleculeError::empty_sequence or
          loaded.error == MoleculeError::protein or loaded.error == MoleculeError::illegal_format);
    return 0;
}
