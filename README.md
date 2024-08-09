This is a microtonal scale research app chiefly dedicated to ternary scales (scales with three distinct step sizes). It focuses on aberrismic theory, developed by groundfault, inthar, and others.

# How to build and run:
On Linux, install the Rust toolchain and run `cargo run --release` in this directory. You should see a link. Ctrl-click on it to go to the main page.

Windows is not supported at this time.

# Completed features:
* Get the set of all scales (up to mode) with a certain step signature.
* Given a step signature, it gives you tuples of JI steps with bounded complexity for the scale (assuming octave equivalence).
* Given step signature, it displays the edo tunings.
* When you select a tuning on the results page, the SonicWeave code is displayed.
* JI-agnostic 2D lattice view for every scale.
* Every scale comes with a Scale Profile that shows properties of the scale selected or queried
* * guide frame (guided generator sequence (changing with tuning); multiplicity or interleaving polyoffset; complexity)
* * monotone MOS properties satisfied (L=M, M=s, s=0)
* * maximum variety
* Filter for
* * whether the scale is a MOS substitution scale
* * length of the guided generator sequence
* * guide frame complexity
* * monotone MOS properties
* * maximum variety

# Planned features:
* Support non-octave equaves
* Remove the need for a server timeout by migrating computational logic to client side, possibly using Rust+WASM
