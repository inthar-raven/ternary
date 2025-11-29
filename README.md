A microtonal scale research app chiefly dedicated to ternary scales (scales with three distinct step sizes). It focuses on [aberrismic theory](https://en.xen.wiki/w/Aberrismic_theory), developed by groundfault, inthar, and others.

![Front page screenshot](https://raw.githubusercontent.com/inthar-raven/ternary/main/static/images/front.png)

# How to build and run

1. [Install `npm`](https://nodejs.org/en/download/package-manager).
1. Run `npm install` in the project directory.
1. Run `npm run serve` in the project directory.
1. Visit `http://localhost:8080/` with your browser. Your browser must support WebAssembly; all major browsers should.

By default the code is compiled in release mode. To compile it in development mode, make sure that the `mode` value is set to `"development"` in the project directory's `webpack.config.js` file.

# Dev scripts

1. `npm run build`: Build the app
2. `npm run serve`: Deploy the app on a development server.
3. `npm run format`: Run `prettier` on all HTML, CSS, JavaScript, and TypeScript files.

# Completed features

- Get the set of all scales (up to mode) with a certain step signature.
- Given a step signature, it gives you tuples of JI steps with bounded complexity for the scale (assuming octave equivalence).
- Given step signature, it displays the edo tunings.
- When you select a tuning on the results page, the SonicWeave code is displayed.
- JI-agnostic 2D lattice view for every scale.
- Every scale comes with a Scale Profile that shows properties of the scale selected or queried
- - guide frame (guided generator sequence (changing with tuning); multiplicity or interleaving polyoffset; complexity)
- - monotone MOS properties satisfied (L=M, M=s, s=0)
- - maximum variety
- Filter for
- - whether the scale is a MOS substitution scale
- - length of the guided generator sequence
- - guide frame complexity
- - monotone MOS properties
- - maximum variety
