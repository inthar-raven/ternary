A microtonal scale research app chiefly dedicated to ternary scales (scales with three distinct step sizes). It focuses on [aberrismic theory](https://en.xen.wiki/w/Aberrismic_theory), developed by groundfault, inthar, and others.

![Front page screenshot](https://raw.githubusercontent.com/inthar-raven/ternary/main/static/images/front.png)

# How to run

Serve the project directory with any static file server. For example:

```bash
python3 -m http.server 8080
```

Then visit `http://localhost:8080/` in your browser.

# Features

- Get the set of all scales (up to mode) with a certain step signature.
- ~~Given a step signature, it gives you tuples of JI steps with bounded complexity for the scale (assuming octave equivalence).~~ (JI tunings are disabled for now)
- Given step signature, it displays the edo tunings.
- When you select a tuning on the results page, the SonicWeave code is displayed.
- JI-agnostic 2D lattice view for every scale.
- Non-octave equaves are supported (enter as a JI ratio like "3/1").
- Configurable tuning bounds:
  - Max ED size (default 111)
  - Min/max smallest step size in cents (default 20–200)
- Every scale comes with a Scale Profile that shows properties of the scale selected or queried
  - guide frame (guided generator sequence; multiplicity or interleaving polyoffset; complexity)
  - monotone MOS properties satisfied (L=M, M=s, s=0)
  - maximum variety
- Filter for
  - whether the scale is a MOS substitution scale
  - length of the guided generator sequence
  - guide frame complexity
  - monotone MOS properties
  - maximum variety
