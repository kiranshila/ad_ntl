# Code Repository for "Design of Wideband Matching Structures using Automatic Differentiation"

## Organization

Simulated data is found in the `data` folder.
Implementations of the various formulas are found in `src`.
The main scripts are found in `scripts`.

The design for the paper was computed with `scripts/wblna.jl`.
Two small examples for resistive matching with an ideal transmission line and microstrip with simulated RLGC data are in `scripts/simple_ntl_matching.jl` and `scripts/microstrip_matching.jl`, respectively.
Post-processing of measured data was performed with `scripts/post_processing.jl`.
