# Bin-count extraction for single-particle ICP-TOF-MS data
Bin-count extraction and split-event correction for single-particle ICP-TOF-MS ice-core data (v1.0.0)
Merges split acquisition events into single particle events and reports
the number of acquisitions each particle spans — its *bin count*.

## Background

A particle entering the plasma produces an ion cloud that is recorded
over one or more consecutive acquisitions. The number of acquisitions a
particle spans — its *bin count* — is not constant, and earlier work has
related this quantity to particle size (Details in Lee et al. (2026)).
The distribution of bin counts within a sample may therefore carry size
information that is complementary to the mass-based particle size derived
from the calibrated signal.

Recovering that information requires split acquisitions to be identified
first, so that a particle is counted once rather than several times and
its bin count is recorded. That is what this module does: it merges
consecutive acquisitions belonging to the same particle, sums their
mass, and reports the number of merged acquisitions. The bin-count
profile it returns is the input to the bin-count analysis.

Lee et al. (2026) apply the bin-count analysis to ice-core
single-particle ICP-TOF-MS measurements and obtain results consistent
with the <independent evidence> presented there. That application is
exploratory: it shows that the approach is usable and informative on
this type of data, and it is the reason this code is released, but it is
not a full methodological characterisation. See *Scope and limitations*
below.

This README documents how to run the code. The measurement setup, the
choice of parameters and the interpretation of the results are described
in the paper.


## Scope and limitations

- The correction itself — merging consecutive acquisitions into single
  events — is a generic operation and should transfer to any comparable
  single-particle time or depth series.
- Bin counts are reported as a measured quantity and are not converted
  to an absolute size. The relation between bin count and particle size
  has not been calibrated against a size standard in this work.
- That relation is statistical, not deterministic. Whether a short
  particle signal falls within one acquisition or straddles two depends
  on when the particle arrives relative to the acquisition boundary, so
  an individual two-bin event is not necessarily a larger particle than
  a one-bin event. Larger particles are only *more likely* to be split
  across several acquisitions. Bin counts are therefore interpretable in
  distribution (a.k.a. bin-count profiles), not particle by particle.
- Bin-count profiles depend on the instrument configuration and the
  measurement conditions — acquisition rate, plasma conditions, ion
  optics settings — and on the sample matrix. The parameter values used
  in the paper apply to that setup and should be re-derived, not copied,
  for other measurements.
- Establishing the bin-count analysis as a quantitative sizing method
  would require validation against particles of known size and a
  characterization of its dependence on the detection threshold. That is
  outside the scope of this release.


## Citation

If you use this code, please cite both the software and the paper:

> Lee, G., Erhardt, T., Larkman, P., Zeppenfeld, C., Jackson, S.,
> Schmitt, J., Delmonte, B., Baccolo, G., Ritz, C., Wilhelms, F.,
> Dahl-Jensen, D., Nikolaus, K. M., Bohleber, P., and Fischer, H.:
> Post-depositional geochemical transformations of aerosol impurities in
> the EPICA Dome C ice core: dissolution, mineral neoformation, and
> immobilization revealed by CFA-sp-ICP-TOFMS, EGUsphere [preprint],
> https://doi.org/10.5194/egusphere-2026-4574, 2026.
