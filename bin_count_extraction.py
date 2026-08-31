"""Bin-count extraction for single-particle ICP-TOF-MS data.

A single particle (SP) entering the plasma produces an ion cloud that is recorded over
one or more consecutive acquisitions.  The number of acquisitions a
particle spans -- its *bin count* -- is not constant, and earlier work has
related this quantity to particle size.  Recovering that information
requires split acquisitions to be identified first, so that a particle is
counted once rather than several times and its bin count is recorded.
That is what this module does: it merges consecutive acquisitions
belonging to the same particle, sums their mass, and reports the number of
merged acquisitions. It also returns the bin-count profiles.

Multiple acquisitions are treated as belonging to the same particle when

1. they are adjacent rows of the input matrix (determined by "consecutive_threshold"), and
2. all contains a SP detection (i.e. neither is NaN).

A run of such acquisitions is capped at "max_bins" acquisitions: a run
longer than the cap is divided into consecutive events of at most
"max_bins" acquisitions each, starting from its leading edge. The cap
is a count of acquisitions

Units
-----
The index (depth or time) and "consecutive_threshold" share the same
unit; metres are assumed here.  Matrix values are element masses per
acquisition (e.g. fg) before correction and element masses per particle
after correction.  Bin counts are acquisitions per particle.

Scope and limitations
---------------------
The merging of consecutive acquisitions implemented here is a generic
operation and should transfer to any comparable single-particle time or
depth series. The interpretation of the resulting bin counts in terms of
particle size is exploratory and is subject to the limitations set out in
the README, in particular:

* bin counts are reported as a measured quantity and are not converted to
  an absolute size;
* the relation between bin count and particle size is statistical, not
  deterministic, so bin counts are interpretable in distribution rather
  than particle by particle;
* both parameters depend on the instrument configuration, the measurement
  conditions and the sample matrix, and should be re-derived rather than
  copied for other measurements.  "consecutive_threshold" therefore has
  no default value; see :func:`consecutive_threshold_from_melt_rate`.

The measurement setup, the choice of parameters and the interpretation of
the results are described in:

    Lee, G., Erhardt, T., Larkman, P., Zeppenfeld, C., Jackson, S.,
    Schmitt, J., Delmonte, B., Baccolo, G., Ritz, C., Wilhelms, F.,
    Dahl-Jensen, D., Nikolaus, K. M., Bohleber, P., and Fischer, H.:
    Post-depositional geochemical transformations of aerosol impurities in
    the EPICA Dome C ice core: dissolution, mineral neoformation, and
    immobilization revealed by CFA-sp-ICP-TOFMS, EGUsphere [preprint],
    https://doi.org/10.5194/egusphere-2026-4574, 2026.

Example
-------
::

    import pandas as pd
    from bin_count_extraction import (
        split_event_correct, bin_count_distribution)

    data = pd.read_hdf("run01_SPProcessed.h5", "SPElemData")

    corrected, mask, bin_counts = split_event_correct(
        data,
        consecutive_threshold=1e-6,   # m, from acquisition rate x melt rate
        max_bins=5,                   # acquisitions per event
        assignment_method="max",
    )

    # Check that the upper tail is not truncated at max_bins
    print(bin_count_distribution(bin_counts))

A runnable example on a small array is given in the docstring of
:func:`split_event_correct`.

Author: Geunwoo Lee
Licence: MIT
"""

import warnings
from typing import Optional, Tuple

import numpy as np
import pandas as pd

__all__ = [
    "split_event_correct",
    "consecutive_threshold_from_melt_rate",
    "bin_count_distribution",
]
__version__ = "1.0.0"

# Fraction of events at max_bins above which the upper tail of the
# bin-count distribution is considered noticeably truncated.
TRUNCATION_WARNING_FRACTION = 0.01


def consecutive_threshold_from_melt_rate(
    acquisition_time: float,
    melt_rate: float,
    melt_rate_unit: str = "cm/min",
    factor: float = 1.0,
) -> float:
    """Depth travelled during one acquisition, for use as the threshold.

    Provided so that `consecutive_threshold` can be re-derived for a given
    measurement rather than copied from another one.

    Parameters
    ----------
    acquisition_time : float
        Duration of one acquisition, in seconds.
    melt_rate : float
        Melt speed of the continuous-flow-analysis melt head.
    melt_rate_unit : {'cm/min', 'mm/min', 'm/min', 'm/s'}
        Unit of `melt_rate`.
    factor : float
        Safety factor on the result. A value slightly above 1 tolerates
        jitter in the depth axis; exactly 1 gives the nominal spacing.

    Returns
    -------
    float
        Depth increment in metres.

    Examples
    --------
    >>> round(consecutive_threshold_from_melt_rate(1.4847e-3, 2.8) * 1e7, 3)
    6.929
    """
    per_second = {
        "cm/min": 0.01 / 60.0,
        "mm/min": 0.001 / 60.0,
        "m/min": 1.0 / 60.0,
        "m/s": 1.0,
    }
    if melt_rate_unit not in per_second:
        raise ValueError(
            f"Unknown melt_rate_unit: '{melt_rate_unit}'. "
            f"Choose one of {sorted(per_second)}."
        )
    if acquisition_time <= 0 or melt_rate <= 0 or factor <= 0:
        raise ValueError("acquisition_time, melt_rate and factor must be positive.")
    return acquisition_time * melt_rate * per_second[melt_rate_unit] * factor


def bin_count_distribution(bin_counts: pd.DataFrame) -> pd.DataFrame:
    """Number of events per bin count, per mass channel.

    A convenience summary for the check recommended in the README: if a
    notable fraction of events sits at `max_bins`, the upper tail of the
    distribution is truncated and should not be interpreted as it stands.

    Parameters
    ----------
    bin_counts : pd.DataFrame
        The third output of :func:`split_event_correct`.

    Returns
    -------
    pd.DataFrame
        Event counts with the bin count as index and one column per mass
        channel.
    """
    counts = {
        column: bin_counts[column].dropna().astype(int).value_counts()
        for column in bin_counts.columns
    }
    table = pd.DataFrame(counts).fillna(0).astype(int).sort_index()
    table.index.name = "bin_count"
    return table


def split_event_correct(
    sp_processed_icp_data: pd.DataFrame,
    consecutive_threshold: float,
    max_bins: Optional[int] = 5,
    assignment_method: str = "max",
    warn_on_truncation: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Integrate split single-particle acquisition events into single bins.

    Parameters
    ----------
    sp_processed_icp_data : pd.DataFrame
        Depth (or time) index by elemental mass columns. Acquisitions
        without a detection must be NaN. The index must increase
        monotonically.
    consecutive_threshold : float
        Maximum gap along the index between two acquisitions for them to
        be treated as the same particle, in the unit of the index. There
        is deliberately no default: the appropriate value follows from the
        acquisition rate and the melt rate of the measurement, see
        :func:`consecutive_threshold_from_melt_rate`.
    max_bins : int or None, default 5
        Maximum number of acquisitions merged into one event. Runs longer
        than this are split from the leading edge, which truncates the
        upper tail of the bin-count distribution. `None` merges runs of
        any length.
    assignment_method : {'max', 'first'}, default 'max'
        Bin that carries the integrated mass of a merged event:
        'max' uses the bin with the highest single-acquisition signal
        (ties resolved towards the leading edge), 'first' uses the
        leading-edge bin. This affects only the reported position of the
        event, never its mass or its bin count.
    warn_on_truncation : bool, default True
        Emit a warning when more than 1 % of the events reach `max_bins`,
        i.e. when the bin-count distribution is noticeably truncated.

    Returns
    -------
    sp_splitcorrected : pd.DataFrame
        Integrated particle masses. All merged bins other than the target
        bin are set to NaN.
    processed_mask : pd.DataFrame
        Boolean, True for every acquisition assigned to an event.
    bin_counts : pd.DataFrame
        Number of acquisitions merged into the event, stored in the same
        cell that carries the integrated mass (1 for unsplit events, NaN
        for background).

    Raises
    ------
    ValueError
        If the index is not monotonically increasing, if it contains NaN,
        or if an argument is out of range.

    See Also
    --------
    bin_count_distribution : Summary of the returned bin counts.

    Notes
    -----
    Total mass is conserved: the sum over `sp_splitcorrected` equals the
    sum over the input. The sum over `bin_counts` equals the number of
    detections in the input.

    Examples
    --------
    >>> import numpy as np, pandas as pd
    >>> data = pd.DataFrame(
    ...     {"Al27": [np.nan, 2.0, 5.0, 1.0, np.nan, 3.0]},
    ...     index=pd.Index(np.arange(6) * 1e-7, name="Depth_top [m]"))
    >>> corrected, mask, bin_counts = split_event_correct(data, 1e-6)
    >>> corrected["Al27"].tolist()
    [nan, nan, 8.0, nan, nan, 3.0]
    >>> bin_counts["Al27"].tolist()
    [nan, nan, 3.0, nan, nan, 1.0]
    """
    if not isinstance(sp_processed_icp_data, pd.DataFrame):
        raise TypeError("sp_processed_icp_data must be a pandas DataFrame.")
    if consecutive_threshold < 0:
        raise ValueError("consecutive_threshold must be non-negative.")
    if max_bins is not None and max_bins < 1:
        raise ValueError("max_bins must be at least 1, or None for no limit.")
    if assignment_method not in ("max", "first"):
        raise ValueError(
            f"Invalid assignment_method: '{assignment_method}'. "
            "Choose 'max' or 'first'."
        )

    index_values = np.asarray(sp_processed_icp_data.index.values, dtype=float)
    if np.isnan(index_values).any():
        raise ValueError("The index contains NaN.")
    if np.any(np.diff(index_values) < 0):
        raise ValueError(
            "The index must increase monotonically; sort the input first. "
            "A decreasing index would silently merge acquisitions across the "
            "reversal, because a negative gap never exceeds the threshold."
        )

    columns = sp_processed_icp_data.columns
    data_matrix = np.asarray(sp_processed_icp_data.values, dtype=float).copy()
    n_rows, n_cols = data_matrix.shape

    processed_mask = np.zeros((n_rows, n_cols), dtype=bool)
    bin_counts = np.full((n_rows, n_cols), np.nan)

    # No limit is expressed as a cap that can never be reached.
    bin_limit = n_rows if max_bins is None else max_bins

    # Process each mass channel independently
    for col_idx in range(n_cols):
        # Copy, because data_matrix is modified in place below and the
        # look-ahead must see the original detections.
        data = data_matrix[:, col_idx].copy()
        i = 0

        while i < n_rows:
            # Skip background
            if np.isnan(data[i]):
                i += 1
                continue

            # Extend the event while the next acquisition is adjacent,
            # carries a detection and lies within the threshold, up to
            # bin_limit acquisitions.
            j = i + 1
            while (j < n_rows) and ((j - i) < bin_limit):
                if np.isnan(data[j]) or (
                    (index_values[j] - index_values[j - 1]) > consecutive_threshold
                ):
                    break
                j += 1

            seq_indices = np.arange(i, j)
            seq_len = len(seq_indices)

            if seq_len == 1:
                # Unsplit event: leave the mass where it is
                bin_counts[i, col_idx] = 1
                processed_mask[i, col_idx] = True
            else:
                # Split event spanning several consecutive acquisitions
                seq_data = data[seq_indices]
                seq_sum = np.sum(seq_data)

                if assignment_method == "max":
                    target_idx = seq_indices[np.argmax(seq_data)]
                else:  # 'first'
                    target_idx = i

                # Clear the merged bins, then place the integrated mass
                data_matrix[seq_indices, col_idx] = np.nan
                data_matrix[target_idx, col_idx] = seq_sum
                bin_counts[target_idx, col_idx] = seq_len
                processed_mask[seq_indices, col_idx] = True

            # Continue after the end of the current event
            i = j

    if warn_on_truncation and max_bins is not None:
        n_events = int(np.isfinite(bin_counts).sum())
        if n_events:
            at_limit = int((bin_counts == max_bins).sum())
            fraction = at_limit / n_events
            if fraction > TRUNCATION_WARNING_FRACTION:
                warnings.warn(
                    f"{at_limit} of {n_events} events ({100 * fraction:.1f} %) "
                    f"reach max_bins={max_bins}, so the upper tail of the "
                    "bin-count distribution is truncated. Increase max_bins "
                    "until the distribution falls off on its own before "
                    "interpreting bin counts as size information.",
                    RuntimeWarning,
                    stacklevel=2,
                )

    index = sp_processed_icp_data.index
    return (
        pd.DataFrame(data_matrix, index=index, columns=columns),
        pd.DataFrame(processed_mask, index=index, columns=columns),
        pd.DataFrame(bin_counts, index=index, columns=columns),
    )
