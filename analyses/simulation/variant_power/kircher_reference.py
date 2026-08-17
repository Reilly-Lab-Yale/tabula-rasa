"""Published-effect reference lines for the pairwise power heatmaps.

Saturation mutagenesis of disease-associated regulatory elements reports a
median significant effect of about a fifth of the wild-type level, which is
+20% for activating variants and -24% for repressing ones (Kircher et al.
2019, doi:10.1038/s41467-019-11526-w). In the |log2 FC| the heatmaps are
gridded on, those are 0.263 and 0.398.

Both are drawn, since the median depends on the direction of effect and the
grid is unsigned. They are medians of two different distributions, so the
interval between them carries no meaning and is left unshaded.

The two lines are encoded redundantly in colour and in shape, so they stay
distinguishable in greyscale and under colour-vision deficiency:

    activating (0.263)  black, long dashes, circle markers
    repressing (0.398)  blue, dotted, square markers

The lines are drawn plain, with no white casing or underlay. Casing a dashed
line fills its gaps as well as its edges, so the rule reads as a row of
white-bordered blocks rather than as one line.
"""

import numpy as np

ACTIVATING_LOG2FC = 0.263   # log2(1.20)
REPRESSING_LOG2FC = 0.398   # |log2(0.759)|

_ACTIVATING_STYLE = dict(color="#000000", linestyle=(0, (7, 3)), marker="o")
_REPRESSING_STYLE = dict(color="#2a78d6", linestyle=(0, (1, 2.2)), marker="s")

_LW = 2.0
_MARKER_SIZE = 4.0
_N_MARKERS = 5


def _row_position(fc_vals, target):
    """Interpolated y for `target` on a categorical heatmap axis.

    `fc_vals` is the pivot index in plot order (descending), each row
    occupying one unit with its centre at i + 0.5. Returns None if the
    target falls outside the plotted range, so a grid that does not reach
    it simply gets no line rather than one pinned to the edge.
    """
    fc_vals = list(fc_vals)
    assert fc_vals == sorted(fc_vals, reverse=True), (
        f"expected a descending fold-change index, got {fc_vals}")
    if target > fc_vals[0] or target < fc_vals[-1]:
        return None
    for i in range(len(fc_vals) - 1):
        hi, lo = fc_vals[i], fc_vals[i + 1]
        if lo <= target <= hi:
            frac = 0.0 if hi == lo else (hi - target) / (hi - lo)
            return (i + 0.5) + frac
    raise AssertionError(f"{target} bracketed by no pair of {fc_vals}")


def _draw_line(ax, y, style):
    """One reference line: the dashed rule, then its markers."""
    ax.axhline(y=y, color=style["color"], linestyle=style["linestyle"],
               lw=_LW, zorder=5)
    marker_x = np.linspace(0.10, 0.90, _N_MARKERS)
    ax.plot(marker_x, np.full(_N_MARKERS, y),
            linestyle="none", marker=style["marker"], ms=_MARKER_SIZE,
            color=style["color"], markeredgecolor="none",
            transform=ax.get_yaxis_transform(), zorder=6)


def annotate(ax, fc_vals):
    """Draw the two medians as distinguishable reference lines.

    Returns the number of lines actually drawn, so a caller can assert the
    grid reaches them.
    """
    drawn = 0
    for target, style in ((ACTIVATING_LOG2FC, _ACTIVATING_STYLE),
                          (REPRESSING_LOG2FC, _REPRESSING_STYLE)):
        y = _row_position(fc_vals, target)
        if y is None:
            continue
        _draw_line(ax, y, style)
        drawn += 1
    return drawn


def annotate_continuous(ax):
    """Same two medians, for axes whose y is log2 FC on a continuous scale."""
    for y, style in ((ACTIVATING_LOG2FC, _ACTIVATING_STYLE),
                     (REPRESSING_LOG2FC, _REPRESSING_STYLE)):
        _draw_line(ax, y, style)
    return 2
