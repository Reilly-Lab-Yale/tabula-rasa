"""Published-effect reference lines for the pairwise power heatmaps.

Saturation mutagenesis of disease-associated regulatory elements reports a
median significant effect of about a fifth of the wild-type level, which is
+20% for activating variants and -24% for repressing ones (Kircher et al.
2019, doi:10.1038/s41467-019-11526-w). In the |log2 FC| the heatmaps are
gridded on, those are 0.263 and 0.398.

Both are drawn, since the median depends on the direction of effect and the
grid is unsigned. The band between them is the honest answer to "where does a
typical published variant effect fall".
"""

ACTIVATING_LOG2FC = 0.263   # log2(1.20)
REPRESSING_LOG2FC = 0.398   # |log2(0.759)|

_COLOR = "#0033cc"


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


def annotate(ax, fc_vals, label=True):
    """Draw the two medians as dashed lines with the band between them shaded.

    Returns the number of lines actually drawn, so a caller can assert the
    grid reaches them.
    """
    y_act = _row_position(fc_vals, ACTIVATING_LOG2FC)
    y_rep = _row_position(fc_vals, REPRESSING_LOG2FC)

    if y_act is not None and y_rep is not None:
        ax.axhspan(min(y_act, y_rep), max(y_act, y_rep),
                   color=_COLOR, alpha=0.08, zorder=3)

    drawn = 0
    for y in (y_act, y_rep):
        if y is None:
            continue
        ax.axhline(y=y, color=_COLOR, linestyle="--", lw=1.3, alpha=0.75,
                   zorder=4)
        drawn += 1

    if label and drawn:
        y_top = min(y for y in (y_act, y_rep) if y is not None)
        ax.text(0.02, y_top - 0.25,
                "median published\nvariant effect",
                transform=ax.get_yaxis_transform(), fontsize=7,
                color=_COLOR, va="bottom", zorder=5)
    return drawn


def annotate_continuous(ax, label=True):
    """Same two medians, for axes whose y is log2 FC on a continuous scale."""
    ax.axhspan(ACTIVATING_LOG2FC, REPRESSING_LOG2FC,
               color=_COLOR, alpha=0.08, zorder=0)
    for y in (ACTIVATING_LOG2FC, REPRESSING_LOG2FC):
        ax.axhline(y=y, color=_COLOR, linestyle="--", lw=1.3, alpha=0.75)
    if label:
        ax.text(ax.get_xlim()[1] * 0.98, REPRESSING_LOG2FC,
                "median published variant effect",
                ha="right", va="bottom", fontsize=8, color=_COLOR)
