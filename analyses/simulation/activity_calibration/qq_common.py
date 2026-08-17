"""Shared machinery for the null-calibration QQ figures.

Two figures draw the same curve: `shendure/plot_calibration_qq_detail.py`
(one panel per cell type, one dataset) and `qq_overlay.py` (every cell type
overlaid, all three datasets). They agree on the form, the reference band and
the downsampling, so those live here.

The form is `observed - expected` against `expected` rather than a plain
[0, 1] QQ. Departures from uniform in this data top out at KS D = 0.036, which
on an undifferenced axis is a vertical excursion of 3.6% of the panel height at
worst -- a line lying on the diagonal, with the whole result below the
resolution of the page. Subtracting the diagonal keeps exactly the same
information and spends the whole y-axis on the part that varies. Ideal
calibration is then the horizontal line at zero.

Below zero the observed p-values are smaller than uniform, i.e. the test is
anti-conservative and the false-positive rate runs above nominal; above zero it
is conservative.

The grey band is the two-sided Kolmogorov 95% acceptance region,
D_crit = 1.358 / sqrt(n), so a curve that stays inside the band is exactly one
whose KS p-value is above 0.05. It is a band for the whole curve, not a
pointwise interval: the deviation curve is a Brownian bridge and is strongly
autocorrelated, so a pointwise interval would invite reading a single excursion
as a result.
"""
import numpy as np

KS_CRIT_95 = 1.358  # Kolmogorov two-sided 95% point; D_crit = KS_CRIT_95 / sqrt(n)

INK, MUTED, HAIR = "#1a1a1a", "#666666", "#c9c9c9"
BAND = "#e4e4e4"

# Okabe-Ito, the palette scripts/figure_styling.py in the manuscript repo
# remaps everything else to. Same key as fpr_dumbbell.py, so a colour means
# the same test in both calibration figures. The pair passes every categorical
# check outright, including contrast against the page.
TEST_COLOR = {"ttest": "#d55e00", "mwu": "#0072b2"}
TEST_LABEL = {"ttest": "t-test", "mwu": "MWU"}
CONDITION_LABEL = {"reporter": "reporter used", "deflated": "reporter withheld"}

# rcParams every figure here needs regardless of its own type scale.
RC_BASE = {
    "svg.fonttype": "none",       # editable <text> downstream
    "pdf.fonttype": 42,
    "axes.unicode_minus": False,  # the repo is plain ASCII
}

N_PRETHIN = 40000     # rank grid before simplification
RDP_FRACTION = 0.30   # simplification tolerance, as a fraction of the linewidth


def thin_ranks(n, n_keep=N_PRETHIN):
    """Rank indices to draw, dense at both tails and even through the bulk.

    One vector element per p-value is what makes the notebook's version of this
    figure a 74 MB SVG. The deviation curve is a Brownian bridge: smooth over
    the bulk, so an even grid reproduces it, and steepest as it leaves zero at
    either end, which the log-spaced ends cover. This is only the first pass;
    `simplify` then drops whatever the drawn stroke would cover anyway.
    """
    assert n_keep >= 8, f"n_keep={n_keep} is too few to shape a curve"
    if n <= n_keep:
        return np.arange(n)
    tail = n_keep // 4
    lo = np.unique(np.round(np.geomspace(1, n / 2, tail)).astype(int)) - 1
    hi = n - 1 - lo
    mid = np.linspace(0, n - 1, n_keep - 2 * tail).astype(int)
    return np.unique(np.concatenate([lo, mid, hi]))


def simplify(x, y, tol_y):
    """Douglas-Peucker keep-mask, distance measured vertically rather than
    perpendicular to the chord.

    A QQ curve is single-valued and strictly increasing in x, and the renderer
    joins the kept points with straight segments, so the vertical gap to the
    chord is exactly the error a reader sees. Bounding that directly makes
    tol_y a guarantee; the usual perpendicular distance is smaller than the
    vertical one on a sloped chord and would let the error past it. Iterative
    rather than recursive: the input is tens of thousands of points and
    Python's stack is not.
    """
    assert len(x) == len(y), f"x has {len(x)} points, y has {len(y)}"
    assert np.all(np.diff(x) > 0), "x must be strictly increasing for a vertical metric"
    keep = np.zeros(len(x), dtype=bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(x) - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        seg = slice(i + 1, j)
        chord = y[i] + (y[j] - y[i]) * (x[seg] - x[i]) / (x[j] - x[i])
        dist = np.abs(y[seg] - chord)
        k = int(np.argmax(dist))
        if dist[k] > tol_y:
            k += i + 1
            keep[k] = True
            stack.append((i, k))
            stack.append((k, j))
    return keep


def expected(n):
    """Expected uniform quantile for each rank, i.e. the QQ x-coordinate."""
    return (np.arange(1, n + 1) - 0.5) / n


def curve(v, form):
    """(x, y) of the full QQ curve for one cell type, every p-value included."""
    exp = expected(len(v))
    return exp, (v if form == "classic" else v - exp)


def reduce_curve(v, form, tol_y):
    """The drawn subset of the curve: rank-thin, then simplify."""
    x, y = curve(v, form)
    keep = thin_ranks(len(v))
    x, y = x[keep], y[keep]
    mask = simplify(x, y, tol_y)
    return x[mask], y[mask]


def check_reduction(pvals, form, tol_y, half_stroke_y, linewidth_pt=1.3, dpi=300,
                    tag=""):
    """Assert the drawn stroke covers the true curve everywhere.

    The criterion is the one that decides whether a reader can see a difference:
    the full-resolution curve must lie within half a linewidth of the polyline
    actually drawn, so the stroke on the page paints over it. Reported in
    rendered pixels as well, since that is the intuitive unit.
    """
    worst, worst_ct, drawn, total = 0.0, None, 0, 0
    for ct, v in pvals.items():
        x, y = curve(v, form)
        xr, yr = reduce_curve(v, form, tol_y)
        err = float(np.abs(np.interp(x, xr, yr) - y).max())
        drawn += len(xr)
        total += len(x)
        if err > worst:
            worst, worst_ct = err, ct
    assert worst < half_stroke_y, (
        f"{tag}reduction is visible: worst error {worst:.3e} ({worst_ct}) exceeds "
        f"half a linewidth {half_stroke_y:.3e} in data units"
    )
    px = half_stroke_y * 2 / (linewidth_pt / 72 * dpi)  # data units per rendered pixel
    print(f"{tag}downsample: {drawn} of {total} points drawn "
          f"({100 * drawn / total:.3f}%); worst departure {worst:.3e} data units "
          f"= {worst / px:.2f} px at {dpi} dpi, half a linewidth is "
          f"{half_stroke_y / px:.2f} px")
