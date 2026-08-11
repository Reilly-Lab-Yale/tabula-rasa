"""How many cell types does a single MPRA study assay?
"""
import pathlib
import sqlite3

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
DB = HERE / "data" / "mprabase_v4_9.3.db"

# Okabe-Ito blue, matching the manuscript figure palette (scripts/figure_styling.py).
BLUE = "#0072b2"
INK = "#1a1a1a"
MUTED = "#6b6b6b"

# A sample's study is the prefix of its sample_id ("DS0030-SID02" -> "DS0030").
# The documented route is sample.library_id -> designed_library.datasets_id, but
# designed_library is missing two libraries, which would silently drop their
# studies. The prefix agrees with the join wherever both resolve (asserted below).
QUERY = """
SELECT SUBSTR(sample_id, 1, 6) AS study,
       COUNT(DISTINCT Cell_line_tissue) AS n_cell_types,
       COUNT(*) AS n_samples
FROM sample
GROUP BY study
ORDER BY n_cell_types DESC, study
"""

VALIDATE = """
SELECT SUM(d.datasets_id = SUBSTR(s.sample_id, 1, 6)) AS agree,
       SUM(d.datasets_id <> SUBSTR(s.sample_id, 1, 6)) AS disagree
FROM sample s JOIN designed_library d ON s.library_id = d.library_id
"""


def load():
    assert DB.exists(), f"missing {DB}; see analyses/mprabase/README.md for the download"
    con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    try:
        agree, disagree = con.execute(VALIDATE).fetchone()
        assert disagree == 0, (
            f"sample_id prefix disagrees with the library join on {disagree} rows")
        rows = con.execute(QUERY).fetchall()
        n_samples = con.execute("SELECT COUNT(*) FROM sample").fetchone()[0]
        n_cells = con.execute(
            "SELECT COUNT(DISTINCT Cell_line_tissue) FROM sample").fetchone()[0]
    finally:
        con.close()

    assert rows, "no studies returned"
    assert sum(r[2] for r in rows) == n_samples, (
        f"grouping lost samples: {sum(r[2] for r in rows)} != {n_samples}")
    assert min(r[1] for r in rows) >= 1, "a study with no cell type"
    # The paper reports 130 experiments across 35 cell types; if these drift, the
    # downloaded database is not the version this analysis was written against.
    assert (n_samples, n_cells) == (130, 35), (
        f"unexpected database contents: {n_samples} samples, {n_cells} cell types")
    print(f"{len(rows)} studies, {n_samples} samples, {n_cells} distinct cell types")
    return rows


def plot(rows, out):
    counts = {}
    for _, n, _ in rows:
        counts[n] = counts.get(n, 0) + 1
    xs = sorted(counts)
    ys = [counts[x] for x in xs]
    total = sum(ys)

    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    ax.bar(xs, ys, width=0.62, color=BLUE, zorder=3)

    for x, y in zip(xs, ys):
        ax.text(x, y + total * 0.015, f"{y}\n{100 * y / total:.0f}%",
                ha="center", va="bottom", fontsize=9, color=INK, linespacing=1.35)

    ax.set_xlabel("Distinct cell types assayed", fontsize=10, color=INK)
    ax.set_ylabel("Studies", fontsize=10, color=INK)
    ax.set_title("Cell types per MPRA study", fontsize=11, color=INK, pad=10)
    ax.set_xticks(xs)
    ax.set_ylim(0, max(ys) * 1.28)
    ax.yaxis.grid(True, color="#e6e6e6", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#cccccc")
    ax.tick_params(colors=MUTED, labelsize=9, length=0)
    fig.text(0.01, 0.015, f"MPRAbase v4.9.3 (2009-2023), n = {total} studies",
             fontsize=7.5, color=MUTED)

    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    print(f"wrote {out.name} and {out.with_suffix('.pdf').name}")
    return xs, ys


def main():
    rows = load()
    xs, ys = plot(rows, HERE / "cell_types_per_study.png")

    csv = HERE / "cell_types_per_study.csv"
    csv.write_text("study,n_cell_types,n_samples\n"
                   + "".join(f"{s},{n},{m}\n" for s, n, m in rows))
    print(f"wrote {csv.name}")

    total = sum(ys)
    print("\ndistribution:")
    for x, y in zip(xs, ys):
        print(f"  {x} cell type{'s' if x > 1 else ' '}: {y:2d} studies "
              f"({100 * y / total:4.1f}%)")
    print(f"  max: {max(xs)}   mean: {sum(x * y for x, y in zip(xs, ys)) / total:.2f}")


if __name__ == "__main__":
    main()
