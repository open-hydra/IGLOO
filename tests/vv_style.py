"""Shared V&V plot style -- Computer Modern (LaTeX) typography for all overlays.

Single source of truth for every case's verify.py and tools/plot_curves.py: each
script puts the tests root on sys.path and calls vv_style.apply() before pyplot
import. Edit here and rerun tests/test.sh to restyle every overlay, current and
future.

Typography is Computer Modern (mathtext.fontset "cm"), baked to vector paths via
`svg.fonttype = "path"`. This is LOAD-BEARING: with `svg.fonttype = "none"` the
math glyphs stay as <text> in the TeX fonts cmmi10/cmsy10, whose cmap is NOT
Unicode -- browsers then render wrong glyphs (comma -> colon, broken operators),
even though matplotlib's own raster is correct (it resolves glyphs by index).
Baking to paths makes the SVG render identically everywhere with no font
dependency. matplotlib wraps each label in `<g id="text_...">`, so the docs site
still recolors plot text with the theme via
`docs/stylesheets/extra.css` (`svg g[id^="text_"] { fill: currentColor }`) without
touching the data lines/markers.

Output is BYTE-REPRODUCIBLE: the 34 overlays are tracked in git, so a re-render
must not show up as a diff unless a curve actually moved. Two sources of per-run
churn are pinned here -- `svg.hashsalt` (matplotlib otherwise seeds clip-path and
glyph-def ids from a random uuid, ~180 changed lines per figure) and the embedded
`<dc:date>` (suppressed by defaulting `metadata={'Date': None}` on savefig, see
`_patch_savefig`). Without both, every `tests/test.sh all` dirtied all 34 files
and a real change was invisible in the noise.
"""
import matplotlib as mpl

STYLE = {
    # transparent, theme-aware (as MOSE)
    "figure.facecolor": "none", "axes.facecolor": "none",
    "savefig.facecolor": "none",
    # bake glyphs to paths -- correct CM in every browser (see module docstring)
    "svg.fonttype": "path",
    # Fixed salt for the ids matplotlib hashes into clip-paths and glyph defs;
    # default None means a fresh uuid per process => byte-churn on every render.
    "svg.hashsalt": "igloo-vv",
    # LaTeX / Computer Modern typography
    "mathtext.fontset": "cm",
    "font.family": "serif",
    "font.serif": ["cmr10", "Computer Modern Roman", "DejaVu Serif"],
    "axes.formatter.use_mathtext": True,
    # Use the ASCII minus on axes: Computer Modern (cmr10) has NO Unicode minus
    # (U+2212), so the default would bake a "missing glyph" box on every negative
    # tick label. U+002D renders correctly and is baked to a path like the rest.
    "axes.unicode_minus": False,
    # Font sizes -- single source of truth: keep per-verify.py legend()/set_title()
    # calls free of explicit fontsize= so these govern every overlay uniformly.
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "figure.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    # Opaque WHITE legend box (was transparent: legends inherited axes.facecolor
    # "none"). The docs recolor plot text to the theme, so extra.css pins legend
    # text dark to stay readable on the white fill in dark mode.
    "legend.frameon": True,
    "legend.facecolor": "white",
    "legend.edgecolor": "0.7",
    "legend.framealpha": 1.0,
}


_savefig_patched = False


def _patch_savefig():
    """Default `metadata={'Date': None}` on SVG saves -- kills the `<dc:date>` line.

    There is no rcParam for it, and patching once here beats threading the kwarg
    through 22 identical `fig.savefig(out, ...)` call sites (and every future one).
    Only SVG targets are touched: PNG metadata takes different keys.
    """
    global _savefig_patched
    if _savefig_patched:
        return
    from matplotlib.figure import Figure
    original = Figure.savefig

    def savefig(self, fname, *args, **kwargs):
        if "metadata" not in kwargs and str(fname).lower().endswith(".svg"):
            kwargs["metadata"] = {"Date": None}
        return original(self, fname, *args, **kwargs)

    Figure.savefig = savefig
    _savefig_patched = True


def apply():
    mpl.rcParams.update(STYLE)
    _patch_savefig()
