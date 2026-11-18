"""Style shared by the figure scripts: one IEEEtran conference column wide, vector PDF output."""
import pathlib

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt

COL = 252.0 / 72.27   # IEEEtran conference column width, 252 pt, in inches
ROOT = pathlib.Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
RESULTS = ROOT / "results"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "Times New Roman", "Nimbus Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.1,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.01,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "pdf.fonttype": 42,
})


def frame(*axes, lw=0.25):
    """Thin box on all four sides, ticks pointing inwards."""
    for ax in axes:
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(True)
        for spine in ax.spines.values():
            spine.set_linewidth(lw)
        ax.tick_params(direction="in", width=lw)


def save(fig, name):
    """Write figures/<name>.pdf, the file name that the paper source uses."""
    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / f"{name}.pdf")
    plt.close(fig)
    print(f"wrote figures/{name}.pdf")
