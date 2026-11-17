"""Output helpers shared by the scripts: a run header and a comparison with the printed values."""
import datetime
import platform

import numpy
import scipy


def header(title):
    now = datetime.datetime.now(datetime.timezone.utc)
    print(title)
    print(f"run {now:%Y-%m-%d %H:%M} UTC | Python {platform.python_version()} | numpy {numpy.__version__} | "
          f"scipy {scipy.__version__} | {platform.system()} {platform.machine()}")
    print()


class Published:
    """Collects computed values next to the values printed in the paper."""

    def __init__(self):
        self.rows = []

    def add(self, label, computed, printed, decimals):
        gap = abs(computed - printed)
        unit = 10.0 ** (-decimals)
        if gap <= 0.5 * unit * (1 + 1e-9):
            verdict = "same"
        elif gap <= 1.5 * unit:
            verdict = "last digit off by one"
        else:
            verdict = "DIFFERS"
        self.rows.append((label, computed, printed, decimals, verdict))

    def report(self):
        """Print the comparison. Returns the number of values that differ by more than one last digit."""
        print()
        print("comparison with the values printed in the paper")
        width = max(len(r[0]) for r in self.rows)
        for label, computed, printed, decimals, verdict in self.rows:
            print(f"  {label:<{width}}  computed {computed:.{decimals + 2}f}  printed {printed:.{decimals}f}  {verdict}")
        same = sum(r[4] == "same" for r in self.rows)
        near = sum(r[4] == "last digit off by one" for r in self.rows)
        bad = len(self.rows) - same - near
        print(f"  {same} of {len(self.rows)} round to the printed value, {near} differ by one in the last digit, {bad} differ by more")
        return bad
