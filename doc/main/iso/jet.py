import numpy as np
import matplotlib
from typing import Tuple


# simply implement LinearSegmentedColormap
# https://matplotlib.org/stable/api/_as_gen/matplotlib.colors.LinearSegmentedColormap.html
def getJet(val: float) -> Tuple[float, float, float]:
    # https://github.com/matplotlib/matplotlib/blob/main/lib/matplotlib/_cm.py#L243
    _jet_data = {
        "red": (
            (0.00, 0, 0),
            (0.35, 0, 0),
            (0.66, 1, 1),
            (0.89, 1, 1),
            (1.00, 0.5, 0.5),
        ),
        "green": (
            (0.000, 0, 0),
            (0.125, 0, 0),
            (0.375, 1, 1),
            (0.640, 1, 1),
            (0.910, 0, 0),
            (1.000, 0, 0),
        ),
        "blue": (
            (0.00, 0.5, 0.5),
            (0.11, 1, 1),
            (0.34, 1, 1),
            (0.65, 0, 0),
            (1.00, 0, 0),
        ),
    }

    floats = []
    for data in _jet_data.values():
        for i in range(len(data) - 1):
            if data[i][0] <= val <= data[i + 1][0]:
                ratio = (val - data[i][0]) / (data[i + 1][0] - data[i][0])
                floats.append(data[i][2] + (data[i + 1][1] - data[i][2]) * ratio)
                break
        else:
            assert False

    return tuple(floats)


def main():
    for i, val in enumerate(np.linspace(0, 1, 10)):
        myColor = getJet(val)
        print(
            f"\\definecolor{{c{i}}}{{rgb}}{{{myColor[0]:.2f},{myColor[1]:.2f},{myColor[2]:.2f}}}"
        )


if __name__ == "__main__":
    main()
