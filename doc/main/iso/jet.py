from typing import Tuple

import matplotlib
import numpy as np


# simply implement LinearSegmentedColormap
# https://matplotlib.org/stable/api/_as_gen/matplotlib.colors.LinearSegmentedColormap.html
def getJet(val: float) -> Tuple[float, float, float]:
    val = float(np.clip(val, 0.0, 1.0))
    cmap = matplotlib.colormaps["viridis"]
    color = cmap(val)
    return (color[0], color[1], color[2])


def main():
    for i, val in enumerate(np.linspace(0, 1, 10)):
        myColor = getJet(val)
        print(
            f"\\definecolor{{c{i}}}{{rgb}}{{{myColor[0]:.2f},{myColor[1]:.2f},{myColor[2]:.2f}}}"
        )


if __name__ == "__main__":
    main()
