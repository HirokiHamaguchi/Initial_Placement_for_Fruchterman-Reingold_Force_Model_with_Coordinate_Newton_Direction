import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def create_legend() -> None:
    plt.rcParams["font.family"] = "Serif"
    plt.rcParams["font.size"] = 70
    plt.rcParams["text.usetex"] = True

    cmap = matplotlib.colormaps.get_cmap("viridis")

    fig, ax = plt.subplots(figsize=(15, 0.8))

    gradient = np.linspace(0, 1, 256).reshape(1, -1)

    ax.imshow(
        gradient,
        aspect="auto",
        cmap=cmap,
        extent=(1, 2, 0, 1),
    )

    # 軸を完全に消す
    ax.set_xticks([])
    ax.set_yticks([])

    # 左ラベル
    ax.text(
        0.95,
        0.45,
        "$1$",
        ha="right",
        va="center",
        transform=ax.transData,
    )

    # 右ラベル
    ax.text(
        2.05,
        0.45,
        "$n$",
        ha="left",
        va="center",
        transform=ax.transData,
    )

    # 説明テキスト
    ax.text(
        0.7,
        0.45,
        "vertex indices:",
        ha="right",
        va="center",
        transform=ax.transData,
    )

    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.savefig(
        "legend.pdf",
        bbox_inches="tight",
        pad_inches=0.02,
    )
    plt.close(fig)


if __name__ == "__main__":
    os.chdir(os.path.dirname(__file__))
    create_legend()
