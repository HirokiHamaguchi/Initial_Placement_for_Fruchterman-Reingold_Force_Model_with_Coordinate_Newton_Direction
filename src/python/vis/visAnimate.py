import glob
import os

import networkx as nx
import numpy as np
import PIL
from packaging import version
from tqdm.auto import tqdm

from src.python.vis.visGraph import visGraph


def visAnimate(G: nx.Graph, optResults: np.ndarray, matrixName: str, methodName: str):
    folderName = "FruchtermanReingoldByRandomSubspace"
    assert os.getcwd().find(folderName) != -1
    cwd = os.getcwd()
    pathToFolder = os.getcwd().split(folderName)[0] + folderName
    os.chdir(pathToFolder)

    assert os.path.exists("temp"), "temp directory does not exist"
    for f in glob.glob("temp/*"):
        os.remove(f)

    for i, pos_and_grad in tqdm(
        enumerate(optResults), total=len(optResults), leave=False
    ):
        visGraph(
            G,
            *pos_and_grad,
            title=f"{matrixName}_{methodName}_{i}",
            savePath=f"temp/{i:03}.png",
        )

    if version.parse(PIL.__version__) < version.parse("3.4"):
        # https://stackoverflow.com/questions/24688802/saving-an-animated-gif-in-pillow
        raise Exception(
            "Pillow in version not supporting making animated GIFs",
            "you need to upgrade library version",
            "see release notes in",
            "https://pillow.readthedocs.io/en/latest/releasenotes/3.4.0.html#append-images-to-gif",
        )
    else:
        files = sorted(glob.glob(f"temp/*.png"))
        images = [PIL.Image.open(file) for file in files]
        images[0].save(
            f"movie/{matrixName}_{methodName}.gif",
            save_all=True,
            append_images=images[1:],
            duration=100,
            loop=0,
        )

    os.chdir(cwd)
    print(os.path.join(pathToFolder, "movie", f"{matrixName}_{methodName}.gif"))
