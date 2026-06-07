import os
from pathlib import Path

import numpy as np
import scipy.io

its = {
    "FR": [150, 50, 100, 150, 150, 150, 150],
    "HC": [150, 150, 150, 150, 150, 150, 150],
    "Eades": [150, 150, 150, 150, 150, 150, 150],
}


DATA_DIR = Path(__file__).parent.parent.parent.parent / "data"


# get preamble from individual.tex
def get_preamble():
    with open(Path(__file__).parent / "individual.tex", "r") as f:
        lines = f.readlines()
    preamble = []
    for line in lines:
        if line.strip() == r"\begin{document}":
            preamble.append(line)
            preamble.append("\n")
            break
        preamble.append(line)
    return "".join(preamble)


def generate_latex_code(matrixNames, mode):
    latex_code = r"\begin{tabular}{cccccc}"

    for matrixNameIdx, matrixName in enumerate(matrixNames):
        mat = scipy.io.mmread(DATA_DIR / f"{matrixName}.mtx")
        n = mat.shape[0] if mat.shape is not None else -1
        m = (mat.nnz - np.count_nonzero(mat.diagonal())) // 2  # type: ignore
        assert n > 0 and m > 0
        density = 100 * m / (n * (n - 1) / 2)
        it_FR = its["FR"][matrixNameIdx]
        it_HC = its["HC"][matrixNameIdx]
        it_Eades = its["Eades"][matrixNameIdx]
        assert it_HC == it_Eades
        matrixNameRep = matrixName.replace("_", "\\_")
        extension = "png"
        if mode == "FR":
            algName = "FR"
            latex_code += f"""
  \\multicolumn{{6}}{{c}}{{\\large \\textbf{{\\texttt{{{matrixNameRep}}}}} $(\\abs{{V}}={n}, \\abs{{E}}={m}, \\text{{density}}={density:.3f}\\text{{\\%}})$ \\quad Figures are at {it_FR} iterations.}} \\\\
  \\raisebox{{-.5\\height}}{{\\includegraphics[width=0.55\\columnwidth]{{plot/{matrixName}_{algName}.pdf}}}} &
  \\makecell{{\\textsf{{RAND-FR}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_FR_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}}-FR (proposed)}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN-FR_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{RAND-L-BFGS}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_L-BFGS_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}}-L-BFGS (proposed)}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN-L-BFGS_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{BEST}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/opt_{matrixName}_{algName}.{extension}}}}} \\\\
"""
        elif mode == "HC_and_Eades":
            latex_code += f"""
  \\multicolumn{{3}}{{c}}{{\\large \\textbf{{\\texttt{{{matrixNameRep}}} \\quad HC}} model \\quad Figures are at {it_HC} iterations.}} & \\multicolumn{{3}}{{c}}{{\\large \\textbf{{\\texttt{{{matrixNameRep}}} \\quad Eades}} model \\quad Figures are at {it_HC} iterations.}} \\\\
  \\raisebox{{-.5\\height}}{{\\includegraphics[width=0.55\\columnwidth]{{plot/{matrixName}_HC.pdf}}}} &
  \\makecell{{\\textsf{{RAND-L-BFGS}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_L-BFGS_HC.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}}-L-BFGS (proposed)}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN-L-BFGS_HC.{extension}}}}} &
  \\raisebox{{-.5\\height}}{{\\includegraphics[width=0.55\\columnwidth]{{plot/{matrixName}_Eades.pdf}}}} &
  \\makecell{{\\textsf{{RAND-L-BFGS}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_L-BFGS_Eades.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}}-L-BFGS (proposed)}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN-L-BFGS_Eades.{extension}}}}} \\\\
"""
        elif mode == "2":
            algName = "2"
            latex_code += f"""
  \\multicolumn{{6}}{{c}}{{\\large \\textbf{{\\texttt{{{matrixNameRep}}}}} $(\\abs{{V}}={n}, \\abs{{E}}={m}, \\text{{density}}={density:.3f}\\text{{\\%}})$ \\quad Figures are at {it_FR} iterations.}} \\\\
  \\raisebox{{-.5\\height}}{{\\includegraphics[width=0.55\\columnwidth]{{plot/{matrixName}_{algName}.pdf}}}} &
  \\makecell{{\\textsf{{RAND-FR}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_FR_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}}-FR (proposed)}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN-FR_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{RAND-L-BFGS}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_L-BFGS_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}}-L-BFGS (proposed)}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN-L-BFGS_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{BEST}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/opt_{matrixName}_{algName}.{extension}}}}} \\\\
"""

    latex_code += """\\end{tabular}"""

    preamble = get_preamble()
    full_latex_code = preamble + latex_code + "\n\\end{document}"

    with open(Path(__file__).parent / f"table_{mode}.tex", "w") as f:
        f.write(full_latex_code)

    for i in range(2):
        os.chdir(Path(__file__).parent)
        os.system(
            f"pdflatex -interaction=nonstopmode -synctex=1 -file-line-error --shell-escape {Path(__file__).parent / f'table_{mode}.tex'}"
        )


def main():
    matrixNames1 = [
        "cycle300",
        "jagmesh1",
        "btree9",
        "1138_bus",
        "dwt_1005",
        "dwt_2680",
        "3elt",
    ]
    matrixNames2 = [
        "jagmesh8",
        "bcsstk14",
        "bcsstk15",
        "bcsstk16",
        "USpowerGrid",
        "bcspwr10",
        "wiki-Vote",
    ]

    generate_latex_code(matrixNames1, mode="FR")
    generate_latex_code(matrixNames1, mode="HC_and_Eades")
    generate_latex_code(matrixNames2, mode="2")


if __name__ == "__main__":
    main()
