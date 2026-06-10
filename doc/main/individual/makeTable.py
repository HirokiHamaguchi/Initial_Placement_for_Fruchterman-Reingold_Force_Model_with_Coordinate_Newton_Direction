import os
from pathlib import Path

import numpy as np
import scipy.io

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
    latex_code = (
        r"\begin{tabular}{cc|cc|cc}"
        if mode == "FR" or mode == "2"
        else r"\begin{tabular}{ccc|ccc}"
    )

    for matrixNameIdx, matrixName in enumerate(matrixNames):
        mat = scipy.io.mmread(DATA_DIR / f"{matrixName}.mtx")
        n = mat.shape[0] if mat.shape is not None else -1
        m = (mat.nnz - np.count_nonzero(mat.diagonal())) // 2  # type: ignore
        assert n > 0 and m > 0
        matrixNameRep = matrixName.replace("_", "\\_")
        extension = "png"
        if mode == "FR" or mode == "2":
            algName = mode
            modePrefix = "" if mode == "FR" else "_2"
            latex_code += f"""
  \\multicolumn{{6}}{{c}}{{\\large \\textbf{{\\texttt{{{matrixNameRep}}}}} $(\\abs{{V}}={n}, \\abs{{E}}={m})$}} \\\\
  \\raisebox{{-.5\\height}}{{\\includegraphics[width=0.55\\columnwidth]{{plot/{matrixName}_{algName}.pdf}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}} (proposed)}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN{modePrefix}.{extension}}}}} &
  \\makecell{{\\textsf{{RAND-SIM}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_FR_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}}-SIM}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN-FR_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{RAND-L-BFGS}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_L-BFGS_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}}-L-BFGS}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN-L-BFGS_{algName}.{extension}}}}} \\\\
"""
        elif mode == "HC_and_Eades":
            latex_code += f"""
  \\multicolumn{{3}}{{c}}{{\\large \\textbf{{\\texttt{{{matrixNameRep}}} \\quad HC}} model}} & \\multicolumn{{3}}{{c}}{{\\large \\textbf{{\\texttt{{{matrixNameRep}}} \\quad Eades}} model}} \\\\
  \\makecell{{\\textsf{{\\textbf{{CN}} (proposed)}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN_HC.{extension}}}}} &
  \\makecell{{\\textsf{{RAND-L-BFGS}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_L-BFGS_HC.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}}-L-BFGS}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN-L-BFGS_HC.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}} (proposed)}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN_Eades.{extension}}}}} &
  \\makecell{{\\textsf{{RAND-L-BFGS}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_L-BFGS_Eades.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}}-L-BFGS}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN-L-BFGS_Eades.{extension}}}}} \\\\
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
        "cylinder_30_30",  # 900
        "gr_30_30",  # 900
        "sierpinski_06",  # 1095
        "jagmesh8",  # 1141
        "USpowerGrid",  # 4941
        "wiki-Vote",  # 7066 (original: 8297)
        "crack",  # 10240
    ]

    generate_latex_code(matrixNames1, mode="FR")
    generate_latex_code(matrixNames1, mode="HC_and_Eades")
    generate_latex_code(matrixNames2, mode="2")


if __name__ == "__main__":
    main()
