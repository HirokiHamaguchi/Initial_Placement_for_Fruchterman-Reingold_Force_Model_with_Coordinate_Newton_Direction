import glob

# n, m, density(%)
info = {
    "cycle300": [300, 300, 0.669],
    "jagmesh1": [936, 2664, 0.609],
    "dwt_1005": [1005, 3808, 0.755],
    "btree9": [1023, 1022, 0.196],
    "1138_bus": [1138, 1458, 0.225],
    "dwt_2680": [2680, 11173, 0.311],
    "3elt": [4720, 13722, 0.123],
}


def generate_latex_code(matrixNames):
    latex_code = f"""
\\begin{{figure*}}[btp]
  \\centering
  \\addtolength{{\\tabcolsep}}{{-0.5em}}
  \\begin{{tabular}}{{ccccc}}"""

    for matrixName in matrixNames:
        n, m, density = info[matrixName]
        matrixNameRep = matrixName.replace("_", "\\_")
        extension = (
            "png"
            if any(name in matrixName for name in ["jagmesh1", "dwt_1005", "btree9"])
            else "pdf"
        )
        latex_code += f"""
    \\multicolumn{{5}}{{c}}{{\\textbf{{\\texttt{{{matrixNameRep}}}}} $(\\abs{{V}}={n}, \\abs{{E}}={m}, \\text{{sparsity}}={density:.3f}\\text{{\\%}})$}} \\\\
    \\raisebox{{-.5\\height}}{{\\includegraphics[width=0.55\\columnwidth]{{individual/plot/{matrixName}.pdf}}}} &
    \\makecell{{\\small{{\\textsf{{FR}}}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{individual/viz/{matrixName}_FR.{extension}}}}} &
    \\makecell{{\\small{{\\textsf{{L-BFGS}}}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{individual/viz/{matrixName}_L_BFGS.{extension}}}}} &
    \\makecell{{\\small{{\\textsf{{CN-FR}}}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{individual/viz/{matrixName}_CN_FR.{extension}}}}} &
    \\makecell{{\\small{{\\textsf{{CN-L-BFGS}}}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{individual/viz/{matrixName}_CN_L_BFGS.{extension}}}}} \\\\"""

    latex_code += f"""
  \\end{{tabular}}
  \\caption{{Numerical experiment results for various graphs.  Please refer Sec.~\\ref{{ssec:exprDetail}} for details.}}
  \\label{{fig:individual}}
\\end{{figure*}}
"""
    return latex_code


def main():
    matrixNames = [
        f.split("/")[-1].split(".")[0]
        for f in glob.glob("doc/main/individual/plot/*.pdf")
    ]
    assert all(matrixName in info for matrixName in matrixNames)
    matrixNames.sort(key=lambda x: info[x][0])

    print(generate_latex_code(matrixNames))


if __name__ == "__main__":
    main()
