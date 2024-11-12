import glob

# n, m, density(%)
info = {
    "cycle300": [300, 300, 0.669, 150],
    "jagmesh1": [936, 2664, 0.609, 50],
    "dwt_1005": [1005, 3808, 0.755, 100],
    "btree9": [1023, 1022, 0.196, 150],
    "1138_bus": [1138, 1458, 0.225, 150],
    "dwt_2680": [2680, 11173, 0.311, 150],
    "3elt": [4720, 13722, 0.123, 150],
}


def generate_latex_code(matrixNames):
    latex_code = f"""
\\begin{{figure*}}[btp]
  \\centering
  \\addtolength{{\\tabcolsep}}{{-0.5em}}
  \\begin{{tabular}}{{cccccc}}"""

    for matrixName in matrixNames:
        n, m, density, it = info[matrixName]
        matrixNameRep = matrixName.replace("_", "\\_")
        extension = "png"
        latex_code += f"""
    \\multicolumn{{6}}{{c}}{{\\textbf{{\\texttt{{{matrixNameRep}}}}} $(\\abs{{V}}={n}, \\abs{{E}}={m}, \\text{{sparsity}}={density:.3f}\\text{{\\%}}) \quad Figures are at {it} iterations.$}} \\\\
    \\raisebox{{-.5\\height}}{{\\includegraphics[width=0.55\\columnwidth]{{individual/plot/{matrixName}.pdf}}}} &
    \\makecell{{\\small{{\\textsf{{FR}}}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{individual/vis/{matrixName}_FR.{extension}}}}} &
    \\makecell{{\\small{{\\textsf{{L-BFGS}}}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{individual/vis/{matrixName}_L-BFGS.{extension}}}}} &
    \\makecell{{\\small{{\\textsf{{CN-FR}}}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{individual/vis/{matrixName}_CN-FR.{extension}}}}} &
    \\makecell{{\\small{{\\textsf{{CN-L-BFGS}}}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{individual/vis/{matrixName}_CN-L-BFGS.{extension}}}}} &
    \\makecell{{\\small{{\\textsf{{BEST}}}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{individual/vis/opt_{matrixName}.{extension}}}}} \\\\
"""

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
