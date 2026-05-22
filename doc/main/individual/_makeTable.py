import pyperclip

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

its = {
    "FR": [150, 50, 100, 150, 150, 150, 150],
    "HC": [150, 150, 150, 150, 150, 150, 150],
    "Eades": [150, 150, 150, 150, 150, 150, 150],
}


def generate_latex_code(matrixNames, mode="FR"):
    latex_code = r"\begin{tabular}{cccccc}"

    for matrixNameIdx, matrixName in enumerate(matrixNames):
        n, m, density = info[matrixName]
        it_FR = its["FR"][matrixNameIdx]
        it_HC = its["HC"][matrixNameIdx]
        it_Eades = its["Eades"][matrixNameIdx]
        assert it_HC == it_Eades
        matrixNameRep = matrixName.replace("_", "\\_")
        extension = "png"
        if mode == "FR":
            algName = "FR"
            latex_code += f"""
  \\multicolumn{{6}}{{c}}{{\\large \\textbf{{\\texttt{{{matrixNameRep}}}}} $(\\abs{{V}}={n}, \\abs{{E}}={m}, \\text{{sparsity}}={density:.3f}\\text{{\\%}})$ \\quad Figures are at {it_FR} iterations.}} \\\\
  \\raisebox{{-.5\\height}}{{\\includegraphics[width=0.55\\columnwidth]{{plot/{matrixName}_{algName}.pdf}}}} &
  \\makecell{{\\textsf{{RAND-FR}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_FR_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}}-FR (proposed)}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN-FR_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{RAND-L-BFGS}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_L-BFGS_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}}-L-BFGS (proposed)}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN-L-BFGS_{algName}.{extension}}}}} &
  \\makecell{{\\textsf{{BEST}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/opt_{matrixName}_{algName}.{extension}}}}} \\\\
"""
        else:
            latex_code += f"""
  \\multicolumn{{6}}{{c}}{{\\large \\textbf{{\\texttt{{{matrixNameRep}}}}} $(\\abs{{V}}={n}, \\abs{{E}}={m}, \\text{{sparsity}}={density:.3f}\\text{{\\%}})$ \\quad Figures are at {it_HC} iterations.}} \\\\
  \\raisebox{{-.5\\height}}{{\\includegraphics[width=0.55\\columnwidth]{{plot/{matrixName}_Eades.pdf}}}} &
  \\makecell{{\\textsf{{RAND-L-BFGS}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_L-BFGS_Eades.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}}-L-BFGS (proposed)}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN-L-BFGS_Eades.{extension}}}}} &
  \\raisebox{{-.5\\height}}{{\\includegraphics[width=0.55\\columnwidth]{{plot/{matrixName}_HC.pdf}}}} &
  \\makecell{{\\textsf{{RAND-L-BFGS}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_L-BFGS_HC.{extension}}}}} &
  \\makecell{{\\textsf{{\\textbf{{CN}}-L-BFGS (proposed)}}\\\\[-0.2em]\\includegraphics[width=0.27\\columnwidth]{{vis/{matrixName}_CN-L-BFGS_HC.{extension}}}}} \\\\
"""
    latex_code += """\\end{tabular}"""
    return latex_code


def main():
    matrixNames = list(info.keys())
    matrixNames.sort(key=lambda x: info[x][0])
    code = generate_latex_code(matrixNames, mode="FR")
    pyperclip.copy(code)
    code = generate_latex_code(matrixNames, mode="HC_and_Eades")
    pyperclip.copy(code)
    print("LaTeX code has been copied to the clipboard.")


if __name__ == "__main__":
    main()
