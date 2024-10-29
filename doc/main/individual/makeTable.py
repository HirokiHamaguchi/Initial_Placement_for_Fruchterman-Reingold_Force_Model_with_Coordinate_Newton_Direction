def generate_latex_code(matrixNames):
    latex_code = f"""
\\begin{{table*}}
  \\centering
  \\addtolength{{\\tabcolsep}}{{-0.5em}}
  \\begin{{tabular}}{{ccccc}}"""

    for matrixName in matrixNames:
        matrixNameRep = matrixName.replace("_", "\\_")
        latex_code += f"""
    \\multicolumn{{4}}{{c}}{{\\textbf{{{matrixNameRep}}}}}\\\\
    \\raisebox{{-.5\\height}}{{\\includegraphics[width=0.58\\columnwidth]{{individual/plot/{matrixName}.pdf}}}} &
    \\makecell{{\\small{{FR}}\\\\[-0.2em]\\includegraphics[width=0.29\\columnwidth]{{individual/viz/{matrixName}_FR.png}}}} &
    \\makecell{{\\small{{L\\_BFGS}}\\\\[-0.2em]\\includegraphics[width=0.29\\columnwidth]{{individual/viz/{matrixName}_L_BFGS.png}}}} &
    \\makecell{{\\small{{RS\\_FR}}\\\\[-0.2em]\\includegraphics[width=0.29\\columnwidth]{{individual/viz/{matrixName}_RS_FR.png}}}} &
    \\makecell{{\\small{{RS\\_L\\_BFGS}}\\\\[-0.2em]\\includegraphics[width=0.29\\columnwidth]{{individual/viz/{matrixName}_RS_L_BFGS.png}}}} \\\\"""

    latex_code += f"""
  \\end{{tabular}}
  \\caption{{graphs.}}
\\end{{table*}}
"""
    return latex_code


import glob

matrixNames = [
    f.split("/")[-1].split(".")[0] for f in glob.glob("doc/main/individual/plot/*.pdf")
]

print(generate_latex_code(matrixNames))
