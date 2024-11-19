import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

with open("main.tex", mode="r") as f:
    lines = f.readlines()
    for i in range(len(lines)):
        if (
            lines[i].startswith("\\section{")
            or lines[i].startswith("\\subsection{")
            or lines[i].startswith("\\subsubsection{")
        ):
            bra = lines[i].find("{")
            ket = lines[i].find("}")
            print(
                f"line {i}: {lines[i][bra + 1:ket]} -> {lines[i][bra + 1:ket].title()}"
            )
            lines[i] = (
                lines[i][: bra + 1] + lines[i][bra + 1 : ket].title() + lines[i][ket:]
            )
