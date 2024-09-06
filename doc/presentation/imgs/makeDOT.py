import os
import numpy as np
import matplotlib.pyplot as plt

num_nodes = 300

colors = plt.cm.jet(np.linspace(0, 1, num_nodes))
hex_colors = [
    "#%02x%02x%02x" % (int(c[0] * 255), int(c[1] * 255), int(c[2] * 255))
    for c in colors
]
dot_content = 'graph G {\ngraph [dpi=30];\n    node [shape=circle, style=filled, width=2, fixedsize=true, label="", penwidth=0];\n'
for i in range(1, num_nodes):
    dot_content += f"    v{i} -- v{i+1};\n"
dot_content += f"    v{num_nodes} -- v1;\n"
for i in range(1, num_nodes + 1):
    dot_content += f'    v{i} [fillcolor="{hex_colors[i-1]}"];\n'
dot_content += "}"

dot_file_path = __file__.replace("makeDOT.py", "circle.dot")
with open(dot_file_path, "w") as f:
    f.write(dot_content)


print(f"DOT file generated at {dot_file_path}")
print("now run the following command to generate the image:")
print(f"fdp -Tsvg {dot_file_path} -o {dot_file_path.replace('.dot', '_fdp.svg')}")
os.system(f"fdp -Tsvg {dot_file_path} -o {dot_file_path.replace('.dot', '_fdp.svg')}")

print("now run the following command to generate the image:")
print(f"sfdp -Tsvg {dot_file_path} -o {dot_file_path.replace('.dot', '_sfdp.svg')}")
os.system(f"sfdp -Tsvg {dot_file_path} -o {dot_file_path.replace('.dot', '_sfdp.svg')}")
