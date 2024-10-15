# convert mtx file to dot file

import sys
import os


def mtx2dot(mtx_file, dot_file):
    with open(mtx_file) as f:
        lines = f.readlines()
        while lines[0].startswith("%"):
            lines.pop(0)
        n = int(lines[0].split()[0])
        m = int(lines[0].split()[1])
        with open(dot_file, "w") as f:
            f.write("graph {\n")
            for i in range(1, len(lines)):
                line = lines[i].split()
                f.write("    %s -- %s;\n" % (line[0], line[1]))
            f.write("}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python mtx2dot.py mtx_file dot_file")
        sys.exit(1)
    mtx_file = sys.argv[1]
    dot_file = sys.argv[2]
    mtx2dot(mtx_file, dot_file)
    print("Convert %s to %s" % (mtx_file, dot_file))
