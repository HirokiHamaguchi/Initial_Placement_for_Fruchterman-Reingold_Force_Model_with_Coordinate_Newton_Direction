# Initial Placement for Fruchterman–Reingold Force Model with Coordinate Newton Direction

This repository contains the implementation of the "Initial Placement for Fruchterman–Reingold Force Model with Coordinate Newton Direction" project.

## Requirements

We use Python and C++ to develop this project.

### Python

You can install the required packages by running the following command:

```bash
pip install -r requirements.txt
```

You can download the matrices used in this project by running the following command:

```bash
uv run src/python/get_matrix.py
uv run src/python/make_matrix.py
```

### C++

You can compile the C++ code by running the following command:

```bash
g++ src/cpp/main.cpp -o src/cpp/main.out -std=c++17 -O3
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
