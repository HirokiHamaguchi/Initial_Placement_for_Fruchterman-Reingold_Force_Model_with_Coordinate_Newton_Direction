import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import matplotlib

matplotlib.use("agg")

plt.rc("text", usetex=True)
plt.rc("font", family="serif")


def f(x, y):
    r = np.sqrt(x**2 + y**2)
    return (r**3) / 3 - np.log(r)


x = np.linspace(-1.5, 1.5, 300)
y = np.linspace(-1.5, 1.5, 300)
x, y = np.meshgrid(x, y)

z = f(x, y)

# Create the 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

# Plot the surface
ax.plot_surface(
    x,
    y,
    z,
    cmap="viridis",
    edgecolor="none",
    antialiased=False,
    rcount=300,
    ccount=300,
    vmin=0,
    vmax=3,
)


# Set plot labels and title
ax.set_xlabel("$x_{i,1}$", fontsize=20)
ax.set_ylabel("$x_{i,2}$", fontsize=20)
ax.set_zlabel("$E_{i,j}(\\| x_i - x_j \\|)$", fontsize=20)
ax.set_zlim(0, 4)
ax.set_zticks([1, 2, 3, 4])

# plt.show()
plt.tight_layout()
plt.subplots_adjust(left=-0.1, right=1.05, top=1.07, bottom=0.02)

# モアレが出るので、pdfではなくpngで保存
plt.savefig("doc/main/energy_3d/energy_3d.png", dpi=300)
