import networkx as nx
from networkx.utils import np_random_state


# Copied from networkx.drawing.layout.py
def rescale_layout(pos, scale=1):
    """Returns scaled position array to (-scale, scale) in all axes.

    The function acts on NumPy arrays which hold position information.
    Each position is one row of the array. The dimension of the space
    equals the number of columns. Each coordinate in one column.

    To rescale, the mean (center) is subtracted from each axis separately.
    Then all values are scaled so that the largest magnitude value
    from all axes equals `scale` (thus, the aspect ratio is preserved).
    The resulting NumPy Array is returned (order of rows unchanged).

    Parameters
    ----------
    pos : numpy array
        positions to be scaled. Each row is a position.

    scale : number (default: 1)
        The size of the resulting extent in all directions.

    Returns
    -------
    pos : numpy array
        scaled positions. Each row is a position.

    See Also
    --------
    rescale_layout_dict
    """
    import numpy as np

    # Find max length over all dimensions
    pos -= pos.mean(axis=0)
    lim = np.abs(pos).max()  # max coordinate for all axes
    # rescale to (-scale, scale) in all directions, preserves aspect
    if lim > 0:
        pos *= scale / lim
    return pos
