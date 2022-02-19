import numpy  as np
from scipy.sparse import diags

def normalize(mx):
    rowsum = np.array(mx.sum(1), dtype=np.float32)
    r_inv = (rowsum ** -1).flatten()
    r_inv[np.isinf(r_inv)] = 0.
    r_mat_inv = diags(r_inv)
    mx = r_mat_inv.dot(mx)
    return mx