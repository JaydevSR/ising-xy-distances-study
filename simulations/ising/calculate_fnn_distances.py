"""
Calculate fnn distances for the Ising model.
"""

import os
import time
import numpy as np
from sklearn.neighbors import NearestNeighbors

# read using np.loadtxt
datapath = os.path.join("d:\\", "Projects", "Dr. Heyl Group", "cluster_data")
print(f"Using data from: {datapath}")

lattice_sizes = [32, 40, 48, 56, 64]
temperatures = [1.800, 1.900, 2.000, 2.100, 2.200, 2.210, 2.220, 2.230,
                2.240, 2.250, 2.260, 2.270, 2.274, 2.278, 2.282, 2.286,
                2.290, 2.294, 2.298, 2.302, 2.306, 2.310, 2.314, 2.318,
                2.322, 2.330, 2.340, 2.350, 2.360, 2.370, 2.380, 2.390,
                2.400, 2.500, 2.600, 2.700]

model = NearestNeighbors(n_neighbors=2, algorithm="auto", n_jobs=-1)

for L in lattice_sizes:
    print(".==================================")
    print(f"| Lattice Size: {L} x {L}        ")
    print(".==================================")
    print("|  ")

    for stepT, T in enumerate(temperatures):
        print(f"| Process strarted for (T = {T}).")
        datafile = os.path.join(datapath, "Ising", f"Size{L}", "uncorr_configs",
                                f"ising_uncorr_configs_temp{T}_L{L}.txt")

        if not os.path.isfile(datafile):
            print(f"| File {datafile} not found.")
            continue

        start = time.perf_counter()
        configs_vecs = np.loadtxt(datafile, delimiter=",", dtype=np.int64).T
        model.fit(configs_vecs)
        dists, idxs = model.kneighbors(configs_vecs)
        fnn = dists[:, 1]
        el = time.perf_counter() - start

        szpath = os.path.join(datapath, "Ising", f"Size{L}", "fnn_dists")
        os.makedirs(szpath, exist_ok=True)
        np.savetxt(os.path.join(szpath, f"ising_fnn_dists_temp{T}_L{L}.txt"), fnn, delimiter=",")
        print(f"| Process completed for (T = {T}) in {el} seconds.")

    print("| Done.")
    print(".==================================")
