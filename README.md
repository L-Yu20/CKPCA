## Directory Structure

### `change/`
This directory contains the following components:

- **choose/**  
  Utilities for parameter selection in CKPCA.

- **generate/**  
  Code for generating simulated change-point data.

- **real_cp/**  
  Scripts for analyzing change points in real datasets.

- **simu_data/**  
  Storage for simulated datasets.

- **simulation_cp/**  
  Scripts for change-point detection experiments on simulated data.

- **pca_python/**  
  Python implementations of the KCP and Multirank change-point methods.

For the change-point simulation experiments, you should first generate synthetic data using the code in the `generate/` directory, and then run the corresponding change-point detection procedures.

---

### `cluster/`
This directory contains:

- **real_clus/**  
  Scripts for clustering analyses on real datasets.

- **simulation_clus/**  
  Scripts for clustering analyses using simulated data.

---

### `plot/`
This directory contains the code used to generate the figures in the paper.

---

### `dataset/`
This directory contains the real datasets used in the study.

Note: The MNIST dataset is not included due to its large size. You can download it from:  
http://yann.lecun.com/exdb/mnist/

---

### `randi/`
Code for computing the Rand Index.
