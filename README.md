# Usage Guide

## 1. Change-Point Detection

You will need the files in `change/`.

### 1.1 Simulations

To run the simulation studies, you first need to generate the data.

You can use the scripts in `generate/`. For example, to run Experiment 1, you can first run:

`generate_e1.R`

In this script:

- `times` represents the number of simulation repetitions;
- `px` represents the data dimension;
- `case` corresponds to the case settings in Experiment 1 of the paper;
- `type` can be set to `ba` or `im`, representing the balanced and imbalanced settings, respectively.

The generated data are stored in `big_array`, and the command

`npySave(file_name, big_array)`

can be used to save them in `.npy` format.

We recommend that you save these files in the folder:

`simu_data/`

The purpose of storing the data in `.npy` format is to make it convenient for you to load the same simulated datasets in both R and Python.

If the generated data are large, saving them in `.npy` format may require substantial memory. In that case, you may reduce the value of `times`. If you only use the R implementations, saving the `.npy` files is optional, and you may skip this step.

To generate the data for Experiments 2 and 3, you can run:

- `generate_e2.R`
- `generate_e3.R`

To generate the outlier setting in Experiment 1, you can run:

- `generate_e1out.R`

After generating the data, you can run the scripts in `simulation_cp/`. For example, for Experiment 1, you may run:

`simue1c1ba100.R`

You need to make sure that the basic settings in `simue1c1ba100.R`, such as `px`, `case`, `type`, and `times`, are consistent with those in `generate_e1.R`.

Similarly:

- for Experiment 2, you may run `simue2a4ba100.R`;
- for Experiment 3, you may run `simue3ba100.R`;
- for the outlier setting in Experiment 1, you may run `simue1c1out100.R`.

It is also worth mentioning that the codes for **KCP** and **Multirank** are located in `pca_python/`, and you need to run these methods in Python.

In addition, you may use the scripts in `choose/` to conduct the parameter sensitivity analysis reported in the paper.

---

### 1.2 Real Data

For the real-data experiments on change-point detection, you can use the scripts in `real_cp/`.

---

## 2. Clustering Analysis

You will need the files in `cluster/`.

### 2.1 Simulations

For the simulation studies on clustering analysis, you can use the scripts in `simulation_clus/`.

Each file corresponds to a clustering method. For example:

- `bull_dbscan.R` corresponds to the clustering method **DBSCAN**.

In these scripts, you can also adjust:

- the data dimension `px`,
- the number of simulation repetitions `times`,
- whether the setting is balanced or imbalanced.

---

### 2.2 Real Data

For the real-data experiments on clustering analysis, you can use the scripts in `real_clus/`.

Please note that the MNIST dataset is not included in this repository because of its large size. You can download it from:

`http://yann.lecun.com/exdb/mnist/`

In `real_clus/`, similarly, each file corresponds to a clustering method. For example:

- `mndbscan.R` corresponds to the clustering method **DBSCAN**.

---

## 3. Plotting

To generate the figures in the paper, you can use the scripts in `plot/`.

---

## Directory Structure

### `change/`

This directory contains the following components:

- **`choose/`**  
  Utilities for parameter selection in CKPCA.
- **`generate/`**  
  Scripts for generating simulated change-point data.
- **`real_cp/`**  
  Scripts for real-data change-point analysis.
- **`simu_data/`**  
  Storage for simulated datasets.
- **`simulation_cp/`**  
  Scripts for simulated change-point detection experiments.
- **`pca_python/`**  
  Python implementations of the KCP and Multirank change-point detection methods.

---

### `cluster/`

This directory contains:

- **`real_clus/`**  
  Scripts for real-data clustering analysis.
- **`simulation_clus/`**  
  Scripts for clustering analysis on simulated data.

---

### `plot/`

This directory contains the scripts used to generate the figures in the paper.

---

### `dataset/`

This directory contains the real datasets used in the study.

**Note:** The MNIST dataset is not included in this repository because of its large size. You can download it from:

`http://yann.lecun.com/exdb/mnist/`

---

### `randi/`

Code for computing the Rand Index.
