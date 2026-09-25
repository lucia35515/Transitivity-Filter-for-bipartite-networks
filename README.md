# Transitivity Filter for Bipartite Networks

This repository contains the code associated with the transitivity filter introduced in:

**L. S. Ramirez et al., "Extracting the geometric backbone of bipartite networks", Nature Communications (2026).**

The filter identifies nodes whose bipartite clustering is significantly larger than expected in a degree-preserving randomized network. In this way, it separates structured connectivity from clustering that can be explained by the degree sequence alone.

## Which files should I use?

If you want to **apply the filter to a new bipartite network**, use the following two files, in this order:

1. `FILTER_PART1_clustering.cpp`
2. `FILTER_PARTII_bipartite_transitivity_filter.ipynb`

Files whose names begin with `NC_` are research files that were fully or partially used in the analyses associated with the paper. They are kept in the repository for reproducibility and reference, but **they are not required to apply the filter to a new network**.

---

## Filter workflow

The filter is applied to the nodes in the second layer of the bipartite network, i.e., the nodes listed in the second column of the input edge list. For a detailed description of the required bipartite edge-list format, see `README_input_format.md`.

The filtering procedure consists of two parts.

### Part I — Bipartite clustering and null-model generation

**File:** `FILTER_PART1_clustering.cpp`

This C++ program takes the original bipartite network as an edge list and performs two tasks.

First, it computes the bipartite clustering coefficient (C_3) of the nodes in the second layer of the bipartite network.

Second, it generates degree-preserving randomized versions of the network using bipartite edge swaps. The degree sequence of both layers is therefore preserved during the randomization. For each randomized realization, the program computes the same bipartite clustering coefficient.

Repeating this procedure over many randomized networks produces a null distribution of clustering values for nodes of a given degree.

This distribution is the reference against which the clustering observed in the original network is evaluated in Part II.

The program produces two main output files:

```text
<dataset>_clustering.dat
<dataset>_rewired_clustering.dat
```

The first contains the degree and clustering coefficient of the nodes in the original network. The second contains clustering coefficients obtained from the randomized networks.

Before running the program, the user should specify the dataset name and the desired number of randomized realizations in `FILTER_PART1_clustering.cpp`.



---

### Part II — Statistical filtering

**File:** `FILTER_PARTII_bipartite_transitivity_filter.ipynb`

The Jupyter notebook uses the files generated in Part I to determine whether the clustering of each node is larger than expected from the degree-preserving null model.

At the beginning of the notebook, the user specifies:

* the original bipartite edge-list file;
* the clustering file generated in Part I;
* the rewired-clustering file generated in Part I;
* the statistical thresholds to use for filtering;
* a degree to inspect when visualizing the null distribution.

For every degree (k), the notebook constructs the empirical null distribution from the randomized networks and calculates the clustering value corresponding to the selected percentile or significance threshold.

A node is identified as structurally significant when its clustering coefficient in the original network exceeds the corresponding degree-dependent threshold of the randomized ensemble.

The notebook also provides two useful visualizations:

* a histogram of the randomized clustering distribution for a user-selected degree, together with the observed clustering values and the selected statistical thresholds;
* the clustering profile as a function of degree, showing the clustering of the original network together with the degree-dependent threshold curves.

Finally, the notebook applies the selected threshold to the original bipartite network and generates the corresponding filtered network.

---

## In short

The complete workflow for a new network is

```text
bipartite edge list
        |
        v
FILTER_PART1_clustering.cpp
        |
        |-- original clustering
        |
        `-- clustering from degree-preserving randomized networks
                         |
                         v
FILTER_PARTII_bipartite_transitivity_filter.ipynb
                         |
                         |-- null distributions P(C3 | k)
                         |-- statistical thresholds
                         |-- diagnostic plots
                         `-- filtered bipartite network
```

Users interested only in applying the method can therefore ignore the `NC_` files and follow **Filter Part I → Filter Part II**.
