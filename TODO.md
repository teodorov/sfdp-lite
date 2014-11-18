sfdp-lite
=========

sfdp-lite is a simple graph layout tool that relies on the Graphviz sfdp algorithm to compute the layout of large graphs in Matrix Market format.

As opposed to Graphviz, sfdp-lite creates the Sparse matrix needed by the algorithm directly from a MatrixMarket representation (triplet representation) of the graph without passing through the Graphviz datastructures. Thus it effectivelly reduces the memory overhead.

The rendering is done using cairo.

It depends on a small patch to Graphviz (https://github.com/teodorov/graphviz) that exposes sfdp algorithm as a linkable library
