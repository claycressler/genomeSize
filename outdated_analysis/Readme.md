## README

This is the repository for Mueller et al. Metamorphosis imposes variable constraints on genome expansion.

To recreate the analysis:

1. Run get_consensus_from_jp_trees.R . The VertLife nexus file was obtained by downloading 1000 random trees from the VertLife website containing the 122 salamander species with genome size information. The full reference tree containing species for which there was molecular data from the Jetz and Pyron paper was obtained by downloading their data from Dryad (https://datadryad.org/stash/dataset/doi:10.5061/dryad.cc3n6j5). Of the 122, 118 had molecular data for tree building, so the tree distribution was filtered to these 118 to ensure a consistent topology, and branch lengths for the 1000 time trees were averaged.