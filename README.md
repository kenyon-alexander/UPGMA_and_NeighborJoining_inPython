# UPGMA_and_NeighborJoining_inPython
My version of the UPGMA and Neighbor Joining phylogenetic algorithms, in Python. Detailed descriptions of these algorithms can be found on Wikipedia.


In order to run each of these algorithms, it is necessary to download the .py and .sh files for each, navigate to the directory of the .sh files, and run the following commands in the terminal:

> sh UPGMA.sh [[0,2,4,6,6,8],[2,0,4,6,6,8],[4,4,0,6,6,8],[6,6,6,0,4,8],[6,6,6,4,0,8],[8,8,8,8,8,0]]  None

or

> sh NeighborJoining.sh [[0,2,4,6,6,8],[2,0,4,6,6,8],[4,4,0,6,6,8],[6,6,6,0,4,8],[6,6,6,4,0,8],[8,8,8,8,8,0]]  None

The "None" at the end of each line is an optional argument for a list of species names. If no list is provided, each algorithm will provide letter names for the species. Furthermore, any distance matrix could be used for these algorithms--the matrices in the example commands above are just one possible matrix. 
