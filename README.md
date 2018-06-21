# grbeefcontam
computing the probability of contamination by E. coli in ground beef batches produced in a large production facility
 
 Copyright 2015 Allan Willms.

This C software source code is distributed under the GNU General Public License, Version 3.

Bug reports and comments should be sent to Allan Willms.
## Download
There is just one C file, instructions on use are in the comments.

## Description
GRBEEFCONTAM takes as input information about each different raw sources contributing to the ground beef production in a sequence of batches including fat percentage and susceptibility of contamination, mass contributed to each batch, average size of pieces from each carcass, number of pieces from each carcass, carcass spread information within the source, and carcass overlap probabilities across sources. Also, a list of hot (contaminated) batches is given. For each hot batch in this list in turn, the program computes the probability that each other batch is also contaminated. The output is a matrix of probabilities (in percentage), each row representing a particular batch, and each column corresponding to the listed hot batches.

The file example_syn provides the input file for the synthetic data example used in the manuscript below.

The algorithm is described in
<ul>
<li>    Petko M. Kitanov and Allan R. Willms, Probability of Escherichia coli Contamination Spread in Ground Beef Production, to appear in Mathematical Biosciences and Engineering, 2017. 
  </ul>
