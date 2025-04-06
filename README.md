# Dimension Sequences

In the paper [Dimension Sequences of Modular Forms](), we investigating when certain families of modular forms spaces can acheive all possible dimensions. This was inspired by a conjecture of Greg Martin in 2004 that `dim S2^new(N)` takes on all possible natural numbers.
This repository contains code to check certain cases of this type of problem.


To run the computations for the full space, one can execute the following. 
To run these for the new space, just replace `FS` with `NS` in the following commnds.
```
$ # Search the level-indexed sequences
$ sage compute.sage search-level FS   
$ # Search the weight-indexed sequences
$ sage compute.sage search-weight FS   
$ # Search the weight-indexed sequences for the sign pattern spaces
$ sage compute.sage search-sgnpatt-weight FS 
$ # Show the N for which the density of S2k^sigma(N) is >= 1
$ sage compute.sage show-densities FS  
$ # Show the sequences of S2k^sigma(N) which have density 1
$ sage compute.sage show-density1-seqs FS  
```
