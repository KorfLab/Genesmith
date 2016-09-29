Genesmith
=========

Gene annotation software that can incorporate protein alignment information and use an alternative decoding algorithm (The Stochastic Viterbi) to bolster genomic predictions 

### Software Dependencies
* StochHMM
* HMMER 3.0
* ik
* HMMstar
* FAlite.pm
* DataBrowser.pm

### Installation Notes

1) Within the Genesmith main directory create three soft links to the main directories of StochHMM, HMMER and ik
```
ln -s <PATH>/StochHMM/ ./StochHMM
ln -s <PATH>/hmmer-3.0/ ./hmmer
ln -s <PATH>/ik/ ./ik
```

2) Once all soft links are setup, compile the code using 'make' 
```
make
```

* **IMPORTANT:** In order to get HMMER to compile with our software you have to comment out line 138 of esl_msa.h because of a 'new' collision. There must be a better way...

To remove program binaries type 'make clean' in the Genesmith main directory 
```
make clean
```




