GeneSmith
=========

Gene finder mixing genomic and protein information


**Installation**

OMG, you have to comment out line 138 of esl_msa.h because of a 'new' collision. There must be a better way...


**Software Dependencies**
* StochHMM
* HMMER 3.0
* HMMstar
* FAlite.pm
* DataBrowser.pm

**Code**
<test_train_sets.pl>
    * Splits inputted data paired combinations of test and training sets for gene modeling an prediction
    * Default:  5 sets of test and training sets
<hmmgen.pl>
<genesmith.cpp>
<evaluator.pl>
<run_genesmith.pl>


