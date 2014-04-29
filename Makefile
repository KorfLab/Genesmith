##########################
# Makefile for GeneSmith #
##########################

APP1 = genesmith
SRC1 = genesmith.cpp

STOCH = -IStochHMM/src -LStochHMM/src -lstochhmm
HMMER = -Ihmmer/src    -Lhmmer/src    -lhmmer
EASEL = -Ihmmer/easel  -Lhmmer/easel  -leasel

###########
# Targets #
###########

default:
	make $(APP1)

$(APP1): $(SRC1)
	g++ -w -o  $(APP1) $(SRC1) -lm $(STOCH) $(HMMER) $(EASEL)

clean:
	rm -f *.o $(APP1)

