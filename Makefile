##########################
# Makefile for GeneSmith #
##########################

APP1 = genesmith
SRC1 = genesmith.cpp

STOCH = -IStochHMM/src -LStochHMM/src -lstochhmm
HMMER = -Ihmmer/src -Lhmmer/src -lhmmer

###########
# Targets #
###########

default:
	make $(APP1)

$(APP1): $(SRC1)
	g++ -o $(APP1) $(SRC1) -lm $(STOCH) $(HMMER) -Werror

clean:
	rm -f *.o $(APP1)


