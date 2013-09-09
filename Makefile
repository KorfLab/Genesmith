##########################
# Makefile for GeneSmith #
##########################

APP1 = genesmith
SRC1 = genesmith.cpp

###########
# Targets #
###########

default:
	make $(APP1)

$(APP1): $(SRC1)
	g++ -o $(APP1) $(SRC1) -IStochHMM/src -LStochHMM/src -lm -lstochhmm

clean:
	rm -f *.o $(APP1)


