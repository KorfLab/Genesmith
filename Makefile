##########################
# Makefile for GeneSmith #
##########################

APP1 = genesmith
SRC1 = genesmith.cpp

APP2 = nucleus
SRC2 = nucleus.cpp

STOCH = -IStochHMM/src -LStochHMM/src -lstochhmm
HMMER = -Ihmmer/src -Lhmmer/src -lhmmer
EASEL = -Ihmmer/easel -Lhmmer/easel -leasel
IK    = -Iik -Lik -lik
LIBS  = StochHMM/src/libstochhmm.a hmmer/src/libhmmer.a hmmer/easel/libeasel.a ik/libik.a

###########
# Targets #
###########

default:
	make $(APP1)
	make $(APP2)

$(APP1): $(SRC1) $(LIBS)
	g++ -w -O2 -o  $(APP1) $(SRC1) -lm $(STOCH) $(HMMER) $(EASEL) $(IK)

$(APP2): $(SRC2) $(LIBS)
	g++ -w -O2 -o  $(APP2) $(SRC2) -lm $(STOCH) $(HMMER) $(EASEL) $(IK)

clean:
	rm -f *.o $(APP1) $(APP2)

