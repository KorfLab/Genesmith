##########################
# Makefile for GeneSmith #
##########################

APP1 = genesmith
SRC1 = genesmith.cpp

STOCH = -IStochHMM/src -LStochHMM/src -lstochhmm
#HMMER = -Ihmmer/src -Lhmmer/src -Ihmmer/easel -Lhmmer/easel -lhmmer
#HMMER = -Ihmmer/src -Ihmmer/src/impl_sse -Ihmmer/easel hmmer/src/libhmmer.a hmmer/src/impl_sse/libhmmerimpl.a hmmer/easel/libeasel.a
HMMER = -Ihmmer-3.1b1/src -Ihmmer-3.1b1/src/impl_sse -Ihmmer-3.1b1/easel hmmer-3.1b1/src/libhmmer.a hmmer-3.1b1/src/impl_sse/libhmmer-impl.stamp hmmer-3.1b1/easel/libeasel.a

###########
# Targets #
###########

default:
	make $(APP1)

$(APP1): $(SRC1)
	g++ -w -o  $(APP1) $(SRC1) -lm $(STOCH) $(HMMER)

clean:
	rm -f *.o $(APP1)

