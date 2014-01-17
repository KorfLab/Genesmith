/*****************************************************************************\
 genesmith.cpp
 
 Copyright (C) 2013 Ravi Dandekar & Ian Korf. All rights reserved.
\*****************************************************************************/

#include <StochHMMlib.h>
#include <hmmer.h>
#include <iostream>
#include <map>


/* GLOBAL variables */
std::map<std::string, int> BP_KEY;        // Dictionary of Base Pair Positions within Translation Table

// Standard Translation Table Hard Coded
static char Codons[5][5][5] = {{{'K','N','K','N', 'X'},{'T','T','T','T', 'X'},{'R','S','R','S', 'X'},{'I','I','M','I', 'X'}, {'A','C','G','T', 'X'}},
						       {{'Q','H','Q','H', 'X'},{'P','P','P','P', 'X'},{'R','R','R','R', 'X'},{'L','L','L','L', 'X'}, {'A','C','G','T', 'X'}},
						       {{'E','D','E','D', 'X'},{'A','A','A','A', 'X'},{'G','G','G','G', 'X'},{'V','V','V','V', 'X'}, {'A','C','G','T', 'X'}},
						       {{'*','Y','*','Y', 'X'},{'S','S','S','S', 'X'},{'*','C','W','C', 'X'},{'L','F','L','F', 'X'}, {'A','C','G','T', 'X'}},
						       {{'A','C','G','T', 'X'},{'A','C','G','T', 'X'},{'A','C','G','T', 'X'},{'A','C','G','T', 'X'}, {'A','C','G','T', 'X'}}};

/* FUNCTIONS */
std::string translate (const std::string *mrna); 
double tb_eval (const std::string *genome, size_t pos, const std::string *mrna, size_t tb) {
// 	std::string aa_seq = translate(mrna);
// 	std::cout << "\n>MRNA\n" << *mrna << "\n>PROTEIN\n" << aa_seq << std::endl;
	return 0.01;
}

static void usage () {
	std::cout
		<< "usage: genesmith [options] <hmm file> <seq file>\n"
		<< "options:\n"
		<< "  -p <file> profile HMM (requires -g)\n"
		<< "  -f <file> protein fasta file (requires -g)\n"
		<< "  -g <file> alternate genetic code\n"
		<< std::endl;
	exit(1);
}

/*========*/
/*  MAIN  */
/*========*/
int main (int argc, char ** argv) {
	char * fasta_file   = NULL;
	char * profile_file = NULL;
	const char* genetic_code = "codon_tbl.txt";
	int c;
	extern char *optarg;
	extern int optind;

	
	/* process command line */
	while ((c = getopt(argc, argv, "g:h:p:")) != -1) {
		switch (c) {
			case 'g': genetic_code = optarg; break;
			case 'p': profile_file = optarg; break;
			case 'f': fasta_file = optarg;   break;
			case 'h': usage();
			default:  usage();
		}
	}
	
	if (argc - optind != 2) usage();
	std::string hmm_file  =	  argv[1];
	std::string seq_file  =	  argv[2];
	
	/* create alphabet */
	std::vector<std::string> dna;
	dna.push_back("A");
	dna.push_back("C");
	dna.push_back("G");
	dna.push_back("T");
	StochHMM::track tr = dna;
	
	/* Insert Position Values for BP_KEY */
	BP_KEY.insert(std::pair<std::string, int>("A", 0));
	BP_KEY.insert(std::pair<std::string, int>("C", 1));
	BP_KEY.insert(std::pair<std::string, int>("G", 2));
	BP_KEY.insert(std::pair<std::string, int>("T", 3));
	
	
	/* decode & output */
	StochHMM::StateFuncs my_trans;
	StochHMM::model hmm;
	StochHMM::seqTracks jobs;
	my_trans.assignTransitionFunction("HMMER", *tb_eval);
	hmm.import(hmm_file, &my_trans);
	jobs.loadSeqs(hmm, seq_file);
	StochHMM::seqJob *job=jobs.getJob();
	
	while (job != NULL) {
		StochHMM::trellis trell(&hmm, job->getSeqs());
		trell.viterbi();
		StochHMM::traceback_path tb(&hmm);
		trell.traceback(tb);
		tb.print_gff(job->getHeader());
		job = jobs.getJob();
	}
	return 0;
}


/*=============*/
/* SUBROUTINES */
/*=============*/

/* Generates protein sequence using inputted DNA sequence and translation table */
std::string translate (const std::string *mrna){
	std::string aa_seq("");
	for (int i=0; i < mrna->length(); i += 3) {
		int pos_one;
		int pos_two;
		int pos_three;
		std::string codon    = mrna->substr(i,3);
		for (size_t s=0; s < codon.length(); s++) {
			std::string bp = codon.substr(s,1);
			std::map<std::string, int>::iterator bp_to_pos = BP_KEY.find(bp);
			if (bp_to_pos != BP_KEY.end()) {
				if      (s == 0) {pos_one   = (*bp_to_pos).second;}
				else if (s == 1) {pos_two   = (*bp_to_pos).second;}
				else if (s == 2) {pos_three = (*bp_to_pos).second;}
			}
		}
		aa_seq.append(1, Codons[pos_one][pos_two][pos_three]);
	}
	return aa_seq;
}