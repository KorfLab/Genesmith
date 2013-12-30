/*****************************************************************************\
 genesmith.cpp
 
 Copyright (C) 2013 Ravi Dandekar & Ian Korf. All rights reserved.
\*****************************************************************************/

#include <StochHMMlib.h>
#include <unistd.h>
#include <hmmer.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>

/* GLOBAL variables */
static std::string PROTEINS("");
static std::vector<std::string> CODONS;

double tb_eval (const std::string *genome, size_t pos, const std::string *mrna, size_t tb) {
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
	
	/* parse translation table */
	std::string codon_tbl;
	codon_tbl = genetic_code;
	
	std::ifstream file;
	file.open(codon_tbl.c_str());
	if (!file.is_open()) {
		std::cerr << "Error reading Codon Table\n";
	}
	
	std::string line;
	std::string proteins("");
	std::vector<std::string> codons;
	
	while (!file.eof()) {
		getline(file, line);
		std::stringstream ss(line);
		std::string elem;
		while(std::getline(ss, elem, '\t')) {
			int slen = elem.length();
			if (slen == 1) {proteins += elem;}
			else           {codons.push_back(elem);}
		}
	}
	PROTEINS = proteins;
	CODONS   = codons;
	
// 	std::cout << PROTEINS << std::endl;
// 	for(size_t i=0; i < CODONS.size(); i++) {
// 		std::cout << CODONS[i] << std::endl;
// 	}
	
	
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

