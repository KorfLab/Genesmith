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

#include <cstring>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <iterator>
#include <string>
#include <vector>

enum BP {
	A = 0,
	C = 1,
	G = 2,
	T = 3,
	N = 4,
};

/* GLOBAL variables */
static std::vector<std::string> PROTEINS;
static std::vector<std::string> CODONS;


std::string translate (const std::string *mrna, std::vector<std::string> &PROTEINS, std::vector<std::string> &CODONS);

double tb_eval (const std::string *genome, size_t pos, const std::string *mrna, size_t tb) {
	//std::string aa_seq = translate(mrna, PROTEINS, CODONS);
	//std::cout << "\n>MRNA\n" << *mrna << "\n>PROTEIN\n" << aa_seq << std::endl;
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
	
	/* parse translation table */
	std::string codon_tbl;
	codon_tbl = genetic_code;
	
	std::ifstream file;
	file.open(codon_tbl.c_str());
	if (!file.is_open()) {
		std::cerr << "Error reading Codon Table\n";
	}
	
	std::string line;
	static const char bps[] = "ACGTN";
	std::vector<std::string> proteins;
	std::vector<std::string> codons;
	
	while (!file.eof()) {
		getline(file, line);
		std::stringstream ss(line);
		std::string elem;
		std::string aa;
		while(std::getline(ss, elem, '\t')) {
			int slen = elem.length();
			if (slen == 1) {
				proteins.push_back(elem);
				aa = elem;
			} else {
				codons.push_back(elem);
				//std::cout << aa << "\t" << elem << std::endl;
			}
		}
	}
	PROTEINS = proteins;
	CODONS   = codons;
	
// 	for(size_t i=0; i < PROTEINS.size(); i++) {
// 		std::cout << PROTEINS[i] << "\t" << CODONS[i] << std::endl;
// 	}
	
	/* open output file */
// 	std::ofstream out;
// 	out.open("test_ann.gff");
// 	if(out.fail()) {
// 		std::cout << "Error writing into <test_ann.gff>\n" << std::endl;
// 		exit(1);
// 	}
	
	/* decode & output */
	StochHMM::StateFuncs my_trans;
	StochHMM::model hmm;
	StochHMM::seqTracks jobs;
// 	std::vector<StochHMM::gff_feature> feat;
	my_trans.assignTransitionFunction("HMMER", *tb_eval);
	hmm.import(hmm_file, &my_trans);
	jobs.loadSeqs(hmm, seq_file);
	StochHMM::seqJob *job=jobs.getJob();
	
	while (job != NULL) {
// 		std::stringstream results;
		StochHMM::trellis trell(&hmm, job->getSeqs());
		trell.viterbi();
		StochHMM::traceback_path tb(&hmm);
		trell.traceback(tb);
		tb.print_gff(job->getHeader());
		
		/* print to output file */
// 		tb.gff(feat, job->getSeqs()->getHeader());
// 		for (size_t i=0; i < feat.size(); i++) {
// 			std::string id    = feat[i].seqname;
// 			std::string struc = feat[i].feature;
// 			if (id.substr(0,1) == ">") {id.erase(0,1);}
// 				if (struc.substr(0,3) != "cds") {
// 					results << id              << "\t" 
// 							<< feat[i].source  << "\t" 
// 							<< feat[i].feature << "\t"
// 							<< feat[i].start   << "\t"
// 							<< feat[i].end     << "\t"
// 							<< feat[i].score   << "\t"
// 							<< feat[i].strand  << "\t"
// 							<< id              << "\n";
// 				}
// 		}
// 		out << results.str();
		job = jobs.getJob();
	}
	
// 	out.close();
	return 0;
}


/*=============*/
/* SUBROUTINES */
/*=============*/

/* Generates protein sequence using inputted DNA sequence and translation table */
std::string translate (const std::string *mrna, std::vector<std::string> &PROTEINS, std::vector<std::string> &CODONS){
	std::string aa_seq("");
	for (int i=0; i < mrna->length(); i += 3) {
		std::string codon = mrna->substr(i,3);
		for (size_t s=0; s < CODONS.size(); s++) {
			if (codon.compare(CODONS[s]) == 0) {
				std::string aa = PROTEINS[s];
				aa_seq += aa;
			}
		}
	}
	return aa_seq;
}
