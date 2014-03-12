/*****************************************************************************\
 genesmith.cpp
 
 Copyright (C) 2013 Ravi Dandekar & Ian Korf. All rights reserved.
\*****************************************************************************/

extern "C" {
#include "p7_config.h"
#include "esl_config.h"
#include "easel.h"
#include "hmmer.h"

#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"
}
#include <StochHMMlib.h>
#include <iostream>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 

/* GLOBAL variables */
std::map<std::string, int> BP_KEY;                                  // Dictionary of Base Pair Positions within Translation Table
static ESL_ALPHABET    *ALPHABET  = NULL;                           // Genomic Alphabet for StochHMM
static P7_HMM          *PROFILE   = NULL;                           // Profile for HMMER
static ESL_ALPHABET    *AA_ALPH   = esl_alphabet_Create(eslAMINO);  // AA Alphabet for S-W Alignment
static char            *KOG_SEQ   = NULL;                           // KOG AA sequence

// Standard Translation Table Hard Coded
static char CODONS[5][5][5] = {{{'K','N','K','N', 'X'},{'T','T','T','T', 'X'},{'R','S','R','S', 'X'},{'I','I','M','I', 'X'}, {'A','C','G','T', 'X'}},
                               {{'Q','H','Q','H', 'X'},{'P','P','P','P', 'X'},{'R','R','R','R', 'X'},{'L','L','L','L', 'X'}, {'A','C','G','T', 'X'}},
                               {{'E','D','E','D', 'X'},{'A','A','A','A', 'X'},{'G','G','G','G', 'X'},{'V','V','V','V', 'X'}, {'A','C','G','T', 'X'}},
                               {{'*','Y','*','Y', 'X'},{'S','S','S','S', 'X'},{'*','C','W','C', 'X'},{'L','F','L','F', 'X'}, {'A','C','G','T', 'X'}},
                               {{'A','C','G','T', 'X'},{'A','C','G','T', 'X'},{'A','C','G','T', 'X'},{'A','C','G','T', 'X'}, {'A','C','G','T', 'X'}}};

/* FUNCTIONS */
std::string translate (const std::string *mrna);
float hmmer_score(const ESL_ALPHABET *ALPHABET, const P7_HMM *PROFILE, std::string aa_seq);  
int pval_bitscore(const ESL_ALPHABET *ALPHABET, const P7_HMM *PROFILE, std::string aa_seq); //read <p7_domaindef.c> code in HMMER for p-val/bitscore object 

double tb_eval (const std::string *genome, size_t pos, const std::string *mrna, size_t tb) {
// 	std::string aa_seq = translate(mrna);
// 	size_t prot_len    = aa_seq.length();
// 	size_t kog_len     = strlen(KOG_SEQ);
// 	double score;
// 	
// 	if (prot_len >= (kog_len*0.60)) {
// 		score = hmmer_score(ALPHABET, PROFILE, aa_seq);
// 		if (score < 0) {
// 			score = 0.01;
// 		} else {
// 			score = 1.0;
// 		}
// 	} else if (prot_len <= (kog_len*0.10)) {
// 		score = 0.00001;
// 	} else {
// 		score = 0.001;
// 	}
// 	double score       = pval_bitscore(ALPHABET, PROFILE, aa_seq);
// 	std::cout << "\n>MRNA:\n"        << *mrna  << std::endl;
// 	std::cout << "\n>PROTEIN:\n"     << aa_seq << std::endl;
// 	std::cout << "\n>HMMER Score:\t" << score  << std::endl;
// 	std::cout << ">Traceback:\t"     << tb     << std::endl;
// 	std::cout << ">Position:\t"      << pos    << std::endl;
	return 0.01;	
// 	return score;
}

static void usage () {
	std::cout
		<< "usage: genesmith [options] <hmm file> <seq file> <profile file> <kog aa seq>\n"
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
// int main(int argc, char* const argv[]){
int main (int argc, char ** argv) {
	char * fasta_file   = NULL;
	char * profile_fh   = NULL;
	const char* genetic_code = "codon_tbl.txt";
	int c;
	extern char *optarg;
	extern int optind;

	
	/* process command line */
	while ((c = getopt(argc, argv, "g:h:p:")) != -1) {
		switch (c) {
			case 'g': genetic_code = optarg; break;
			case 'p': profile_fh   = optarg; break;
			case 'f': fasta_file   = optarg; break;
			case 'h': usage();
			default:  usage();
		}
	}
	
	if (argc - optind != 4) usage();
	std::string hmm_file     =   argv[1];  // Gene Model for StochHMM
	std::string seq_file     =   argv[2];  // KOG genomic DNA sequence
	char       *profile_file =   argv[3];  // Convert this to an cmd line option once optimized
	std::string protein_seq  =   argv[4];  // KOG protein sequence
	KOG_SEQ                  = (char*)protein_seq.c_str();
	
	/* Setup HMMER Profile */
	P7_HMMFILE *hfp;
	if (p7_hmmfile_Open(profile_file, NULL, &hfp) != eslOK)
	    p7_Fail("Failed to open HMM file %s", profile_file);
	if (p7_hmmfile_Read(hfp, &ALPHABET, &PROFILE) != eslOK)
	    p7_Fail("Failed to read HMM");
	p7_hmmfile_Close(hfp); 
	
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
	std::vector<StochHMM::gff_feature> feat;
	
	std::string prev_st = "none";
	std::string start_id;
	std::string end_id;
	size_t cds_start(0);
	size_t cds_end(0);
	
	while (job != NULL) {
		StochHMM::trellis trell(&hmm, job->getSeqs());
		trell.viterbi();
		StochHMM::traceback_path tb(&hmm);
		trell.traceback(tb);
		tb.gff(feat, job->getSeqs()->getHeader());
		for (size_t i=0; i < feat.size(); i++) {
			std::string st = feat[i].feature;
			std::string id  = feat[i].seqname;
			if (id.substr(0,1) == ">") {id.erase(0,1);}
			if (prev_st == "GU0" and st == "start0") {
				start_id  = id;
				cds_start = feat[i].start;
			}
			if (prev_st.substr(0,3) == "cds" and st.substr(0,1) == "D") {
				cds_end = feat[i].end;
				std::cout   << id              << "\t" 
						    << feat[i].source  << "\t" 
						    << "CDS"           << "\t"
						    << cds_start       << "\t"
						    << cds_end         << "\t"
						    << feat[i].score   << "\t"
						    << feat[i].strand  << "\t"
						    << id              << "\n";
			}
			if (prev_st.substr(0,1) == "A" and st.substr(0,3) == "cds") {
				cds_start = feat[i].start;
			}
			if (prev_st == "stop1" and st == "stop2") {
				end_id  = id;
				cds_end = feat[i].end;
				if (start_id == end_id) {
					std::cout << id              << "\t" 
						      << feat[i].source  << "\t" 
						      << "CDS"           << "\t"
						      << cds_start       << "\t"
						      << cds_end         << "\t"
						      << feat[i].score   << "\t"
						      << feat[i].strand  << "\t"
						      << id              << "\n";
				}
			}
			prev_st = st;
		}
		//tb.print_gff(job->getHeader());
		job = jobs.getJob();
		feat.clear();
	}
	/* Global Cleanup */
	p7_hmm_Destroy(PROFILE);
	esl_alphabet_Destroy(ALPHABET);
	esl_alphabet_Destroy(AA_ALPH);
	return 0;
}


/*=============*/
/* SUBROUTINES */
/*=============*/

/* Generates protein sequence using inputted DNA sequence and Global char array of AAs */
std::string translate (const std::string *mrna) {
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
		aa_seq.append(1, CODONS[pos_one][pos_two][pos_three]);
	}
	return aa_seq;
}

/* Calculates HMMER score */
float hmmer_score(const ESL_ALPHABET *ALPHABET, const P7_HMM *PROFILE, std::string aa_seq) {
	char *seq = (char*)aa_seq.c_str();   // Convert string into char array
		
	P7_BG        *bg     = NULL; /* background */
	P7_PROFILE   *gm     = NULL; /* generic model */
	P7_GMX       *mx     = NULL; /* viterbi matrix */
	ESL_DSQ      *dsq    = NULL; /* digital sequence */
	int           L      = 0;    /* length of sequence */
	float         score;
	
	/* digitize sequence */
	L = strlen(seq);
	esl_abc_CreateDsq(ALPHABET, seq, &dsq);
	
	/* background */
	bg = p7_bg_Create(ALPHABET);
	p7_bg_SetLength(bg, L);
	
	/* profile */
	gm = p7_profile_Create(PROFILE->M, ALPHABET);
	p7_ProfileConfig(PROFILE, bg, gm, L, p7_GLOCAL);
	
	/* viterbi */
	mx = p7_gmx_Create(gm->M, L);
    try{
        p7_GViterbi(dsq, L, gm, mx, &score);
    }
    catch(...){
        std::cerr << "Error\n";
    }
    
    /* local clean up */
	p7_gmx_Destroy(mx);
	p7_profile_Destroy(gm);
	p7_bg_Destroy(bg);
	free(dsq);
	
	return score;
}

/* Calculates P-value and Bitscore */
int pval_bitscore(const ESL_ALPHABET *ALPHABET, const P7_HMM *PROFILE, std::string aa_seq) {
// 	char *seq = (char*)aa_seq.c_str();   // Convert string into char array
// 	
// 	P7_DOMAINDEF   *ddef    = NULL;
// 	ESL_RANDOMNESS *r       = esl_randomness_CreateFast(42);
// 	ESL_SQ         *sq      = NULL;
// 	P7_TRACE       *tr      = NULL;
// 	P7_GMX         *fwd     = NULL;
// 	P7_GMX         *bck     = NULL;
// 	float           overall_sc;
// 	float           sc;
// 	int             d;           /* iterator */
// 	int             tot_true;
// 	int             tot_found;
// 		
// 	P7_BG        *bg     = NULL; /* background */
// 	P7_PROFILE   *gm     = NULL; /* generic model */
// 	ESL_DSQ      *dsq    = NULL; /* digital sequence */
// 	int           L      = 0;    /* length of sequence */
// 	
// 	/* digitize sequence */
// 	L = strlen(seq);
// 	esl_abc_CreateDsq(ALPHABET, seq, &dsq);
// 	
// 	/* background */
// 	bg = p7_bg_Create(ALPHABET);
// 	p7_bg_SetLength(bg, L);
// 	
// 	/* profile */
// 	gm = p7_profile_Create(PROFILE->M, ALPHABET);
// 	p7_ProfileConfig(PROFILE, bg, gm, L, p7_LOCAL);
// 	
// 	/* viterbi */
// 	sq = esl_sq_CreateDigital(ALPHABET);
// 	ddef = p7_domaindef_Create(r);
// 	
// 	fwd = p7_gmx_Create(gm->M, L);
// 	bck = p7_gmx_Create(gm->M, L);
// 	tr  = p7_trace_Create();
// 	p7_FLogsumInit();
// 	
// 	/* First configure to actual protein length */
// 	tot_true = tot_found = 0;
// 	p7_ProfileEmit(r, PROFILE, bg, sq, tr);
// 	p7_trace_Index(tr);
// 	
// 	/* Reconfigure to length of viterbi seq */
// 	
//     try{
//         p7_GViterbi(dsq, L, gm, fwd, &overall_sc);
//         p7_domaindef_ByViterbi(gm, seq, fwd, bck, ddef);
//         
// 		for (d=0; d < ddef->ndom; d++) {
// 			std::cout << "P-value: "    << ddef->dcl[d].pvalue 
// 					  << "\tBitscore: " << ddef->dcl[d].bitscore << std::endl;
// 		} 
//     }
//     catch(...){
//         std::cerr << "Error\n";
//     }
//     
//     
//     /* local clean up */
// 	p7_gmx_Destroy(fwd);
// 	p7_gmx_Destroy(bck);
// 	p7_profile_Destroy(gm);
// 	p7_bg_Destroy(bg);
// 	esl_randomness_Destroy(r);
// 	free(dsq);
// 	
// 	return overall_sc;
// 
}
