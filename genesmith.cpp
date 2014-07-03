/*****************************************************************************\
 genesmith.cpp
 
 Copyright (C) 2013 Ravi Dandekar & Ian Korf. All rights reserved.
\*****************************************************************************/

extern "C" {
#include "p7_config.h"
#include "esl_config.h"
#include "easel.h"
#include "hmmer.h"
#include "ik/toolbox.h"
#include "ik/sequence.h"
#include "ik/align.h"

}
#include <StochHMMlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <getopt.h>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 

/* GLOBAL variables */
static ESL_ALPHABET    *ALPHABET = NULL; // Genomic Alphabet for StochHMM
static P7_HMM          *PROFILE  = NULL; // Profile for HMMER
static char            *KOG_SEQ  = NULL; // KOG AA sequence
static double           ASCALE   = 1.0;  // alignment scaling factor
static double           PSCALE   = 1.0;  // profile scaling factor

double tb_eval (const std::string *genome, size_t pos, const std::string *mrna, size_t tb) {
// 	std::string aa_seq = translate(mrna);
// 	char *seq          = (char*)aa_seq.c_str();
// 	double score       = hmmer_score(ALPHABET, PROFILE, seq);
// 	double loc_score   = sw_score(seq, KOG_SEQ);
	return 0.01;	
}

static float hmmer_score(const ESL_ALPHABET *ALPHABET, const P7_HMM *PROFILE, char *seq) {
// 	char *seq = (char*)aa_seq.c_str();   // Convert string into char array
		
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

static void usage () {
	std::cout
		<< "usage: genesmith [options] <hmm file> <seq file>\n"
		<< "options:\n"
		<< "  -p <file> profile HMM\n"
		<< "  -P <float> weight for profile HMM\n"
		<< "  -s <file> protein for smith-waterman alignment\n"
		<< "  -S <float> weight for smith-waterman\n" 
		<< "  -g <file> alternate genetic code\n"
		<< std::endl;
	exit(1);
}

/*========*/
/*  MAIN  */
/*========*/
int main (int argc, char ** argv) {
	extern char *optarg;
	extern int optind;
	int c;
	char  * fasta_file   = NULL;
	char  * profile_file = NULL;
	char  * protein_file = NULL;
	char  * phmm_coeff   = NULL;
	char  * sw_coeff     = NULL;
	char  * genetic_code = NULL;
		
	/* process command line */
	while ((c = getopt(argc, argv, "g:p:P:s:S:h:")) != -1) {
		switch (c) {
			case 'g': genetic_code = optarg; break;
			case 'p': profile_file = optarg; break;
			case 'P': phmm_coeff   = optarg; break;
			case 's': protein_file = optarg; break;
			case 'S': sw_coeff     = optarg; break;
			case 'h': usage();
			default:  usage();
		}
	}
	if (sw_coeff   != NULL) {ASCALE = atof(sw_coeff);}
	if (phmm_coeff != NULL) {PSCALE = atof(phmm_coeff);}
	
	
	if (argc - optind != 2) usage();
	std::string hmm_file = argv[optind];   // Gene Model for StochHMM
	std::string seq_file = argv[optind+1]; // KOG genomic DNA sequence
		
	if (protein_file != NULL) {
		std::ifstream aa_file;
		std::string protein_seq;
		std::string line;
		aa_file.open(protein_file);
		while (!aa_file.eof()) {
			aa_file >> line;
			if (line.substr(0,1) != ">") {
				protein_seq = line;
			}
		}
		aa_file.close();
		KOG_SEQ = (char*)protein_seq.c_str();
	}
	
	/* Setup HMMER Profile */
	if (profile_file != NULL) {
		P7_HMMFILE *hfp;
		if (p7_hmmfile_Open(profile_file, NULL, &hfp) != eslOK)
			p7_Fail("Failed to open HMM file %s", profile_file);
		if (p7_hmmfile_Read(hfp, &ALPHABET, &PROFILE) != eslOK)
			p7_Fail("Failed to read HMM");
		p7_hmmfile_Close(hfp);
	} 
	
	/* create alphabet */
	std::vector<std::string> dna;
	dna.push_back("A");
	dna.push_back("C");
	dna.push_back("G");
	dna.push_back("T");
	StochHMM::track tr = dna;	
		
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
			if (prev_st.substr(0,3) == "cds" and st.substr(0,3) == "don") {
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
			if (prev_st.substr(0,5) == "accep" and st.substr(0,3) == "cds") {
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
	if (PROFILE != NULL) {
		p7_hmm_Destroy(PROFILE);
	}
	esl_alphabet_Destroy(ALPHABET);
	return 0;
}



