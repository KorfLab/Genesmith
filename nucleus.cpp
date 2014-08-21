/*****************************************************************************\
 nucleus.cpp
 
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
static char            *PROTEIN  = NULL; // KOG AA sequence

/* FUNCTIONS */
char* translate (const std::string *mrna) {
	char *seq = (char*)mrna->c_str();
	return ik_translate(seq, 1);
}

static float hmmer_score(const ESL_ALPHABET *ALPHABET, const P7_HMM *PROFILE, const char *aa_seq) {		
	P7_BG        *bg     = NULL; /* background */
	P7_PROFILE   *gm     = NULL; /* generic model */
	P7_GMX       *mx     = NULL; /* viterbi matrix */
	ESL_DSQ      *dsq    = NULL; /* digital sequence */
	int           L      = 0;    /* length of sequence */
	float         score;
	
	/* digitize sequence */
	L = strlen(aa_seq);
	esl_abc_CreateDsq(ALPHABET, aa_seq, &dsq);
	
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
		<< "usage: nucleus [options] <hmm file> <seq file>\n"
		<< "options:\n"
		<< "  -t <int>   transcripts\n"
		<< "  -p <file>  profile HMM\n"
		<< "  -P <float> threshold for profile HMM\n"
		<< "  -s <file>  protein for alignment\n"
		<< "  -S <float> threshold for alignment\n" 
		<< "  -g <file>  alternate genetic code\n"
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
	char  * hmmer_file   = NULL;
	char  * align_file   = NULL;
	char  * hmmer_cutoff = NULL;
	char  * align_cutoff = NULL;
	char  * genetic_code = NULL;
	char  * transcripts  = NULL;
	double t_hmmer, t_align;
	int    tracebacks = 1000;
	
	/* process command line */
	while ((c = getopt(argc, argv, "g:p:P:s:S:h:")) != -1) {
		switch (c) {
			case 't': transcripts  = optarg; break;
			case 'g': genetic_code = optarg; break;
			case 'p': hmmer_file   = optarg; break;
			case 'P': hmmer_cutoff = optarg; break;
			case 's': align_file   = optarg; break;
			case 'S': align_cutoff = optarg; break;
			case 'h': usage();
			default:  usage();
		}
	}
	if (transcripts != NULL) tracebacks = atoi(transcripts);
	if (hmmer_file != NULL)  t_hmmer = atof(hmmer_file);
	if (align_file != NULL)  t_align = atof(align_file);
	
	
	if (argc - optind != 2) usage();
	std::string hmm_file = argv[optind];
	std::string seq_file = argv[optind+1];
	
	/* protein file? */
	if (align_file != NULL) {
		std::ifstream aa_file;
		std::string protein_seq;
		std::string line;
		aa_file.open(align_file);
		while (!aa_file.eof()) {
			aa_file >> line;
			if (line.substr(0,1) != ">") {
				protein_seq = line;
			}
		}
		aa_file.close();
		PROTEIN = (char*)protein_seq.c_str();
	}
	
	/* HMMER profile? */
	if (hmmer_file != NULL) {
		P7_HMMFILE *hfp;
		if (p7_hmmfile_Open(hmmer_file, NULL, &hfp) != eslOK)
			p7_Fail("Failed to open HMM file %s", hmmer_file);
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
	hmm.import(hmm_file);
	jobs.loadSeqs(hmm, seq_file);
	StochHMM::seqJob *job = jobs.getJob();
	std::vector<StochHMM::gff_feature> feat;

	/* main loop */
	
	while (job != NULL) {
		StochHMM::trellis trell(&hmm, job->getSeqs()); // maybe use 'simple' versions?
		trell.stochastic_viterbi();
		
		for (int i = 0; tracebacks; i++) {
			StochHMM::traceback_path path(&hmm);
			trell.stochastic_traceback(path);
			/* now do something with the paths, like catalog them
			
				see path documentation 
				
				may want to hash the trace-back string or some facsimile
			
			*/
		}
		

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



