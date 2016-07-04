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
#include <iomanip>

/* GLOBAL variables */
static ESL_ALPHABET    *ALPHABET = NULL; // Genomic Alphabet for StochHMM
static P7_HMM          *PROFILE  = NULL; // Profile for HMMER
static char            *KOG_SEQ  = NULL; // KOG AA sequence
static double           ASCALE   = 1.0;  // alignment scaling factor
static double           PSCALE   = 1.0;  // profile scaling factor
//static long           FCOUNT   = 0;    // counts number of iterations through 'tb_eval' function

/* FUNCTIONS */
char* translate (const std::string *mrna) {
	char *seq = (char*)mrna->c_str();
	return ik_translate(seq, 1);
}

static float hmmer_score(const ESL_ALPHABET *ALPHABET, const P7_HMM *PROFILE, const char *aa_seq, size_t pos) {
// 	char *seq = (char*)aa_seq.c_str();   // Convert string into char array
		
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
	p7_ProfileConfig(PROFILE, bg, gm, L, p7_LOCAL);  // p7_GLOCAL used previously and had errors
	
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
//	fprintf(stderr, "freeing..");
	free(dsq);
//	fprintf(stderr, "freed\n");
	
	return score;
}

double tb_eval (const std::string *genome, size_t pos, const std::string *mrna, size_t tb) {
	ik_align align;
	char *sub_kogseq;
	
	if (PROFILE != NULL or KOG_SEQ != NULL) {
		char *aa_seq = translate(mrna);
		double score = 0;
		
		if (PROFILE != NULL) {
			score += PSCALE * hmmer_score(ALPHABET, PROFILE, aa_seq, pos);
		}
		/*
		if (KOG_SEQ != NULL) {
// 			printf("KOG: %s\n", KOG_SEQ);
// 			printf("AA: %s\n", aa_seq);
			int seqlen = strlen(aa_seq);
			int koglen = strlen(KOG_SEQ);
			sub_kogseq = (char*)std::malloc(strlen(aa_seq) + 1);
			std::strncpy(sub_kogseq, KOG_SEQ, seqlen);
			sub_kogseq[seqlen] = '\0';
// 			printf ("SUB: %s\n\n", sub_kogseq);			
			if (seqlen + 1 <= koglen) {			
				align = sw_mat(aa_seq, (char*)sub_kogseq, 62);
				score -= 0.01 * align->score;
				//printf("%f %s\n", score, align->alignment);
			} else {
				align = sw_mat(aa_seq, KOG_SEQ, 62);
				score -= 0.01 * align->score;
				//printf("%f %s\n", score, align->alignment);
			}
			//Original Method
// 			align = sw_mat(aa_seq, KOG_SEQ, 62);
// 			score -= 0.01 * align->score;
			
			ik_align_free(align);
			
			free(sub_kogseq);
		}
		*/
		
		free(aa_seq);
		return score;
	} else {
		return 0;
	}
}


static void usage () {
	std::cout
		<< "usage: genesmith [options] <hmm file> <seq file>\n"
		<< "options:\n"
		<< "  -p <file> profile HMM\n"
		<< "  -P <float> weight for profile HMM\n"
		<< "  -s <file> protein for smith-waterman alignment\n"
		<< "  -S <float> weight for smith-waterman\n" 
		<< "  -V <size_t> number of repetitions used for Stochastic Viterbi, silence Viterbi"
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
    char  * repetitions  = NULL;
	size_t reps(1000);   // for Stochastic Viterbi
	ik_pipe pipe;
	ik_fasta fasta;
    
	/* process command line */
	while ((c = getopt(argc, argv, "g:p:P:s:S:V:h")) != -1) {
		switch (c) {
			case 'g': genetic_code = optarg; break;
			case 'p': profile_file = optarg; break;
			case 'P': phmm_coeff   = optarg; break;
			case 's': protein_file = optarg; break;
			case 'S': sw_coeff     = optarg; break;
            case 'V': repetitions  = optarg; break;
			case 'h': usage();
			default:  usage();
		}
	}
	if (sw_coeff    != NULL) {ASCALE = atof(sw_coeff);}
	if (phmm_coeff  != NULL) {PSCALE = atof(phmm_coeff);}
    if (repetitions != NULL) {reps  = atoi(repetitions);}
	
	if (argc - optind != 2) usage();
	std::string hmm_file = argv[optind];   // Gene Model for StochHMM
	std::string seq_file = argv[optind+1]; // KOG genomic DNA sequence
	
	
		
	if (protein_file != NULL) {
		pipe = ik_pipe_open(protein_file, "r");
		fasta = ik_fasta_read(pipe->stream);
		KOG_SEQ = fasta->seq;
	
	/*
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
	*/
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
	my_trans.assignTransitionFunction("TRANSLATE", *tb_eval);
	
	hmm.import(hmm_file, &my_trans);
	jobs.loadSeqs(hmm, seq_file);
	StochHMM::seqJob *job=jobs.getJob();
	std::vector<StochHMM::gff_feature> feat;  // Vector for processing GFF output (Viterbi only)
	
	
	std::string prev_st = "none";
	std::string start_id;
	std::string end_id;
	size_t cds_start(0);
	size_t cds_end(0);
    double perc_reps(0);  // for Stochastic Viterbi
    
	while (job != NULL) {
		StochHMM::trellis trell(&hmm, job->getSeqs());
		
		if (repetitions != NULL) {
			// Stochastic Viterbi ---------------------------------------------
			trell.simple_stochastic_viterbi();
			std::map<StochHMM::traceback_path, int > paths;
			for(size_t i=0;i<reps;i++){
				StochHMM::traceback_path pth(&hmm);
				trell.stochastic_traceback(pth);
				paths[pth]++;
			}
		
			size_t path_count(0);   // Helps keeps track of which gene isoforms you are on
			std::map<StochHMM::traceback_path,int>::iterator it;
			for( it =paths.begin(); it!=paths.end();it++){
				int count = (*it).second;
				StochHMM::traceback_path tb = (*it).first;
				//(*it).first.print_gff(job->getSeqs()->getHeader());
			
				std::vector<StochHMM::gff_feature> svfeat;
				tb.gff(svfeat, job->getSeqs()->getHeader());
				for (size_t i=0; i < svfeat.size(); i++) {
					std::string st = svfeat[i].feature;
					std::string id = svfeat[i].seqname;
					if (id.substr(0,1) == ">") {id.erase(0,1);}
					perc_reps = ((double)count / reps) * 100;
				
					if (prev_st == "GU0" and st == "start0") {
						start_id  = id;
						cds_start = svfeat[i].start;
					}
					if (prev_st.substr(0,3) == "cds" and st.substr(0,3) == "don") {
						cds_end = svfeat[i].end -1;
						std::cout << id             << "\t"
								<< svfeat[i].source << "\t"
								<< "CDS"            << "\t"
								<< cds_start        << "\t"
								<< cds_end          << "\t"
								<< svfeat[i].score  << "\t"
								<< svfeat[i].strand << "\t"
								<< id               << "." << path_count << ":percentage_reps=" <<  std::fixed << std::setprecision(3) << perc_reps << "%\n";
					}
					if (prev_st.substr(0,5) == "accep" and st.substr(0,3) == "cds") {
						cds_start = svfeat[i].start;
					}
					if (prev_st == "stop1" and st == "stop2") {
						end_id  = id;
						cds_end = svfeat[i].end;
						if (start_id == end_id) {
							std::cout << id     << "\t"
							<< svfeat[i].source << "\t"
							<< "CDS"            << "\t"
							<< cds_start        << "\t"
							<< cds_end          << "\t"
							<< svfeat[i].score  << "\t"
							<< svfeat[i].strand << "\t"
							<< id               << "." << path_count << ":percentage_reps="  << std::fixed << std::setprecision(3) << perc_reps << "%\n";
						}
					}
					// Basic Model 1 CDS state
					if (prev_st.substr(0,3) == "CDS" and st.substr(0,3) == "don") {
						cds_end = svfeat[i].end -1;
						std::cout << id     << "\t"
						<< svfeat[i].source << "\t"
						<< "CDS"            << "\t"
						<< cds_start        << "\t"
						<< cds_end          << "\t"
						<< svfeat[i].score  << "\t"
						<< svfeat[i].strand << "\t"
						<< id               << "." << path_count << ":percentage_reps=" <<  std::fixed << std::setprecision(3) << perc_reps << "%\n";
					}
					if (prev_st.substr(0,5) == "accep" and st.substr(0,3) == "CDS") {
						cds_start = svfeat[i].start;
					}
				
					// Basic Model 1 CDS state with no Start/Stop States
					if (prev_st == "GU0" and st == "CDSfull0") {
						start_id  = id;
						cds_start = svfeat[i].start;
					}
					if (prev_st.substr(0,5) == "accep" and st == "CDSfull0") {
						cds_start = svfeat[i].start;
					}
					if (prev_st == "CDSfull0" and st == "GD0") {
						end_id  = id;
						cds_end = svfeat[i].end -1;
						if (start_id == end_id) {
							std::cout << id     << "\t"
							<< svfeat[i].source << "\t"
							<< "CDS"            << "\t"
							<< cds_start        << "\t"
							<< cds_end          << "\t"
							<< svfeat[i].score  << "\t"
							<< svfeat[i].strand << "\t"
							<< id               << "." << path_count << ":percentage_reps=" << std::fixed << std::setprecision(3) << perc_reps << "%\n";
						}
					}
				
					prev_st = st;
				}
				svfeat.clear();
				path_count++;
			}
        } else {
			// Viterbi ------------------------------------------------------
			trell.viterbi();
			StochHMM::traceback_path tb(&hmm);
			trell.traceback(tb);
			tb.gff(feat, job->getSeqs()->getHeader());
		
			for (size_t i=0; i < feat.size(); i++) {
				std::string st = feat[i].feature;
				std::string id  = feat[i].seqname;
				if (id.substr(0,1) == ">") {id.erase(0,1);}
				
				// DEBUG
				//std::cout << i << "\t" << st << "\t" << feat[i].start << "\t" << feat[i].end << "\n";
				//======
				
				// Standard Model 3 CDS states
				if (prev_st == "GU0" and st == "start0") {
					start_id  = id;
					cds_start = feat[i].start;
				}
				if (prev_st.substr(0,3) == "cds" and st.substr(0,3) == "don") {
					cds_end = feat[i].end -1;
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
		
				// Basic Model 1 CDS state
				if (prev_st.substr(0,3) == "CDS" and st.substr(0,3) == "don") {
					cds_end = feat[i].end -1;
					std::cout   << id              << "\t" 
								<< feat[i].source  << "\t" 
								<< "CDS"           << "\t"
								<< cds_start       << "\t"
								<< cds_end         << "\t"
								<< feat[i].score   << "\t"
								<< feat[i].strand  << "\t"
								<< id              << "\n";
				}
				if (prev_st.substr(0,5) == "accep" and st.substr(0,3) == "CDS") {
					cds_start = feat[i].start;
				}

				// Basic Model 1 CDS state with no Start/Stop States
				if (prev_st == "GU0" and st == "CDSfull0") {
					start_id  = id;
					cds_start = feat[i].start;
				}
				if (prev_st.substr(0,5) == "accep" and st == "CDSfull0") {
					cds_start = feat[i].start;
				}
				if (prev_st == "CDSfull0" and st == "GD0") {
					end_id  = id;
					cds_end = feat[i].end -1;
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
        }
// 		tb.print_gff(job->getHeader());               // Viterbi: Single entry print
//        paths.print_gff(job->getSeqs()->getHeader()); // Stochastic Viterbi: Single entry print
		job = jobs.getJob();
		feat.clear();
	}	
	
	/* Global Cleanup */
	if (PROFILE != NULL) {
		p7_hmm_Destroy(PROFILE);
	}
	esl_alphabet_Destroy(ALPHABET);
	return 0;
	
// 	printf("%d tb_eval functions called\n", FCOUNT);
}

