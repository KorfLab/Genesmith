/*****************************************************************************\
 genesmith.cpp
 
 Copyright (C) 2013 Ravi Dandekar & Ian Korf. All rights reserved.
\*****************************************************************************/


#include <StochHMMlib.h>
#include <unistd.h>

static void usage () {
	std::cout
		<< "usage: genesmith [options] <hmm file> <seq file>\n"
		<< "options:\n"
		<< "  -g <file> genetic code (assumes -p and or -f)\n"
		<< "  -p <file> profile HMM (requires -g)\n"
		<< "  -f <file> protein fasta file (requires -g)\n"
		<< std::endl;
	exit(1);
}

int main (int argc, char ** argv) {
	char * fasta_file   = NULL;
	char * profile_file = NULL;
	char * genetic_code = NULL;
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
	
	/* decode & output */
	StochHMM::model hmm;
	StochHMM::seqTracks jobs;
	jobs.loadSeqs(hmm, seq_file);
	StochHMM::seqJob *job=jobs.getJob();
	while (job != NULL) {
		StochHMM::trellis trell(&hmm, job->getSeqs());
		trell.simple_viterbi();
		StochHMM::traceback_path tb(&hmm);
		trell.traceback(tb);
		tb.print_gff(job->getHeader());
		job = jobs.getJob();
	}
	
	return 0;
}

