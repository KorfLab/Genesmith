/*****************************************************************************\
 genesmith.cpp
 
 Copyright (C) 2013 Ravi Dandekar & Ian Korf. All rights reserved.
\*****************************************************************************/

#include <StochHMMlib.h>



int main (int argc, const char * argv[]) {

	/* command line */
	std::string usage = "usage: genesmith <hmm file> <seq file>";
	
	if (argc != 2) {
		std::cout << usage << std::endl;
		exit(2);
	}
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
