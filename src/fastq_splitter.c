#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>

// Using code from Heng Li's kseq example - http://lh3lh3.users.sourceforge.net/parsefastq.shtml
// modified from https://github.com/vasisht/fastq_splitter/blob/master/fastq_splitter.c

// declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

int is_valid_fd(int fd) { // Function to check if file handle is open
	return fcntl(fd,F_GETFD) != 1 || errno != EBADF;
}

int main(int argc, char *argv[])
{
	gzFile fp;
	gzFile out;
	kseq_t *seq;
	int l,a,i=0,c=0,numreads;
	char *prefix, *fname;
	
	while ((a = getopt(argc,argv, "p:n:")) >= 0) {
		switch(a) {
			case 'p': prefix = optarg; break;
			case 'n': numreads = atoi(optarg); break;
		}
	} 

	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: fastq_splitter [options] <in.fq.gz>\n");
		fprintf(stderr, "Options -p STR Output prefix \n");
		fprintf(stderr, "        -n INT Number of reads per file \n");
		fprintf(stderr,"\n");
		return 1;
	}

	fp = strcmp(argv[optind],"-")? gzopen(argv[optind],"r") : gzdopen(fileno(stdin),"r"); // Handle both stdin and file 
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file.\n", __func__);
		return 1;
	}
	seq = kseq_init(fp); //  initialize seq
	while ((l = kseq_read(seq)) >= 0) { // read sequence
		i += 1;
		if (i % numreads == 1) {
			c += 1;
			if (is_valid_fd((int) out)) gzclose(out);
			fname = (char*) malloc(strlen(prefix) + 18);
			sprintf(fname,"%s-S%06d.fa.gz",prefix,c); //Set output filename based on the prefix 
			out = gzopen(fname,"w");
			free(fname);
		}
		gzprintf(out,">%s\n",seq->name.s);
		gzputs(out, "\n");
		gzputs(out, seq->seq.s);
		gzputs(out, "\n");
	}
	kseq_destroy(seq); 
	gzclose(fp);  // close input file
	if (is_valid_fd((int) out)) gzclose(out);
	return 0;
}

