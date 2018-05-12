#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char *argv[]) {
	char path[1000], seqp[100], seq[4];
	int n = 0, t = 0, i, prob = 0, help = 1, cnt = 0;
	double pa = 0.0, pc = 0.0, pg = 0.0, pt = 0.0, count = 0.0;
	FILE* fa;
	//parsing command line
	if (argc < 1) help = 1;
	for (i = 1; i < argc; i++)
    	{
		if (strcmp(argv[i], "-genf") == 0) { strcpy(&path[0], argv[++i]); help = 0; continue; }		//destination file
		if (strcmp(argv[i], "-n") == 0) { n = atoi(argv[++i]); help = 0; continue; }         		//length of reads
		if (strcmp(argv[i], "-t") == 0) { t = atoi(argv[++i]); help = 0; continue; }  				//number of reads
		if (strcmp(argv[i], "-pa") == 0) { pa = atof(argv[++i]); help = 0; continue; }     			//probability A
		if (strcmp(argv[i], "-pc") == 0) { pc = atof(argv[++i]); help = 0; continue; }				//probability C
		if (strcmp(argv[i], "-pg") == 0) { pg = atof(argv[++i]); help = 0; continue; }          	//probability G     
		if (strcmp(argv[i], "-pt") == 0) { pt = atof(argv[++i]); help = 0; continue; }          	//probability t     
		if (strcmp(argv[i], "-h") == 0) { help = 1; continue; }                         			//help           
    	}
	//probability check
	if((pa == 0.0) && (pc == 0.0) && (pc == 0.0) && (pt == 0.0)) {
		prob = 0;
	}
	else {
		if(pa+pc+pg+pt < 1.0) {
			count = 0.0;
			if(pa == 0.0) count++;
			if(pc == 0.0) count++;
			if(pg == 0.0) count++;
			if(pt == 0.0) count++;
			if(pa == 0.0) pa = (1.0 - (pc+pg+pt))/count;
			if(pc == 0.0) pc = (1.0 - (pa+pg+pt))/count;
			if(pg == 0.0) pg = (1.0 - (pc+pa+pt))/count;
			if(pt == 0.0) pt = (1.0 - (pc+pg+pa))/count;
			prob = 1;
		}
		if(pa+pc+pg+pt > 1.0) {
			prob = 0;
			printf("Sum of probability A, C, G, T greater than 1.0\n\n");
			exit(1);
		}
	}
	//help print
	if ((help == 1) || (strcmp(path,"") == 0) || (n < 1) || (t < 1)) {
		printf("Solution developed by Martini Davide and Nigro Simone\n\n");
		printf("-genf <pathfile> //insert the file which the positions of turbines\n");
		printf("-n <value> //insert the length of reads (required integer non-negative values)\n");
		printf("-t <value> //insert the number of reads (required integer non-negative values)\n");
		printf("-pa <value> //insert the probability of nucleotide A (required double non-negative values < 1.0)\n");
		printf("-pc <value> //insert the probability of nucleotide C (required double non-negative values < 1.0)\n");
		printf("-pg <value> //insert the probability of nucleotide G (required double non-negative values < 1.0)\n");
		printf("-pt <value> //insert the probability of nucleotide T (required double non-negative values < 1.0)\n");
		printf("-h //help\n");
		exit(1);
	}
	//ok now i have all i need
	else {
		//open file
		fa = fopen(&path[0], "w");
		if(fa ==NULL) { perror("Error creating file\n"); exit(1); }
		fprintf(fa, ">FASTA Dataset\n");
		if(prob == 0) {
			seq[0] = 'A';
			seq[1] = 'C';
			seq[2] = 'G';
			seq[3] = 'T';
			for(int j = 0; j < t; j++) {
				for(int k = 0; k < n; k++) {
					i = rand() % 4;
					if(seq[i] == 'A') fprintf(fa,"A");
					if(seq[i] == 'C') fprintf(fa,"C");
					if(seq[i] == 'G') fprintf(fa,"G");
					if(seq[i] == 'T') fprintf(fa,"T");
				}
				fprintf(fa,"\n");
			}
			fprintf(fa,"\n");	
		}
		else {
			cnt = 0; 
			while(cnt < 100) {
				if(cnt < pa*100.0)					//insert A
					seqp[cnt] = 'A';
				else if(cnt < (pa+pc)*100.0)		//insert C
					seqp[cnt] = 'C';
				else if(cnt < (pa+pc+pg)*100.0)		//insert G
					seqp[cnt] = 'G';
				else if(cnt < (pa+pc+pg+pt)*100.0)	//insert T
					seqp[cnt] = 'T';
				cnt++;
			}
			for(int j = 0; j < t; j++) {
				for(int k = 0; k < n; k++) {
					i = rand() % 100;
					if(seqp[i] == 'A') fprintf(fa,"A");
					if(seqp[i] == 'C') fprintf(fa,"C");
					if(seqp[i] == 'G') fprintf(fa,"G");
					if(seqp[i] == 'T') fprintf(fa,"T");
				}
				fprintf(fa,"\n");
			}
			fprintf(fa,"\n");			
		}
		fclose(fa);	
	}
}