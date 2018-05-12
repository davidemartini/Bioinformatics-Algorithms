#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TMAX 201
#define KBIG 64
#define KSMALL 30

int FASTAgen(char *path, int n, int t, double pa, double pc, double pg, double pt) {
	int prob;
	double count;
	FILE *fa;
	char seq[4], seqp[100];
	if ((pa == 0.0) && (pc == 0.0) && (pc == 0.0) && (pt == 0.0)) {
		prob = 0;
	}
	else {
		if (pa + pc + pg + pt < 1.0) {
			count = 0.0;
			if (pa == 0.0) count++;
			if (pc == 0.0) count++;
			if (pg == 0.0) count++;
			if (pt == 0.0) count++;
			if (pa == 0.0) pa = (1.0 - (pc + pg + pt)) / count;
			if (pc == 0.0) pc = (1.0 - (pa + pg + pt)) / count;
			if (pg == 0.0) pg = (1.0 - (pc + pa + pt)) / count;
			if (pt == 0.0) pt = (1.0 - (pc + pg + pa)) / count;
			prob = 1;
		}
		if (pa + pc + pg + pt > 1.0) {
			prob = 0;
			printf("Sum of probability A, C, G, T greater than 1.0\n\n");
			return 1;
		}
	}
	fa = fopen(&path[0], "w");
	if (fa == NULL) { perror("Error creating file\n"); exit(1); }
	fprintf(fa, ">FASTA Dataset\n");
	if (prob == 0) {
		seq[0] = 'A';
		seq[1] = 'C';
		seq[2] = 'G';
		seq[3] = 'T';
		for (int j = 0; j < t; j++) {
			for (int k = 0; k < n; k++) {
				int i = rand() % 4;
				if (seq[i] == 'A') fprintf(fa, "A");
				if (seq[i] == 'C') fprintf(fa, "C");
				if (seq[i] == 'G') fprintf(fa, "G");
				if (seq[i] == 'T') fprintf(fa, "T");
			}
			fprintf(fa, "\n");
		}
		fprintf(fa, "\n");
	}
	else {
		int cnt = 0;
		while (cnt < 100) {
			if (cnt < pa*100.0)							//insert A
				seqp[cnt] = 'A';
			else if (cnt < (pa + pc)*100.0)				//insert C
				seqp[cnt] = 'C';
			else if (cnt < (pa + pc + pg)*100.0)		//insert G
				seqp[cnt] = 'G';
			else if (cnt < (pa + pc + pg + pt)*100.0)	//insert T
				seqp[cnt] = 'T';
			cnt++;
		}
		for (int j = 0; j < t; j++) {
			for (int k = 0; k < n; k++) {
				int i = rand() % 100;
				if (seqp[i] == 'A') fprintf(fa, "A");
				if (seqp[i] == 'C') fprintf(fa, "C");
				if (seqp[i] == 'G') fprintf(fa, "G");
				if (seqp[i] == 'T') fprintf(fa, "T");
			}
			fprintf(fa, "\n");
		}
		fprintf(fa, "\n");
	}
	fclose(fa);
	return 0;
}

int FASTAgenbig(char *path, int n, double t, double pa, double pc, double pg, double pt) {
	int prob;
	double count;
	FILE *fa;
	char seq[4], seqp[100];
	if ((pa == 0.0) && (pc == 0.0) && (pc == 0.0) && (pt == 0.0)) {
		prob = 0;
	}
	else {
		if (pa + pc + pg + pt < 1.0) {
			count = 0.0;
			if (pa == 0.0) count++;
			if (pc == 0.0) count++;
			if (pg == 0.0) count++;
			if (pt == 0.0) count++;
			if (pa == 0.0) pa = (1.0 - (pc + pg + pt)) / count;
			if (pc == 0.0) pc = (1.0 - (pa + pg + pt)) / count;
			if (pg == 0.0) pg = (1.0 - (pc + pa + pt)) / count;
			if (pt == 0.0) pt = (1.0 - (pc + pg + pa)) / count;
			prob = 1;
		}
		if (pa + pc + pg + pt > 1.0) {
			prob = 0;
			printf("Sum of probability A, C, G, T greater than 1.0\n\n");
			return 1;
		}
	}
	fa = fopen(&path[0], "w");
	if (fa == NULL) { perror("Error creating file\n"); exit(1); }
	fprintf(fa, ">FASTA Dataset\n");
	if (prob == 0) {
		seq[0] = 'A';
		seq[1] = 'C';
		seq[2] = 'G';
		seq[3] = 'T';
		for (int j = 0; j < t; j++) {
			for (int k = 0; k < n; k++) {
				int i = rand() % 4;
				if (seq[i] == 'A') fprintf(fa, "A");
				if (seq[i] == 'C') fprintf(fa, "C");
				if (seq[i] == 'G') fprintf(fa, "G");
				if (seq[i] == 'T') fprintf(fa, "T");
			}
			fprintf(fa, "\n");
		}
		fprintf(fa, "\n");
	}
	else {
		int cnt = 0;
		while (cnt < 100) {
			if (cnt < pa*100.0)							//insert A
				seqp[cnt] = 'A';
			else if (cnt < (pa + pc)*100.0)				//insert C
				seqp[cnt] = 'C';
			else if (cnt < (pa + pc + pg)*100.0)		//insert G
				seqp[cnt] = 'G';
			else if (cnt < (pa + pc + pg + pt)*100.0)	//insert T
				seqp[cnt] = 'T';
			cnt++;
		}
		for (int j = 0; j < t; j++) {
			for (int k = 0; k < n; k++) {
				int i = rand() % 100;
				if (seqp[i] == 'A') fprintf(fa, "A");
				if (seqp[i] == 'C') fprintf(fa, "C");
				if (seqp[i] == 'G') fprintf(fa, "G");
				if (seqp[i] == 'T') fprintf(fa, "T");
			}
			fprintf(fa, "\n");
		}
		fprintf(fa, "\n");
	}
	fclose(fa);
	return 0;
}

void set_timeint(int n, int t, int k, char *time) {
	sprintf(time, "../../../../../../usr/bin/time --output=./Prestation/n%dt%dk%d.txt --verbose", n, t, k);
}

void set_timedbl(int n, double t, int k, char *time) {
	sprintf(time, "../../../../../../usr/bin/time --output=./Prestation/n%dt%.fk%d.txt --verbose", n, t, k);
}

char* parse_prestation(char *pre, char* print) {
	FILE *fp;
	char x, buffer[100000];
	int res, index = 0, row = 0, column = 0;
	fp = fopen(pre, "r");
	if (fp == NULL) { perror("Error open file prestation"); exit(1); }
	res=fscanf(fp,"%c",&x);
	while(index < 2) {
		if(x == '"') index++;
		res=fscanf(fp,"%c",&x);
	}
	index = 0;
	column = 0;
	while(res!=EOF) {
		if(x!='\t') {
			if((x>=48)&&(x<=57)) {
				while(x!='\n') {
					buffer[index]=x;
					index++;
					buffer[index]='\0';
					res=fscanf(fp,"%c",&x);
				}
				if ((column == 0) || (column == 2) || (column == 3) || (column == 8)) {
					if (row == 0) sprintf(print, "%s", buffer);
					else sprintf(print, "%s;%s", print, buffer);
					row++;
				}
				column++;
				index = 0;
			}
		}
		res=fscanf(fp,"%c",&x);
	}
	fclose(fp);
	return &print[0];
}

char* first_line(char* print) {
	sprintf(print, "User time (seconds);Percent of CPU;Wall clock time (h:mm:ss or m:ss);Maximum resident set size (kbytes);n;t;k");
	return print;
}

char* first_linefin(char* print) {
	sprintf(print, "User time (seconds);Percent of CPU;Wall clock time (h:mm:ss or m:ss);Maximum resident set size (kbytes);");
	sprintf(print, "%sUser time (seconds);Percent of CPU;Wall clock time (h:mm:ss or m:ss);Maximum resident set size (kbytes);", print);
	sprintf(print, "%sUser time (seconds);Percent of CPU;Wall clock time (h:mm:ss or m:ss);Maximum resident set size (kbytes);n;t;k\n", print);
	return print;
}

void rm_csv() {
	system("rm BFKC0.csv");
	system("rm BFKC1.csv");
	system("rm BFKC2.csv");
	system("rm DSK0.csv");
	system("rm DSK1.csv");
	system("rm DSK2.csv");
	system("rm TURTLE0.csv");
	system("rm TURTLE1.csv");
	system("rm TURTLE2.csv");
	system("rm KMC0.csv");
	system("rm KMC1.csv");
	system("rm KMC2.csv");
	//big
	system("rm DSKbig0.csv");
	system("rm DSKbig1.csv");
	system("rm DSKbig2.csv");
	system("rm TURTLEbig0.csv");
	system("rm TURTLEbig1.csv");
	system("rm TURTLEbig2.csv");
	system("rm KMCbig0.csv");
	system("rm KMCbig1.csv");
	system("rm KMCbig2.csv");
}

void make_csvfin(FILE* ftime, FILE* f1, FILE* f2, FILE* f3) {
	int res, column = 0, flag = 0;
	char x;
	for (int i = 0; i < 111; i++) {
		res = fscanf(f1, "%c", &x);
		res = fscanf(f2, "%c", &x);
		res = fscanf(f3, "%c", &x);
	}
	column = 0;
	flag = 0;
	res = fscanf(f1, "%c", &x);
	while (res != EOF) {
		if (x == ';') column++;
		if (column <= 3) {
			if (flag == 0) {
				if (x != '%') fprintf(ftime, "%c", x);
				res = fscanf(f1, "%c", &x);
			}
			else if (flag == 1) {
				if (x != '%') fprintf(ftime, "%c", x);
				res = fscanf(f2, "%c", &x);
			}
			else if (flag == 2) {
				if (x != '%') fprintf(ftime, "%c", x);
				res = fscanf(f3, "%c", &x);
			}
		}
		else {
			if (flag == 0) {
				if (x != '%') fprintf(ftime, "%c", x);
				flag++;
				column = 0;
				while (x != '\n') res = fscanf(f1, "%c", &x);
				//res = fscanf(f1, "%c", &x);
				res = fscanf(f2, "%c", &x);
			}
			else if (flag == 1) {
				if (x != '%') fprintf(ftime, "%c", x);
				flag++;
				column = 0;
				while (x != '\n') res = fscanf(f2, "%c", &x);
				//res = fscanf(f2, "%c", &x);
				res = fscanf(f3, "%c", &x);
			}
			else if (flag == 2) {
				flag = 0;
				column = 0;
				while (x != '\n') {
					if (x != '%') fprintf(ftime, "%c", x);
					res = fscanf(f3, "%c", &x);
				}
				//res = fscanf(f3, "%c", &x);
			}
		}
	}
}

void close_file(FILE *f1, FILE *f2, FILE *f3, FILE *f4) {
	fclose(f1);
	fclose(f2);
	fclose(f3);
	fclose(f4);
}

void init_filecsv(FILE *ftime) {
	char path[1000];
	for(int i = 0; i < 3; i++) {
		//BFKC
		sprintf(path, "BFKC%d.csv", i);
		ftime = fopen(path, "w");
		if (ftime == NULL) { perror("Error open file BFKC.csv"); exit(1); }
		fprintf(ftime, "%s\n", first_line(&path[0]));
		fclose(ftime);
		//DSK
		sprintf(path, "DSK%d.csv", i);
		ftime = fopen(path, "w");
		if (ftime == NULL) { perror("Error open file DSK.csv"); exit(1); }
		fprintf(ftime, "%s\n", first_line(&path[0]));
		fclose(ftime);
		//TURTLE
		sprintf(path, "TURTLE%d.csv", i);
		ftime = fopen(path, "w");
		if (ftime == NULL) { perror("Error open file TURTLE.csv"); exit(1); }
		fprintf(ftime, "%s\n", first_line(&path[0]));
		fclose(ftime);
		//KMC
		sprintf(path, "KMC%d.csv", i);
		ftime = fopen(path, "w");
		if (ftime == NULL) { perror("Error open file KMC.csv"); exit(1); }
		fprintf(ftime, "%s\n", first_line(&path[0]));
		fclose(ftime);
	}
}

void init_big_filecsv(FILE *ftime) {
	char path[1000];
	for (int i = 0; i < 3; i++) {
		//DSKbig
		sprintf(path, "DSKbig%d.csv", i);
		ftime = fopen(path, "w");
		if (ftime == NULL) { perror("Error open file DSKbig.csv"); exit(1); }
		fprintf(ftime, "%s\n", first_line(&path[0]));
		fclose(ftime);
		//TURTLEbig
		sprintf(path, "TURTLEbig%d.csv", i);
		ftime = fopen(path, "w");
		if (ftime == NULL) { perror("Error open file TURTLEbig.csv"); exit(1); }
		fprintf(ftime, "%s\n", first_line(&path[0]));
		fclose(ftime);
		//KMCbig
		sprintf(path, "KMCbig%d.csv", i);
		ftime = fopen(path, "w");
		if (ftime == NULL) { perror("Error open file KMCbig.csv"); exit(1); }
		fprintf(ftime, "%s\n", first_line(&path[0]));
		fclose(ftime);
	}
}

void mkdir() {
	system("mkdir Kmers");
	system("mkdir Dataset");
	system("mkdir DSKoutput");
	system("mkdir Turtleoutput");
	system("mkdir KMCoutput");
	system("mkdir Prestation");
}

void rmdir() {
	system("rm -rf Kmers");
	//system("rm -rf Dataset");
	system("rm -rf DSKoutput");
	system("rm -rf Turtleoutput");
	system("rm -rf KMCoutput");
	system("rm -rf Prestation");
	system("rm trash*");
}

int main(int argc, char *argv[]) {
	//variables
	int flag = 0, n, t;
	FILE *ftime, *f1, *f2, *f3;
	char command[10000], path[1000], time[100], print[10000];
	mkdir();
	init_filecsv(ftime);
	init_big_filecsv(ftime);
	for(int i = 0; i < 3; i++) {
		//****************************************************************
		//FIRST TEST: BFKC, DSK, TURTLE, KMC for small dataset
		//****************************************************************
		flag = 0;
		//start counting
		while (flag <= 3) {
			n = 40;
			t = 40;
			while (t < TMAX) {
				//generate FASTA dataset
				sprintf(path, "./Dataset/datasetn%dt%d%d.fa", n, t, i);
				FASTAgen(&path[0], n, t, 0.0, 0.0, 0.0, 0.0);
				if (flag == 0) {
					//BRUTEFORCE TEST
					sprintf(path, "BFKC%d.csv", i);
					ftime = fopen(path,"a");
					if (ftime == NULL) { perror("Error open file BFKC.csv"); exit(1); }
					for (int k = 5; k <= KSMALL; k=k+5) {
						//sprintf(command, "../BFKC/BFKC -fa ./Dataset/datasetn%dt%d.fa -k %d", n, t, k);
						set_timeint(n, t, k, &time[0]);
						sprintf(command, "%s ../BFKC/BFKC -fa ./Dataset/datasetn%dt%d%d.fa -k %d", &time[0], n, t, i, k);
						system(command);
						printf("%s\n", command);
						sprintf(path, "./Prestation/n%dt%dk%d.txt", n, t, k);
						fprintf(ftime, "%s;%d;%d;%d\n", parse_prestation(&path[0], &print[0]), n, t, k);
					}
					fclose(ftime);
				}
				if (flag == 1) {
					//DSK TEST
					sprintf(path, "DSK%d.csv", i);
					ftime = fopen(path,"a");
					if (ftime == NULL) { perror("Error open file DSK.csv"); exit(1); }
					for (int k = 5; k <= KSMALL; k=k+5) {
						set_timeint(n, t, k, &time[0]);
						sprintf(command, "%s ../../Tools/dsk-2.1.0-Linux/bin/dsk -file ./Dataset/datasetn%dt%d%d.fa -verbose 0 -kmer-size %d -out ./Kmers/kmersn%dt%d > ./DSKoutput/n%dt%dk%d.txt", &time[0], n, t, i, k, n, t, n, t, k);
						system(command);
						printf("%s\n", command);
						sprintf(path, "./Prestation/n%dt%dk%d.txt", n, t, k);
						fprintf(ftime, "%s;%d;%d;%d\n", parse_prestation(&path[0], &print[0]), n, t, k);
					}
					fclose(ftime);
				}
				if (flag == 2) {
					//TURTLE TEST
					sprintf(path, "TURTLE%d.csv", i);
					ftime = fopen(path,"a");
					if (ftime == NULL) { perror("Error open file TURTLE.csv"); exit(1); }
					for (int k = 5; k <= KSMALL; k=k+5) {
						set_timeint(n, t, k, &time[0]);
						sprintf(command, "%s ../../Tools/Turtle-0.3/scTurtle32 -i ./Dataset/datasetn%dt%d%d.fa -o ./Kmers/kmersn%dt%d -k %d -t 1 -n 5000000 > ./Turtleoutput/n%dt%dk%d.txt", &time[0], n, t, i, n, t, k, n, t, k);
						system(command);
						printf("%s\n", command);
						sprintf(path, "./Prestation/n%dt%dk%d.txt", n, t, k);
						fprintf(ftime, "%s;%d;%d;%d\n", parse_prestation(&path[0], &print[0]), n, t, k);
					}
					fclose(ftime);
				}
				if (flag == 3) {
					//KMC TEST
					sprintf(path, "KMC%d.csv", i);
					ftime = fopen(path,"a");
					if (ftime == NULL) { perror("Error open file KMC.csv"); exit(1); }
					for (int k = 5; k <= KSMALL; k=k+5) {
						set_timeint(n, t, k, &time[0]);
						sprintf(command, "%s ../../Tools/KMC-master/bin/kmc -k%d -m8 -fa ./Dataset/datasetn%dt%d%d.fa ./KMCoutput/n%dt%dk%d.txt ./KMCoutput", &time[0], k, n, t, i, n, t, k);
						system(command);
						printf("%s\n", command);
						sprintf(path, "./Prestation/n%dt%dk%d.txt", n, t, k);
						fprintf(ftime, "%s;%d;%d;%d\n", parse_prestation(&path[0], &print[0]), n, t, k);
					}
					fclose(ftime);
				}
				if (n < 80) n=n+5;
				t=t+5;
			}
			flag++;
		}
		//****************************************************************
		//SECOND TEST: DSK, TURTLE, KMC for BIG dataset
		//***************************************************************
		flag = 0;
		n = 80;
		double trow = 1e+4;
		while (flag <= 2) {
			sprintf(path, "./Dataset/datasetn%dt%.f%d.fa", n, trow, i);
			printf("%s %d\n", path, flag);
			FASTAgenbig(&path[0], n, trow, 0.0, 0.0, 0.0, 0.0);
			//generate FASTA dataset
			if (flag == 0) {
				//DSK TEST
				sprintf(path, "DSKbig%d.csv", i);
				ftime = fopen(path, "a");
				if (ftime == NULL) { perror("Error open file DSKbig.csv"); exit(1); }
				for (int k = 50; k <= KBIG; k=k+5) {
					if (k == KBIG+1) k = KBIG;
					set_timedbl(n, trow, k, &time[0]);
					sprintf(command, "%s ../../Tools/dsk-2.1.0-Linux/bin/dsk -file ./Dataset/datasetn%dt%.f%d.fa -verbose 0 -kmer-size %d -out ./Kmers/kmersn%dt%.f > ./DSKoutput/n%dt%.fk%d.txt", &time[0], n, trow, i, k, n, trow, n, trow, k);
					system(command);
					printf("%s\n", command);
					sprintf(path, "./Prestation/n%dt%.fk%d.txt", n, trow, k);
					fprintf(ftime, "%s;%d;%.f;%d\n", parse_prestation(&path[0], &print[0]), n, trow, k);
				}
				fclose(ftime);
			}
			if (flag == 1) {
				//TURTLE TEST
				sprintf(path, "TURTLEbig%d.csv", i);
				ftime = fopen(path, "a");
				if (ftime == NULL) { perror("Error open file TURTLEbig.csv"); exit(1); }
				for (int k = 50; k <= KBIG; k=k+5) {
					if (k == KBIG+1) k = KBIG;
					set_timedbl(n, trow, k, &time[0]);
					sprintf(command, "%s ../../Tools/Turtle-0.3/scTurtle64 -i ./Dataset/datasetn%dt%.f%d.fa -o ./Kmers/kmersn%dt%.f -k %d -t 1 -n 5000000 > ./Turtleoutput/n%dt%.fk%d.txt", &time[0], n, trow, i, n, trow, k, n, trow, k);
					system(command);
					printf("%s\n", command);
					sprintf(path, "./Prestation/n%dt%.fk%d.txt", n, trow, k);
					fprintf(ftime, "%s;%d;%.f;%d\n", parse_prestation(&path[0], &print[0]), n, trow, k);
				}
				fclose(ftime);
			}
			if (flag == 2) {
				//KMC TEST
				sprintf(path, "KMCbig%d.csv", i);
				ftime = fopen(path, "a");
				if (ftime == NULL) { perror("Error open file KMCbig.csv"); exit(1); }
				for (int k = 50; k <= KBIG; k=k+5) {
					if (k == KBIG+1) k = KBIG;
					set_timedbl(n, trow, k, &time[0]);
					sprintf(command, "%s ../../Tools/KMC-master/bin/kmc -k%d -m8 -fa ./Dataset/datasetn%dt%.f%d.fa ./KMCoutput/n%dt%.fk%d.txt ./KMCoutput", &time[0], k, n, trow, i, n, trow, k);
					system(command);
					printf("%s\n", command);
					sprintf(path, "./Prestation/n%dt%.fk%d.txt", n, trow, k);
					fprintf(ftime, "%s;%d;%.f;%d\n", parse_prestation(&path[0], &print[0]), n, trow, k);
				}
				fclose(ftime);
			}
			trow = trow * 10;
			if (trow > 1e+8) {
				flag++;
				trow = 1e+4;
			}
		}		
	}
	rmdir();
	//open files for BFKC
	ftime = fopen("BFKCfin.csv", "w");
	if (ftime == NULL) { perror("Error open file BFKCfin.csv"); exit(1); }
	f1 = fopen("BFKC0.csv", "r");
	if (f1 == NULL) { perror("Error open file BFKC0.csv"); exit(1); }
	f2 = fopen("BFKC1.csv", "r");
	if (f2 == NULL) { perror("Error open file BFKC1.csv"); exit(1); }
	f3 = fopen("BFKC2.csv", "r");
	if (f3 == NULL) { perror("Error open file BFKC2.csv"); exit(1); }
	fprintf(ftime, "%s", first_linefin(print));
	make_csvfin(ftime, f1, f2, f3);
	close_file(ftime, f1, f2, f3);
	//open files for DSK
	ftime = fopen("DSKfin.csv", "w");
	if (ftime == NULL) { perror("Error open file DSKfin.csv"); exit(1); }
	f1 = fopen("DSK0.csv", "r");
	if (f1 == NULL) { perror("Error open file DSK0.csv"); exit(1); }
	f2 = fopen("DSK1.csv", "r");
	if (f2 == NULL) { perror("Error open file DSK1.csv"); exit(1); }
	f3 = fopen("DSK2.csv", "r");
	if (f3 == NULL) { perror("Error open file DSK2.csv"); exit(1); }
	fprintf(ftime, "%s", first_linefin(print));
	make_csvfin(ftime, f1, f2, f3);
	close_file(ftime, f1, f2, f3);
	//open files for Turtle
	ftime = fopen("TURTLEfin.csv", "w");
	if (ftime == NULL) { perror("Error open file TURTLEfin.csv"); exit(1); }
	f1 = fopen("TURTLE0.csv", "r");
	if (f1 == NULL) { perror("Error open file TURTLE0.csv"); exit(1); }
	f2 = fopen("TURTLE1.csv", "r");
	if (f2 == NULL) { perror("Error open file TURTLE1.csv"); exit(1); }
	f3 = fopen("TURTLE2.csv", "r");
	if (f3 == NULL) { perror("Error open file TURTLE2.csv"); exit(1); }
	fprintf(ftime, "%s", first_linefin(print));
	make_csvfin(ftime, f1, f2, f3);
	close_file(ftime, f1, f2, f3);
	//open files for KMC
	ftime = fopen("KMCfin.csv", "w");
	if (ftime == NULL) { perror("Error open file KMCfin.csv"); exit(1); }
	f1 = fopen("KMC0.csv", "r");
	if (f1 == NULL) { perror("Error open file KMC0.csv"); exit(1); }
	f2 = fopen("KMC1.csv", "r");
	if (f2 == NULL) { perror("Error open file KMC1.csv"); exit(1); }
	f3 = fopen("KMC2.csv", "r");
	if (f3 == NULL) { perror("Error open file KMC2.csv"); exit(1); }
	fprintf(ftime, "%s", first_linefin(print));
	make_csvfin(ftime, f1, f2, f3);
	close_file(ftime, f1, f2, f3);
	//open files for DSKbig
	ftime = fopen("DSKbigfin.csv", "w");
	if (ftime == NULL) { perror("Error open file DSKbigfin.csv"); exit(1); }
	f1 = fopen("DSKbig0.csv", "r");
	if (f1 == NULL) { perror("Error open file DSKbig0.csv"); exit(1); }
	f2 = fopen("DSKbig1.csv", "r");
	if (f2 == NULL) { perror("Error open file DSKbig1.csv"); exit(1); }
	f3 = fopen("DSKbig2.csv", "r");
	if (f3 == NULL) { perror("Error open file DSKbig2.csv"); exit(1); }
	fprintf(ftime, "%s", first_linefin(print));
	make_csvfin(ftime, f1, f2, f3);
	close_file(ftime, f1, f2, f3);
	//open files for TURTLEbig
	ftime = fopen("TURTLEbigfin.csv", "w");
	if (ftime == NULL) { perror("Error open file TURTLEbigfin.csv"); exit(1); }
	f1 = fopen("TURTLEbig0.csv", "r");
	if (f1 == NULL) { perror("Error open file TURTLEbig0.csv"); exit(1); }
	f2 = fopen("TURTLEbig1.csv", "r");
	if (f2 == NULL) { perror("Error open file TURTLEbig1.csv"); exit(1); }
	f3 = fopen("TURTLEbig2.csv", "r");
	if (f3 == NULL) { perror("Error open file TURTLEbig2.csv"); exit(1); }
	fprintf(ftime, "%s", first_linefin(print));
	make_csvfin(ftime, f1, f2, f3);
	close_file(ftime, f1, f2, f3);
	//open files for KMCbig
	ftime = fopen("KMCbigfin.csv", "w");
	if (ftime == NULL) { perror("Error open file KMCbigfin.csv"); exit(1); }
	f1 = fopen("KMCbig0.csv", "r");
	if (f1 == NULL) { perror("Error open file KMCbig0.csv"); exit(1); }
	f2 = fopen("KMCbig1.csv", "r");
	if (f2 == NULL) { perror("Error open file KMCbig1.csv"); exit(1); }
	f3 = fopen("KMCbig2.csv", "r");
	if (f3 == NULL) { perror("Error open file KMCbig2.csv"); exit(1); }
	fprintf(ftime, "%s", first_linefin(print));
	make_csvfin(ftime, f1, f2, f3);
	close_file(ftime, f1, f2, f3);
	rm_csv();
	system("shutdown now");
}