#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char *argv[]) {
	int help = 1, k = 0, n = 0, t = 0, read, flag, *s;
	char path[1000], buffer[120], file[1000];
	char  x, *res, **dna;
	FILE* fa;
	if (argc < 1) help = 1;
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-fa") == 0) { strcpy(&path[0], argv[++i]); help = 0; continue; }		//destination file    	
		if (strcmp(argv[i], "-k") == 0) { k = atoi(argv[++i]); help = 0; continue; }						//k
		if (strcmp(argv[i], "-h") == 0) { help = 1; continue; }																	//help  
    }
	if (help == 1) {
		//help print
		printf("Solution developed by Davide Martini and Simone Nigro\n");
		printf("-fa <filename> //FASTA dataset\n");
		printf("-k <value> //length of k-mers (integer non-negative values required\n");
		printf("-h //help\n");
	}
	else {
		//open file
		fa = fopen(path, "r");
		if (fa == NULL) {
			perror("Error opening file\n");
			exit(1);
		}
		t = -2;
		n = 0;
		while (1) {
			res = fgets(buffer, 120, fa);
			if (t == 0) {
				for (int i = 0; i < 120; i++) {
					if ((buffer[i] == 'A') || (buffer[i] == 'C') || (buffer[i] == 'G') || (buffer[i] == 'T'))
						n++;
					else
						i = 120;
				}
			}
			if (res == NULL)
				break;
			t++;
		}
		fclose(fa);
		//s positions vector
		s = (int*)malloc(t * sizeof(int*));
		for (int i = 0; i < t; i++)
			s[i] = 0;
		//dna sequences
		dna = (char**)malloc(t * sizeof(char*));
		for (int i = 0; i < t; i++) {
			dna[i] = (char*)malloc(n * sizeof(char));
		}
		//open file
		fa = fopen(path, "r");
		if (fa == NULL) {
			perror("Error opening file\n");
			exit(1);
		}
		for (int i = 0; i < 15; i++) {
			read = fscanf(fa, "%c", &x);
		}
		for (int i = 0; i < t; i++) {
			for (int j = 0; j < n; j++) {
				read = fscanf(fa, "%c", &x);
				dna[i][j] = x;
			}
			read = fscanf(fa, "%c", &x);
		}
		fclose(fa);
		int count = 0;
		for (int j = 0; j < t; j++) {
			for (int i = 0; i < n - k+1; i++) {
				count++;
			}
		}
		//list of kmers
		char **list = (char**)malloc(count * sizeof(char*));
		for (int i = 0; i < count; i++) {
			list[i] = (char*)malloc(k * sizeof(char));
		}
		int indexr = 0;
		int indexc = 0;
		for (int j = 0; j < t; j++) {
			for (int i = 0; i < n - k+1; i++) {
				for (int z = i; z < i+k; z++) {
					list[indexr][indexc] = dna[j][z];
					indexc++;
				}
				indexc = 0;
				indexr++;
			}
		}
		sprintf(file, "./Kmers/kmersn%dt%dk%d.txt", n, t, k);
		//open file
		fa = fopen(file, "w");
		if(fa == NULL) {
			perror("Error open file kmers\n");
			exit(1);
		}
		int flag = 0;
		int lim = 0;
		int exclude[count];
		for(int i = 0; i < count; i++) {
			exclude[i] = 1;
		}
		//count kmers
		for(int i = 0; i < count-1; i++) {
			lim = 1;
			for(int j = i+1; j < count; j++) {
				if(exclude[j]!=0) {
					flag = 0;
					for(int z = 0; z < k; z++) {
						if(list[i][z] == list[j][z]) {
							flag++;
						}
					}
					if(flag == k) {
						lim++;
						exclude[j]=0;
					}
				}
			}
			if((exclude[i]!=0)&&(lim>1)) {
				for(int z = 0; z < k; z++) {
					fprintf(fa, "%c",list[i][z]);
				}
				fprintf(fa," %d\n", lim);
			}
		}
		fclose(fa);
	}
}
