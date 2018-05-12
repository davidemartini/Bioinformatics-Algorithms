/*
Turtle: Identifying frequent k-mers with cache-efficient algorithms
Copyright 2014 Rajat Shuvro Roy, Alexander Schliep

This file is part of Turtle-0.3.

    Turtle-0.3 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    Turtle-0.3 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Turtle-0.3.  If not, see <http://www.gnu.org/licenses/>.
*/

//#define T __uint128_t
#define T uint64_t
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <template_helpers.hpp>
#include <batch_pbf.hpp>
#include <time.h>

using namespace std;
bool quake;

void usage(){
	cout<<"aTurtle32 Usage:"<<endl;
	cout <<"aTurtle32 [arguments]"<<endl;
	cout<<"example: ./aTurtle32 -f 1Mreads.fq -o kmer_counts -k 31 -n 35000000"<<endl;
	cout <<"-i \t input reads file in fasta format."<<endl;
	cout <<"-f \t input reads file in fastq format. This is mutually exclusive with -i."<<endl;
	cout <<"-o \t output file name. k-mers and their counts are stored in fasta format (headers indicating frequency)."<<endl;
	cout <<"-q \t output file name. k-mers and their counts are stored in tab delimited fromat (quake compatible)."<<endl;
	cout <<"-k \t k-mer length."<<endl;
	cout <<"-n \t Expected number of k-mers. "<<endl;
	cout <<"-s \t The approximate amount of space (in GB) to be used. It is used to indirectly compute -n and is mutually exclusive with -n. When both -n and -s are specified, the one that appears last is used."<<endl;
	cout <<"-h \t Print this help menu."<<endl;
	cout <<"-v \t Print software version."<<endl;
	cout <<""<<endl;
	}
int main(int argc, char *argv[]){
	char input_file_name[1000], output_file_name[1000];
	input_file_name[0]=0;
	output_file_name[0]=0;
	int kmer_length=0 , no_of_threads=0;
	long gen_length=0;
	cout<<"Turtle Copyright (C) 2014 Rajat Shuvro Roy, Alexander Schliep.\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, and you are welcome to redistribute it  under certain conditions. For details see the document COPYING.\n"<<endl;
	cout<<"Parameters received:"<<endl;
	void (*get_read)(ifstream &, char*);
	  while ((argc > 1) )
		{
			if (argv[1][0] != '-'){
				printf("Wrong Argument: %s\n", argv[1]);
				usage();
				return 0;
				}
				
			switch (argv[1][1])
			{
				case 'i':
					cout<< "fasta input \t"<< argv[2];
					strcpy(input_file_name, argv[2]);
					get_read=&getNextRead;
					break;

				case 'f':
					cout <<"fastq input \t"<<argv[2];
					//fastq_input=true;
					strcpy(input_file_name, argv[2]);
					get_read=&getNextRead_fq;
					break;
				case 'o':
					cout <<"ouput \t"<<argv[2];
					strcpy(output_file_name, argv[2]);
					quake=false;
					break;
				case 'q':
					cout <<"ouput \t"<<argv[2];
					strcpy(output_file_name, argv[2]);
					quake=true;
					break;
					
				case 'k':
					kmer_length=atoi(argv[2]);
					cout <<"k-mer length \t"<<kmer_length;
					if (kmer_length>32){
						cout<<" ERROR: aTurtle32 can handle k-mers of length up to 32. Please try scTurtle64."<<endl;
						return 0;
						}
					break;
				case 'n':
					gen_length=atol(argv[2]);
					cout <<"no of k-mers \t"<<gen_length;
					break;	
				case 's':
					gen_length=(atof(argv[2])*(1<<30))/41;
					cout <<"no of k-mers \t"<<gen_length;
					break;	
				case 'h':
					
					usage();
					return 0;
					break;
				case 'v':
					cout<<"Turtle frequent k-mer counter version 0.3"<<endl;
					return 0;
					break;
				
				default:
					printf("Wrong Argument: %s\n", argv[1]);
					usage();
					return 0;
			}
			cout<<endl;
			argv+=2;
			argc-=2;
		}
	if (input_file_name[0]==0) {
		cout<<"Please specify an input file."<<endl;
		usage();
		return 0;
		}
	if (output_file_name[0]==0) {
		cout<<"Please specify an output file."<<endl;
		usage();
		return 0;
		}
	if (kmer_length==0) {
		cout<<"Please specify a k-mer length."<<endl;
		usage();
		return 0;
		}
	if (gen_length==0) {
		cout<<"Please specify the number of expected frequent kmers."<<endl;
		usage();
		return 0;
		}
	long no_of_kmers=gen_length*2;
	clock_t start, end;	
	start=clock();
	
	ifstream reads_file(input_file_name);
	vector<kmer> buffer(no_of_kmers,0); 
	unsigned long kmer_i=0, sorted_i=0;
	if (!reads_file.good()) {
		cout<<"file not found"<<endl;
		return 0;
		}
	char  rd[MAX_RD_LEN];
	get_read(reads_file,rd );
	if (rd==NULL) {cout<< "null rd pointer"<<endl; return 0;}
	int tok_i=0, tok_j;
	long kmer_seen=0, no_of_reads=0, no_of_N=0;
	long no_of_items=0;
	while (rd[0]){//while read not empty
		//cout << rd<<endl;
		no_of_reads++;
		tok_i=0;
		int read_len=strlen(rd);
		T bit_read, bit_read_rc, bit_rd, tmp;
			
		//if (frag==NULL) cout << "null"<<endl;
		while (tok_i<=(read_len-kmer_length)){
			while(rd[tok_i]=='N') {tok_i++;no_of_N++;}
			bool kmer_found=false;
			while (tok_i<=(read_len-kmer_length) && !kmer_found ){
				for (tok_j=tok_i+1; tok_j< tok_i+kmer_length;tok_j++)
					if(rd[tok_j]=='N') {no_of_N++; break;}
				if (tok_j==	tok_i+kmer_length) kmer_found=true;
				else tok_i=tok_j+1;
				}	
			if (!kmer_found) break;
			//we will calculate this only once and then extend per base
			//bit_read = get_bit_read_128(rd+tok_i,kmer_length);
			bit_read=0;
			int i0;
			for(i0=0; i0<(kmer_length-16); i0+=16){
				bit_read<<=32;
				//cout<<rd+tok_i+i0<<endl;
				bit_read|=sse_get_bit_read(rd+tok_i+i0, 16);
				
				}
			//cout<<rd+tok_i+i0<<endl;
			bit_read<<=(2*(kmer_length-i0));
			bit_read|=sse_get_bit_read(rd+tok_i+i0, kmer_length-i0);
			
			/*if (bit_read!=get_bit_read_128(rd+tok_i,kmer_length)){
				cout<<"wrong bit_read extraction"<<endl;
				return 0;
				}*/
			//else cout<<"correct bit-read extraction"<<endl;
			bit_read_rc=get_rev_comp(bit_read,kmer_length);
			unsigned int  bit_seq_i;
			unsigned char *bit_seq_buff, *bit_seq_buff_rc;
			bit_seq_i=16;
			for (; tok_i<= (read_len-kmer_length)&& rd[tok_i+kmer_length-1]!='N' ;tok_i++){
				//for(int i=tok_i; i<tok_i+kmer_length; i++) cout<<rd[i];
				//cout<<endl;
			
				if (bit_read<bit_read_rc) bit_rd=bit_read;
				else bit_rd=bit_read_rc;
				kmer_seen++;
				if(__builtin_expect(kmer_i>=no_of_kmers*.85, 0)){
					compress(buffer, sorted_i, kmer_i);
					cout<<"compressing .."<<endl;
					}
					
				if(__builtin_expect(kmer_i>=no_of_kmers*.85, 0)){
					sorted_i= compress_kmers(buffer, kmer_i);
					kmer_i=sorted_i;
					cout<<"full compress.."<<endl;
					}
				if(__builtin_expect(kmer_i>=no_of_kmers*.85, 0)){
					cout<<"please increase memory"<<endl;
					return 0;
					}
				buffer[kmer_i].bit_kmer=bit_rd;
				buffer[kmer_i++].count=1;
				no_of_items++;
				
				
				
				__m128i m, m_rc;
				if (bit_seq_i==16){
					m=sse_next_16_bit_rd(rd+tok_i+kmer_length);
					//inlining this function below
					bit_seq_buff=(unsigned char *)&m;
					
					m_rc=sse_next_16_bit_rc(rd+tok_i+kmer_length);
					bit_seq_buff_rc=(unsigned char *)&m_rc;
					}
				next_bit_read(bit_read, bit_read_rc, rd[tok_i+kmer_length], kmer_length);
				
				}
			tok_i+=kmer_length;
			}
		//begin=clock();
		get_read(reads_file,rd );
		
		if (strlen(rd)==0) break;
		//strcpy(rd,read.c_str());
	}
	//cout<<"final compression.."<<endl;
	//compress(buffer, sorted_i, kmer_i);
	kmer_i=compress_kmers(buffer, kmer_i);
	
	//cout<<"STATs:"<<"\nNo of reads:\t "<<no_of_reads<<"\nNo of N's:\t "<<no_of_N<<"\nNo of k-mers:\t "<<kmer_seen<<endl;
	cout<<"STATs:"<<"\nNo of reads:\t "<<no_of_reads<<"\nNo of k-mers:\t "<<kmer_seen<<endl;
	cout<< " no of unique kmers found :"<<kmer_i<< endl;
	ofstream *kmer_file=new ofstream(output_file_name);
	char tmp_rd[1000];
	for (long i=0; i<kmer_i ; i++){
		get_string_read(buffer[i].bit_kmer, tmp_rd,kmer_length);
		if (quake)
			(*kmer_file) << tmp_rd<<'\t'<<buffer[i].count<< '\n';	
		else
			(*kmer_file) << '>'<<buffer[i].count<<'\n'<<tmp_rd<< '\n';	
		}	
	(*kmer_file).close();
	/*ofstream *kmer_file=new ofstream(output_file_name, ios::out |ios::binary);
	(*kmer_file).write((char *)(&buffer[0]), kmer_i*sizeof(kmer));
	(*kmer_file).close();
	long freq_item=0;
	for (long i=0; i<kmer_i; i++)
		if (buffer[i].count>1)
			freq_item++;
	cout<< " Frequent items "<<freq_item<<" infrequent item "<< kmer_i-freq_item<< ","<<100.0*(kmer_i-freq_item)/kmer_i<< endl;		
	cout<<"number of items:"<<no_of_items<<endl;*/
	//end=clock();
	//cout <<"total time "<<(end-start)/CLOCKS_PER_SEC<<endl;
	
	return 0;
}
