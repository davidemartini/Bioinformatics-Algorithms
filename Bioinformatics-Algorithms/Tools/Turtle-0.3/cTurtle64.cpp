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

#define T __uint128_t
//#define T uint64_t
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <random>
#include <template_helpers.hpp>
#include <patt2_bf.hpp>
//#include <pthread.h>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <atomic>
#include <queue>
#include <time.h>

using namespace std;
bool producer_alive;

mutex *mutices; 
condition_variable is_not_full;
condition_variable *is_not_empty;
mutex *consumer_mutices, prod_mutex; 
condition_variable is_prodQ_not_empty;
condition_variable *is_consumerQ_not_empty;
int no_of_threads,no_of_smb;//no of small buffers
queue<int> producer_Q;

struct small_buff{
	atomic<int> thread_i;
	atomic<bool> full;
	long i;
	T* buffer;
	small_buff(){thread_i=-1;i=0; full=false;};
	~small_buff(){delete buffer;};
	void initialize(long n){buffer=new T[n];};
	} ;


small_buff *sm_buffs;
	
struct comp_thread_data{
	//int job_queue[100];
	//atomic<int> job_q_i;
	queue<int> consumer_Q;
	uint64_t no_of_kmers;
	unsigned long kmer_i, buffer_i, sorted_i;
	//FILE *kmer_file;
	pattern_bf *pbf;
	ofstream *kmer_file;
	T * buffer; 
	long buf_size;
	int kmer_length, thread_i;
	comp_thread_data( uint64_t N, unsigned char* , int, int, int, long, char*);
	comp_thread_data(){return ;}
	~comp_thread_data();
	void insert();
	void transfer();
	void full_transfer();
	void initialize(uint64_t N, unsigned char* , int, int , int,long, char*);
	};
	
void comp_thread_data::initialize ( uint64_t N, unsigned char* pattern, int ran_bits_size, int km_l, int th_i, long bsize, char* file_name){
	no_of_kmers=N;
	buf_size=bsize;
	//cout<<"initializing.."<<endl;
	pbf=new pattern_bf(no_of_kmers*13, 4, pattern, ran_bits_size);
	buffer= new T[buf_size];
	//km_list=new vector<kmer> (no_of_kmers,0); 
	if (pbf==NULL || buffer==NULL )
		cout <<"unallocated stuff"<<endl;
	kmer_length=km_l;
	kmer_i=0;
	buffer_i=0;
	sorted_i=0;
	thread_i=th_i;
	
	string s=file_name;
	s.append(to_string(thread_i));
	kmer_file=new ofstream(s);
	
	}
comp_thread_data::comp_thread_data ( uint64_t N, unsigned char* pattern,int ran_bits_size, int km_l, int i, long bsize, char* file_name){
	initialize(N , pattern, ran_bits_size,km_l, i, bsize, file_name);
	}
comp_thread_data::~comp_thread_data	(){
	if (pbf!=NULL)delete pbf;
	//cout <<"deleted pbf"<<endl;
	if (buffer!=NULL) delete buffer;
	//cout <<"deleted km_list"<<endl;
	//fclose(kmer_file);
	}

void comp_thread_data::full_transfer(){
	long s_i;
	int q_size;
	while(true){
		
			unique_lock<mutex> lock(consumer_mutices[thread_i]);
			//if (consumer_Q.empty() && producer_alive){
			if (consumer_Q.empty()){
				//cout<<"consumer Q empty:"<<thread_i<< endl;
				if(producer_alive)
					is_not_empty[thread_i].wait(lock);//wait if queue is empty
				else {
					//cout<<"found producer dead..exiting"<<endl;
					is_not_full.notify_all();
					break;
					}
				if (consumer_Q.empty() ){ //this should only happen when the producer is dead
					//cout<<"queue empty after wakeup"<<endl;
					is_not_full.notify_all();
					break;
					}
				else continue;
				}
			//lock.lock();//lock the q for transferring the content	
			/*if(!producer_alive) {
				cout<<"found producer dead..exiting"<<endl;
				is_not_full.notify_all();
				break;
				}*/
			s_i=consumer_Q.front();
			consumer_Q.pop();
			//q_size=consumer_Q.size();
			lock.unlock();
			//if ((sm_buffs[s_i].i+buffer_i)>= no_of_kmers) cout<<"Buffer overflow!!"<<endl;
			for(; sm_buffs[s_i].i>0 && buffer_i< buf_size; ){
				buffer[buffer_i++]=sm_buffs[s_i].buffer[--sm_buffs[s_i].i];
				}
			if(sm_buffs[s_i].i<=0){
				sm_buffs[s_i].i=0;
				}
			else{cout<<"items left out "<<sm_buffs[s_i].i<<endl;}
			//now this has to be returned to the producer queue
			unique_lock<mutex> lock_p(prod_mutex);
			producer_Q.push(s_i);
			lock_p.unlock();
			is_not_full.notify_all();
			//nanosleep(0, NULL); //this_thread::yield()
			insert();
			//if(!producer_alive && !q_size) break;
		}
	//cout<<"seen 2 times "<<t<<endl;
	//cout<<"thread exiting "<<thread_i<<" q size"<<q_size<<endl;
	}


void comp_thread_data::insert(){
	long comp_i, t=0;
	//cout <<"inserting..."<<buffer_i<< "items"<<endl;
	char tmp_rd[1000];
			
	for (long i =0 ; i<buffer_i; i++){
		//uint64_t rc=get_rev_comp(buffer[i],kmer_length);
		
		if (pbf->insert(buffer[i])==1) {//return of 1 means seen for the second time
			//write it to disk
			kmer_i++;
			get_string_read(buffer[i], tmp_rd,kmer_length);
			(*kmer_file) << '>'<<kmer_i<<'\n'<<tmp_rd<< '\n';
			//fwrite(&buffer[i], sizeof(buffer[i]),1, kmer_file);
			}
		}
	buffer_i=0;
	}
	
void * thread_insert(void* obj){
	comp_thread_data *comp=(comp_thread_data *) obj;
	nanosleep(0, NULL); //this_thread::yield()
	comp->full_transfer();
	return NULL;
	}

void usage(){
	cout<<"cTurtle64 Usage:"<<endl;
	cout <<"cTurtle64 [arguments]"<<endl;
	cout<<"example: ./cTurtle64 -f 1Mreads.fq -o kmer_counts -k 51 -n 6000000 -t 9"<<endl;
	cout <<"-i \t input reads file in fasta format."<<endl;
	cout <<"-f \t input reads file in fastq format. This is mutually exclusive with -i."<<endl;
	cout <<"-o \t ouput files prefix. k-mers and their counts are stored in fasta format (headers indicating frequency) in multiple files named prefix0, prefix1... which the user can concatenate if desired."<<endl;
	cout <<"-k \t k-mer length."<<endl;
	cout <<"-t \t Number of threads."<<endl;
	cout <<"-n \t Expected number of frequent k-mers. For uniform coverage libraries this is usually close to genome length. For single-cell libraries, 2-3 times the gemome length is recommended."<<endl;
	cout <<"-h \t Print this help menu."<<endl;
	cout <<"-v \t Print software version."<<endl;
	cout <<""<<endl;
	}	
queue<small_buff*> empty_buff_queue;
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
					strcpy(input_file_name, argv[2]);
					get_read=&getNextRead_fq;
					break;
				case 'o':
					cout <<"ouput prefix\t"<<argv[2];
					strcpy(output_file_name, argv[2]);
					break;
			
				case 't':
					no_of_threads=atoi(argv[2]);
					if (no_of_threads<2) no_of_threads=2;
					cout <<"no of threads \t"<<no_of_threads;
					
					no_of_threads--;//one for the producer
					if (no_of_threads%2==0) no_of_threads--;//even number of workers is not good, so we force it to be odd	
					break;
				case 'k':
					kmer_length=atoi(argv[2]);
					cout <<"k-mer length \t"<<kmer_length;
					if (kmer_length>64){
						cout<<"cTurtle64 can handle k-mers of length up to 64."<<endl;
						return 0;
						}
					if (kmer_length<=32){
						cout<<". Please consider using cTurtle32 for k-mers of length <=32. cTurtle32 is more memory and time efficient for smaller k-mers."<<endl;
						}
					break;
				case 'n':
					gen_length=atol(argv[2]);
					cout <<"Freq. k-mers\t"<<gen_length;
					break;	
				case 's':
					gen_length=(atof(argv[2])*(1<<30))/41;
					cout <<"Freq. k-mers\t"<<gen_length;
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
	if (no_of_threads==0) {
		cout<<"Please specify the number of worker threads."<<endl;
		usage();
		return 0;
		}
	if (gen_length==0) {
		cout<<"Please specify the number of expected frequent kmers."<<endl;
		usage();
		return 0;
		}
	long no_of_kmers=gen_length*2, kmer_i=0;
	
	ifstream reads_file(input_file_name);
	
	if (!reads_file.good()) {
		cout<<"file not found"<<endl;
		return 0;
		}
	//cout<<"generating rands..."<<endl;
	const int ran_bits_size = 8388608;
	//unsigned char random_bits[ran_bits_size];
	unsigned char *random_bits=new unsigned char [ran_bits_size];
    typedef std::minstd_rand G;
    G g;
    typedef std::uniform_int_distribution<> D;
    D d(0, 63);
    for (int i = 0; i < ran_bits_size; i++) {
        random_bits[i]=(unsigned char)d(g);
        
		}
	comp_thread_data *ctds =new comp_thread_data[no_of_threads];
	no_of_smb=no_of_threads*3;
	sm_buffs=new small_buff[no_of_smb];
	long buf_size= 100*1000000;
	long small_buff_size=buf_size/(4*no_of_threads);
	
	for (int j=0; j< no_of_smb; j++){
		sm_buffs[j].initialize(small_buff_size);
		}
	//vector<comp_thread_data> ctds(2);
	//cout<<"point 1"<<endl;
	int *curr_buff=new int[no_of_threads];	
	/*
	char command[1000]="rm -f ";
	strcat(command, output_file_name);
	
	system(command);
	*/
	for (int j=0; j< no_of_threads; j++){
		ctds[j].initialize(no_of_kmers/no_of_threads, random_bits, ran_bits_size, kmer_length,j, buf_size/no_of_threads, output_file_name);
		sm_buffs[j].thread_i=j;
		curr_buff[j]=j;//thread j gets the j-th small buffer at first
		}
	for(int j=no_of_threads; j<no_of_smb;j++)	
		producer_Q.push(j);//queue of empty buffers.
		
	char  rd[MAX_RD_LEN];
	get_read(reads_file,rd );
	if (rd==NULL) {cout<< "null rd pointer"<<endl; return 0;}
	
	int counter=0, len=0;
	int count_ins=0;
	thread *threads=new thread[no_of_threads];
	mutices = new mutex [no_of_threads]; 
	consumer_mutices=new mutex [no_of_threads]; 
	is_not_empty=new condition_variable [no_of_threads]; 
	
	int thread_i=0;
	producer_alive=true;
	int tok_i=0, tok_j;
	
	for (int j=0; j< no_of_threads; j++){
		threads[j]=thread(thread_insert, &ctds[j]);
		}
	long kmer_seen=0, no_of_reads=0, no_of_N=0;
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
				//thread_i=bit_rd%no_of_threads;//which bloom filter it should go to.
				thread_i=((unsigned int)bit_rd)%no_of_threads;//which bloom filter it should go to.
				if(__builtin_expect(sm_buffs[curr_buff[thread_i]].i>= small_buff_size, 0)	){
				//if (sm_buffs[curr_buff[thread_i]].i>= small_buff_size) {
					unique_lock<mutex> lock_c(consumer_mutices[thread_i]);
					ctds[thread_i].consumer_Q.push(curr_buff[thread_i]);
					lock_c.unlock();
					is_not_empty[thread_i].notify_all();
					
					unique_lock<mutex> lock_p(prod_mutex);
					if(producer_Q.empty()) {
						//cout<<"producer found Q empty:"<<thread_i<<endl;
						is_not_full.wait(lock_p);//wait if queue is empty
						//cout<<"producer unblocked "<<prod_block_total<<endl;
						}
					curr_buff[thread_i]=producer_Q.front();
					producer_Q.pop();
					lock_p.unlock();
					
					}
				sm_buffs[curr_buff[thread_i]].buffer[sm_buffs[curr_buff[thread_i]].i++]=bit_rd;
				
				
				__m128i m, m_rc;
				if (bit_seq_i==16){
					m=sse_next_16_bit_rd(rd+tok_i+kmer_length);
					//inlining this function below
					bit_seq_buff=(unsigned char *)&m;
					m_rc=sse_next_16_bit_rc(rd+tok_i+kmer_length);
					bit_seq_buff_rc=(unsigned char *)&m_rc;
					}
				next_bit_read(bit_read, bit_read_rc, rd[tok_i+kmer_length], kmer_length);
				/*
				bit_read<<=2;
				tmp=3;
				tmp<<=(kmer_length*2);
				bit_read&=(~tmp);
				bit_read|=bit_seq_buff[bit_seq_i];
				bit_read_rc>>=2;
				tmp=bit_seq_buff_rc[bit_seq_i];
				tmp<<=(kmer_length*2-8);
				bit_read_rc|=tmp;
				bit_seq_i++;*/
				}
			tok_i+=kmer_length;
			}
		//begin=clock();
		get_read(reads_file,rd );
		//end=clock();
		//read_time+=end-begin;
		//read=get_read(reads_file);
		if (strlen(rd)==0) break;
		//strcpy(rd,read.c_str());
	}
	
	for (int i0=0; i0<no_of_threads; i0++){
		unique_lock<mutex> lock_c(consumer_mutices[i0]);
		ctds[i0].consumer_Q.push(curr_buff[i0]);
		lock_c.unlock();
		is_not_empty[i0].notify_all();
		//cout<<"notifing "<<i0<<endl;
		nanosleep(0, NULL); //this_thread::yield()
		}
	producer_alive=false;
	
	for(int i=0; i<no_of_threads; i++){
		is_not_empty[i].notify_all();
		}
	for (int j=0; j< no_of_threads; j++)
		threads[j].join();
	//for (int j=0; j< no_of_threads; j++)
		//ctds[j].kmer_i=compress_kmers((*ctds[j].km_list), ctds[j].kmer_i);
	//cout<<"STATs:"<<"\nNo of reads:\t "<<no_of_reads<<"\nNo of N's:\t "<<no_of_N<<"\nNo of k-mers:\t "<<kmer_seen<<endl;
	cout<<"STATs:"<<"\nNo of reads:\t "<<no_of_reads<<"\nNo of k-mers:\t "<<kmer_seen<<endl;
	cout<< "no of frequent k-mers found :";
	long t=0;
	for (int j=0; j< no_of_threads; j++)
		{//cout<< ctds[j].kmer_i<< ',';
		t+=ctds[j].kmer_i;
		}
	cout <<t<<endl;
	/*
	command[0]=0;//blank
	strcat(command, "cat ");
	strcat(command, output_file_name);
	strcat(command, "* > ");
	strcat(command, output_file_name);	
	system(command);
	command[0]=0;//blank
	strcat(command, "rm ");
	strcat(command, output_file_name);
	strcat(command, "?*");
	system(command);
	* */
	//return 0;
	delete []ctds;//the com_thread_datas
	delete []threads;
	delete random_bits;
	delete [] mutices;
	delete [] consumer_mutices;
	//delete [] is_not_full;
	delete [] is_not_empty;
	delete [] sm_buffs;
	delete curr_buff;
	//delete km_list;
	return 0;
}
