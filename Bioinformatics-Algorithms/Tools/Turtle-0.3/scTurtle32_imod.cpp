//#define T __uint128_t
#define T uint64_t
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <random>
#include <vector>
#include <template_helpers.hpp>
#include <batch_pbf.hpp>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <time.h>

using namespace std;
bool producer_alive, quake;
int mod_mask;

mutex *mutices; 
condition_variable is_not_full;
condition_variable *is_not_empty;
mutex *consumer_mutices, prod_mutex, out_file_mutex; 
condition_variable is_prodQ_not_empty;
condition_variable *is_consumerQ_not_empty;
int no_of_threads,no_of_smb;//no of small buffers
queue<int> producer_Q;
//ofstream *kmer_file;

struct small_buff{
	//atomic<int> thread_i;
	int thread_i;
	//atomic<bool> full;
	bool full;
	long i;
	T* buffer;
	small_buff(){thread_i=-1;i=0; full=false;};
	~small_buff(){delete buffer;};
	void initialize(long n){buffer=new T[n];};
	} ;


small_buff *sm_buffs;

	
struct comp_thread_data{
	queue<int> consumer_Q;
	uint64_t no_of_kmers;
	uint64_t kmer_i, buffer_i, sorted_i;
	pattern_bf *pbf;
	vector<kmer> *km_list; 
	long bkt_i;
	ofstream *kmer_file;
	int kmer_length, thread_i;
	comp_thread_data( uint64_t N, unsigned char* , int, int, int, char*);
	comp_thread_data(){return ;}
	~comp_thread_data();
	void insert();
	void transfer();
	void full_transfer();
	void initialize(uint64_t N, unsigned char* , int, int , int, char*);
	};
	
void comp_thread_data::initialize ( uint64_t N, unsigned char* pattern, int ran_bits_size, int km_l, int th_i, char* file_name){
	no_of_kmers=N;
	//cout<<"initializing.."<<endl;
	pbf=new pattern_bf(no_of_kmers*7, 4, pattern, ran_bits_size);
	km_list=new vector<kmer> (no_of_kmers,0); 
	if (pbf==NULL || km_list==NULL )
		cout <<"unallocated stuff"<<endl;
	kmer_length=km_l;
	kmer_i=0;
	buffer_i=0;
	sorted_i=0;
	bkt_i=0;
	thread_i=th_i;
	
	string s=file_name;
	s.append(to_string(thread_i));
	kmer_file=new ofstream(s);
	//if (kmer_file==NULL) cout<<"could not open kmer file " << thread_i<<endl;
	/*
	char tmp_fn[200];
	strcpy(tmp_fn, file_name);
	char thread_str[10];
	sprintf(thread_str, "%d", thread_i);
	strcat(tmp_fn, thread_str);
	kmer_file=new ofstream(tmp_fn);
	*/
	}
comp_thread_data::comp_thread_data ( uint64_t N, unsigned char* pattern,int ran_bits_size, int km_l, int i, char* file_name){
	initialize(N , pattern, ran_bits_size,km_l, i, file_name);
	}
comp_thread_data::~comp_thread_data	(){
	if (pbf!=NULL)delete pbf;
	//cout <<"deleted pbf"<<endl;
	if (km_list!=NULL) delete km_list;
	//cout <<"deleted km_list"<<endl;
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
			for(; sm_buffs[s_i].i>0 && buffer_i< no_of_kmers; ){
				(*km_list)[buffer_i++].bit_kmer=sm_buffs[s_i].buffer[--sm_buffs[s_i].i];
				}
			/*if(sm_buffs[s_i].i<=0){
				sm_buffs[s_i].i=0;
				}
			else{cout<<"items left out "<<sm_buffs[s_i].i<<endl;}
			*/
			
			//now this has to be returned to the producer queue
			unique_lock<mutex> lock_p(prod_mutex);
			producer_Q.push(s_i);
			lock_p.unlock();
			is_not_full.notify_all();
			//nanosleep(0, NULL); //this_thread::yield()
			insert();
			if(__builtin_expect(kmer_i>=no_of_kmers*.85, 0)){
				cout<<"The number of frequent k-mers is much more than expected. Please force quit (ctrl+z) and rerun using a higher value for -n."<<endl;
				//is_not_full.notify_all();
				return;
				}
			//if(!producer_alive && !q_size) break;
		}
	//cout<<"seen 2 times "<<t<<endl;
	//cout<<"within thread final compress"<<endl;
	//compress((*km_list), sorted_i, kmer_i);
	kmer_i=compress_kmers(*km_list, kmer_i);
	//all counts needs to be increased by 1 to counter the effect of counting 1 when they were first seen
	for (long i=0; i<kmer_i; i++)
		(*km_list)[i].count++;
	//cout<<"thread exiting "<<thread_i<<" q size"<<q_size<<endl;
	//unique_lock<mutex> lock(out_file_mutex);
	//cout<<"thread "<<thread_i<< " writing "<<kmer_i<<" k-mers"<<endl;
	//cout <<fwrite(km_list, sizeof(kmer),kmer_i, kmer_file);
	
	char tmp_rd[1000];
	for (long i=0; i<kmer_i ; i++){
		get_string_read((*km_list)[i].bit_kmer, tmp_rd,kmer_length);
		if (quake)
			(*kmer_file) << tmp_rd<<'\t'<<(*km_list)[i].count<< '\n';	
		else
			(*kmer_file) << '>'<<(*km_list)[i].count<<'\n'<<tmp_rd<< '\n';	
		}	
	(*kmer_file).close();	
	//lock.unlock();
	}


void comp_thread_data::insert(){
	pbf->insert(*km_list,kmer_i, buffer_i);
	
	for (long i =kmer_i ; i<buffer_i; i++){
		if ((*km_list)[i].count==1){
			(*km_list)[kmer_i].bit_kmer=(*km_list)[i].bit_kmer;
			(*km_list)[kmer_i].count=1;
			kmer_i++;
			}
		}
	if(__builtin_expect(kmer_i>=no_of_kmers*.85, 0))
	//if (kmer_i>=no_of_kmers*.85) 
		compress((*km_list), sorted_i, kmer_i);
	if(__builtin_expect(kmer_i>=no_of_kmers*.85, 0)){
	//if ( kmer_i >=no_of_kmers*.85){
		//cout<<"sort the whole thing"<<endl;
		sorted_i= compress_kmers((*km_list), kmer_i);
		kmer_i=sorted_i;
		}
	//kmer_i= compress_kmers((*km_list), kmer_i);	
	buffer_i=kmer_i;
		
}	
void * thread_insert(void* obj){
	comp_thread_data *comp=(comp_thread_data *) obj;
	nanosleep(0, NULL); //this_thread::yield()
	comp->full_transfer();
	return NULL;
	}
void * thread_compress(void* obj){
	comp_thread_data *comp=(comp_thread_data *) obj;
	nanosleep(0, NULL); //this_thread::yield()
	comp->kmer_i=compress_kmers(*(comp->km_list), comp->kmer_i);
	return NULL;
	}
	
queue<small_buff*> empty_buff_queue;

void usage(){
	cout<<"scTurtle32 Usage:"<<endl;
	cout <<"scTurtle32 [arguments]"<<endl;
	cout<<"example: ./scTurtle32 -f 1Mreads.fq -o kmer_counts -k 31 -n 6000000 -t 9"<<endl;
	cout <<"-i \t input reads file in fasta format."<<endl;
	cout <<"-f \t input reads file in fastq format. This is mutually exclusive with -i."<<endl;
	cout <<"-o \t output file name. k-mers and their counts are stored in fasta format (headers indicating frequency)."<<endl;
	cout <<"-q \t output file name. k-mers and their counts are stored in tab delimited fromat (quake compatible)."<<endl;
	cout <<"-k \t k-mer length."<<endl;
	cout <<"-t \t Number of workers threads. It is recommended that a prime (or at least an odd) number of workers are used and it is less than the number of cores on the machine."<<endl;
	cout <<"-n \t Expected number of frequent k-mers. For uniform coverage libraries this is usually close to genome length. For single-cell libraries, 2-3 times the gemome length is recommended."<<endl;
	cout <<"-s \t The approximate amount of space (in GB) to be used. It is used to indirectly compute -n and is mutually exclusive with -n. When both -n and -s are specified, the one that appears last is used."<<endl;
	cout <<"-h \t Print this help menu."<<endl;
	cout <<""<<endl;
	}
int main(int argc, char *argv[]){
	char input_file_name[1000], output_file_name[1000];
	input_file_name[0]=0;
	output_file_name[0]=0;
	int kmer_length=0 , no_of_threads=0;
	long gen_length=0;
	
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
					
				case 't':
					no_of_threads=atoi(argv[2]);
					if (no_of_threads<1) no_of_threads=1;
					cout <<"worker threads \t"<<no_of_threads;
					break;
				case 'k':
					kmer_length=atoi(argv[2]);
					cout <<"k-mer length \t"<<kmer_length;
					if (kmer_length>32){
						cout<<" ERROR: scTurtle32 can handle k-mers of length up to 32. Please try scTurtle64."<<endl;
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
	//return 0;
	long no_of_kmers=gen_length*2, kmer_i=0;
	mod_mask=1;
	//for (int tmp_threads=no_of_threads; tmp_threads && mod_mask!=no_of_threads; tmp_thread>>=1)
		//mod_mask<<=1
	mod_mask<<= (int)ceil(log(no_of_threads)/log(2));
	mod_mask-=1;
	cout<<"mod mask:"<< mod_mask<<endl;
	
	ifstream reads_file(input_file_name);
	
	if (!reads_file.good()) {
		cout<<"file not found"<<endl;
		return 0;
		}
	//cout<<"generating rands..."<<endl;
	const int ran_bits_size = 8388608L;
	//unsigned char random_bits[ran_bits_size];
	unsigned char *random_bits=new unsigned char [ran_bits_size];
    typedef std::minstd_rand G;
    G g;
    typedef std::uniform_int_distribution <>D;
    D d(0, 63);
    for (int i = 0; i < ran_bits_size; i++) {
        random_bits[i]=(unsigned char)d(g);
        
		}
	//cout<<"finished."<<endl;
	
	comp_thread_data *ctds =new comp_thread_data[no_of_threads];
	no_of_smb=no_of_threads*4;
	sm_buffs=new small_buff[no_of_smb];
	long small_buff_size=no_of_kmers/(8*no_of_threads);
	
	for (int j=0; j< no_of_smb; j++){
		sm_buffs[j].initialize(small_buff_size);
		}
	
	int *curr_buff=new int[no_of_threads];	
	//removing any existing file with name output_file_name
	char command[1000]="rm -f ";
	strcat(command, output_file_name);
	
	system(command);
	
	for (int j=0; j< no_of_threads; j++){
		ctds[j].initialize(no_of_kmers/no_of_threads, random_bits, ran_bits_size, kmer_length,j, output_file_name);
		sm_buffs[j].thread_i=j;
		curr_buff[j]=j;//thread j gets the j-th small buffer at first
		}
	for(int j=no_of_threads; j<no_of_smb;j++)	producer_Q.push(j);//queue of empty buffers.
		
	char  rd[1000];
	//getNextRead(reads_file,rd );
	get_read(reads_file, rd);
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
	//clock_t prod_block_start, prod_block_end, prod_block_total=0;
	//kmer_file=fopen(output_file_name,"wb");
	//kmer_file=new ofstream(output_file_name);
	long kmer_seen=0;
	long counts[20];
	fill_n(counts,20,0);
	while (rd[0]){//while read not empty
		//cout << rd<<endl;
		tok_i=0;
		int read_len=strlen(rd);
		T bit_read, bit_read_rc, bit_rd, tmp;
			
		//if (frag==NULL) cout << "null"<<endl;
		while (tok_i<=(read_len-kmer_length)){
			while(rd[tok_i]=='N') tok_i++;
			bool kmer_found=false;
			while (tok_i<=(read_len-kmer_length) && !kmer_found ){
				for (tok_j=tok_i+1; tok_j< tok_i+kmer_length;tok_j++)
					if(rd[tok_j]=='N') break;
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
				//kmer_seen++;
				//thread_i=bit_rd%no_of_threads;//which bloom filter it should go to.
				//thread_i=((unsigned int)bit_rd)%no_of_threads;//which bloom filter it should go to.
				/*
				thread_i=(bit_rd& mod_mask);
				T tmp_bit_rd=bit_rd;
				while (thread_i >=no_of_threads && tmp_bit_rd){
					tmp_bit_rd>>=4;
					thread_i^=(tmp_bit_rd& mod_mask);
					}
				if (thread_i >=no_of_threads)
					thread_i=bit_rd%no_of_threads;*/
				thread_i=(bit_rd>>4)%no_of_threads;	
				counts[thread_i]++;
				
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
	cout<<"total load :";
	for (int j=0; j< no_of_threads; j++)
		cout<< counts[j]<<',';
	cout<<endl;
	cout<< "no of frequent k-mers found :";
	long t=0;
	
	for (int j=0; j< no_of_threads; j++)
		{
		cout<< ctds[j].kmer_i<< ',';
		t+=ctds[j].kmer_i;
		}
	cout <<t<<endl;
	//concatenating output files and removing intermediate ones
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
