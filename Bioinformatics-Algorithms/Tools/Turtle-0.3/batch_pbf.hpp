//#include <boost/dynamic_bitset.hpp>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstdlib>
using namespace std;
class pattern_bf{
	unsigned long n,  ran_bits_size, *bit_array,size; // expected no of items
	int k;//optimal number of hash functions
	unsigned char * random_bits;
	void get_indices(unsigned long n, unsigned long rc);
	unsigned long hash64shift(unsigned long);
	unsigned long hash128shift(__uint128_t key);
	public:
		pattern_bf( unsigned long N, int K, unsigned char*, int);
		~pattern_bf();
		int insert(unsigned long n, unsigned long rc);
		int insert(unsigned long n);
		int insert(__uint128_t n);
		int insert(void *, int);
		void insert(vector<kmer> &km_list, long start, long end);
		void insert_only(vector<kmer> &km_list, long start, long end);
		void insert_only(T *km_list, long start, long end);
		void check(vector<kmer> &km_list, long start, long end);
		bool check(T bit_kmer);
		void clear(){for (long i=0; i<size;i++) bit_array[i]=0;};
		//int insert (kmer&);
		
	};

unsigned long pattern_bf:: hash64shift(unsigned long key)
{
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}

unsigned long pattern_bf:: hash128shift(__uint128_t key)
{
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return (unsigned long)key;
}

unsigned long hash_bytes(unsigned char *str ,int length)
    {
        unsigned long hash = 5381;
        int c;

        for(int i=0; i<length; i++){
			c = str[i];
			//str++;
            hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

			} 
        return hash;
    }

pattern_bf:: pattern_bf(unsigned long N, int K,unsigned char * pattern, int i){
	//m=M+(8-M%8);//we want m to be an exact multiple of 8
	n=N;
	k=K;
	size=n/8;
	ran_bits_size=i;
	random_bits=pattern;
	
	bit_array=new unsigned long[size+k];
	for (long i=0; i<size;i++) bit_array[i]=0;//empty bloom filter
	//cout<<"created the bf "<<endl;
	}
pattern_bf::~pattern_bf(){
	delete bit_array;
	//delete random_bits;
	}

int pattern_bf::insert(unsigned long n, unsigned long rc){
	unsigned long rd=hash64shift(n+rc),bit_i, pat_i;
	bit_i=rd%(size-k);//the start of the bit_array where it will set bits
	pat_i=rd%(ran_bits_size-k);//start of the random bit pattern
	//bit_i=rd&(size-1);
	//pat_i=rd&(ran_bits_size-1);
	bool present=true;
	unsigned int mask=63,t;
	unsigned long t1;
	for(int i=0 ; i<k; i++, bit_i++, pat_i++){
		t1=1;
		t=random_bits[pat_i]& mask;//which bit in the bit array is set for this element
		t1<<=t;
		if (!(t1& bit_array[bit_i])){//this bit is not set
			present= false;
			bit_array[bit_i]|=t1;//set the bit
			}	
		}
	if(!present) return 0;
	else return 1;
	}

int pattern_bf::insert(unsigned long n){
	unsigned long rd=hash64shift(n),bit_i, pat_i;
	//rd=hash_bytes((unsigned char*)&n, sizeof(n));
	//srand(rd);
	//bit_i=rand()%(size-k);
	bit_i=rd%(size-k);//the start of the bit_array where it will set bits
	pat_i=rd%(ran_bits_size-k);//start of the random bit pattern
	bool present=true;
	unsigned int mask=63,t;
	unsigned long t1;
	for(int i=0 ; i<k; i++, bit_i++, pat_i++){
		t1=1;
		t=random_bits[pat_i]& mask;//which bit in the bit array is set for this element
		t1<<=t;
		if (!(t1& bit_array[bit_i])){//this bit is not set
			present= false;
			bit_array[bit_i]|=t1;//set the bit
			}	
		}
	if(!present) return 0;
	else return 1;
	}
int pattern_bf::insert(__uint128_t n){
	unsigned long rd,bit_i, pat_i;
	//rd=hash128shift(n);
	rd=hash_bytes((unsigned char*)&n, sizeof(n));
	//unsigned long rd=hash64shift(n),bit_i, pat_i;
	//srand(rd);
	//bit_i=rand()%(size-k);
	bit_i=rd%(size-k);//the start of the bit_array where it will set bits
	pat_i=rd%(ran_bits_size-k);//start of the random bit pattern
	bool present=true;
	unsigned int mask=63,t;
	unsigned long t1;
	for(int i=0 ; i<k; i++, bit_i++, pat_i++){
		t1=1;
		t=random_bits[pat_i]& mask;//which bit in the bit array is set for this element
		t1<<=t;
		if (!(t1& bit_array[bit_i])){//this bit is not set
			present= false;
			bit_array[bit_i]|=t1;//set the bit
			}	
		}
	if(!present) return 0;
	else return 1;
	}
struct kmer_count {
  bool operator() (kmer i,kmer j) { return (i.count<j.count);}
} kmer_count_obj;

void pattern_bf::insert(vector<kmer> &km_list, long start, long end){
	unsigned long rd,bit_i, pat_i;
	int x=0;
		
	for(long j =start; j<end; j++){
		//bit_i=km_list[j].count;//the start of the bit_array where it will set bits
		//pat_i=km_list[j].count%(ran_bits_size-k);//start of the random bit pattern
		rd=hash_bytes((unsigned char*)&(km_list[j].bit_kmer), sizeof(km_list[j].bit_kmer));
		bit_i=rd%(size-k);
		pat_i=rd%(ran_bits_size-k);
		bool present;
		present=true;
		unsigned int mask=63,t;
		unsigned long t1;
		for(int i=0 ; i<k; i++, bit_i++, pat_i++){
			t1=1;
			t=random_bits[pat_i]& mask;//which bit in the bit array is set for this element
			t1<<=t;
			if (!(t1& bit_array[bit_i])){//this bit is not set
				present= false;
				bit_array[bit_i]|=t1;//set the bit
				}	
			}
		if(!present) km_list[j].count=0;
		else {km_list[j].count=1;x++;}
		}
	//cout<< "found duplicates "<< x<<endl;
	}	
void pattern_bf::insert_only(vector<kmer> &km_list, long start, long end){
	unsigned long rd,bit_i, pat_i;
	
	int x=0;
		
	for(long j =start; j<end; j++){
		//bit_i=km_list[j].count;//the start of the bit_array where it will set bits
		//pat_i=km_list[j].count%(ran_bits_size-k);//start of the random bit pattern
		rd=hash_bytes((unsigned char*)&(km_list[j].bit_kmer), sizeof(km_list[j].bit_kmer));
		bit_i=rd%(size-k);
		pat_i=rd%(ran_bits_size-k);
		bool present;
		present=true;
		unsigned int mask=63,t;
		unsigned long t1;
		for(int i=0 ; i<k; i++, bit_i++, pat_i++){
			t1=1;
			t=random_bits[pat_i]& mask;//which bit in the bit array is set for this element
			t1<<=t;
			if (!(t1& bit_array[bit_i])){//this bit is not set
				present= false;
				bit_array[bit_i]|=t1;//set the bit
				}	
			}
		}
	//cout<< "found duplicates "<< x<<endl;
	}	
void pattern_bf::insert_only(T *km_list, long start, long end){
	unsigned long rd,bit_i, pat_i;
	
	int x=0;
		
	for(long j =start; j<end; j++){
		//bit_i=km_list[j].count;//the start of the bit_array where it will set bits
		//pat_i=km_list[j].count%(ran_bits_size-k);//start of the random bit pattern
		rd=hash_bytes((unsigned char*)(km_list+j), sizeof(T));
		bit_i=rd%(size-k);
		pat_i=rd%(ran_bits_size-k);
		bool present;
		present=true;
		unsigned int mask=63,t;
		unsigned long t1;
		for(int i=0 ; i<k; i++, bit_i++, pat_i++){
			t1=1;
			t=random_bits[pat_i]& mask;//which bit in the bit array is set for this element
			t1<<=t;
			if (!(t1& bit_array[bit_i])){//this bit is not set
				//present= false;
				bit_array[bit_i]|=t1;//set the bit
				}	
			}
		}
	}	

void pattern_bf::check(vector<kmer> &km_list, long start, long end){
	unsigned long rd,bit_i, pat_i;
	
	int x=0;
		
	for(long j =start; j<end; j++){
		//bit_i=km_list[j].count;//the start of the bit_array where it will set bits
		//pat_i=km_list[j].count%(ran_bits_size-k);//start of the random bit pattern
		rd=hash_bytes((unsigned char*)&(km_list[j].bit_kmer), sizeof(km_list[j].bit_kmer));
		bit_i=rd%(size-k);
		pat_i=rd%(ran_bits_size-k);
		bool present;
		present=true;
		unsigned int mask=63,t;
		unsigned long t1;
		for(int i=0 ; i<k; i++, bit_i++, pat_i++){
			t1=1;
			t=random_bits[pat_i]& mask;//which bit in the bit array is set for this element
			t1<<=t;
			if (!(t1& bit_array[bit_i])){//this bit is not set
				present= false;
				break;
				}	
			}
		if(!present) km_list[j].count=0;
		else {km_list[j].count=1;x++;}
		}
	//cout<< "found duplicates "<< x<<endl;
	}	
bool pattern_bf::check(T bit_kmer){
	unsigned long rd,bit_i, pat_i;
	
	int x=0;
		
	
		//bit_i=km_list[j].count;//the start of the bit_array where it will set bits
		//pat_i=km_list[j].count%(ran_bits_size-k);//start of the random bit pattern
	rd=hash_bytes((unsigned char*)&bit_kmer, sizeof(bit_kmer));
	bit_i=rd%(size-k);
	pat_i=rd%(ran_bits_size-k);
	bool present;
	present=true;
	unsigned int mask=63,t;
	unsigned long t1;
	for(int i=0 ; i<k; i++, bit_i++, pat_i++){
		t1=1;
		t=random_bits[pat_i]& mask;//which bit in the bit array is set for this element
		t1<<=t;
		if (!(t1& bit_array[bit_i])){//this bit is not set
			present= false;
			break;
			}	
		}
	if(!present) return false;
	else return true;
	
	//cout<< "found duplicates "<< x<<endl;
	}	

int pattern_bf::insert(void *obj, int length){
	
	unsigned long rd=hash_bytes((unsigned char*)obj, length),bit_i, pat_i;
	//srand(rd);
	//bit_i=rand()%(size-k);
	bit_i=rd%(size-k);//the start of the bit_array where it will set bits
	pat_i=rd%(ran_bits_size-k);//start of the random bit pattern
	bool present=true;
	unsigned int mask=63,t;
	unsigned long t1;
	for(int i=0 ; i<k; i++, bit_i++, pat_i++){
		t1=1;
		t=random_bits[pat_i]& mask;//which bit in the bit array is set for this element
		t1<<=t;
		if (!(t1& bit_array[bit_i])){//this bit is not set
			present= false;
			bit_array[bit_i]|=t1;//set the bit
			}	
		}
	if(!present) return 0;
	else return 1;
	}
