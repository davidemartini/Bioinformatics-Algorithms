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
		int insert(T n);
		int insert(void *, int);
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
	random_bits=pattern;
	ran_bits_size=i;
	/*ran_bits_size=(1<<24);
	random_bits=new unsigned char[ran_bits_size];
	for(int i=0; i<ran_bits_size; i++)
		random_bits[i]=(unsigned char)rand();*/
	bit_array=new unsigned long[size];
	for (long i=0; i<size;i++) bit_array[i]=0;//empty bloom filter
	//cout<<"created the bf "<<endl;
	}
pattern_bf::~pattern_bf(){
	delete bit_array;
	//delete random_bits;
	}


int pattern_bf::insert(T n){
	unsigned long rd,bit_i, pat_i;
	rd=hash_bytes((unsigned char*)&n, sizeof(n));
	bit_i=rd%(size-2*k);//the start of the bit_array where it will set bits
	pat_i=rd%(ran_bits_size-2*k);//start of the random bit pattern
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
	else{
		//it has been seen at least once we want to know if it has been seen twice
		bool present2=true;
		mask=63,t;
		
		for(int i=k ; i<2*k; i++, bit_i++, pat_i++){
			t1=1;
			t=random_bits[pat_i]& mask;//which bit in the bit array is set for this element
			t1<<=t;
			if (!(t1& bit_array[bit_i])){//this bit is not set
				present2= false;
				bit_array[bit_i]|=t1;//set the bit
				}	
			}
		if (!present2) 	return 1;
		else return 2;
		} 
	}

