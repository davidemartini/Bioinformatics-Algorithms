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

#include <iostream>
#include <x86intrin.h>
#define MAX_RD_LEN 15000
using namespace std;

struct kmer{
	uint16_t count;
	T bit_kmer;
	
	kmer(int) ;
	kmer();
	friend bool operator< ( kmer one, kmer other);
	void operator= ( kmer other);
	void operator+= ( kmer other);
	} __attribute__ ((packed));

kmer::kmer(int x){return;}
kmer::kmer(){return;}

bool operator< (kmer one, kmer other){
	return (one.bit_kmer < other.bit_kmer);
	}

void kmer::operator= (kmer other){
	bit_kmer= other.bit_kmer;
	count= other.count;
	
	}

void kmer::operator+= (kmer other){
	count+= other.count;
	
	}

void getNextRead(ifstream & fp, char* read){
	//in this version, the read pointer is passed and the read is returned in that char array
	char line[MAX_RD_LEN];
	read[0]=0;//read should be empty
	while (!fp.eof()) {
		fp.getline(line,MAX_RD_LEN);
		if (strlen(line)>0 && line[0]=='>') fp.getline(line,MAX_RD_LEN);
		
		while (strlen(line)>0 && line[0]!='>') {
			strcat(read,line);
			fp.getline(line,MAX_RD_LEN);
			}
		if(strlen(read)) return;
		}
	}
void getNextRead_fq(ifstream & fp, char* read){
	char line[MAX_RD_LEN];
	read[0]=0;//read should be empty
	//char buf[1000];
/*@1_1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+1_1
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@1_2*/

	while (!fp.eof()) {
		fp.getline(line,MAX_RD_LEN);
		//cout<<"line:"<<line<<endl;
		if (line[0]=='@' && strlen(line)>0 ) fp.getline(line,MAX_RD_LEN);
		
		
		while (line[0]!='+' && strlen(line)>0) {
			strcat(read,line);
			
			fp.getline(line,MAX_RD_LEN);
			}
		while (line[0]!='@' && strlen(line)>0) {
			fp.getline(line,MAX_RD_LEN);
			}
		//cout<<read<<endl;
		if(strlen(read)) return;
		}
	}
string getNextRead(ifstream & fp){
	string line, read;
	//ifstream fp1=*fp;
	char buf[100];
	while (fp.good()) {
		getline(fp,line);
		if (!line.empty() && line.at(0)=='>') getline(fp,line);
		
		while (!line.empty() && line.at(0)!='>') {
			read.append(line);
		
			getline(fp,line);
			}
		return read;
		}
	}
T get_rev_comp(T l, int kmer_length){
	T comp_l=~l, rev_comp=0,temp=3;
	for (int i=0; i<kmer_length; i++){
		temp=3;
		temp&=comp_l;
		rev_comp<<=2;
		comp_l>>=2;
		rev_comp|=temp;
		
		}
	return rev_comp;
	}

T get_bit_read_128(char* s, int len){
	T l=0;
	//a=1, c=0, g=3 , t=2
	for (int i=0; i<len;i++){
		l<<=2;
		switch(s[i]){
			case 'A':{
				l|=1;
				break;
				}
				
			case 'C':{
				l|=0;
				break;
				}
			
			case 'G':{
				l|=3;
				break;
				}
			case 'T':{
				l|=2;
				break;
				}
			
			}
		
		}
	return l;
	}
unsigned int inline sse_get_bit_read( char * s, int len){
	//__m128i m0=_mm_set_epi8(s[15],  s[14],  s[13],  s[12],  s[11],  s[10],  s[9],  s[8],  s[7],  s[6],  s[5],  s[4],  s[3],  s[2],  s[1],  s[0]);
	__m128i m0=_mm_loadu_si128((const __m128i *)s);
	__m128i I_a, I_c, I_g, I_t;
	/*
	I_c=_mm_set_epi8(0,  0,  0, 0,0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	I_a=_mm_set_epi8(1,  1,  1, 1,1 , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
	I_g=_mm_set_epi8(3,  3,  3, 3,3 , 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3);
	I_t=_mm_set_epi8(2,  2,  2, 2,2 , 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2);
	*/
	I_c=_mm_setzero_si128();
	I_a=_mm_set1_epi8(1);
	I_g=_mm_set1_epi8(3);
	I_t=_mm_set1_epi8(2);
	/*
	__m128i comp_a=	_mm_cmpeq_epi8(m0, _mm_set_epi8('A',  'A',  'A', 'A','A' , 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'));
	__m128i comp_c=	_mm_cmpeq_epi8(m0,  _mm_set_epi8('C',  'C',  'C', 'C','C' , 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'));
	__m128i comp_g=	_mm_cmpeq_epi8(m0, _mm_set_epi8('G',  'G',  'G', 'G','G' , 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G'));
	__m128i comp_t=	_mm_cmpeq_epi8(m0,  _mm_set_epi8('T',  'T',  'T', 'T','T' , 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T'));
	*/
	__m128i comp_a=	_mm_cmpeq_epi8(m0, _mm_set1_epi8('A'));
	__m128i comp_c=	_mm_cmpeq_epi8(m0, _mm_set1_epi8('C'));
	__m128i comp_g=	_mm_cmpeq_epi8(m0, _mm_set1_epi8('G'));
	__m128i comp_t=	_mm_cmpeq_epi8(m0, _mm_set1_epi8('T'));
	
	comp_a= _mm_and_si128(I_a, comp_a);
	comp_g= _mm_and_si128(I_g, comp_g);
	comp_t= _mm_and_si128(I_t, comp_t);
	
	m0=_mm_or_si128(comp_a, comp_g);
	m0=_mm_or_si128(m0, comp_t);
	
	unsigned char *ptr= (unsigned char*) &m0;
	unsigned int l=0;
	for (int i=0; (i<len);i++){
		l<<=2;
		l|=ptr[i];
		}
	return l;
	}
__m128i inline sse_next_16_bit_rd( char * s){
	//__m128i m0=_mm_set_epi8(s[0],  s[1],  s[2],  s[3],  s[4],  s[5],  s[6],  s[7],  s[8],  s[9],  s[10],  s[11],  s[12],  s[13],  s[14],  s[15]);//this is places in the reverse order so that while reading them out of the __m128i, we can just progress linearly
	//__m128i m0=_mm_set_epi8(s[15],  s[14],  s[13],  s[12],  s[11],  s[10],  s[9],  s[8],  s[7],  s[6],  s[5],  s[4],  s[3],  s[2],  s[1],  s[0]);
	__m128i m0=_mm_loadu_si128((const __m128i *)s);
	__m128i I_a, I_c, I_g, I_t;
	/*
	I_c=_mm_set_epi8(0,  0,  0, 0,0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	I_a=_mm_set_epi8(1,  1,  1, 1,1 , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
	I_g=_mm_set_epi8(3,  3,  3, 3,3 , 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3);
	I_t=_mm_set_epi8(2,  2,  2, 2,2 , 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2);
	*/
	I_c=_mm_setzero_si128();;
	I_a=_mm_set1_epi8(1);
	I_g=_mm_set1_epi8(3);
	I_t=_mm_set1_epi8(2);
	/*
	__m128i comp_a=	_mm_cmpeq_epi8(m0, _mm_set_epi8('A',  'A',  'A', 'A','A' , 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'));
	__m128i comp_c=	_mm_cmpeq_epi8(m0,  _mm_set_epi8('C',  'C',  'C', 'C','C' , 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'));
	__m128i comp_g=	_mm_cmpeq_epi8(m0, _mm_set_epi8('G',  'G',  'G', 'G','G' , 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G'));
	__m128i comp_t=	_mm_cmpeq_epi8(m0,  _mm_set_epi8('T',  'T',  'T', 'T','T' , 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T'));
	*/
	__m128i comp_a=	_mm_cmpeq_epi8(m0, _mm_set1_epi8('A'));
	__m128i comp_c=	_mm_cmpeq_epi8(m0, _mm_set1_epi8('C'));
	__m128i comp_g=	_mm_cmpeq_epi8(m0, _mm_set1_epi8('G'));
	__m128i comp_t=	_mm_cmpeq_epi8(m0, _mm_set1_epi8('T'));
	
	comp_a= _mm_and_si128(I_a, comp_a);
	comp_g= _mm_and_si128(I_g, comp_g);
	comp_t= _mm_and_si128(I_t, comp_t);
	
	m0=_mm_or_si128(comp_a, comp_g);
	m0=_mm_or_si128(m0, comp_t);
	
	return m0;
	}	
__m128i inline sse_next_16_bit_rc( char * s){
	//__m128i m0=_mm_set_epi8(s[15],  s[14],  s[13],  s[12],  s[11],  s[10],  s[9],  s[8],  s[7],  s[6],  s[5],  s[4],  s[3],  s[2],  s[1],  s[0]);
	__m128i m0=_mm_loadu_si128((const __m128i *)s);
	//__m128i m0=_mm_set_epi8(s[0],  s[1],  s[2],  s[3],  s[4],  s[5],  s[6],  s[7],  s[8],  s[9],  s[10],  s[11],  s[12],  s[13],  s[14],  s[15]);//this is places in the reverse order so that while reading them out of the __m128i, we can just progress linearly
	__m128i I_a, I_c, I_g, I_t;
	I_c=_mm_set1_epi8(0xC0);
	I_a=_mm_set1_epi8(0x80);
	I_g=_mm_setzero_si128();
	I_t=_mm_set1_epi8(0x40);
	
	__m128i comp_a=	_mm_cmpeq_epi8(m0, _mm_set1_epi8('A'));
	__m128i comp_c=	_mm_cmpeq_epi8(m0, _mm_set1_epi8('C'));
	__m128i comp_g=	_mm_cmpeq_epi8(m0, _mm_set1_epi8('G'));
	__m128i comp_t=	_mm_cmpeq_epi8(m0, _mm_set1_epi8('T'));
	
	comp_a= _mm_and_si128(I_a, comp_a);
	comp_c= _mm_and_si128(I_c, comp_c);
	comp_t= _mm_and_si128(I_t, comp_t);
	
	m0=_mm_or_si128(comp_a, comp_c);
	m0=_mm_or_si128(m0, comp_t);
	
	return m0;
	}	
void next_bit_read(T &bit_read,T &bit_read_rc,char c,int kmer_length){
	T tmp=1;
	
	T tmp1=1,tmp2=3;
	tmp<<=(kmer_length*2-1);
	tmp1<<=(kmer_length*2-2);
	tmp2<<=(kmer_length*2);
	//tmp<<=63;
	//tmp1<<=62;
	
	bit_read<<=2;
	bit_read&=(~tmp2);
	bit_read_rc>>=2;
		//a=1, c=0, g=3 , t=2
	//cout <<frag[i+k1mer_length]<<endl;
	switch(c){
		case 'A':{
			bit_read|=1;
			//bit_read_rc|=(1<<63);//9223372036854775808;
			bit_read_rc|=tmp;
			break;
			//13835058055282163712;
			}
		case 'T':{
			bit_read|=2;
			//bit_read_rc|=(1<<62);//4611686018427387904;
			bit_read_rc|=tmp1;
			break;
			//13835058055282163712;
			}
		case 'C':{
			//bit_read|=3; c is 0
			bit_read_rc|=tmp;
			bit_read_rc|=tmp1;
			
			//bit_read_rc|=(1<<63|1<<62);//13835058055282163712;
			break;
			
			}
		case 'G':{
			bit_read|=3; 
			//bit_read_rc|=13835058055282163712;
			break;
			
			}
		default:
			break;
		}
	}

char * get_string_read(T n, char * read,int len){
	int tmp,i;
	//char *read=new char[len+1];
		//a=1, c=0, g=3 , t=2
	for (i=0 ; i<len; i++,  n>>=2){
		tmp=n&3;
		switch(tmp){
			case 1:{
				read[len-i-1]='A';
				break;
				}
			case 2:{
				read[len-i-1]='T';
				break;
				}	
			case 3:{
				read[len-i-1]='G';
				break;
				}
			case 0:{
				read[len-i-1]='C';
				break;
				}
			}
		}
	read[i]=0; // null terminator for strings
	return read;
	}	
void print_binary(T n, int kmer_length) {
	char str[200];
    for (int i=0 ; i<kmer_length*2; i++, n>>=1){
		if (n&1) str[kmer_length*2-1-i]='1';
		else str[kmer_length*2-1-i]='0';
		}
	str[kmer_length*2]=0;
	cout<<str<<endl;
    
}
long compress_kmers(vector<kmer> &kmers, long length){
	
	sort (kmers.begin(), kmers.begin()+length);
	
	long i=0,j,k;
	
	int count=0;
	
	for(j=0, k=1; j<length-1;){
		count=kmers[j].count;
		k=j+1;
		while(kmers[j].bit_kmer==kmers[k].bit_kmer){
			count+=kmers[k].count;
			k++;
			}
		kmers[i].bit_kmer=kmers[j].bit_kmer;
		kmers[i].count=count;
		i++;
		j=k;
		}
	if (j==length-1)//the one item at the end was unique
		kmers[i++]=kmers[j];
	
	return i;
	
	}
	
long compress(vector<kmer> &kmers, uint64_t &sorted_i, uint64_t &kmer_i){
	if(__builtin_expect(sorted_i==0, 0)){
		sorted_i=compress_kmers(kmers, kmer_i);
		kmer_i=sorted_i;
	
		return kmer_i;
		}
	
	// sorte the unsorted part
	sort (kmers.begin()+sorted_i, kmers.begin()+kmer_i);
	
	long i=sorted_i;
	
	//if therea are elements in the new part thats presented in the sorted part,
	// we want to merge them.
	i=0;
	long new_i=sorted_i;
	for(long j=sorted_i, k; (j<kmer_i-1)&& (i< sorted_i);){
		//find the element thats >= kmers[j]
		k=j;
		while(kmers[i].bit_kmer <kmers[j].bit_kmer && (i< sorted_i)) i++;
		if (i>= sorted_i) {
			//we need to compact the rest of the newly sorted part
			//cout <<"reached end of sorted array"<<endl;
			while(j<kmer_i-1){
				k=j+1;
				while((kmers[j].bit_kmer==kmers[k].bit_kmer) && (k<kmer_i)){
					kmers[j].count+=kmers[k].count;
					k++;
					}
				kmers[new_i].bit_kmer=kmers[j].bit_kmer;
				kmers[new_i].count=kmers[j].count;
				new_i++;
				j=k;
				}
			//cout<<"compacted newly sorted part to "<<(new_i-sorted_i) <<" elements"<<endl;
			//cout<< sorted_i << " " << kmer_i<<"," <<new_i<<endl;
			sorted_i=kmer_i=new_i;
			
			return new_i;
			}
		while ((kmers[i].bit_kmer == kmers[k].bit_kmer) && (k<kmer_i)) 
			kmers[i].count+=kmers[k++].count;
		if (j==k){//the element was not found in the sorted part;
			k=j+1;
			while((kmers[j].bit_kmer==kmers[k].bit_kmer) && (k<kmer_i)){
				kmers[j].count+=kmers[k].count;
				k++;
				}
			kmers[new_i].bit_kmer=kmers[j].bit_kmer;
			kmers[new_i].count=kmers[j].count;
			new_i++;
			j=k;
			}
		else {j=k;}
		}
	//cout<<"compacted newly sorted part to "<<(new_i-sorted_i) <<" elements"<<endl;
	//cout<< sorted_i << " " << kmer_i<<"," <<new_i<<endl;
	
	sorted_i=kmer_i=new_i;
	return new_i;
	}
int remu63(unsigned int n) { 
	n = (0x04104104*n + (n >> 4) + (n >> 10)) >> 26; 
	return n&((n-63)>>6); //Change63to0.
}
/*
int remu31(unsigned int n) { 
	n = (0x04210842*n + (n >> 4) + (n >> 9)) >> 27; 
	//return n&((n-63)>>6); //Change63to0.
	if (n==0) return 16;
	if (n==31) return 0;
	else return n;
}
*/
int remu31(unsigned int n) { 
	n = (0x04210842*n + (n >> 4) + (n >> 9)) >> 27; 
	//return n&((n-63)>>6); //Change63to0.
	//if ( n) return n&((int)(n-31)>>6);
	if(__builtin_expect(n, 1)) return n&((int)(n-31)>>6);
	else return 16;
	//if (n==31) return 0;
	//else return n;
}
int remu32(unsigned int n) { 
	n = (0x04210842*n + (n >> 4) + (n >> 9)) >> 27; 
	return n;
}
int remu32_0(unsigned int n) { 
	n = (0x04104104*n + (n >> 4) + (n >> 10)) >> 26; 
	return n>>1;
}
