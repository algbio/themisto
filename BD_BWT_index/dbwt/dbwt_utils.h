/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef dbwt_utils_h
#define dbwt_utils_h
#ifndef uchar
typedef unsigned char uchar;
#endif

#ifndef uint
typedef unsigned int uint;
#endif

#ifndef ulong
typedef unsigned long ulong;
#endif

#ifndef min
#define min(x,y) (((x)<(y))?(x):(y))
#endif
#ifndef max
#define max(x,y) (((x)>(y))?(x):(y))
#endif

#ifndef byte
typedef unsigned char byte;
#endif
#ifndef word
typedef unsigned short word;
#endif
#ifndef dword
typedef unsigned int dword;
#endif

#if 1
typedef unsigned int qword;
typedef word pb;
#define logD 4
#else
typedef unsigned long long qword;
typedef dword pb;
#define logD 5
#endif
#define PBS (sizeof(pb)*8)
#define D (1<<logD)

void * dbwt_mymalloc(size_t n);
void * dbwt_myrealloc(void *ptr, size_t new, size_t old);
void  dbwt_myfree(void *p, size_t s);
void dbwt_report_mem(char *s);
extern size_t dbwt_cur_alloc, dbwt_max_alloc;

unsigned int dbwt_getbits(unsigned short *B, unsigned long i, int d);

typedef struct {
  ulong n;
  int w;
  pb *b;
} packed_array;

int dbwt_blog(ulong x);
packed_array * dbwt_allocate_packed_array(ulong n, int w);
void dbwt_free_packed_array(packed_array *p);
ulong dbwt_pa_get(packed_array *p, ulong i);
void dbwt_pa_set(packed_array *p, ulong i, long x);

pb * dbwt_allocate_vector(ulong n);
int dbwt_getbit(pb *B, ulong i);
int dbwt_setbit(pb *B, ulong i,int x);
dword dbwt_getbits(pb *B, ulong i, int d);
int dbwt_setbits(pb *B, ulong i, int d, ulong x);
#endif
