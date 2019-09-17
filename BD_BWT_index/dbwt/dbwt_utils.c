/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include "dbwt_utils.h"

int dbwt_blog(ulong x)
{
int l;
  l = -1;
  while (x>0) {
    x>>=1;
    l++;
  }
  return l;
}

size_t dbwt_cur_alloc=0, dbwt_max_alloc=0;
void * dbwt_mymalloc(size_t n)
{
  void *p;

  p = malloc(n);
  if (p == NULL) {
    printf("malloc failed.\n");
    exit(1);
  }
  dbwt_cur_alloc += n;
  if (dbwt_cur_alloc > dbwt_max_alloc) {
    dbwt_max_alloc = dbwt_cur_alloc;
    //printf("allocated %ld\n",max_alloc);
  }

  if (n == 0) {
    printf("warning: 0 bytes allocated p=%p\n",p);
  }
//  printf("alloc_pointer is %p with size %lu\n",p,(long unsigned int )n);
  return p;
}

void *dbwt_myrealloc(void *ptr, size_t new, size_t old)
{
  void *p;

  p = realloc(ptr, new);
  if (new > 0 && p == NULL) {
    printf("realloc failed. ptr=%p new=%d old=%d\n",ptr, (int) new, (int) old);
    exit(1);
  }
  dbwt_cur_alloc += new - old;
  if (dbwt_cur_alloc > dbwt_max_alloc) {
    dbwt_max_alloc = dbwt_cur_alloc;
    //printf("allocated %ld\n",max_alloc);
  }
//  printf("free alloc_pointer %p with size %lu\n",ptr,(long unsigned int )old);
//  printf("alloc_pointer is %p with size %lu\n",p,(long unsigned int )new);
  return p;
}

void dbwt_report_mem(char *s)
{
  puts(s);
  printf("allocated total %lu   max %lu\n", 
	(unsigned long) dbwt_cur_alloc,(unsigned long) dbwt_max_alloc);
}


void dbwt_myfree(void *p, size_t s)
{
  free(p);
  dbwt_cur_alloc -= s;
//  printf("free alloc_pointer %p with size %lu\n",p,(long unsigned int )s);

}

int dbwt_setbit(pb *B, ulong i,int x)
{
  ulong j,l;

  j = i / D;
  l = i % D;
  if (x==0) B[j] &= (~(1<<(D-1-l)));
  else if (x==1) B[j] |= (1<<(D-1-l));
  else {
    printf("error setbit x=%d\n",x);
    exit(1);
  }
  return x;
}

int dbwt_setbits0(pb *B, ulong i, int d, ulong x)
{
  ulong j;

  for (j=0; j<d; j++) {
    dbwt_setbit(B,i+j,(x>>(d-j-1))&1);
  }
  return x;
}

int dbwt_getbit(pb *B, ulong i)
{
  ulong j,l;

  //j = i / D;
  //l = i % D;
  j = i >> logD;
  l = i & (D-1);
  return (B[j] >> (D-1-l)) & 1;
}

#if 1
dword dbwt_getbits(pb *B, ulong i, int d)
{
  qword x,z;

  B += (i >> logD);
  i &= (D-1);
  if (i+d <= 2*D) {
    x = (((qword)B[0]) << D) + B[1];
    x <<= i;
    x >>= (D*2-1-d);
    x >>= 1;
  } else {
    x = (((qword)B[0])<<D)+B[1];
    z = (x<<D)+B[2];
    x <<= i;
    x &= ((1L<<D)-1)<<D;
    z <<= i;
    z >>= D;
    x += z;
    x >>= (2*D-d);
  }
  
  return x;
}
#else
dword dbwt_getbits(pb *B, ulong i, int d)
{
  dword j,x;

  x = 0;
  for (j=0; j<d; j++) {
    x <<= 1;
    x += dbwt_getbit(B,i+j);
  }
  return x;
}
#endif

int dbwt_setbits(pb *B, ulong i, int d, ulong x)
{
  ulong y,m;
  int d2;


  //  BB = B;  ii = i;  dd = d;  xx = x;

  B += (i>>logD);
  i %= D;

  while (i+d > D) {
    d2 = D-i; // x の上位 d2 ビットを格納
    y = x >> (d-d2);
    m = (1<<d2)-1;
    *B = (*B & (~m)) | y;
    B++;  i=0;
    d -= d2;
    x &= (1<<d)-1; // x の上位ビットを消去
  }
  m = (1<<d)-1;
  y = x << (D-i-d);
  m <<= (D-i-d);
  *B = (*B & (~m)) | y;

#if 0
  if (getbits(BB,ii,dd) != xx) {
    printf("setbits ??? x=%ld %ld\n",xx,getbits(BB,ii,dd));
  }
#endif
  return x;
}

pb * dbwt_allocate_vector(ulong n)
{
  ulong i,x;
  pb *b;

  x = (n+PBS-1)/PBS;
  b = dbwt_mymalloc(x*sizeof(pb));
  for (i=0; i<x; i++) b[i] = 0;
  return b;
}

packed_array * dbwt_allocate_packed_array(ulong n, int w)
{
  ulong i,x;
  packed_array *p;

  if (w >= 32) {
    printf("warning: w=%d\n",w);
  }

  p = dbwt_mymalloc(sizeof(packed_array));
  p->n = n;  p->w = w;
  x = (n / PBS)*w + ((n % PBS)*w + PBS-1) / PBS+1;
//  x= (n * w + PBS - 1) / PBS;
  p->b = dbwt_mymalloc(x*sizeof(pb));
  for (i=0; i<x; i++) p->b[i] = 0;
  return p;
}

void dbwt_free_packed_array(packed_array *p)
{
  ulong x;
//  x = (p->n/PBS+1) * p->w;
  x = (p->n / PBS)*p->w + ((p->n % PBS)*p->w + PBS-1) / PBS+1;
//  x= (p->n * p->w + PBS - 1) / PBS;
  dbwt_myfree(p->b,x*sizeof(pb));
  dbwt_myfree(p,sizeof(packed_array));
}

ulong dbwt_pa_get(packed_array *p, ulong i)
{
  int w;
  pb *b;

  w = p->w;
  b = p->b + (i>>logD)*w;
  i = (i % D)*w;

  return (ulong) dbwt_getbits(b,i,p->w);
}

void dbwt_pa_set(packed_array *p, ulong i, long x)
{
  int w;
  pb *b;
#if 0
  if (x < 0 || x > (1<<p->w)) {
    printf("pa_set: x=%ld w=%d\n",x,p->w);
  }
  if (i < 0 || i >= p->n) {
    printf("pa_set: i=%ld n=%d\n",i,p->n);
  }
#endif
  w = p->w;
  b = p->b + (i>>logD)*w;
  i = (i % D)*w;

  dbwt_setbits(b,i,p->w,x);
}

