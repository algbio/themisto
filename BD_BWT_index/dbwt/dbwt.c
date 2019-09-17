/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
//#include "mman.h"
#include "dbwt_utils.h"
#include "dbwt_queue.h"
#include "dbwt.h"

#define HSIZ 67777
//#define HSIZ 375559
#define HBSIZ (256*2)

typedef struct {
  int rest[HSIZ];
  int head[HSIZ];
  uchar *buf;
  int bufsiz;
} htbl;

htbl *h1;

static long collision;

static htbl *init_hashtable(void)
{
  htbl *t;
  int i;
  t = (htbl *) dbwt_mymalloc(sizeof(htbl));
  if (t == NULL) {
    printf("init_hashtable: malloc (1) failed.\n");
    exit(1);
  }
  for (i=0; i<HSIZ; i++) {
    t->rest[i] = 0;
    t->head[i] = 0;
  }
  t->bufsiz = 0;
  t->buf = NULL;
  collision = 0;
  return t;
}

static int hfunc(int m, uchar *p)
{
  ulong x;
  int i;

  x = 0;
  for (i=0; i<m; i++) {
    x *= 101;
    x += *p++;
  }
  x %= HSIZ;

  return x;
}

static ulong getpointer(uchar *p, int k)
{
  ulong x;
  int i;
  
  x = 0;
  for (i=0; i<k; i++) {
    x <<= 8;  x += *p++;
  }

  return x;
}

static void setpointer(uchar *p, ulong x, int k)
{
  int i;

  p += k;

  for (i=k; i>0; i--) {
    *--p = x & 0xff;  x >>= 8;
  }

}

static void setlen(uchar **p, ulong x)
{
  uchar *q,c;
  int i,w;
  ulong y;

  q = *p;
  y = x;
  w = 1;
  while (y > 0x7f) {
    y >>= 7;
    w++;
  }
  for (i=w-1; i>0; i--) {
    c = (x >> (i*7)) | 0x80;
    *q++ = c;
  }
  c = x & 0x7f;
  *q++ = c;

  *p = q;
}

ulong getlen(uchar **p)
{
  uchar *q,c;
  ulong x;

  q = *p;
  x = 0;
  while (1) {
    c = *q++;
    x += (c & 0x7f);
    if ((c & 0x80) == 0) break;
    x <<= 7;
  }
  *p = q;

  return x;
}

static int lenlen(ulong x)
{
  int w;

  w = 1;
  while (x > 0x7f) {
    x >>= 7;
    w++;
  }
  return w;
}

#define CONT 2

// format: length (1 byte), string (m bytes)
// if length==CONT, the following 4 bytes are the pointer to next element
// if length==0, end of list

static int insert_string(htbl *t, int m, uchar *p)
{
  int i,l,l2;
  int h;
  uchar *r,*r2;
  uchar *buf;
  int q, pp, q2;
  int w = sizeof(int); // length of name

  buf = t->buf;

  h = hfunc(m,p);

  q = t->head[h];
  while (q != 0) {
    pp = q;
    r = &buf[q];
    l = getlen(&r); // length
    if (l == 0) { // not found
      r = &buf[pp];
      break;
    }
    if (l == CONT) {// end of block
      q = getpointer(r,w);
      collision++;
      continue;
    }
    if (l != m) {
      q += lenlen(l);
      q += l;
      q += 1; // for sentinel
      q += w; // for name
      collision++;
      continue;
    }

    l2 = lenlen(l);
    for (i=0; i<m; i++) {
      if (p[i] != buf[q+1+l2+i]) break;
    }
    if (i == m) { // match
      //printf(" match\n");
      return 0; // existing string
    }
    q += l2;
    q += l;
    q += 1; // for sentinel
    q += w; // for name
    collision++;
  }

  if (m+lenlen(m)+1+w >= t->rest[h]-10) { // +1 stands for the space to store len
    r2 = (uchar *) dbwt_myrealloc(buf, t->bufsiz + HBSIZ + m+lenlen(m)+1+w, t->bufsiz);
    if (r2 != buf) {
//      printf("buf has moved from %p to %p\n", buf, r2);
      t->buf = buf = r2;
    }
    q2 = t->bufsiz+1;
    t->bufsiz += HBSIZ + m+lenlen(m)+1+w;
    t->rest[h] = HBSIZ + m+lenlen(m)+1+w;
    if (q == 0) { // first block
      t->head[h] = q2;
    } else {
      r2 = &buf[pp];
      setlen(&r2,CONT);
      setpointer(r2,q2,sizeof(w));
    }
    r = &buf[q2];
  }

  setlen(&r,m);
  *r++ = p[0]+1; // sentinel
  for (i=0; i<m; i++) *r++ = p[i];
  setpointer(r,0x99999999,w);
  r += w;

  r[0] = 0;

  t->rest[h] -= m+lenlen(m)+1+w;

  return 1; // new string
}

static int getname(htbl *t, int m, uchar *p)
{
  int i,l,l2;
  int h;
  uchar *r;
  uchar *buf;
  int q;
  int w = sizeof(int);

  buf = t->buf;

  h = hfunc(m,p);

  q = t->head[h];
  while (q != 0) {
    r = &buf[q];
    l = getlen(&r); // length
    if (l == 0) { // not found
      break;
    }
    if (l == CONT) {// end of block
      q = getpointer(r,w);
      continue;
    }
    if (l != m) {
      q += lenlen(l);
      q += l;
      q += 1;
      q += w; // for name
      continue;
    }

    l2 = lenlen(l);
    for (i=0; i<m; i++) {
      if (p[i] != buf[q+1+l2+i]) break;
    }
    if (i == m) { // match
      return getpointer(&buf[q+1+l2+m],w);
    }
    q += l2;
    q += l;
    q += 1;
    q += w; // for name
  }
  printf("??? not found\n");
  return -1;
}


#ifndef TYPE_S
#define TYPE_S 0
#define TYPE_L 1
#define TYPE_LMS 2
#endif

#define gc(i) (T[i])

static int LMS_compare(const void *p1, const void *p2)
{
  int c1,c2;
  int l1,l2;
  uchar *q1, *q2;

  q1 = *(uchar **)p1;  q2 = *(uchar **)p2;
  l1 = getlen(&q1);  l2 = getlen(&q2);
  
  q1++;  q2++;

  while (l1 > 0 && l2 > 0) {
    c1 = *q1++;  c2 = *q2++;
    if (c1 != c2) break;
    l1--;  l2--;
  }
  if (l1 == 0) return 1;
  if (l2 == 0) return -1;
  return c1 - c2;
}
/*
static void printlstr(uchar *p)
{
  int i,l;
  l = getlen(&p);
  //printf("(%d) ",l);
  //for (i=0; i<l; i++) printf("%02x ",p[i]);
  for (i=0; i<l-1; i++) printf("%c",p[i+1]);
}
*/
static uchar *min_ptr, *max_ptr;

static uchar **sort_LMS(int n, htbl *h)
{
  uchar **s, *r, *q;
  int p;
  int i,j,l,m;

  s = (uchar **) dbwt_mymalloc((n+2) * sizeof(*s));
//  dbwt_report_mem("allocate S* ptr");
  // s は S*-substring へのポインタを適当な順番に並べたもの
  // s[0] は文字列先頭のsubstringを格納するためで，S* ではない
  // s[1] から s[n-1] には S* が入る
  // s[n] には $ を入れる (実際には入っていない)
  
  j = 1; // j=0 is for the head string
  for (i = 0; i < HSIZ; i++) {
    p = h->head[i];
    while (p != 0) {
  // r からのメモリには [文字列長][番兵][文字列][name] が入る
  // 文字列長は可変長符号
  // 番兵は1バイト．文字列の先頭の文字よりも大きい文字
  // (文字列の先頭はTYPE_Sなので 0xff になることはない)
  // nameは4バイト
      r = &h->buf[p];
      q = r;
      l = getlen(&r);
      if (l == 0) break;
      if (l == CONT) {
        p = getpointer(r,sizeof(int));
        continue;
      }
      if (j > n) {
        printf("j = %d\n",j);
      }
      s[j++] = q;
      p += lenlen(l);
      p += 1; // for sentinel;
      p += l;
      p += sizeof(int);
    }
  }
  m = j-1;

  qsort(s+1,m,sizeof(uchar *),LMS_compare);

  for (i=1; i<=m; i++) {
    q = s[i];
//    printf("s[%d] = %p ",i,q);  printlstr(q);  printf("\n");
    l = getlen(&q);
    q += 1;
    q += l;
    setpointer(q,(ulong)i,sizeof(int));
  }
  // ソート後は，S* に 1 から n-1 の数字がつく
  // $ は 0 にする (ここでは代入していない)

  return s;
}

#define SIGMA (256+1)



uchar * dbwt_bwt(uchar * T,long n,unsigned int *_last,unsigned int free_text)
{
  long i,j;
  int t,tt;
  long p,q;
  long n1;
  long s1; // alphabet size
  long m1; // total length of substrings
  int c;
  ulong last;
  uchar *lastptr;

  uchar **S;
  long C[SIGMA+2]; // 文字 c から始まるS*文字列の数
  long M[SIGMA+2]; // 文字 c の頻度 (c は -1 から SIGMA-1)
  long NL[SIGMA+2]; // TYPE_Lの文字 c の数
  long M2[SIGMA+2], M3[SIGMA+2];
  long C2[SIGMA+2], cc;

  dbwt_queue *Q[3][SIGMA+2];

  uint *sa; // uint??
  long sa_size;
  uchar *bw;
  packed_array *T2;
  uchar *bwp_base;
  int bwp_w;
//  FILE *fp;
//  MMAP *map;

///////////////////////////////////////////////
// ファイルの読み込み
/*
#if 1
  fp = fopen(fname,"rb");
  if (fp == NULL) {
    printf("??? %s\n",fname);
    exit(1);
  }
  fseek(fp,0,SEEK_END);
  n = ftell(fp);
  fseek(fp,0,SEEK_SET);

  T = mymalloc((n+1)*sizeof(uchar));
  report_mem("read string");
  
  for (i=0; i<n; i += (1<<20)) {
    long size;
    printf("%ld \r",i/(1<<20));
    fflush(stdout);
    size = min(n-i,1<<20);
    fread(&T[i],1,size,fp);
  }
  fclose(fp);
#else
  map = mymmap(fname);
  n = map->len;
  T = map->addr;
#endif */
///////////////////////////////////////////////

  h1 = init_hashtable();

  s1 = 0; // S*の数
  m1 = 0; // S*の長さの合計
  n1 = 0; // 次のレベルの文字列長 (S*の数)

  for (i=0; i<=SIGMA+1; i++) {
    C[i] = C2[i] = 0;
    M[i] = M2[i] = 0;
    NL[i] = 0; // TYPE_L の数
  }

  // C[c] : c で始まるS*
  // M[c]: c の数


//////////////////////////
// T[n] = $ が1つのS*-string
  t = TYPE_S; // i = n;
  M[(-1)+1]++; // 文字の頻度

//////////////////////////
// T[n-1] は必ずTYPE_L
  tt = t; // 1つ前のTYPE
  i = n-1;
  t = TYPE_L;
  M[gc(i)+1]++;
  NL[gc(i)+1]++; // TYPE_Lの文字の頻度

//////////////////////////
// T[p..q] = T[n..n] が S*-string
  p = i+1;
  q = p;
  C[(-1)+1]++; // T[p] = $ から始まるS*の数
  s1++; // S*の種類の数
  m1 += 1; // S*の長さの合計
  n1++; // 次のレベルの文字列長 (S*の数)
  tt = t;

  for (i=n-2; i>=0; i--) {
/*    if (i % (1<<20) == 0) {
      printf("%ld \r",i>>20);
      fflush(stdout);
    }*/
    M[gc(i)+1]++;
    if (gc(i) < gc(i+1)) { // TYPE_S
      t = TYPE_S;
    } else if (gc(i) > gc(i+1)) { // TYPE_L
      if (tt == TYPE_S) { // 直前の文字がTYPE_Sならば S* T[i+1..q] を作る
        p = i+1;
        C[gc(p)+1]++;
        if (insert_string(h1,q-p+1,&T[p])) { // T[p..q] を挿入
          s1++; // 新しいS*ならば種類の数を増やす
          m1 += q-p+1; // S*の長さを足す
        }
        n1++;
        q = p;
      }
      t = TYPE_L; // T[i] はTYPE_L
      NL[gc(i)+1]++;
    } else { // TYPEは前の文字と同じ
      if (t == TYPE_L) NL[gc(i)+1]++;
    }
    tt = t;
  }
//  printf("n=%ld n1=%ld\n",n,n1);
//  printf("s1=%ld m1=%ld\n",s1,m1);
  //printf("space for level 1 %ld bytes\n",n1/8*blog(s1)+m1);
//  printf("collision1 %ld\n",collision);

//  dbwt_report_mem("compute S*");

  { // T[0..p] は S* ではないがBW変換には必要なので保存
    uchar *r;
    r = dbwt_myrealloc(h1->buf, h1->bufsiz + p+1+1+lenlen(p+1) + 4, h1->bufsiz);
    if (r != h1->buf) {
//      printf("buf has moved from %p to %p\n", h1->buf, r);
      h1->buf = r;
    }
    r = &h1->buf[h1->bufsiz];
    h1->bufsiz += p+1+1+lenlen(p+1) + 4;

    S = sort_LMS(s1, h1);

    S[s1] = r;
    setlen(&r,p+1);
    *r++ = T[0]+1; // sentinel
    lastptr = r; // T[0] へのポインタ
    for (i=0; i<p+1; i++) *r++ = T[i]; // 先頭の文字列 T[0..p]

    S[0] = r;
    setlen(&r,1);
    *r = 0; // 最後の文字列 T[n..n]．本当は -1
  }

  min_ptr = max_ptr = S[0];
  for (i=1; i<=s1; i++) {
    if (S[i] < min_ptr) min_ptr = S[i];
    if (S[i] > max_ptr) max_ptr = S[i];
  } 

  T2 = dbwt_allocate_packed_array(n1+1,dbwt_blog(s1+2)+1);
//  dbwt_report_mem("create new string T2");

//////////////////////////
// T[n] = $ が1つのS*-string
  t = TYPE_S; // i = n;
  tt = t;
  q = n;
  j = n1;

//////////////////////////
// T[n-1] は必ずTYPE_L
  i = n-1;
  t = TYPE_L;
  p = i+1;
  dbwt_pa_set(T2,j--,0); // T2[n1] = 0 ($のname)

  q = p;
  tt = t;

  for (i=n-2; i>=0; i--) {
//    if (i % (1<<20) == 0) {
  //    printf("%ld \r",i>>20);
//      fflush(stdout);
  //  }
    if (gc(i) < gc(i+1)) {
      t = TYPE_S;
    } else if (gc(i) > gc(i+1)) {
      if (tt == TYPE_S) {
        p = i+1;
        dbwt_pa_set(T2,j--,getname(h1,q-p+1,&T[p])); // T[p..q] のnameを T2 に書く
        q = p;
      }
      t = TYPE_L;
    }
    tt = t;
  }
  dbwt_pa_set(T2,0,s1); // T2[0] = s1 (Tの先頭の文字列のname)

#if 1
  if(free_text)
    free(T);
/* {
    dbwt_myfree(T,(n+1)*sizeof(uchar));
    dbwt_report_mem("free T");
  };*/
#else
  mymunmap(map);
#endif

  sa_size=(n1+1)*sizeof(*sa);
  sa = (uint *) dbwt_mymalloc(sa_size);
//  dbwt_report_mem("allocate sa");
////////////////////////////////////
// T2[1..n1] の文字列の接尾辞配列を作成
//  printf("sorting level-1 suffixes using IS...\n");
  dbwt_sais_main((const unsigned char *) T2->b, (int *)sa, 0, n1, s1+1, -T2->w);

  for (i=0; i<n1; i++) {
    sa[i]++; // sa[i] は 1 から s1
  }
//  printf("done\n");
////////////////////////////////////

  bwp_base = min_ptr-1;
  bwp_w = dbwt_blog(max_ptr - min_ptr + 2)+1;

  for (i=0; i<n1; i++) {
    int l;
    uchar *q;
    p = sa[i];
    q = S[dbwt_pa_get(T2,p-1)]; // bw[i] = T2[sa[i]-1] の文字列
    l = getlen(&q);
    q += 1; // for sentinel
    q += l-1; // 文字列の最後の文字を指す (後ろからBWに書いていくから)

    sa[i] = q - bwp_base;
  }
  dbwt_free_packed_array(T2);
//  dbwt_report_mem("free T2");

  for (i=0; i<=SIGMA; i++) {
    Q[TYPE_LMS][i] = dbwt_init_queue(bwp_w);
    Q[TYPE_L][i] = dbwt_init_queue(bwp_w);
    Q[TYPE_S][i] = dbwt_init_queue(bwp_w);
  }

  j = n1+1; // sa のサイズ
  for (i=n1-1; i>=0; i--) {
    uchar *q;
//    uint *sa2;
    q = sa[i] + bwp_base;
    if (i == 0) {
      c = -1;
    } else {
      c = q[0];
    }
    dbwt_enqueue_l(Q[TYPE_LMS][c+1], q-bwp_base);
    if (i % 1024 == 0) {
//      sa2 = sa;
      sa_size=i*sizeof(*sa);
      sa = dbwt_myrealloc(sa,sa_size,j*sizeof(*sa));
      //if (sa != sa2) printf("sa = %p\n",sa);
      j = i;
    }
  }
//  dbwt_report_mem("allocate queue");

  bw = (uchar *) dbwt_mymalloc(n+1);
    
  cc = 0;
  for (i=0; i<=SIGMA; i++) {
    M3[i] = M2[i] = cc; // M2[i] is the first index for i
    cc += M[i];
  }

//  dbwt_report_mem("allocate tmp");
    
//  printf("induced-sorting-L n=%ld\n",n);

  // C2 は TYPE_S バケットの先頭
  for (i=0; i<=SIGMA; i++) C2[i] = M2[i] + NL[i];

  // 使う配列
  // M2: コピー先のアドレス, 初期値はバケットの先頭
  // C2: bw を一時的に格納するスタックのアドレス, 初期値はSバケットの位置
  // M3: バケットの位置 (変化しない)
    
  {
    ////////////////////////////////////////////
    // $ のバケットの文字列 (1つだけ) を取り出す
    // q はS*の最後の文字を指す
    // c2 = q[0] = $ (-1)
    // c1 = q[-1] がBW変換の文字
    // c1 のバケットに q-1 を書く
    uchar *q;
    int c1,c2;
    i = 0;
    q = bwp_base + dbwt_dequeue(Q[TYPE_LMS][i]);
    c1 = q[-1];  c2 = -1;
    bw[0] = c1;
    dbwt_enqueue(Q[TYPE_L][c1+1],(q-1)-bwp_base);
    bw[M2[c1+1]++] = q[-2];
    bw[C2[c2+1]++] = c1; // 一時的に格納
  }

  for (c = 1; c <= SIGMA; c++) {
    uchar *q;
    int c1;
    long m;
    int t;
//    printf("%d \r",c);
//    fflush(stdout);
    for (t=1; t<=2; t++) { // TYPE_L と TYPE_LMS を順に処理
      m = M3[c]; // Lバケットの先頭
      if (t == TYPE_LMS) m += M[c] - C[c]; // LMSバケットの先頭
      while (!dbwt_emptyqueue(Q[t][c])) {
        q = bwp_base + dbwt_dequeue(Q[t][c]);
        if (q == lastptr) { // T[0] ならば，その前は $ なのでその辞書順をlastに格納
          last = m;
        } else {
          c1 = q[-1];  // c2 = q[0]; // c2+1 = c
          if (c1 >= (c-1)) { // TYPE_L
            dbwt_enqueue(Q[TYPE_L][c1+1], (q-1)-bwp_base);
            bw[M2[c1+1]++] = q[-2];
            if (t == TYPE_LMS) {
              bw[C2[c]++] = c1; // 一時的に格納
            }
          } else {
            dbwt_enqueue_l(Q[TYPE_S][c], (q)-bwp_base);
          }
        }
        m++;
      }
    }
  }

  for (c=0; c<=SIGMA; c++) {
    dbwt_free_queue(Q[TYPE_LMS][c]);
    dbwt_free_queue(Q[TYPE_L][c]);
  }
    
  for (c=0; c<=SIGMA; c++) {
    Q[TYPE_L][c] = dbwt_init_queue(bwp_w);
  }
    
//  dbwt_report_mem("allocate queue");
    
//  printf("induced-sorting-S n=%ld\n",n);

  // 使う配列: M2, C2
  // M2: コピー先のアドレス, 初期値はバケットの最後
  // C2: bw を一時的に格納するスタックのアドレス, 初期値は前の続き

  M2[0] = 0;  for (c=1; c<=SIGMA; c++) M2[c] = M2[c-1] + M[c];

  for (c = SIGMA; c>=0; c--) {
    uchar *q;
    int c1,c0;
    int t;
//    printf("%d \r",c);
//    fflush(stdout);
    for (t = 1; t >= 0; t--) { // TYPE_L, TYPE_S の順に処理
      while (!dbwt_emptyqueue(Q[t][c])) {
        q = bwp_base + dbwt_dequeue(Q[t][c]);
        c1 = q[-1];  // c2 = q[0]; // c2+1 = c
        if (c1 <= (c-1)) { // TYPE_S
          dbwt_enqueue(Q[TYPE_L][c1+1],(q-1) - bwp_base); // 入れるのは常に L
          if (q-1 == lastptr) {
            last = M2[c1+1]--;
          } else {
            c0 = q[-2];
            bw[M2[c1+1]--] = (c0 <= c1) ? c0 : bw[--C2[c1+1]];
          }
        }
      }
    }
  }

/*  printf("writing...\n");
  fp = fopen("output.bw","wb");
  if (fp == NULL) {
    printf("fopen\n");
    exit(1);
  }

  {
    uchar *p;
    long s,t,w;
    s = last;
    p = bw;
    w = 1<<20;
    while (s > 0) {
      printf("%ld \r",(p-bw)/w);
      fflush(stdout);
      t = min(s,w);
      fwrite(p,1,t,fp);
      s -= t;
      p += t;
    }
    s = n-last;
    p = &bw[last+1];
    while (s > 0) {
      printf("%ld \r",(p-bw)/w);
      fflush(stdout);
      t = min(s,w);
      fwrite(p,1,t,fp);
      s -= t;
      p += t;
    }
  }

  fclose(fp);

  fp = fopen("output.lst","w");
  if (fp == NULL) {
    printf("fopen2\n");
    exit(1);
  }
  fprintf(fp,"%lu",last);
  fclose(fp);*/
  

 // dbwt_report_mem("done");
//  printf("%.2f bpc\n",(double) dbwt_max_alloc/n);
  //printf("end");  getchar();
  for (c=0; c<=SIGMA; c++) {
    dbwt_free_queue(Q[TYPE_L][c]);
    dbwt_free_queue(Q[TYPE_S][c]);
  }
  dbwt_myfree(S,(s1+2)*sizeof(*S));
  dbwt_myfree(sa,sa_size);
//  myfree(bw,n+1);
  dbwt_myfree(h1->buf,h1->bufsiz);
  dbwt_myfree(h1,sizeof(htbl));
  (*_last)=last;
//  dbwt_report_mem("Freed all memory");
  return bw;
}
/*
#if 1
int main(int argc, char *argv[])
{

  printf("sizeof(uchar *)=%ld sizeof(uint)=%ld sizeof(ulong)=%ld\n",sizeof(uchar *),sizeof(uint),sizeof(ulong));

  bwt(argv[1]);
  
  return 0;
}
#endif*/
