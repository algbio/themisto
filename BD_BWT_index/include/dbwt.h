#ifndef dbwt_h
#define dbwt_h
#ifndef uchar
typedef unsigned char uchar;
#endif
int dbwt_sais_main(const unsigned char *T, int *SA, int fs, int n, int k, int cs);
int dbwt_sais_int(const int *T, int *SA, int n, int k);
int dbwt_sais(const unsigned char *T, int *SA, int n);
/**
 * Performs bwt on a given unsigned char array.
 *
 * @param T input string
 * @param n length of input string
 * @param *_last index of the end of the string
 * @param free_text a boolean that determines some option, leave at 0
 *
 * @return Burrows-Wheeler transformed unsigned char array
 */
uchar * dbwt_bwt(uchar * T,long n,unsigned int *_last,unsigned int free_text);

#endif
