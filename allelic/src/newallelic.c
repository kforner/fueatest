#include <stdio.h>

#include <math.h>
#include <stdlib.h>

/*
This is the only setting that you may want to tune : this is the max size of 
the precomputed log factorials array
It means that the sum of the cells in your table must be lower than that maximum
otherwise an error message will be printed on standard error, and the function
will return -1 
*/
#define TABLE_OF_LOG_FACTORIALS_SIZE 50000

#define MYMIN(x,y) ((x) < (y) ? (x) : (y))
#define MYMAX(x,y) ((x) > (y) ? (x) : (y))


/**
 init a static precomputed table of log factorials : size is in the defined TABLE_OF_LOG_FACTORIALS_SIZE

 Arg size_needed: the max log factorial needed

 Arg from_where: the name of the user function

 Return: the pointer on the table, or 0 otherwise

*/
double* init_table_of_log_factorial(int size_needed, char* from_where) {
  static double T[TABLE_OF_LOG_FACTORIALS_SIZE];
  static int is_table_init = 0;

  /* == check table size ===*/
  if ( size_needed >= TABLE_OF_LOG_FACTORIALS_SIZE ) {
    fprintf(stderr,"[ called by %s]:log factorial needed too high (n=%i) for precomputed table, increase the TABLE_OF_LOG_FACTORIALS_SIZE #define value !!!!!!!!!\n", from_where, size_needed );
    return 0;
  }

  /* === precompute table if needed ===*/
  if ( ! is_table_init ) {
    int i;
    T[0] = 0;
    for (i = 1; i < TABLE_OF_LOG_FACTORIALS_SIZE; i++) 
      T[i] = T[i-1] + log( i );
    
    is_table_init = 1;
  }

  return (double*)&T;
}

double approximatePearsonAllelicExactTestOnGenotypicData(int aa,int bb,int c,int d,int e,int f)  
{

  double* T = 0;
  int b_min = 0, b_max = 0, ai,a_start,a_end,bi,b_start,b_end,ci,ei,fi,i;
  double pvalue = 0, num = 0, p = 0;
  
  int c1 = aa + d;
  int c2 = bb + e;
  int c3 = c + f;
  int l1 = aa + bb + c;
  int l2 = d + e + f;
  int n = l1 + l2;
  int N = n + n;

  int b0 = c2 + l1 + c1 - n;

  int A1 = aa + aa + bb;
  int N_AS = (c1 + c1 + c2) * (l1 + l1);

  int a_min = MYMAX(0,c1 - l2);
  int a_max = MYMIN(c1,l1);
  
  int observed = abs(N*A1 - N_AS);

  int alpha = (int)floor((float)(N_AS - observed)/(float)N) + 1;
  int beta = (int)ceil((float)(N_AS + observed)/(float)N) - 1;

  int c_alpha_c2 = (int)ceil( (float)(alpha-c2)/2.0);
  int beta_c3_l1 = beta + c3 - l1;
  int f_b2 = (int)floor((float)(beta)/2.0);

  int starts[4] = { 0,0,0,0};
  int ends[4] = { 0,0,0,0};

  if (alpha > beta)
    return 1;

  /* ========== PRECOMPUTED LOG FACTORIALS TABLE ==========*/

  T = init_table_of_log_factorial(n, "approximatePearsonAllelicExactTestOnGenotypicData");
  if ( 0 == T )
    return -1;
  
  num = T[c1] + T[c2] + T[c3] + T[l1] + T[l2] - T[n];
  
  starts[0] = c_alpha_c2;
  ends[0] = MYMIN(MYMIN(l1-c2,l1-c3),beta_c3_l1);

  starts[1] = MYMAX(c_alpha_c2,l1-c3+1);
  ends[1] = MYMIN(l1-c2,f_b2);

  starts[2] = MYMAX(l1-c2+1,alpha-l1);
  ends[2] = MYMIN(l1-c3,beta_c3_l1);

  starts[3] = MYMAX(MYMAX(alpha-l1,l1-c2+1),l1-c3+1);
  ends[3] = f_b2;
  for (i = 0; i < 4; i++) {
    a_start = MYMAX(starts[i],a_min);
    a_end = MYMIN(ends[i],a_max);
    for (ai = a_start; ai <= a_end; ai++) {
      b_min = MYMAX(0, b0 - ai);
      b_max = MYMIN(l1 - ai,c2);
      b_start = MYMAX(alpha-ai-ai,b_min);
      b_end = MYMIN(beta-ai-ai,b_max);
      ci = l1 - ai - b_start;
      ei = c2 - b_start;
      fi = c3 - l1 + ai + b_start;
      p = exp(num -T[ai]-T[b_start]-T[ci]-T[c1-ai]-T[ei]-T[fi] );
      pvalue += p;
      for (bi = b_start + 1; bi <= b_end; bi++) {
	++fi;
	p *= (double)(ei*ci)/(double)(bi*fi);
	pvalue += p;
	--ei;
	--ci;
      }
    }
  }

  pvalue = MYMIN(1 - pvalue, 1);
  pvalue = MYMAX(0,pvalue);
 
  return pvalue;
}

void newallelic(int* a, int* b, int* c,int* d,int *e, int* f, double* pv) {
  *pv = approximatePearsonAllelicExactTestOnGenotypicData(*a,*b,*c,*d,*e,*f);
}

