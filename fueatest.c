/*
Compilation: gcc -O2  -o fast.allelic.exact.test fast.allelic.exact.test.c -lm
*/

/*
 Author : Karl FORNER, Copyright (c) Serono Research Institute, Switzerland
 LICENSE :  GPL version 2 or newer
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MYMIN(x,y) ((x) < (y) ? (x) : (y))
#define MYMAX(x,y) ((x) > (y) ? (x) : (y))
#define TABLE_OF_LOG_FACTORIALS_SIZE 50000

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
    fprintf(stderr,"[ called by %s]:log factorial needed too high (n=%i) for precomputed table, increase the TABLE_OF_LOG_FACTORIALS_SIZE #define value !!!!!!!!!\n",from_where, size_needed);
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

/*
 This version does exactly what is described in the article
" A fast, unbiased and exact allelic test" 
*/
double approximatePearsonAllelicExactTestOnGenotypicData_opt(int aa,int bb,int c,int d,int e,int f)  
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


  T = init_table_of_log_factorial(n, "approximatePearsonAllelicExactTestOnGenotypicData_opt");
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

void usage() {
  printf(
"NAME\n"
"       fast.allelic.exact.test - A fast, unbiased and exact allelic exact test\n"
"\n"
"SYNOPSIS\n"
"       fast.allelic.exact.test [options] [file] | [d0 d1 d2 h0 h1 h2]\n"
"\n"
"        Options:\n"
"          -help            brief help message\n"
"          -man             full documentation\n"
"          -debug           execute in debug mode\n"
"\n"
"       Take the input data filename given as argument\n"
"\n"
"DESCRIPTION\n"
"       This is the implementation in C of a new association test\n"
"       described in \"A fast, unbiased and exact allelic exact test for case-\n"
"       control association studies. Human Heredity (in press)\".\n"
"\n"
"       It appears that in most cases the classical chi-square test used for\n"
"       testing for allelic association on genotype data is biased.  Our test\n"
"       is unbiased, exact but fast through careful optimization.\n\n"
"       See http://stat.genopole.cnrs.fr/software/fueatest\n"
"\n"

"INPUT FILE FORMAT\n"
"       Notations\n"
"\n"
"       Here is the 2x3 contingency table:\n"
"\n"
"                           aa aA AA\n"
"        [case (diseased)]  d0 d1 d2\n"
"        [control(healthy)] h0 h1 h2\n"
"\n"
"       Format\n"
"\n"
"       simple : one table per line, the six counts separated by semicolons.\n"
"\n"
"       example\n"
"\n"
"        226;57;5;249;63;4\n"
"        1;109;191;0;110;221\n"
"        7;110;174;6;132;191\n"
"\n"
"OUTPUT\n"
"       Output the contingency tables and their corresponding p-values in sci‐\n"
"       entific notation separated by tabs, one line per input table.\n"
"\n"
"EXAMPLES\n"
"       Suppose you have a file toto.csv such as the example just before\n"
"\n"
"       fast.allelic.exact.test  toto.csv\n"
"\n"
"\n"
"       Or if you want to experiment with tables, just type the counts on the\n"
"       command line:\n"
"\n"
"        fast.allelic.exact.test 226 57 5 249 63 4\n"
"\n"
"AUTHOR\n"
"       Karl FORNER <karl.forner@gmail.com>\n"
"\n"
"MAINTAINER\n"
"       Mickael Guedj <guedj@genopole.cnrs.fr>\n"
"\n"
"VERSION\n"
"       $Id:  filename  revision  date  time  author  state\n"
"\n"
"       $$Id: fast.allelic.exact.test.c,v 1.3 2006/06/06 07:59:00 kforner Exp $\n"
" \n"
	 );
}

void process_table(int* t) {
  double pv = approximatePearsonAllelicExactTestOnGenotypicData_opt(t[0],t[1],t[2],t[3],t[4],t[5]);
  printf("%i\t%i\t%i\t%i\t%i\t%i\t%e\n",t[0],t[1],t[2],t[3],t[4],t[5],pv);


}

int main(int argc, char** argv) {
  int i = 0;
  int line = 1;
  int t[] = { 0,0,0,0,0,0 };
  double pv = 0;
  FILE* file = 0;

  if (argc < 2 ) {
    usage();
    exit(1);
  }

  if (argc > 6) {
    for (i = 1; i <= 6; i++) {
      t[i-1] = atoi( argv[i] );
    }
    process_table( t );
  } else {
    if ( ! (file = fopen(argv[1],"r")) ) {
      fprintf(stderr,"unable to open file [%s], aborting...",argv[1]);
      exit(1);
    }

    while ( (i = fscanf(file,"%i;%i;%i;%i;%i;%i\n",t,t+1,t+2,t+3,t+4,t+5)) != EOF) {
      if (6 != i) {
	fprintf(stderr,"Error while parsing input file line %i, exiting...\n",line);
      }
      process_table( t );
      line++;
    }
  }
}

