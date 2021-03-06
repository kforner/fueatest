#! /usr/bin/env perl

# Author : Karl FORNER, Copyright (c) Serono Research Institute, Switzerland
# LICENSE :  GPL version 2 or newer

use warnings; 
use strict;  
use Carp qw(carp cluck confess croak);
use POSIX qw(floor ceil);
use Getopt::Long;
use autouse 'Pod::Usage' =>  'pod2usage';    # postpone load of module until pod2usage is used !!

my $SEP = ';'; # INPUT FORMAT SEPARATOR
########## ARGS PROCESSING ##########
my %OPTIONS;
GetOptions(\%OPTIONS,qw(help|? man debug)) or pod2usage(2);
my $DEBUG = $OPTIONS{debug} || 0;
pod2usage(1) if $OPTIONS{help};
pod2usage( -verbose => 2 ) if $OPTIONS{man};

my $file = $ARGV[0];
if ( ! defined($file) && -t STDIN ) {
  pod2usage("No input data\n");
}

if ( @ARGV < 6 ) {
  my @table;
  my $ln = [];
  my $pv;
  while( <> ) {
    chomp;
    @table = split $SEP;
    process_table(\@table);
  }
} else {
  process_table(\@ARGV);
}

sub process_table {
  my ($t) = @_;
  my $pv = approximatePearsonAllelicExactTestOnGenotypicData_opt(@$t);
  print join("\t",@$t,sprintf("%e",$pv)),"\n";
}

=head1 NAME

fast.allelic.exact.test.pl - A fast, unbiased and exact allelic exact test

=head1 SYNOPSIS

fast.allelic.exact.test.pl [options] [<tables] [file1 ... fileN] [d0 d1 d2 h0 h1 h2]

 Options:
   -help            brief help message
   -man             full documentation
   -debug           execute in debug mode

Take the input data either from the filenames given as argument, or from STDIN, or from
the table cell counts given in the command line.



=head1 DESCRIPTION

This is the implementation in perl and B<C> (using the wonderful L<Inline::C> module) of a new association
test described in B<"A fast, unbiased and exact allelic exact test
for case-control association studies. Human Heredity (in press)">.

This version is much faster than the corresponding I<pure perl>
version named B<allelic.exact.test.pl>.

It appears that in most cases the classical chi-square test used
for testing for B<allelic> association on genotype data is biased.
Our test is B<unbiased>, B<exact> but B<fast> through careful optimization.

See http://stat.genopole.cnrs.fr/software/fueatest

=head1 INSTALLATION

This script requires the B<Inline::C> module, available from CPAN (www.cpan.org)
, and a C compiler (only tested with B<gcc> available from B<gcc.gnu.org> ).

=head1 INPUT FILE FORMAT

=head2 Notations

Here is the 2x3 contingency table:

                    aa aA AA
 [case (diseased)]  d0 d1 d2
 [control(healthy)] h0 h1 h2

=head2 Format

simple : one table per line, the six counts separated by semicolons.

=head1 OUTPUT

Output the contingency tables and their corresponding p-values in scientific format
separated by tabs, one line per input table.

=head2 example

 226;57;5;249;63;4
 1;109;191;0;110;221
 7;110;174;6;132;191


=head1 EXAMPLES

Suppose you have a file toto.csv such as the example just before

 fast.allelic.exact.test.pl toto.csv

or 

 fast.allelic.exact.test.pl < toto.csv

Or if you want to experiment with tables, just type the counts on the command line:

 fast.allelic.exact.test.pl 226 57 5 249 63 4


=head1 AUTHOR

Karl FORNER <karl.forner@gmail.com>

=head1 MAINTAINER

Mickael Guedj <guedj@genopole.cnrs.fr>

=head1 VERSION

$Id:  filename  revision  date  time  author  state

$$Id: fast.allelic.exact.test.pl,v 1.3 2006/06/06 07:58:57 kforner Exp $

=cut


use Inline "C" => << 'END_OF_C_CODE';

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

  // == check table size ===
  if ( size_needed >= TABLE_OF_LOG_FACTORIALS_SIZE ) {
    fprintf(stderr,"[ called by %s]:log factorial needed too high (n=%i) for precomputed table, increase the TABLE_OF_LOG_FACTORIALS_SIZE #define value !!!!!!!!!\n",from_where,size_needed);
    return 0;
  }

  // === precompute table if needed ===
  if ( ! is_table_init ) {
    int i;
    T[0] = 0;
    for (i = 1; i < TABLE_OF_LOG_FACTORIALS_SIZE; i++) 
      T[i] = T[i-1] + log( i );
    
    is_table_init = 1;
  }

  return &T;
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
  int N_AS = (c1 + c1 + c2) * (l1 + l1);// genotypes_dest has to be allocated (more than nb_genotypes cells)

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

  // printf("A1=%i N_AS=%i observed=%i\n",A1,N_AS,observed);
  
  if (alpha > beta)
    return 1;

  // ========== PRECOMPUTED LOG FACTORIALS TABLE ==========

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
  //fprintf(stderr,"l1-c3=%i f_b2=%i\n",l1-c3,f_b2);
  for (i = 0; i < 4; i++) {
    a_start = MYMAX(starts[i],a_min);
    a_end = MYMIN(ends[i],a_max);
    //        fprintf(stderr,"a_max=%i end=%i\n",a_max,ends[i]);
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
      //fprintf(stderr,"i=%i ai=%i b_start=%i b_end=%i pvalue=%e\n",i,ai,b_start,b_end,pvalue);
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

END_OF_C_CODE

