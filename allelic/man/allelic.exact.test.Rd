\name{allelic.exact.test}
\alias{allelic.exact.test}
\title{Fast Unbiased Exact Allelic Test}
\description{This is the implementation in R+C of a new association
	test described in "A fast, unbiased and exact allelic exact test
	for case-control association studies" (Submitted).
	It appears that in most cases the classical chi-square test used
	for testing for allelic association on genotype data is biased.
	Our test is unbiased, exact but fast throught careful
	optimization.
}
\usage{
 allelic.exact.test(d0,d1,d2,h0,h1,h2)
}
\arguments{

  Takes the 2x3 contingency table on which to compute the test

  \tabular{rlll}{

    \tab aa \tab aA \tab AA \cr
  [case (diseased)] \tab d0  \tab d1 \tab d2 \cr
  [control(healthy)] \tab h0  \tab h1 \tab h2 \cr
  }

  
  \item{d0}{nb of first homozygous among cases}
  \item{d1}{nb of heterozygous among cases}
  \item{d2}{nb of second homozygous among cases}
  \item{h0}{nb of first homozygous among controls}
  \item{h1}{nb of heterozygous among controls}
  \item{h2}{nb of second homozygous among controls}
  
}

\value{return the p-value of the test, or -1 if the sum of all cells in
  table is greater than TABLE\_OF\_LOG\_FACTORIALS\_SIZE, a C symbol defined
  in src/newallelic.c
}

\examples{
 allelic.exact.test(160,80,60,160,160,30)
}

\seealso{\code{\link{chisq.test}}, \code{\link{fisher.test}}}
\keyword{univar}
