\name{vc.score}
\alias{vc.score}

\title{
Perform Score Test
}
\description{
This function performs marker set based association (score test) based on the variance component model.
}
\usage{
vc.score(Tped, n1, n2, app.method = "davi")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tped}{
Pedigree data matrix. The first six columns are the Family ID (FID), Inidvidual ID (IID),  
Paternal ID (FAT),  Maternal ID (MAT), Sex (SEX), Phenotype (pheno), respectively.
SNP genotypes (column 7 onwards) are coded 0/1/2 (minor allele count). 	
}
  \item{n1}{
Number of nuclear families.	
}
  \item{n2}{
 Nunber of unrelated individuals.
}
  \item{app.method}{
%%     ~~Describe \code{app.method} here~~
}
}

\value{
Returns p-value
}
%% \references{
%% ~put references to the literature/web site here ~
%% }
\author{
Xuefeng Wang
}


\seealso{
fgls
}
\examples{
		data(example)
		vc.score(Tped,100,400)
}

