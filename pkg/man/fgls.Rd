\name{fgls}
\alias{fgls}

\title{
Feasible Generalized Least Squares
}
\description{
The function is based on the same function in "MixABEL" package. 
		It performs marker set based association based on the Feasible Generalized Least Squares (FGLS) model. 
}
\usage{
fgls(Tped, n1, n2)
}

\arguments{
  \item{Tped}{
Pedigree data matrix. The first six columns are the Family ID (FID), Inidvidual ID (IID),  Paternal ID (FAT),  Maternal ID (MAT), Sex (SEX), Phenotype (pheno), respectively.
							SNP genotypes (column 7 onwards) are coded 0/1/2 (minor allele count). 	
}
  \item{n1}{
Number of nuclear families.			
}
  \item{n2}{Nunber of unrelated individuals.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Returns p-value
}
\references{
Li et al. 2011. A rapid generalized least squares model for a genome-wide quantitative trait association analysis in families. Hum Hered. 71:67-82.
}
\author{
Xuefeng Wang
}



\seealso{
vc.score
}
\examples{
		data(example)
		fgls(Tped,100,400)
}
