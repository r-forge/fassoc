\name{ped.simu}
\alias{ped.simu}

\title{
	Simulation Function 
}
\description{
	This function simulates both phenotypes and genotypes of a combination of family and unreated data. 
}
\usage{
ped.simu(hap = hap, pos = pos, nfa = 100, nun = 400, h2 = 0.03,
           rfreq = 0.02, hrisk = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hap}{
	haplotype information
}
  \item{pos}{
	
}
  \item{nfa}{
Number of nulcear families.
}
  \item{nun}{
Number of unrelated individuals
}
  \item{h2}{
The proportion of variation that can be explained by the simualted genetic region.
}
  \item{rfreq}{
The rare variants threshold of MAF used in the simulation
}
  \item{hrisk}{
The proportion of risk haplotypes in the haplotype pool. 
}
}

\value{
It returns a pedigree matrix ready for use by functions vc.score and fgls in this package.
}
\references{
		T Feng, RC Elston and X Zhu. 2011. Detecting rare and common variants for complex traits: sibpair and odds ratio weighted sum statistics (SPWSS, ORWSS). Genet Epidemiol, 35: 398-409.
}
\author{
Xuefeng Wang
}


\seealso{
vc.score
}
\examples{
		data(example)
		Tped<-ped.simu(hap,pos,h2=0.02)
		vc.score(Tped,100,400)
		fgls(Tped,100,400)	
}
