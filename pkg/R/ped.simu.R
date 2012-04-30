ped.simu<- function (hap=hap,pos=pos,nfa=100,nun=400,h2=0.03,rfreq=0.02,hrisk=0.1) {
		#rfreq      #rare variants threshold
		#hrisk		#proportion of risk haplotypes

		pos$SNP=1:nrow(pos)	
		rsnp<-pos[(pos$FREQ1<=rfreq|pos$FREQ2<=rfreq),]$SNP   #get rare SNP no.
		index<-numeric(nrow(hap))  #index of casual haplotype
		snp.index<-numeric(ncol(hap)-2) #index of casual rare snps
		while (sum(index)<(nrow(hap)*hrisk)){
		r1=sample(rsnp,1)  # a random rare variants
		snp.index[r1]=1
		if (pos$FREQ1[r1]<pos$FREQ2[r1]){index[hap[,r1+2]==1] =1} else {
			index[hap[,r1+2]==2] =1}
		}
		riskp=sum(index)/nrow(hap) #haplotype risk

		#sample haplotypes for family
		sample.index=sample(1:nrow(hap),(nfa*4+nun*2),replace=T)
		pheno=index[sample.index] # risk index for the sampled haplo
		hap.sample<-hap[sample.index, c(-1,-2)]

		#transfrom 1,2 coding (1=the derived allele,) to 0,1 (1=Minor allele)
		N<-nrow(hap.sample); halfN=N/2
		maf_one_col<-function (acol) {
		  cval<-numeric(N)
		  if (sum(i<-acol==1)<halfN)  cval[i]<-1
		  else if (sum(i<-acol==2)<halfN) cval[i]<-1 
		  cval
		}
		haps<-apply(hap.sample,2,maf_one_col)
		#haplos for families
		hapsF<-cbind(pheno=pheno[1:(nfa*4)], haps[1:(nfa*4),])
		#haplos for unrelated
		hapsU<-cbind(pheno=pheno[(nfa*4+1):(nfa*4+nun*2)], haps[(nfa*4+1):(nfa*4+nun*2),])	

		##Family genotype
		hapspar<-hapsF[seq(1,nfa*4,2),]+ hapsF[seq(2,nfa*4,2),]  #sum up odds and even rows
		ped.par<-cbind(FID=rep(1:nfa,each=2), IID=rep(c(21,22),nfa), FAT=0, MAT=0, hapspar)
		#get no. of children from zero-truncated poisson
		rtpois <- function(N, lambda) qpois(runif(N, dpois(0, lambda), 1), lambda)
		nkids=rtpois(nfa,2)
		ped.kid<-list()
		for (i in 1:nfa) {
		m<-1+(i-1)*4  # every 4 rows
			#randomly choose one from the first 2 rows
			random.fat=hapsF[sample(c(m,m+1),nkids[i],replace=T),]
			#choose one from second 2 rows
			random.mat=hapsF[sample(c(m+2,m+3),nkids[i],replace=T),]
			if (nkids[i]==1) random.family=matrix(random.fat+random.mat,nrow=1) 
			else random.family=random.fat+random.mat #combind haplo to get geno
		random.family=cbind(FID=i, IID=1:nkids[i], FAT=21, MAT=22, random.family)
		colnames(random.family)=c("FID","IID","FAT","MAT","pheno", paste("s",1:(ncol(hapsF)-1),sep="")) 
		ped.kid[[i]]<-random.family
		}
		ped2=do.call(rbind,ped.kid) #collapse dataframes
		#Merge parents and kids
		colnames(ped.par)=colnames(ped2)
		ped<-rbind(ped.par,ped2)  #merge
		ped<-ped[order(ped[,1]),] #sort according to FID 

		##Unrelated genotype
		hapsUnr<-hapsU[seq(1,nun*2,2),]+ hapsU[seq(2,nun*2,2),]  #sum up odds and even rows
		ped.unr<-cbind(FID=(nfa+1):(nfa+nun), IID=1, FAT=0, MAT=0, hapsUnr)
		colnames(ped.unr)<-colnames(ped)

		##
		Tped<-rbind(ped,ped.unr)
		Tped<-cbind(Tped[,1:4], SEX=0,Tped[,-c(1:4)]) #add sex

		##Generate quantitative phenotype
		Tped<-Tped[order(Tped[,1],Tped[,3]),] #sort to make sure 1st 2 rows are parents 
		Tpedh<-Tped[,1:6] #only include head of Tped file
		h2=h2 #heritability
		beta=1.5 #effect size
		g=beta*(Tpedh[,6]-2*riskp)	 #g=beta*centered.x
		#var.e=var(g)*(1/h2-1)
		var.e=beta^2*2*riskp*(1-riskp)*(1/h2-1) #theoretical value
		vp=var.e*(37/(37+60)); ve=var.e*(60/(37+60)); # 37 polygenic and 60% error
		e=rnorm(length(g),mean=0,sd=sqrt(ve))  #get random error
		p<-do.call(c, lapply(split(as.data.frame(Tpedh),Tpedh[,1]),function (x) {
		  #x=split(as.data.frame(Tpedh),Tpedh[,1])[[1]]
		  if (nrow(x)==1){poly=rnorm(1,mean=0,sd=sqrt(vp)) }
		  else { rn=rnorm(nrow(x),mean=0,sd=sqrt(vp));  #random normals
		  poly=rbind(cbind(diag(2),matrix(0,2,nrow(x)-2)), 
			cbind(matrix(0.5,nrow(x)-2,2), diag(nrow(x)-2)*(1/sqrt(2)))) %*% rn
			#pk=(pm+pf)/2+(1/sqrt(2)*N(0,vp)
		 }; return(poly) }))
		Pheno=g+p+e #u=0 
		Tped[,6]=Pheno
  return(Tped)
}