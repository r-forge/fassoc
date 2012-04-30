vc.score<-function (Tped,n1,n2,app.method="davi") {


		nu.ped=as.matrix(Tped[Tped[,1]<=n1,])        
		#creat kinship matrix of nuclear familes
		K=lapply(split(as.data.frame(nu.ped[,1:6]),nu.ped[,1]),function (adata) {
			if (nrow(adata)>1){
			  k=.kinship(adata[,2],adata[,4],adata[,3])} else{k=as.matrix(0.5)}
			  return(k)
		})
		K=bdiag(K); K2=2*K
		K2=bdiag(K2,diag(n2))

		gc()
		X=matrix(rep(1,nrow(Tped))) #no covariates
		MLE<-emma.REMLE(y=Tped[,6],X=X, K=K2)  #phone~2K
		V0=(MLE$vg)*K2+(MLE$ve)*diag(nrow(K2)) 

		## calculate A matrix (gene similarity matrix)
		#get standardized genotype matrix W
		G=Tped[,-(1:6)]   #genotype matrix
		G=G[,colSums(G)!=0]  #remove non-var snps
		W=apply(G, 2, function (x) {
		  p=sum(x)/(2*length(x))
		  w=(x-2*p)/sqrt(2*p*(1-p))
		  return(w)
		})
		A=tcrossprod(W)/ncol(W)

		##Score test
		y=data.matrix(Tped[,6]) #quantitative phenotype
		y=y-mean(y)
		#beta=0
		iV0=solve(V0)  #inverse of V0
		U=t(y)%*%(0.5* (VS<-iV0 %*% A) %*% iV0 )%*%y  #score statistic
		#
		if (app.method=="davi" ){
		#Davies
		V.eig<-eigen(V0)
		V.sqrt=V.eig$vectors %*% diag(sqrt(V.eig$values)) %*% solve(V.eig$vectors)
		D=V.sqrt%*%(0.5*(VS%*%iV0))%*%V.sqrt
		si=eigen(D,only.values=T)$values
		p_val=davies(as.numeric(U),si,rep(1, length(si)))$Qq

		} else {
		##SATTER
		#EU=0.5*(trace(iV0%*%A))=0.5*(sum(diag(iV0%*%A)))=0.5*(sum(iV0*t(A))) 
		EU=0.5*(sum(diag(VS)))
		#varU=0.5*(sum(diag(iV0%*%A%*%iV0%*%A)))=0.5*(sum(diag(VS%*%VS)))
		Itt=0.5*sum((VS)*t(VS))
		#get partical information
		#Itp=0.5*trace(VS%*%iV0%*%K2)=0.5*sum(diag(VS%*%iV0%*%K2))
		Itp=0.5*sum(VS*t(VP<-iV0%*%K2))
		Ite=0.5*sum(VS*t(iV0))#=0.5*trace(VS%*%iV0) 
		Ipp=0.5*sum(VP*t(VP))        
		Ipe=0.5*sum(VP*t(iV0)) #=0.5*sum(diag(VP%*%iV0))=0.5*sum(diag(iV0%*%K2%*%iV0))
		Iee=0.5*sum(iV0*t(iV0))
		Itn=matrix(c(Itp,Ite),nrow=1)
		Inn=matrix(c(Ipp,Ipe,Ipe,Iee),2)
		It=Itt-Itn%*%solve(Inn)%*%t(Itn)#partical information to replace Itt/VarU
		#Gamma approximation
		a=EU^2/It
		b=It/EU 
		p_val=1-pgamma(as.numeric(U),shape=a, scale=b)
	}	
  return(p_val)

}