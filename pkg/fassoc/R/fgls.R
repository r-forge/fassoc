fgls<-function (Tped,n1,n2) {

        nu.ped=as.matrix(Tped[Tped[,1]<=n1,])        
        #creat kinship matrix of nuclear familes
        K=lapply(split(as.data.frame(nu.ped[,1:6]),nu.ped[,1]),function (adata) {
            if (nrow(adata)>1){
              k=.kinship(adata[,2],adata[,4],adata[,3])} else{k=as.matrix(0.5)}
              return(k)
        })
		K=bdiag(K); K2=2*K
		K2=bdiag(K2,diag(n2))
		X=matrix(rep(1,nrow(Tped))) #no covariates
		MLE<-emma.REMLE(y=Tped[,6],X=X, K=K2)  #phone~2K
		V0=(MLE$vg)*K2+(MLE$ve)*diag(nrow(K2)) 

        G=Tped[,-(1:6)]   #genotype matrix 
        G=G[,colSums(G)!=0]  #remove non-var snps        
        G=apply(G, 2, function (x) {
          p=sum(x)/(2*length(x))
          w=(x-2*p)/sqrt(2*p*(1-p))
          return(w)
        })
        sumG=rowSums(G)

        y=data.matrix(Tped[,6])

        ##FGLS
        W=solve(V0)
        fgls.test=FGLS(Y=y,X=cbind(1,data.matrix(sumG)),as.matrix(W),
			test="score",whichtest=c(FALSE,TRUE))
        pval=1-pchisq(fgls.test$T2,fgls.test$df)  #chi-square

  return(pval)
}