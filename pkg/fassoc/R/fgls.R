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
		MLE<-.emma.REMLE(y=Tped[,6],X=X, K=K2)  #phone~2K
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
        fgls.test=.FGLS(Y=y,X=cbind(1,data.matrix(sumG)),as.matrix(W),
			test="score",whichtest=c(FALSE,TRUE))
        pval=1-pchisq(fgls.test$T2,fgls.test$df)  #chi-square

  return(pval)
}


####################################################
#FGLS function was taken from the MixABEL -Yurii Aulchenko
#    http://www.genabel.org/packages/MixABEL         
#    We have put it here instead of using the MixABEL 
#    package to make it easier to install fassoc.
####################################################
#' Feasible GLS
#' 
#' Feasible Generalised Least Squares
#' 
#' @param Y dependent variable
#' @param X design matrix (including intercept, if necessary)
#' @param test test to be applied, one of 'wald', 'score' or 'robust'
#' @param whichtest which independent variables to be tested (set to 'TRUE')
#' @param W for GLS, inverse variance-covariance matrix, as such returned by 
#' GenABEL's polygenic(...)$InvSigma, or NULL for LS
#' 
#' @return List with elements 'beta' -- estimates fo the regression 
#' coefficients; 'V' -- variance covariance matrix for parameters 
#' estimates; 'T2' -- test statistics (distributed as Chi-squared under 
#' the null) for the testing of whichtest parameters; 'df' -- the number of 
#' degrees of freedom of the T2 test; 'tested' -- which parameters were 
#' tested with T2; 'meanColX' -- mean value of variable in columns of X;
#' 'n' -- length of Y (== height of X)
#'
#' @author Yurii Aulchenko
#'
.FGLS <- function(Y,X,W=NULL,test="wald",whichtest = c(FALSE,rep(TRUE,dim(X)[2]-1))) 
{
	# do checks
	if (!is(whichtest,"logical") || length(whichtest)!=dim(X)[2]) 
		stop("argument whichtest should be logical of length dim(X)[2]")
	if (dim(X)[1] != length(Y))
		stop("dimensions of X and Y do not match")
	if (any(is.na(Y)) || any(is.na(X))) {
		warning("missing data points in Y or X, dropping")
		IScomplete <- !is.na(Y) & complete.cases(X)
		Y <- Y[IScomplete]
		X <- X[IScomplete,]
		if (!is.null(W)) W <- W[IScomplete,IScomplete]
	}
	
	# precompute tXW
	if (is.null(W)) {
		tXW <- t(X)
	} else {
		tXW <- t(X) %*% W
	}
	#print(tXW)
	
	# estimate beta
	XpWXm1 <- ginv(tXW %*% X)
	XpWY <- tXW %*% Y
	betaW <- XpWXm1 %*% XpWY
	#print(XpWXm1)
	
	# estimate V
	if (test=="wald") {
		YmXB <- Y - X %*% betaW
		if (is.null(W)) {
			sigma2 <- as.vector(((t(YmXB) %*% YmXB)/(dim(X)[1] - dim(X)[2])))
		} else {
			sigma2 <- as.vector(((t(YmXB) %*% W %*% YmXB)/(dim(X)[1] - dim(X)[2])))
		}
		#print(sigma2)
		V <- sigma2*XpWXm1
	} else if (test=="score") {
		Xm <- X[,!whichtest,drop=FALSE]
		if (is.null(W)) {
			bt <- ginv(t(Xm) %*% Xm) %*% (t(Xm) %*% Y)
			YmXB <- Y - Xm %*% bt
			sigma2 <- as.vector(((t(YmXB) %*% YmXB)/(dim(Xm)[1] - dim(Xm)[2])))
		} else {
			bt <- ginv(t(Xm) %*% W %*% Xm) %*% (t(Xm) %*% W %*% Y)
			YmXB <- Y - Xm %*% bt
			sigma2 <- as.vector(((t(YmXB) %*% W %*% YmXB)/(dim(Xm)[1] - dim(Xm)[2])))
		}
		V <- sigma2*XpWXm1
	} else if (test=="robust") {
		
		YmXB <- Y - X %*% betaW
		if (is.null(W)) {
			sigma2vec <- as.vector(YmXB * YmXB)
		} else {
			sigma2vec <- as.vector(as.vector(t(YmXB) %*% W) * YmXB)
		}

		V <- XpWXm1 %*% (tXW %*% diag(sigma2vec) %*% X) %*% XpWXm1

	} else {
		stop("test not recognised")
	}
	
	# do the test
	if (sum(whichtest)>1) {
		Vinv <- ginv(V)
		Vbg <- Vinv[whichtest,whichtest]-Vinv[whichtest,!whichtest] %*% ginv(Vinv[!whichtest,!whichtest]) %*% Vinv[!whichtest,whichtest]
		T2 <- t(betaW[whichtest]) %*% Vbg %*% betaW[whichtest]
	} else if (sum(whichtest) == 1) {
		T2 <- betaW[whichtest]^2/diag(V)[whichtest]
	} else {
		stop("unreachable statement")
	}
	
	rownames(betaW) <- colnames(X)
	dimnames(V) <- list(colnames(X),colnames(X))
	out <- list(beta=betaW,V=V,T2=T2,df=sum(whichtest),tested=whichtest,
			meanColX = apply(X,FUN=mean,MARGIN=2), n = dim(X)[1])
	class(out) <- "FGLS"
	out
}