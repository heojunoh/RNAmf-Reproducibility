#----------------------------------------------------------------------------------------------#
#------INTERNAL FUNCTIONS
#----------------------------------------------------------------------------------------------#
voronoi <- function(X,D,y,Loo){
		nX <- dim(as.matrix(X))[1]
		nD <- dim(as.matrix(D))[1]
		d <- dim(as.matrix(D))[2]
	#-- LOO variances
		varloo <- Loo$sd^2
	#-- LOO errors
		e2loo <- (Loo$mean - y)^2
	#-- Voronoi cells tesselation
		grid <- expand.grid(as.matrix(X)[,1],as.matrix(D)[,1])
		h <- (grid[,1]-grid[,2])^2
		if(d > 1){
			for(i in 2:d){
				grid <- expand.grid(X[,i],D[,i])
				h <-  h +(grid[,1]-grid[,2])^2
			}
		}
		matker <- matrix(h,nrow=nX,ncol= nD)
		indice <- max.col(-(matker))	
	#-- adjustment
		aj <- (e2loo/varloo)[indice]
	return(aj)
}

adjustment <- function(var,adj){
	return(var*(1+adj))
}


CVapplyLevel1 <- function(indice,mymodel){
	return((CrossValidationMuFicokmAll(mymodel,indice))$CVerr[[1]])
}
CVvarapplyLevel1 <- function(indice,mymodel){
	return((CrossValidationMuFicokmAll(mymodel,indice))$CVvar[[1]])
}
CVapply1l <- function(indice,mymodel){
	return((CrossValidationMuFicokm(mymodel,indice))$CVerr)
}
CVvarapply1l <- function(indice,mymodel){
	return((CrossValidationMuFicokm(mymodel,indice))$CVvar)
}
CVapply <- function(indice,mymodel){
	return((CrossValidationMuFicokmAll(mymodel,indice))$CVerrall)
}
CVvarapply <- function(indice,mymodel){
	return((CrossValidationMuFicokmAll(mymodel,indice))$CVvarall)
}

varpredx <- function(x,model){
	ypred <- predict(model, newdata=x, type="UK")
	return(ypred$sd)
}

varpredxMH <- function(x,model){
	output <- - Inf
	if(min(x) <0){
		output <- - Inf
	}else{
		if(max(x) > 1){
			output <- - Inf
		}else{
			x <- t(as.matrix(x))
			ypred <- predict(model, newdata=x, type="UK")
			output <- ypred$sd
		}
	}
	return(output)
}
	
clmean <- function(icl,xech, nclass,modelseq){
		(cl <- kmeans(xech, nclass ))	
		ycl <- predict(modelseq, newdata=cl$centers, type="UK",cov.compute=TRUE)
		return(min(ycl$sd^2))
}

clnclass <- function(nclass,xech,modelseq,nclmean){
	ycldetmean <- apply(matrix(1:nclmean),1,clmean,xech, nclass,modelseq)
	return(mean(ycldetmean))
}


varcokpredx <- function(x,model){
	output <- - Inf
	if(max(x) > 1){
		output <- - Inf
	}else{
		if(min(x) < 0){
			output <- - Inf
		}else{
			x <- as.matrix(x)
			if(dim(x)[2]==1){
				x <- t(x)
			}
			ypred <- predict(model, newdata=x, type="UK")
			output <- ypred$sig2
		}
	}
	return(output)
}
	
CVapplyDelta <- function(indice,mymodel,rho){
	CVMFk <- CrossValidationMuFicokmAll(mymodel,indice)
	return(CVMFk$CVerrall-rho*CVMFk$CVerr[[1]])
}
CVvarapplyDelta <- function(indice,mymodel,rho){
	CVMFk <- CrossValidationMuFicokmAll(mymodel,indice)
	return(CVMFk$CVvarall-rho^2*CVMFk$CVvar[[1]])
}

clmeancok <- function(icl,xech, nclass,modelseq){
		(cl <- kmeans(xech, nclass ))	
		ycl <- predict(modelseq, newdata=cl$centers, type="UK",cov.compute=TRUE)
		return(min(ycl$sig2))
}

clnclasscok <- function(nclass,xech,modelseq,nclmean){
	ycldetmean <- apply(matrix(1:nclmean),1,clmeancok,xech, nclass,modelseq)
	return(mean(ycldetmean))
}

FuncLie <- function(x){
	return(x[,1])
}


SubstDesign <- function(PX2,PX1){
	d <- dim(as.matrix(PX2))[2]
	n2 <- dim(as.matrix(PX2))[1]
	n1 <- dim(as.matrix(PX1))[1]
	
	dist <- 0 
	for(i in 1:d){
		grid <- expand.grid(PX2[,i],PX1[,i])
		dist <- dist + (grid[,1]-grid[,2])^2
	}
	Matdist <- matrix(dist,n2,n1)
	indice <- max.col(-(Matdist))
	
	PX1 <- PX1[-indice,]
	PX1 <- rbind(PX1,PX2)
	return(list(PX = PX1,le = length(indice)))
}

NestedLHS <- function(names,design=NULL,n, k,nlevel, pop=100, gen=4, pMut=.1, dup=1, maxSweeps=2, eps=.1){
	
	if(identical(design,NULL)){design<- list() ; for(i in 1:nlevel){design[[i]] <- -1}}
	if(length(names)==1){names <- rep(names,nlevel)}
	if(length(pop)==1){pop <- rep(pop,nlevel)}
	if(length(gen)==1){gen <- rep(gen,nlevel)}
	if(length(pMut)==1){pMut <- rep(pMut,nlevel)}
	if(length(dup)==1){dup <- rep(dup,nlevel)}
	if(length(maxSweeps)==1){maxSweeps<- rep(maxSweeps,nlevel)}
	if(length(eps)==1){eps<- rep(eps,nlevel)}

	PX <- list()
	indices <- list()

	for(i in 1:nlevel){
		if(identical(names[i],"geneticLHS")){
			if(identical(design[[i]],-1)){
				PX[[i]] <- geneticLHS(n = n[i],k = k,pop = pop[i], gen = gen[i],pMut = pMut[i])
			}else{
				PX[[i]] <- design[[i]]
			}
		}
		if(identical(names[i],"randomLHS")){
			if(identical(design[[i]],-1)){
			PX[[i]] <- randomLHS(n = n[i],k = k)
			}else{
				PX[[i]] <- design[[i]]
			}
		}
		if(identical(names[i],"maximinLHS")){
			if(identical(design[[i]],-1)){
			PX[[i]] <- maximinLHS(n = n[i],k = k,dup = dup[i])
			}else{
				PX[[i]] <- design[[i]]
			}
		}
		if(identical(names[i],"optimumLHS")){
			if(identical(design[[i]],-1)){
			PX[[i]] <- optimumLHS(n = n[i],k = k,maxSweeps=maxSweeps[i], eps=eps[i])
			}else{
				PX[[i]] <- design[[i]]
			}
		}
	}
	for(i in (nlevel-1):1){
		SB <- SubstDesign(PX[[i+1]],PX[[i]])
		PX[[i]] <- SB$PX
		n <- dim(SB$PX)[1]
		indices[[i]] <- seq(n-SB$le+1,n,by=1)
	}
	for(i in 1:(nlevel-1)){
		
	}
	PX <- NestedDesign(PX[[1]], nlevel = nlevel , indices = indices)
	return(PX)
}
