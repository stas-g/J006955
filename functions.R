#v1.04 --> 27.10.14
#################MISSING MARKER IMPUTATION##########################
# imputes missing markers according to impirical allele frequency. can handly any type of snp coding (-1, 0, +1) or (0, 1, 2) etc. assumes markers are in columns.(25.04.14)
marker.imp <- function(geno){
	n <- nrow(geno)
	apply(geno, 2, FUN = function(x){
		na.ind <- is.na(x)
			if(sum(na.ind) == 0) return(x) else{
				entries <- unique(x[!na.ind])
				freq <- sapply(entries, FUN = function(a) sum(x == a & !na.ind)/sum(!na.ind))
				x[na.ind] <- sample(entries, sum(na.ind), replace = TRUE, prob = freq)		
				x
			}
		}	
	)
}

#####################REMOVE MONOMORPHIC MARKERS#####################

#removes monomorphic markers. assumes markers are ordered in columns, samples - in rows (08.04.14)
#- note add an option to pass custom missing value(s) 
#v3, 02.10.2015
mono.remove <- function(geno, ind = FALSE){
    ind_ <- apply(geno, 2, FUN = function(x){
      x <- x[!is.na(x)]
      length(unique(x)) == 1
    } )
	geno <- geno[, !ind_]
	if(ind == TRUE) attributes(geno)$index <- ind_
	geno
}



#################ENTRY.TRANSFORM###################################
#heavy duty entry.transform. subsetting xs is faster then using apply for each column (13.03.14)

entry.transform <- function(geno, entries = c(1, 0, -1)){
  aa <- entries[1]; ab <- entries[2]; bb <- entries[3]
  d <- which(dim(geno) == min(dim(geno)))
  ans <-apply(geno, d, FUN = function(x){
	x[x == 'AA'] <- aa
	x[x == 'AB'] <- ab
	x[x == 'BB'] <- bb
	x[x != aa & x != ab & x != bb] <- NA
	as.numeric(x)
	}
		)
	if(d == 1) ans <- t(ans)
		dimnames(ans) <- dimnames(geno)
			return(ans)
}



######################thin.markers##################################
#remove markers with more than a specified proportion of values missing. ex NA.REMOVE
thin.markers <- function(geno, p){
  n <- nrow(geno)
  freq <- apply(geno, 2, FUN = function(x){sum(is.na(x))/n})
  geno[, freq < p]
  
}

#################DRAWREP###################################
#v1.2 --> 22.02.2016

drawrep <- function(dat, avrg = FALSE){
	if(avrg == TRUE){
		ans <- apply(dat[,-1, drop = FALSE], 2, FUN = function(x) tapply(x, dat[,1], mean, na.rm = TRUE ))
		} else {
			ans <- apply(dat[,-1, drop = FALSE], 2, FUN = function(x) tapply(x, dat[,1], FUN = function(y){sample(y[!is.na(y)], 1 )} ))}
	return(ans)
}


#########UNITE TRAITS###############################

# utility function used in both 'unite.traits' and 'drawrep'. if a string has a numerical ending numstrip deletes it or, if num argument is set to TRUE, returns it. so for example numstrip('test20') returns 'test' and numstrip('test20', num = T) returns '20'. ok (01.04.2014)

numstrip <- function(string, num = F){
	r <- regexpr('[0-9]+$', string)
	if(num == T) return(regmatches(string, r))
	else{if(r[[1]] != -1) return(substring(string, 1, r[1] - 1))
		else return(string)
}	
}

# takes a data.frame or a matrix of traits as an argument, returns a data.frame of trait averages. i.e. if there are traits with names only differing in the end number (e.g. 'yield1', 'yield2'), they are averaged for each plant (the new trait name is stripped of the numbers, so the average of 'yield1' and 'yield2' would be called 'yield'). non repeating traits do not have to be numeric (i.e. can be factors or strings) and will be returned unchanged. - v1.1 (12.06.2014)

unite.traits <- function(pdata){
	nonum <- sapply(colnames(pdata), numstrip)
	uname <- unique(nonum)
	ind <- lapply(uname, FUN = function(x) x == nonum)
	pdata.av <- lapply(ind, FUN = function(i){
		if(sum(i) == 1) pdata[, i]
			else rowMeans(pdata[, i], na.rm = TRUE)
	})
		pdata.av <- as.data.frame(pdata.av, stringsAsFactors = FALSE)
		colnames(pdata.av) <- unlist(uname)
		pdata.av
}


#####################REMOVE BLANK/MISSING TRAITS#####################
#should work both on matrices and data.frames of numeric values. removes any column (a trait) that only contains 0's and NAs

blank.remove <- function(pheno){
	ind <- apply(pheno, 2, FUN = function(x){
		   sum(x == 0 | is.na(x)) != nrow(pheno)
  })
  pheno[, ind]
}

##########################AGENT 007#################################

# takes a number or a vector of numbers and an optional argument n (default: n = 3) and adds the required number of 0's in front of each number in a vector to buff them up to length n (v1.1 - 15.08.2014)
agent007 <- function(nums, n = 3){
  len <- nchar(nums)
  subs <- sapply(1 : n, FUN = function(i) paste(rep('0', i - 1), collapse = ''))
  nums[len < n] <- paste(subs[n + 1 - len[len < n]], nums[len < n], sep = '')
  nums
}

# takes a string of numbers of the form 'xxx123..' as an argument, applies agent007 to the end number and returns augmented strings: e.g. agent008('agent7', n = 3) --> 'agent007' (v1.1 - 15.08.2014)
agent008 <- function(x, n = 3){
  nums <- numstrip(x, num = TRUE)
  bodies <- numstrip(x)
  paste(bodies, agent007(nums, n), sep = '')
}

#---------------------------------------------------------------------------------------------------#
#			             BLUP                 					    #
#---------------------------------------------------------------------------------------------------#
#OK! - 28.02.14

#cv.folds is a helper function from the lars package
cv.folds <- function (n, folds = 10) 
{
    split(sample(1:n), rep(1 : folds, length = n))
}


#GEBV predictions using fitted rrBLUP model
#changed geno for new.geno, b for new.b to avoid confusion 09.04.14 - OK!

blup.pred <- function(fit, new.geno, new.b = NULL){
	if(is.null(new.b)){
		new.geno %*% as.numeric(fit$u) + as.numeric(fit$beta)
	} else{
			new.geno %*% as.numeric(fit$u) + as.numeric(fit$beta %*% t(new.b))
		}
	}

#------------------------------cross-validated rel. error for BLUP---------------------------

blup.cv <- function(pheno, geno, b = NULL, n = 10, folds = NULL, err = FALSE, parallel = FALSE){
	if(is.null(folds)){
		folds <- cv.folds(length(pheno), folds = n)
			} else{folds <- folds
			  n <- length(folds)
				}
			if(parallel == FALSE){
	pred <- 
		lapply(1 : n, FUN = function(i){
			test <- folds[[i]]
			fit <- mixed.solve(pheno[-test], X = b[-test, ], Z = geno[-test, ], method = 'REML')
			blup.pred(fit, geno[test, ], b[test, ])	
		})
			} else{
	pred <- 
		foreach(i = 1 : n)  %dopar% {
			test <- folds[[i]]
			fit <- mixed.solve(pheno[-test], X = b[-test, ], Z = geno[-test, ], method = 'REML')
			blup.pred(fit, geno[test, ], b[test, ])
	}
		}
			pred <- unlist(pred)[order(unlist(folds))]
				if(err == TRUE){
					pheno <- pheno[unlist(folds)][order(unlist(folds))]
					return(mean((pred - pheno)^2)/var(pheno))
				} else{
					return(pred) 
					 }
}


################################CONTIG_TO_FASTA#############################################
# 20.05.14
# file_: path to a .csv file with chunks (first column = names, second column = chunks in the format XXXX[X/X]XXXX
# ref: if ref = TRUE, take first snp value (ref), otherwise take second (alt) when gluing stuff together
# destination (optional): if you don't provide destination, then the fasta file (in txt format) will be saved in your work directory under the name 'name_of_your_original_csv_file_fasta.txt'. alternatively you can provde a full path\name of your new fasta file in the format 'path_to_file\\name_of_fasta.txt' or 'name_of_fasta.txt' (then 'name_of_fasta.xtx' is saved to working directory).
#info: if TRUE creates a data frame containint alt, ref values of snp, its position in the chunk and the downstream and upstream regions. if FALSE default behaviour, i.e. create and save .fasta file.

contig_to_fasta <- function(file_, ref = TRUE, destination = NULL, info = FALSE,...){
	snp <- read.csv(file_, stringsAsFactors = FALSE, ...)
	r <- gregexpr('[[:punct:]]', snp[, 2])
	m <- regmatches(snp[, 2], r, invert = TRUE)
	name <- snp[, 1]
		if(info == FALSE){
			if(ref == TRUE) snp.ind <- 3 else snp.ind <- 2
			new_file <- regmatches(file_, regexec('[\\\\|\\/]([^\\\\|\\/]+)\\..+', file_))[[1]][2]
			new_file <- paste(new_file, '_fasta', '.txt', sep = '')
			if(!is.null(destination)) new_file <- destination
			write(paste('>', name[1], sep = ''), file = new_file)
			write(paste(m[[1]][-snp.ind], collapse = ''), file = new_file, append = TRUE)	
			invisible(mapply(m[-1], name[-1], FUN = function(seqn, n){
				write(paste('>', n, sep = ''), file = new_file, append = TRUE)
				write(paste(seqn[-snp.ind], collapse = ''), file = new_file, append = TRUE)
			}))
		} else{
				snp.info <- data.frame(pos = sapply(m, FUN = function(x) nchar(x[[1]]) + 1), ref = sapply(m, '[', 2), alt = sapply(m, '[', 3),  downstream = sapply(m, '[', 1), upstream = sapply(m, '[', 4) )
				rownames(snp.info) <- name
				return(snp.info)
			}
}



#-------------------\\\\\\\\\\\\\\\\\\\\\\\\\\\GENERAL-PURPOSE FUNCTIONS and PLOTTING/////////////////////////-----------------------#

#returns x if x is positive and 0 otherwise -- vectorised (v1 -- 06.10.2014)
maxz <- function(x){sapply(x, FUN = function(z) ifelse(z > 0, z, 0))}

#return indices (default) or a logic vector (if logic = TRUE) of intersecting values between vectors x and y (relative to x). (v2 17.10.2014)
intersect.ind <- function(x, y, logic = FALSE){
  vals <- intersect(x, y)
  ans <- lapply(vals, FUN = function(v) which(x == v))
  ans <- sort(unlist(ans))
		if(logic == FALSE){return(ans)}else{
			ans.l <- rep(F, length(ans))
			ans.l[ans] <- TRUE
			return(ans.l)
  }
}

#display column numbers for variables
trait.index <- function(dat){
data.frame(trait = colnames(dat))
}


# plot columns of dat as lines on one picture frame. colours provided optionally (v1.1 - 06.06.14)
boo.plot <- function(dat.list, names_, col_ = NULL, ylim_ = NULL, legend_ = TRUE, ...){
	n <- length(dat.list)
	if(is.null(col_)) col_ <- rainbow(n)
	if(is.null(ylim_)) ylim_ <- c(min(unlist(dat.list)), max(unlist(dat.list)))
	plot(dat.list[[1]], type = 'l', col = col_[1], ylim = ylim_, ...)
	for(i in 2 : n) points(dat.list[[i]], type = 'l', col = col_[i])
	if(legend_ == TRUE) legend('top', names_, lty = rep(1, n), col = col_)
}    
      
    

# plot a histogram of data contained in each element of list dat.list all on one frame with a legend. colours provided optionally (v1 - 03.06.14)
boo.hist <- function(dat.list, col_ = NULL, names_, ...){
	require(scales)
	n <- length(dat.list)
	if(is.null(col_)) col_ <- rainbow(n)
	den <- unlist(lapply(dat.list, FUN = function(x) hist(x, plot = FALSE)$density))
	plot(c(range(unlist(dat.list))), c(range(den)), type = 'n', ylab = 'Density', ...) 
	for(i in 1 : n) hist(dat.list[[i]], freq = FALSE, col = alpha(col_[i], 0.5), add = TRUE)
	legend('topright', pch = rep(15, n), legend = names_, col = col_)
}	
	

#accepts both vectorised obs and pred as well as matrices, in which case rsq, or mse, is calculated pairwise for columns.
# v3 - 28.10.2014	
calc.mse <- function(obs, pred, rsq = FALSE){
						if(!is.matrix(obs) | !is.matrix(pred)){obs <- as.matrix(obs); pred <- as.matrix(pred)}
						rs <- (obs - pred)^2						
						if(rsq == FALSE){colMeans(rs, na.rm = TRUE)}else{
							meanmat <- matrix(rep(colMeans(obs, na.rm = TRUE), nrow(obs)), nrow = nrow(obs), byrow = TRUE)
							1 - colSums(rs, na.rm = TRUE)/colSums((obs - meanmat)^2, na.rm = TRUE)
						}
							}

#--------------

# v1.2 - 18.11.2014	
fold.mse <- function(obs, pred, folds, rsq = FALSE){
		rsqs <- lapply(folds, FUN = function(x){
				if(is.vector(obs)){
					obs <- obs[x]
					pred <- pred[x]
						}else{
							obs <- obs[x, ]
							pred <- pred[x, ]
					}
			calc.mse(obs, pred, rsq = rsq)})
	rsqs <- do.call(rbind, rsqs)	
	sdev <- apply(rsqs, 2, sd)
	m <- colMeans(rsqs)
	rsq <- calc.mse(obs, pred, rsq = rsq)
	list(rsq = rsq, mean = m, sdev = sdev)
}


# LASSO.BETAS (06.08.2014)
lasso.betas <- function(cvlasso.obj, scale = FALSE, nonzero = TRUE, lambdamin = FALSE, names = FALSE){
	require('glmnet')
		if(lambdamin == TRUE){
			b.parse <- coef(cvlasso.obj, s = 'lambda.min')
			betas <- as.vector(b.parse)[-1]
			names(betas) <- rownames(b.parse)[-1]
		} else{
				b.parse <- coef(cvlasso.obj, s = 'lambda.1se')
				betas <- as.vector(b.parse)[-1] 
				names(betas) <- rownames(b.parse)[-1]
			}
					betas <- betas[order(abs(betas), decreasing = TRUE)]
					if(nonzero == TRUE) betas <- betas[betas != 0]
					if(scale == TRUE) betas <- abs(betas)/sum(abs(betas))
					if(names == TRUE) return(names(betas))
					else return(betas)
	}


#transofrm each column of a matrix/df into a factor; output as a df.
#v1 --> 15.01.2015
	
columns.to.factors <- function(mat, ord = FALSE){
			levs <- unique(as.vector(mat))
			ans <- as.data.frame(lapply(as.data.frame(mat), factor, levels = levs, ordered = ord))
			if(!is.null(dimnames(mat))) dimnames(ans) <- dimnames(mat)
			ans
	}
	


#given a matrix/df with plant names as rownames (userow = TRUE) or as a first column (userow = FALSE), normalise and scale each sub-population separately, then stick back together 
#v1 --> 16.10.2014	

subscale <- function(x, userow = TRUE, ...){	
	if(userow == TRUE){
		#ind <- as.factor(numstrip(rownames(x)))
		ind <- numstrip(rownames(x))
		ind <- factor(ind, levels = unique(ind))
		ans <- as.data.frame(do.call(rbind, by(x, ind, scale, ...)))
		return(ans)
			} else {
				ind <- numstrip(x[, 1])
				ind <- factor(ind, levels = unique(ind))
				ans <- as.data.frame(do.call(rbind, by(x[, -1], ind, scale, ...)))
				ans <- data.frame(name = x[, 1], ans)
				return(ans)
			}
}


#v2 --> 29.10.2015	
rf.importance <- function(x, n = NULL, order_ = TRUE, k = 1, ...){
	require(randomForest)
	imp.score <- importance(x, ...)[, k]
	if(order_) imp.score <- sort(imp.score, decreasing = TRUE) 
	if(is.null(n)) n <- length(imp.score)
	imp.score[1 : n]
}
