### FLXSA
### This file contains code just for the core class FLXSA
# $Id: FLXSA.R,v 1.54 2008/05/07 01:26:34 imosqueira Exp $

### class ######################################################################
validFLXSA.control <- function(object){
	if (object@tol <= 0)
		return("value of tol must be > 0")
	if (object@maxit <= 0)
		return("value of maxit must be > 0")
	if (object@min.nse < 0)
		return("value of min.nse must be > 0")
	if (object@fse < 0)
		return("value of fse must be > 0")
	if (object@rage < -1)
		return("value of rage must be >= -1")
	if (object@qage < 0)
		return("value of qage must be >= 0")
	if (object@shk.yrs < 0)
		return("value of shk.yrs must be >= 0")
	if (object@shk.ages < 0)
		return("value of shk.ages must be >= 0")
	if (object@window <= 5)
		return("value of window must be >= 5")                               

	# Everything is fine
	return(TRUE)
}

setClass("FLXSA.control",
	representation(
		tol          ="numeric",
		maxit        ="integer",
		min.nse      ="numeric",
		fse          ="numeric",
		rage         ="integer",
		qage         ="integer",
		shk.n        ="logical",
		shk.f        ="logical",
		shk.yrs      ="integer",
		shk.ages     ="integer",
		window       ="integer",
  	tsrange      ="integer",
  	tspower      ="integer",
  	vpa          ="logical"),
  prototype=prototype(
    tol          =as.double(10e-10),
  	maxit        =as.integer(30),
  	min.nse      =as.double(0.3),
  	fse          =as.double(0.5),
  	rage         =as.integer(-1),
  	qage         =as.integer(10),
  	shk.n        =TRUE,
  	shk.f        =TRUE,
  	shk.yrs      =as.integer(5),
  	shk.ages     =as.integer(5),
  	window       =as.integer(100),
  	tsrange      =as.integer(20),
  	tspower      =as.integer(3),
  	vpa          =FALSE),
  validity=validFLXSA.control
)

setValidity("FLXSA.control", validFLXSA.control)
remove(validFLXSA.control)	# We do not need this function any more

## FLXSA ######################
validFLXSA <- function(object){
	# All FLQuant objects must have same dimensions
	Dim <- dim(object@stock.n)
	if (!all(dim(object@harvest) == Dim))
		return("n and f arrays must have same dimensions")
	# Everything is fine
	return(TRUE)
}

setClass("FLXSA",
  contains='FLAssess',
	representation(
		survivors="FLQuant",
		se.int   ="FLQuant",
		se.ext   ="FLQuant",
		n.fshk   ="FLQuant",
		n.nshk   ="FLQuant",
		var.fshk ="FLQuant",
		var.nshk ="FLQuant",
		q.hat    ="FLQuants",
		q2.hat   ="FLQuants",
		diagnostics="data.frame",
		control  ="FLXSA.control"),
	prototype=prototype(
		survivors=FLQuant(),
		se.int   =FLQuant(),
		se.ext   =FLQuant(),
		n.fshk   =FLQuant(),
		n.nshk   =FLQuant(),
		var.fshk =FLQuant(),
		var.nshk =FLQuant(),
		q.hat    =FLQuants(),
		q2.hat   =FLQuants(),
		diagnostics=new("data.frame"),
		control  =new("FLXSA.control")),
	validity=validFLXSA
)

setValidity("FLXSA", validFLXSA)
remove(validFLXSA)	# We do not need this function any more


### Methods #############################################################
FLXSA <- function(stock, indices, control=FLXSA.control(), desc, diag.flag=TRUE){
  
  if (!is(diag.flag,"logical")) diag.flag=FALSE
  Call <- match.call()
	if (!inherits(stock, "FLStock"))
		stop("stock must be an 'FLStock' object!")
	if (inherits(indices, "FLIndex"))
    indices<-FLIndices(indices)
	if (!inherits(indices, "FLIndices"))
  	stop("indices must be an 'FLIndices' object!")
  for (i in 1:length(indices))
     {
     if (is.na(indices[[i]]@range["startf"]) || is.na(indices[[i]]@range["endf"]))
  	     stop(paste("Must supply startf & endf for range in FLIndex",i))

      if (!all(names(indices[[i]]@range) == c("min","max","plusgroup","minyear","maxyear","startf","endf")))
         stop("Range must have names 'min','max','plusgroup','minyear','maxyear','startf','endf'")

      indices[[i]]@range["min"]     <- max(indices[[i]]@range["min"], dims(indices[[i]])$min, stock@range["min"])
      indices[[i]]@range["max"]     <- min(indices[[i]]@range["max"], dims(indices[[i]])$max, stock@range["max"])
      indices[[i]]@range["minyear"] <- max(indices[[i]]@range["minyear"], dims(indices[[i]])$minyear, stock@range["minyear"])
      indices[[i]]@range["maxyear"] <- min(indices[[i]]@range["maxyear"], dims(indices[[i]])$maxyear, stock@range["maxyear"])

      age <- indices[[i]]@range["min"]:indices[[i]]@range["max"]
      year<- indices[[i]]@range["minyear"]:indices[[i]]@range["maxyear"]
      
      indices[[i]]<-trim(indices[[i]],age=age,year=year)
      }

  	if (!is(control, "FLXSA.control"))
  		stop("control must be an 'FLXSA.control' object!")
  	if (!validObject(stock))
  		stop("stock is not valid!")
  	if (!validObject(indices))
  		stop("FLIndices is not valid!")
  	if (!validObject(control))
  		stop("control is not valid!")

    ## check range
    if ("minage"    %in% names(stock@range)) minage <- stock@range["minage"]   else
    if ("min"       %in% names(stock@range)) minage <- stock@range["min"]      else
    if ("minquant"  %in% names(stock@range)) minage <- stock@range["minquant"] else  stop("'minage' not found in range")

    if ("maxage"    %in% names(stock@range)) maxage <- stock@range["maxage"]   else
    if ("max"       %in% names(stock@range)) maxage <- stock@range["max"]      else
    if ("maxquant"  %in% names(stock@range)) maxage <- stock@range["maxquant"] else  stop("'maxage' not found in range")
    
    if ("plusgroup" %in% names(stock@range)) stock@range["plusgroup"] <- maxage

    if (maxage<minage | stock@range["maxyear"]<stock@range["minyear"]) stop("Error in range")
    if (is.na(stock@range["plusgroup"])) stop("Plus Group must be specified")
    
    ## trim stock
    stock   <- trim(stock, year=stock@range["minyear"]:stock@range["maxyear"])

    if (!is.na(stock@range["plusgroup"]) & stock@range["plusgroup"] < dims(stock@catch.n)$max)
       stock   <- setPlusGroup(stock, stock@range["plusgroup"])

	  stock@m <- stock@m[as.character(minage:maxage), , , ,]

	  if (all(is.na(stock@catch.n))){
       stop("catch.n is not available")
       }

    stock@stock.n <- new("FLQuant")
    stock@harvest <- new("FLQuant")

#    for (i in 1:length(indices))
#      {
#      start=max(dims(indices[[i]])$minyear,dims(stock@m)$maxyear-control@window+1)
#      end  =min(dims(indices[[i]])$maxyear,dims(stock@m)$maxyear)
#      indices[[i]]<-window(indices[[1]],start=start,end=end)
#      }

     fqs<-function(assess) {
            assess@index <- new("FLQuants", lapply(assess@index,FLQuant))
            assess@index.hat <- new("FLQuants", lapply(assess@index.hat,FLQuant))
            assess@index.var <- new("FLQuants", lapply(assess@index.var,FLQuant))
            assess@index.res <- new("FLQuants", lapply(assess@index.res,FLQuant))
            assess@q.hat <- new("FLQuants", lapply(assess@q.hat,FLQuant))
            assess@q2.hat <- new("FLQuants", lapply(assess@q2.hat,FLQuant))
            if (validObject(assess))
                return(assess)
            else stop("not valid")
        }

    iters.stock  <-max(unlist(qapply(stock, function(x) dims(x)$iter)))
    iters.indices<-max(unlist(lapply(indices@.Data, function(x) max(unlist(qapply(x, function(x2) dims(x2)$iter))))))
    
    if ((iters.stock>1 && iters.indices>1) && diag.flag)
       return("Multiple iters only allowed if diag.flag=FALSE")

    if(!diag.flag) 
       {
       res<-.Call("FLXSA", iter(stock,1), lapply(indices, iter,1), control, FALSE)

       iters<-max(iters.stock,iters.indices)
       if (iters>1)
           {
           res@stock.n<-propagate(FLQuant(res@stock.n@.Data),iters)
           res@harvest<-propagate(FLQuant(res@harvest@.Data),iters)

           for (i in as.character(2:iters))
              {
              res.              <-.Call("FLXSA", iter(stock,i), lapply(indices, iter,i), control, FALSE)
              #res.              <-FLXSA(iter(stock,i), lapply(indices, iter,i), control, diag.flag)
              iter(res@stock.n,i)<-FLQuant(res.@stock.n@.Data)
              iter(res@harvest,i)<-FLQuant(res.@harvest@.Data)
              }
           }       

       res@harvest@units <- "f"

       return(res)
       }
       
    res <-fqs(.Call("FLXSA", stock, indices, control, diag.flag))

    if (class(res) != "FLXSA") return(res)
	     res@call <- as.character(Call)

	  if (!missing(desc)) res@desc <- as.character(desc)

    ## put wts amd nhats into a data frame
    df <- as.data.frame(res@wts)
    df1 <- (df[4])
    df1[df1 >= 1, 1] <- paste("index", df1[df1 >= 1, 1])
    df1[df1 == -1, 1] <- "fshk"
    df1[df1 == -2, 1] <- "nshk"
    df <- cbind(df[-4], df1)

  	names(df) <- c("w", "nhat", "yrcls", "age", "year", "source")
  	for(i in 1:length(indices)){
        v <- paste("index", i)
        df$source[df$source==v] <- slot(indices[[i]],"name")   
    }

    wts.df <-df[,c(4,5,1,6)]
    names(wts.df)[3] <-"data"
    nhat.df<-df[,c(4,5,2,6)]
    names(nhat.df)[3]<-"data"

    index.names<-sort(unique(df[,"source"]))
    index.names<-index.names[substr(index.names,1,5)=="index"]
  
    fill.flq<-function(obj){
        dms <-dims(obj)
        dmns<-dimnames(obj)
        dmns[[1]]<-dms[[2]]:dms[[3]]
        dmns[[2]]<-dms[[5]]:dms[[6]]
        
        res<-as.FLQuant(NA,dimnames=dmns)
        res[dimnames(obj)[[1]],dimnames(obj)[[2]],,,]<-obj[dimnames(obj)[[1]],dimnames(obj)[[2]],,,]
        
        return(res)
        }
          
    res2          <- new("FLXSA")
    res2@index.var<-new('FLQuants')
    res2@index.hat<-new('FLQuants')
    
    j=0
    for (i in index.names) {
       j=j+1
       res2@index.var[[j]]<-1.0/fill.flq(as.FLQuant(wts.df[wts.df[,  "source"]==i,-4]))
       res2@index.hat[[j]]<-    fill.flq(as.FLQuant(nhat.df[nhat.df[,"source"]==i,-4]))
       res2@index.var[[j]]<-1.0/as.FLQuant(wts.df[wts.df[,  "source"]==i,-4])
       res2@index.hat[[j]]<-as.FLQuant(nhat.df[nhat.df[,"source"]==i,-4])

       dmns                        <-dimnames(res2@index.hat[[j]])
       index                       <-trim(indices[[j]]@index,age=dmns$age,year=dmns$year)

       #print(index)
       res2@index[    is.na(index)]<-NA
       res2@index.hat[is.na(index)]<-NA
       res2@index.var[is.na(index)]<-NA
       res2@index.res[is.na(index)]<-NA
       }
    #names(wts) <-index.names
    #names(nhat)<-index.names

#    if (any(unique(df[,"source"])=="fshk")){
#       res2@var.fshk<-1.0/fill.flq(as.FLQuant(wts.df[wts.df[,  "source"]=="fshk",-4]))
#       res2@n.fshk  <-    fill.flq(as.FLQuant(nhat.df[nhat.df[,"source"]=="fshk",-4]))
#       }

#    if (any(unique(df[,"source"])=="nshk")){
#       res2@var.nshk<-1.0/fill.flq(as.FLQuant(wts.df[wts.df[,  "source"]=="nshk",-4]))
#       res2@n.nshk  <-    fill.flq(as.FLQuant(nhat.df[nhat.df[,"source"]=="nshk",-4]))
#       }

## F shrinkage
fshk   <-df[df[,"source"]=="fshk",]
if (length(fshk[,1])>0)
  {
  y.range<-range(fshk[,"year"])
  a.range<-range(fshk[,"age"])
  max.yr <-fshk[fshk[,"age"] ==a.range[2],]
  max.age<-fshk[fshk[,"year"]==y.range[2],]

  res2@n.fshk  <-as.FLQuant(NA,dimnames=list(age=a.range[1]:a.range[2],year=y.range[1]:y.range[2]))
  res2@n.fshk[as.character(max.age[,"age"]),as.character(y.range[2])]<-max.age[,"nhat"]
  res2@n.fshk[as.character(a.range[2]),as.character(max.yr[,"year"])]<-max.yr[,"nhat"]

  res2@var.nshk  <-as.FLQuant(NA,dimnames=list(age=a.range[1]:a.range[2],year=y.range[1]:y.range[2]))
  res2@var.nshk[as.character(max.age[,"age"]),as.character(y.range[2])]<-1/max.age[,"w"]
  res2@var.nshk[ as.character(a.range[2]),as.character(max.yr[,"year"])]<-1/max.yr[,"w"]
  }

## N shrinkage
nshk   <-df[df[,"source"]=="nshk",]
if (length(nshk[,1])>0)
  {
  y.range<-range(nshk[,"year"])
  a.range<-range(nshk[,"age"])
  max.yr <-nshk[nshk[,"age"] ==a.range[2],]
  max.age<-nshk[nshk[,"year"]==y.range[2],]

  res2@n.nshk  <-as.FLQuant(NA,dimnames=list(age=a.range[1]:a.range[2],year=y.range[1]:y.range[2]))
  res2@n.nshk[,as.character(y.range[2])]<-max.age[,"nhat"]
  res2@n.nshk[ as.character(a.range[2])]<-max.yr[,"nhat"]

  res2@var.nshk  <-as.FLQuant(NA,dimnames=list(age=a.range[1]:a.range[2],year=y.range[1]:y.range[2]))
  res2@var.nshk[,as.character(y.range[2])]<-1/max.age[,"w"]
  res2@var.nshk[ as.character(a.range[2])]<-1/max.yr[,"w"]
  }

    res2@diagnostics<- df
    res2@index.hat  <- res@index.hat
    #res2@index.hat  <- nhat
    #res2@index.var  <- wts
    res2@stock.n    <- FLQuant(res@stock.n@.Data)
    res2@harvest    <- FLQuant(res@harvest@.Data)
    res2@survivors  <- FLQuant(res@survivors@.Data)
    res2@se.int     <- FLQuant(res@se.int@.Data)
    res2@se.ext     <- FLQuant(res@se.ext@.Data)
    res2@n.fshk     <- FLQuant(res@n.fshk@.Data)
    res2@n.nshk     <- FLQuant(res@n.nshk@.Data)
    res2@var.fshk   <- FLQuant(res@var.fshk@.Data)
    res2@var.nshk   <- FLQuant(res@var.nshk@.Data)
    res2@harvest@units <- "f"
    res2@q.hat      <- res@q.hat
    res2@q2.hat     <- res@q2.hat
    res2@index.res  <- res@index.res
    res2@index      <- res@index

    res2@index.hat  <-FLQuants(res@index.hat)
    res2@index.res  <-FLQuants(res@index.res)
    res2@index.var  <-FLQuants(res@index.var)
    res2@q.hat      <-FLQuants(res@q.hat)    
    res2@q2.hat     <-FLQuants(res@q2.hat)   

    for (i in 1:length(indices))
       {
       res2@index.range[[i]]<-indices[[i]]@range
       res2@index.name[i]   <-indices[[i]]@name

       res2@index[[i]]      <-res@index[[i]]
       res2@index.hat[[i]]  <-FLQuant(res@index.hat[[i]], dimnames=dimnames(res@index.hat[[i]]))
       res2@index.res[[i]]  <-FLQuant(res@index.res[[i]], dimnames=dimnames(res@index.res[[i]]))
       res2@index.var[[i]]  <-FLQuant(res@index.var[[i]], dimnames=dimnames(res@index.var[[i]]))
       res2@q.hat[[i]]      <-FLQuant(res@q.hat[[i]],     dimnames=dimnames(res@q.hat[[i]]))
       res2@q2.hat[[i]]     <-FLQuant(res@q2.hat[[i]],    dimnames=dimnames(res@q2.hat[[i]]))

       dmns<-dimnames(indices[[i]]@index)
       na. <-is.na(indices[[i]]@index[dmns$age[ dmns$age  %in% dimnames(res2@index[[i]])$age],
                                      dmns$year[dmns$year %in% dimnames(res2@index[[i]])$year]])
       
       res2@index[[i]][    na.]<-NA
       res2@index.hat[[i]][na.]<-NA
       res2@index.res[[i]][na.]<-NA
       res2@index.var[[i]][na.]<-NA

       dimnames(res2@index[[i]]    )$unit   <-dimnames(indices[[i]]@index)$unit
       dimnames(res2@index.hat[[i]])$unit   <-dimnames(indices[[i]]@index)$unit
       dimnames(res2@index.res[[i]])$unit   <-dimnames(indices[[i]]@index)$unit
       dimnames(res2@index.var[[i]])$unit   <-dimnames(indices[[i]]@index)$unit
       dimnames(res2@q.hat[[i]])$unit       <-dimnames(indices[[i]]@index)$unit
       dimnames(res2@q2.hat[[i]])$unit      <-dimnames(indices[[i]]@index)$unit
       dimnames(res2@index[[i]]    )$season <-dimnames(indices[[i]]@index)$season
       dimnames(res2@index.hat[[i]])$season <-dimnames(indices[[i]]@index)$season
       dimnames(res2@index.res[[i]])$season <-dimnames(indices[[i]]@index)$season
       dimnames(res2@index.var[[i]])$season <-dimnames(indices[[i]]@index)$season
       dimnames(res2@q.hat[[i]])$season     <-dimnames(indices[[i]]@index)$season
       dimnames(res2@q2.hat[[i]])$season    <-dimnames(indices[[i]]@index)$season
       dimnames(res2@index[[i]]    )$area   <-dimnames(indices[[i]]@index)$area
       dimnames(res2@index.hat[[i]])$area   <-dimnames(indices[[i]]@index)$area
       dimnames(res2@index.res[[i]])$area   <-dimnames(indices[[i]]@index)$area
       dimnames(res2@index.var[[i]])$area   <-dimnames(indices[[i]]@index)$area
       dimnames(res2@q.hat[[i]])$area       <-dimnames(indices[[i]]@index)$area
       dimnames(res2@q2.hat[[i]])$area      <-dimnames(indices[[i]]@index)$area
       }

    res2@control <- res@control
    res2@call    <- res@call
    res2@desc    <- res@desc
    res2@range   <- stock@range
    units(res2@harvest)<-"f"

	  return(res2)
#	  return(convert6d(res2))
    }
    
FLXSA.control <- function(x=NULL, tol=10e-10, maxit=30, min.nse=0.3, fse=0.5, rage=0, qage=10, shk.n=TRUE,
                        	shk.f=TRUE, shk.yrs=5, shk.ages=5, window=100, tsrange=20, tspower=3, vpa=FALSE){
	if (is.null(x)){
		res <- new("FLXSA.control", tol=tol, maxit=as.integer(maxit), min.nse=min.nse, fse=fse,
		rage=as.integer(rage), qage=as.integer(qage), shk.n=as.logical(shk.n)[1],
		shk.f=as.logical(shk.f)[1], shk.yrs=as.integer(shk.yrs), shk.ages=as.integer(shk.ages),
		window=as.integer(window), tsrange=as.integer(tsrange), tspower=as.integer(tspower),
		vpa=as.logical(vpa)[1])
	  } else {	# We reuse an FLXSA.control object embedded in an FLXSA object
        if (!is.null(x) & !(is(x, "FLXSA") | is(x, "FLXSA.control")))
    		  	stop("FLXSA must be an 'FLXSA' or an 'FLXSA.control' object!")

        if (is(x, "FLXSA"))
    		   res <- x@control
        if (is(x, "FLXSA.control"))
    		   res <- x
        if (is.null(x))
    		   res <- new("FLXSA.control")
        		    
 		   # ... and possibly redefine some of its parameters
       if (!missing(tol))
       		res@tol <- tol
       if (!missing(maxit))
    			res@maxit <- as.integer(maxit)
    	 if (!missing(min.nse))
    			res@min.nse <- min.nse
   		 if (!missing(fse))
    			res@fse <- fse
     	 if (!missing(rage))
    			res@rage <- as.integer(rage)
       if (!missing(qage))
       		res@qage <- as.integer(qage)
       if (!missing(shk.n))
       		res@shk.n <- as.logical(shk.n)[1]
       if (!missing(shk.f))
       		res@shk.f <- as.logical(shk.f)[1]
       if (!missing(shk.yrs))
        	res@shk.yrs <- as.integer(shk.yrs)
       if (!missing(shk.ages))
       		res@shk.ages <- as.integer(shk.ages)
       if (!missing(window))
       		res@window <- as.integer(window)
       if (!missing(tsrange))
       		res@tsrange <- as.integer(tsrange)
       if (!missing(tspower))
      		res@tspower <- as.integer(tspower)
       if (!missing(vpa))
      		res@vpa <- as.logical(vpa)[1]

  		# Verify that this object is valid
 	  	test <- validObject(res)

		  if (!test)
        stop("Invalid object:", test)
	    }

	return(res)
  }

setMethod("assess", signature(control="FLXSA.control"),
   function(control, stock, indices, ...){

   if (!is(stock,   "FLStock"))   stop("stock not of type FLStock")
   if ( is(indices, "FLIndex"))   indices<-FLIndices(list(indices))
   if (!is(indices, "FLIndices")) stop("indices not of type FLIndices")
   
   print("FLXSA")
   FLXSA(stock,indices,control)   
   }
)

# Test if an object is of FLXSA class
is.FLXSA <- function(x)
	return(inherits(x, "FLXSA"))

# Test if an object is of FLXSA.control class
is.FLXSA.control <- function(x)
	return(inherits(x, "FLXSA.control"))

## show (a replacement for print of S3 classes)
setMethod("show", signature(object="FLXSA.control"),
	function(object){
      n.<-slotNames(object)
	   for (i in 1:length(n.))
         cat(n.[i],"\t\t",slot(object,n.[i]),"\n")
	}
)


# DiagsXSA::FLXSA
# Author     : MP, KH & JJP, EJ
# Modified to become summary fof FLXSA
# ---------------------------------------------------------------------------------

setGeneric("diagnostics", function(object, ...){
	standardGeneric("diagnostics")
	}
)

setMethod("diagnostics", signature(object="FLXSA"), function(object, ...){

    indices<-new("FLIndices")
    for (i in 1:length(object@index))
        {
#        indices[[i]]       <-new("FLIndex")
#        indices[[i]]@index <-object@index[[i]]
        indices[[i]] <- FLIndex(index=object@index[[i]])
        indices[[i]]@name  <-object@index.name[i]
        indices[[i]]@range <-object@index.range[[i]]
        }
    control <-object@control #<- eval(parse(text=xsa@call[4]))

    cat("FLR XSA Diagnostics ",  as.character(Sys.time()),"\n\nCPUE data from ", object@call[3], "\n\n",
        "Catch data for ", dims(object@stock.n)$year," years. ", dims(object@stock.n)$minyear," to ",
        dims(object@stock.n)$maxyear, ". Ages ",dims(object@stock.n)$min," to ",dims(object@stock.n)$max,".\n\n", sep="")

    # print general tuning series information
    idx.info  <- NULL
    for (i in 1:length(object@index)) {
       idx.info  <-  rbind(idx.info,c(indices[[i]]@name, (dims(object@index[[i]]))$min,
          (dims(object@index[[i]]))$max,(dims(object@index[[i]]))$minyear,(dims(object@index[[i]]))$maxyear,
          indices[[i]]@range["startf"], indices[[i]]@range["endf"]))
    }

    dimnames(idx.info) <- list(NULL,c("fleet","first age","last age","first year","last year","alpha","beta"))
    print(as.data.frame(idx.info))
    cat("\n\n","Time series weights :\n\n")
    cat(ifelse(control@tsrange==0|control@tspower==0, "   Tapered time weighting not applied\n\n",
        paste("   Tapered time weighting applied\n", "  Power =  ",control@tspower,"over ",control@tsrange,
        "years\n\n", sep=" ")))

    cat("Catchability analysis :\n\n")
    cat(ifelse(as.numeric(control@rage) < dims(object@stock.n)$min, "    Catchability independent of size for all ages\n\n",
        paste("    Catchability independent of size for ages >  ",control@rage,"\n\n",sep=" ")))
    cat(ifelse(as.numeric(control@qage) < dims(object@stock.n)$min, "    Catchability independent of age for all ages\n\n",
        paste("    Catchability independent of age for ages >  ",control@qage,"\n\n",sep=" ")))

    cat("Terminal population estimation :\n\n")
    cat(ifelse(control@shk.f, paste("    Survivor estimates shrunk towards the mean F\n",
        "   of the final  ",control@shk.yrs,"years or the ",control@shk.ages,"oldest ages.\n\n",
        "   S.E. of the mean to which the estimates are shrunk =  ", control@fse,"\n",sep=" "),
        "    Final estimates not shrunk towards mean F\n"))
    cat(ifelse(as.numeric(control@min.nse)==0, "\n", paste("\n", "   Minimum standard error for population\n",
        "   estimates derived from each fleet = ",control@min.nse,"\n\n", sep=" ")))
    cat("   prior weighting not applied\n\n")

    cat("Regression weights\n")
    ### Calculation of time series weighting
    yr.range <- (dims(object@harvest)$maxyear-9):dims(object@harvest)$maxyear
    regWt <- FLQuant(dimnames=list(age = 'all', year = yr.range))
    for(y in yr.range) regWt[,as.character(y)] <- (1-((max(yr.range)-y)/control@tsrange)^control@tspower)^control@tspower
    print(matrix(round(regWt,3),dims(regWt)$age,dimnames=list(age="all",year=yr.range)))

    cat("\n\n Fishing mortalities\n")
    print(matrix(round(trim(object@harvest,year=yr.range),3), dims(object@harvest)$age,
        dimnames=list(age=dims(object@harvest)$min:dims(object@harvest)$max, year=yr.range)))

    cat("\n\n XSA population number (Thousand)\n")
    print(t(matrix(round(trim(object@stock.n,year=yr.range),0), dims(object@stock.n)$age,
        dimnames=list(age=dims(object@stock.n)$min:dims(object@stock.n)$max, year=yr.range))))

    nextyear  <- dims(object@survivors)$maxyear
    cat("\n\n Estimated population abundance at 1st Jan ",nextyear,"\n")
    print(t(matrix(round(object@survivors[,as.character(nextyear)]),
        dimnames=list(age=dims(object@survivors)$min:dims(object@survivors)$max, year=nextyear))))

    ## tuning info
    for (f in 1:length(object@index)) {
        cat("\n\n Fleet: ",indices[[f]]@name,"\n\n","Log catchability residuals.\n\n")

        print(matrix(round(object@index.res[[f]],3), nrow=dims(object@index.res[[f]])$age,
            dimnames=list(age=dimnames(object@index.res[[f]])$age, year=dimnames(object@index.res[[f]])$year)))

        if (control@rage < dims(object@index[[f]])$max){
          cat("\n\n Mean log catchability and standard error of ages with catchability \n",
              "independent of year class strength and constant w.r.t. time \n\n")

          q.tab <- rbind(Mean_Logq=round(log(object@q.hat[[f]]),4), S.E_Logq=round(sd(matrix(object@index.res[[f]],
              dim(object@index.res[[f]])[2],dim(object@index.res[[f]])[1],byrow=T),na.rm=T),4))
          colnames(q.tab) <- dimnames(object@q.hat[[f]])$age

          if (dims(object@index[[f]])$min <= control@rage ) {
              print(q.tab[,as.character((control@rage+1):max(as.numeric(dimnames(object@index[[f]])$age)))])
          } else {print(q.tab)}
        }
        # print reg stats powermodel, note that maximum printed age is min(rage, max(age in tun series))
        if (dims(object@index[[f]])$min <= control@rage ) {
            cat("\n Regression statistics \n", "Ages with q dependent on year class strength \n")
            print(cbind((matrix(object@q2.hat[[f]][as.character(dims(object@index[[f]])$min:(min(control@rage,dims(object@index[[f]])$max)))],dimnames=list(age=paste("Age ",dims(object@index[[f]])$min:(min(control@rage,dims(object@index[[f]])$max)),sep=""),"slope"))),
            (matrix(object@q.hat[[f]][as.character(dims(object@index[[f]])$min:(min(control@rage,dims(object@index[[f]])$max)))],dimnames=list(age=dims(object@index[[f]])$min:(min(control@rage,dims(object@index[[f]])$max)),"intercept")))))
        }
    }

    cat("\n\n Terminal year survivor and F summaries: \n ")
    for ( age in sort(unique(object@diagnostics$age))){
        cat("\n Age ",age, " Year class =",  max(object@diagnostics$year) - age ," \n\n","source \n", sep="")
        weights <- object@diagnostics[(object@diagnostics$age==age) & (object@diagnostics$year== max(object@diagnostics$year)),]
        # calc surivors and scaled wts
        weights$survivors <- round(exp(weights$nhat))
        weights$scaledWts <- round(weights$w / sum(weights$w) ,3)
        row.names(weights) <- weights$source
        print(weights[ ,c("scaledWts","survivors","yrcls") ])
    }
    invisible()
})
