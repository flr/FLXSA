# FLXSA.R - 
# FLXSA/R/FLXSA.R

# Copyright 2003-2017 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Iago Mosqueira, JRC

# FLXSA class {{{


#' Class FLXSA
#' 
#' A class for the results of an XSA analysis.
#' 
#' 
#' @name FLXSA-class
#' @aliases FLXSA-class show,FLXSA-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("FLXSA", ...)} and are output by calls to \code{\link{FLXSA}}.
#' @author Laurence Kell
#' @seealso \code{\link{FLXSA}}
#' @references Darby, C. D., and Flatman, S. 1994. Virtual Population Analysis:
#' version 3.1 (Windows/Dos) user guide. Info. Tech. Ser., MAFF Direct. Fish.
#' Res., Lowestoft, (1): 85pp.
#' 
#' Shepherd, J.G. 1992. Extended survivors analysis: an improved method for the
#' analysis of catch-at-age data and catch-per-unit-effort data. Working paper
#' No. 11 ICES Multi-species Assessment Working Group, June 1992, Copenhagen,
#' Denmark. 22pp. (mimeo).
#' 
#' Shepherd, J.G. 1994. Prediction of yearclass strength by calibration
#' regression analysis of multiple recruit index series. ICES J. Mar. Sci. In
#' Prep.
#' @keywords classes

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
		survivors=FLQuant(quant='age'),
		se.int   =FLQuant(quant='age'),
		se.ext   =FLQuant(quant='age'),
		n.fshk   =FLQuant(quant='age'),
		n.nshk   =FLQuant(quant='age'),
		var.fshk =FLQuant(quant='age'),
		var.nshk =FLQuant(quant='age'),
		catch.n =FLQuant(quant='age'),
		stock.n =FLQuant(quant='age'),
		harvest =FLQuant(quant='age'),
		q.hat    =FLQuants(),
		q2.hat   =FLQuants(),
		diagnostics=new("data.frame"),
		control  =new("FLXSA.control"))
)

setValidity("FLXSA", function(object) {

	# All FLQuant objects must have same dimensions
	return(TRUE)
	
  Dim <- dim(object@stock.n)
	if (!all(dim(object@harvest) == Dim))
		return("n and f arrays must have same dimensions")
	
  # Everything is fine
	return(TRUE)
})
# }}}

# FLXSA() {{{

#' Create a new FLXSA object -run an XSA analysis-
#' 
#' This function runs an XSA (extended survivor analysis) and creates an FLXSA
#' object used to analyse its results.
#' 
#' 
#' Virtual population analysis and cohort analysis are essentially accountancy
#' methods whereby a stock's historical population structure may be
#' reconstructed from total catch data given a particular level of natural
#' mortality. Firstly, however, numbers at age in the last year and age have to
#' be found since both methods iterate backwards down a cohort. The main
#' problem in many sequential age based assessment methods is therefore to
#' estimate these terminal population numbers. In XSA these are found from the
#' relationship between catch per unit effort (CPUE), abundance and year class
#' strength.
#' 
#' Estimates of the catchability for the oldest age in an assessment, tuned by
#' the ad hoc or XSA procedures, are directly dependent on the terminal
#' population or F values used to initialise the underlying VPA. Catchability
#' at the oldest age is therefore under-determined and cannot be utilised
#' without additional information. Within the ad hoc tuning procedures the
#' additional information is obtained by making the assumption that the
#' exploitation pattern on the oldest ages was constant during the assessment
#' time series. F on the oldest age in the final year is estimated as a
#' proportion of an average of the F for preceding ages in the same year. XSA
#' uses an alternative approach by making the assumption that fleet
#' catchability is constant (independent of age) above a certain age. The age
#' (constant for all fleets) is user-defined. For each fleet, the catchability
#' value estimated at the specified age, is used to derive population abundance
#' estimates for all subsequent ages in the fleet data set.
#' 
#' @param stock An FLStock object to be used for the analysis
#' @param indices An FLIndices object holding the indices of abundance to
#' consider in the model
#' @param control An \code{FLXSA.control} object giving parameters of the model
#' (see \code{\link{FLXSA.control}})
#' @param desc A short description of this analysis
#' @param diag.flag If TRUE returns all diagnostics, if FALSEonly returns
#' stock.n, harvest and control
#' @return
#' 
#' An \code{FLXSA} object is returned, whith slots: \item{n }{An FLQuant with
#' the number of individuals at age} \item{f }{An FLQuant with the fishing
#' mortality} \item{swt }{An FLQuant with the stock weight} \item{mat }{An
#' FLQuant with the maturity indices} \item{qres }{A list with residuals for q}
#' \item{cpue }{A list with the various cpues} \item{wts }{A list with the
#' various weights} \item{control}{The \code{FLXSA.control} object that was
#' used for this analysis} \item{call }{A copy of the call to run this
#' analysis} \item{desc }{A description of the analysis}
#' @note See \code{\link{update}} to learn how to update stock data according
#' to an XSA analysis
#' @author Laurence Kell and Philippe Grosjean
#' @seealso \code{\link{FLXSA.control}}, \code{\link[FLCore]{FLStock-class}}
#' @references Darby, C. D., and Flatman, S. 1994. Virtual Population Analysis:
#' version 3.1 (Windows/Dos) user guide. Info. Tech. Ser., MAFF Direct. Fish.
#' Res., Lowestoft, (1): 85pp.
#' 
#' Shepherd, J.G. 1992. Extended survivors analysis: an improved method for the
#' analysis of catch-at-age data and catch-per-unit-effort data. Working paper
#' No. 11 ICES Multi-species Assessment Working Group, June 1992, Copenhagen,
#' Denmark. 22pp. (mimeo).
#' 
#' Shepherd, J.G. 1994. Prediction of yearclass strength by calibration
#' regression analysis of multiple recruit index series. ICES J. Mar. Sci. In
#' Prep.
#' @keywords classes
#' @examples
#'
#' data(ple4)
#' data(ple4.indices)
#' 
#' res <- FLXSA(ple4, ple4.indices)
#' 

setGeneric("FLXSA", function(stock, indices, ...)
	standardGeneric("FLXSA"))

setMethod("FLXSA", signature(stock="FLStock", indices="FLIndex"),
  function(stock, indices, control=FLXSA.control(), diag.flag=TRUE) {
    FLXSA(stock=stock, indices=FLIndices(indices), control=control, diag.flag=diag.flag)
  }
)

setMethod("FLXSA", signature(stock="FLStock", indices="FLIndices"),
  function(stock, indices, control=FLXSA.control(), diag.flag=TRUE) {
    
    Call <- match.call()
    # check FLIndices input
    for (i in 1:length(indices)) {
      # startf & endf present in @range
      if (is.na(indices[[i]]@range["startf"]) || is.na(indices[[i]]@range["endf"]))
  	     stop(paste("Must supply startf & endf for range in FLIndex",i))
      # @range has all elements
      if (!all(c("min","max","plusgroup","minyear","maxyear","startf","endf") %in%
        names(indices[[i]]@range)))
        stop("Range must have names 'min','max','plusgroup','minyear','maxyear',
          'startf','endf'")
      # adjust ranges to dims()
      indices[[i]]@range["min"] <- max(indices[[i]]@range["min"], dims(indices[[i]])$min,
        stock@range["min"])
      indices[[i]]@range["max"] <- min(indices[[i]]@range["max"], dims(indices[[i]])$max,
        stock@range["max"])
      indices[[i]]@range["minyear"] <- max(indices[[i]]@range["minyear"],
        dims(indices[[i]])$minyear, stock@range["minyear"])
      indices[[i]]@range["maxyear"] <- min(indices[[i]]@range["maxyear"],
        dims(indices[[i]])$maxyear, stock@range["maxyear"])

      # trim according to range
      age <- indices[[i]]@range["min"]:indices[[i]]@range["max"]
      year<- indices[[i]]@range["minyear"]:indices[[i]]@range["maxyear"]
      
      indices[[i]] <- trim(indices[[i]], age=age, year=year)
    }
      
    # Double check validity
    if (!validObject(stock))
  	  stop("stock is not valid!")
  	if (!validObject(indices))
  	  stop("FLIndices is not valid!")
  	if (!validObject(control))
  	  stop("control is not valid!")

    # adjust range for minage
    if ("minage" %in% names(stock@range))
      minage <- stock@range["minage"]
    else if ("min" %in% names(stock@range))
      minage <- stock@range["min"]
    else if ("minquant"  %in% names(stock@range))
      minage <- stock@range["minquant"]
    else
      stop("'minage' not found in range")
      
    # adjust range for maxage
    if ("maxage" %in% names(stock@range))
      maxage <- stock@range["maxage"]
    else if ("max" %in% names(stock@range))
      maxage <- stock@range["max"]
    else if ("maxquant" %in% names(stock@range))
      maxage <- stock@range["maxquant"]
    else
      stop("'maxage' not found in range")
    
    # adjust plsugroup
    if ("plusgroup" %in% names(stock@range))
      stock@range["plusgroup"] <- maxage

    if (maxage<minage | stock@range["maxyear"] < stock@range["minyear"])
      stop("Error in range")
    
    if (is.na(stock@range["plusgroup"]))
      stop("Plus Group must be specified")
    
    # trim stock
    stock <- trim(stock, year=stock@range["minyear"]:stock@range["maxyear"])

    if (!is.na(stock@range["plusgroup"]) &
      stock@range["plusgroup"] < dims(stock@catch.n)$max)
      stock <- setPlusGroup(stock, stock@range["plusgroup"])

    stock@m <- stock@m[as.character(minage:maxage),,,,]

    # check catch.n is available
    if (all(is.na(stock@catch.n)))
      stop("catch.n is not available")

    stock@stock.n[] <- new('FLQuant')
    stock@harvest[] <- new('FLQuant')

    # fqs
    fqs <- function(assess) {
      assess@index <- new("FLQuants", lapply(assess@index,FLQuant))
      assess@index.hat <- new("FLQuants", lapply(assess@index.hat,FLQuant))
      assess@index.var <- new("FLQuants", lapply(assess@index.var,FLQuant))
      assess@index.res <- new("FLQuants", lapply(assess@index.res,FLQuant))
      assess@q.hat <- new("FLQuants", lapply(assess@q.hat,FLQuant))
      assess@q2.hat <- new("FLQuants", lapply(assess@q2.hat,FLQuant))
       
      if (validObject(assess))
        return(assess)
      else
        stop("not valid")
    }

    #
    iters.stock  <-max(unlist(qapply(stock, function(x) dims(x)$iter)))
    iters.indices<-max(unlist(lapply(indices@.Data,
      function(x) max(unlist(qapply(x, function(x2) dims(x2)$iter))))))
      
    if ((iters.stock>1 | iters.indices>1) && missing(diag.flag))
      diag.flag <- FALSE

    if ((iters.stock>1 && iters.indices>1) && diag.flag)
      return("Multiple iters only allowed if diag.flag=FALSE")

    if(!diag.flag) {
      
      res<-.Call("runFLXSA", iter(stock, 1), lapply(indices, iter, 1), control, FALSE)
      
      iters <- max(iters.stock,iters.indices)
       
      if (iters>1) {
          
        res@stock.n<-propagate(FLQuant(res@stock.n@.Data),iters)
        res@harvest<-propagate(FLQuant(res@harvest@.Data),iters)
        for (i in as.character(2:iters)) {
          res. <- .Call("runFLXSA", iter(stock,i), lapply(indices, iter,i), control, FALSE)
          iter(res@stock.n,i)<-FLQuant(res.@stock.n@.Data)
          iter(res@harvest,i)<-FLQuant(res.@harvest@.Data)
        }
      }       

      res@harvest@units <- "f"
      res@range   <- stock@range
      
      return(res)
      }

    res <-fqs(.Call("runFLXSA", stock, indices, control, diag.flag))
      
    if (class(res) != "FLXSA")
      return(res)
    res@call <- as.character(Call)

    # put wts amd nhats into a data frame
    df <- as.data.frame(res@wts)
    df1 <- (df[4])
    df1[df1 >= 1, 1] <- paste("index", df1[df1 >= 1, 1])
    df1[df1 == -1, 1] <- "fshk"
    df1[df1 == -2, 1] <- "nshk"
    df <- cbind(df[-4], df1)

    names(df) <- c("w", "nhat", "yrcls", "age", "year", "source")
  	
    for(i in 1:length(indices)) {
      v <- paste("index", i)
      if (length(slot(indices[[i]],"name"))>0)
         df$source[df$source==v] <- slot(indices[[i]],"name")   
    }

    wts.df <-df[,c(4,5,1,6)]
    names(wts.df)[3] <-"data"
    nhat.df<-df[,c(4,5,2,6)]
    names(nhat.df)[3]<-"data"

    index.names<-sort(unique(df[,"source"]))
    index.names<-index.names[substr(index.names,1,5)=="index"]
  
    fill.flq<-function(obj) {
      dms <-dims(obj)
      dmns<-dimnames(obj)
      dmns[[1]]<-dms[[2]]:dms[[3]]
      dmns[[2]]<-dms[[5]]:dms[[6]]
        
      res<-as.FLQuant(NA,dimnames=dmns)
      res[dimnames(obj)[[1]],dimnames(obj)[[2]],,,]<-obj[dimnames(obj)[[1]],
        dimnames(obj)[[2]],,,]
        
      return(res)
    }
          
    res2  <- new("FLXSA")
    res2@index.var<-new('FLQuants')
    res2@index.hat<-new('FLQuants')
    res2@index    <-new('FLQuants')
 
    j=0
    
    for (i in index.names) {
      j=j+1
      res2@index.var[[j]]<-1.0/fill.flq(as.FLQuant(wts.df[wts.df[,  "source"]==i,-4]))
      res2@index.hat[[j]]<-    fill.flq(as.FLQuant(nhat.df[nhat.df[,"source"]==i,-4]))
      res2@index.var[[j]]<-1.0/as.FLQuant(wts.df[wts.df[,  "source"]==i,-4])
      res2@index.hat[[j]]<-as.FLQuant(nhat.df[nhat.df[,"source"]==i,-4])
      dmns <- dimnames(res2@index.hat[[j]])
      index  <- trim(indices[[j]]@index,age=dmns$age,year=dmns$year)

      res2@index[    is.na(index)]<-NA
      res2@index.hat[is.na(index)]<-NA
      res2@index.var[is.na(index)]<-NA
      res2@index.res[is.na(index)]<-NA
    }

    # F shrinkage
    fshk  <- df[df[,"source"]=="fshk",]

    if (length(fshk[,1])>0) {
      y.range <- range(fshk[,"year"])
      a.range<-range(fshk[,"age"])
      max.yr <-fshk[fshk[,"age"] ==a.range[2],]
      max.age<-fshk[fshk[,"year"]==y.range[2],]

      res2@n.fshk  <-as.FLQuant(NA,dimnames=list(age=a.range[1]:a.range[2],
          year=y.range[1]:y.range[2]))
      res2@n.fshk[as.character(max.age[,"age"]),as.character(y.range[2])]<-max.age[,
        "nhat"]
      res2@n.fshk[as.character(a.range[2]),as.character(max.yr[,"year"])]<-max.yr[,"nhat"]

      res2@var.nshk  <-as.FLQuant(NA,dimnames=list(age=a.range[1]:a.range[2],
          year=y.range[1]:y.range[2]))
      res2@var.nshk[as.character(max.age[,"age"]),as.character(y.range[2])]<-
        1/max.age[,"w"]
      res2@var.nshk[ as.character(a.range[2]),as.character(max.yr[,"year"])]<-
        1/max.yr[,"w"]
    }

    # N shrinkage
    nshk   <-df[df[,"source"]=="nshk",]
    if (length(nshk[,1])>0) {
      y.range<-range(nshk[,"year"])
      a.range<-range(nshk[,"age"])
      max.yr <-nshk[nshk[,"age"] ==a.range[2],]
      max.age<-nshk[nshk[,"year"]==y.range[2],]

      res2@n.nshk  <-as.FLQuant(NA,dimnames=list(age=a.range[1]:a.range[2],
          year=y.range[1]:y.range[2]))
      res2@n.nshk[,as.character(y.range[2])] <- max.age[,"nhat"]
      res2@n.nshk[ as.character(a.range[2])] <- max.yr[,"nhat"]

      res2@var.nshk <- as.FLQuant(NA,dimnames=list(age=a.range[1]:a.range[2],
          year=y.range[1]:y.range[2]))
      res2@var.nshk[,as.character(y.range[2])]<-1/max.age[,"w"]
      res2@var.nshk[ as.character(a.range[2])]<-1/max.yr[,"w"]
    }
    
    res2@diagnostics<- df
    res2@index.hat  <- res@index.hat
    res2@stock.n    <- FLQuant(res@stock.n@.Data)
      units(res2@stock.n) <- units(stock@catch.n)
    res2@harvest    <- FLQuant(res@harvest@.Data)
    res2@survivors  <- FLQuant(res@survivors@.Data)
      units(res2@survivors) <- units(stock@catch.n)
    res2@se.int     <- FLQuant(res@se.int@.Data)
      units(res2@se.int) <- ""
    res2@se.ext     <- FLQuant(res@se.ext@.Data)
      units(res2@se.ext) <- ""
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

    for (i in 1:length(indices)) {
       res2@index.range[[i]]<-indices[[i]]@range
       res2@index.name[i]   <-indices[[i]]@name
       res2@index[[i]]      <-res@index[[i]]
       res2@index.hat[[i]]  <-FLQuant(res@index.hat[[i]],
         dimnames=dimnames(res@index.hat[[i]]))
       res2@index.res[[i]]  <-FLQuant(res@index.res[[i]],
         dimnames=dimnames(res@index.res[[i]]))
       res2@index.var[[i]]  <-FLQuant(res@index.var[[i]],
         dimnames=dimnames(res@index.var[[i]]))
       res2@q.hat[[i]]      <-FLQuant(res@q.hat[[i]],
         dimnames=dimnames(res@q.hat[[i]]))
       res2@q2.hat[[i]]     <-FLQuant(res@q2.hat[[i]],
        dimnames=dimnames(res@q2.hat[[i]]))

       dmns<-dimnames(indices[[i]]@index)
       na. <-is.na(indices[[i]]@index[
           dmns$age[dmns$age %in% dimnames(res2@index[[i]])$age],
           dmns$year[dmns$year %in% dimnames(res2@index[[i]])$year]])
       
       res2@index[[i]][na.] <- NA
       res2@index.hat[[i]][na.]<- NA
       res2@index.res[[i]][na.]<- NA
       res2@index.var[[i]][na.]<- NA

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

    # index
    names(res2@index) <- names(res2@index.res) <- names(res2@index.var) <-
     names(res2@index.hat)  <- names(res2@q2.hat) <- names(res2@q.hat) <- names(indices)

   for(i in names(indices)) {
      units(res2@index[[i]]) <- units(indices[[i]]@index)
      units(res2@index.hat[[i]]) <- units(indices[[i]]@index)
      units(res2@index.res[[i]]) <- ""
      units(res2@index.var[[i]]) <- ""
      units(res2@q.hat[[i]]) <- ""
      units(res2@q2.hat[[i]]) <- ""
    }

    res2@control <- res@control
    res2@call    <- res@call
    res2@desc    <- paste("FLXSA run:", res@desc)
    res2@range   <- stock@range
    units(res2@harvest)<-"f"

    res2@control@shk.n   <-as.logical(res2@control@shk.n)
    res2@control@shk.f   <-as.logical(res2@control@shk.f)
    res2@control@vpa     <-as.logical(res2@control@vpa)

	  return(res2)
    }
)
# }}}

# assess {{{
setMethod("assess", signature(control="FLXSA.control"),
   function(control, stock, indices, ...){

   if (!is(stock,   "FLStock"))   stop("stock not of type FLStock")
   if ( is(indices, "FLIndex"))   indices<-FLIndices(list(indices))
   if (!is(indices, "FLIndices")) stop("indices not of type FLIndices")
   
   print("FLXSA")
   FLXSA(stock,indices,control)   
   }
)
# }}}

# is.FLXSA {{{


#' is.FLXSA
#' 
#' These two functions return \code{code} if objects are of class FLXSA and
#' FLXSA.control, respectively.
#' 
#' 
#' @aliases is.FLXSA is.FLXSA.control
#' @param x An object to be tested
#' @return is.FLXSA returns \code{TRUE} if its argument is of class
#' \code{FLXSA} (that is, has "FLXSA" amongst its classes) and FALSE otherwise.
#' is.FLXSA.control returns \code{TRUE} if its argument is of class
#' \code{FLXSA.control} (that is, has "FLXSA.control" amongst its classes) and
#' FALSE otherwise.
#' @keywords attribute
#' @examples
#' 
#' xsa <- FLXSA.control()
#' is.FLXSA.control(xsa)
#' 
is.FLXSA <- function(x)
	return(inherits(x, "FLXSA"))
# }}}

# diagnostics {{{
setGeneric("diagnostics", function(object, ...){
	standardGeneric("diagnostics")
	}
)

#' XSA diagnostics
#' 
#' Provides the diagnostics table used in ICES WG to analyse the results of the
#' XSA run.
#' 
#' 
#' @name diagnostics-method
#' @aliases diagnostics diagnostics-method diagnostics,FLXSA-method
#' @docType methods
#' @section Generic Function: \describe{ \item{call}{diagnostics(object, ...}
#' \item{object}{An object with the results of a VPA method.} }
#' @keywords methods

setMethod("diagnostics", signature(object="FLXSA"), 
  tempfun <- function(object, sections=rep(T, 8), ...){
#print(1)
    indices<-new("FLIndices")
    for (i in 1:length(object@index))
        {
        indices[[i]]       <- FLIndex(index=object@index[[i]])
        indices[[i]]@name  <- object@index.name[i]
        }
    control <-object@control 

    titledat <- paste("FLR XSA Diagnostics ",  as.character(Sys.time()),
                    "\n\nCPUE data from ", object@call[3], "\n\n",
                    "Catch data for ", dims(object@stock.n)$year," years ", 
                     dims(object@stock.n)$minyear," to ",dims(object@stock.n)$maxyear, 
                    ". Ages ",dims(object@stock.n)$min," to ",dims(object@stock.n)$max,".\n\n", sep="")
#   print(titledat)

    # general tuning series information
    idx.info  <- NULL
    for (i in 1:length(object@index)) {
       idx.info  <-  rbind(idx.info,c(indices[[i]]@name, (dims(object@index[[i]]))$min,
          (dims(object@index[[i]]))$max,(dims(object@index[[i]]))$minyear,(dims(object@index[[i]]))$maxyear,
          indices[[i]]@range["startf"], indices[[i]]@range["endf"]))
    }
    dimnames(idx.info) <- list(NULL,c("fleet","first age","last age","first year","last year","alpha","beta"))
#    print(as.data.frame(idx.info))

    set1 <- paste("\n\n","Time series weights :\n\n")
    set2 <- paste(ifelse(control@tsrange==0|control@tspower==0, "   Tapered time weighting not applied\n\n",
        paste("   Tapered time weighting applied\n", "  Power =  ",control@tspower,"over ",control@tsrange,
        "years\n\n", sep=" ")))

    set3 <- "Catchability analysis :\n\n"
    set4 <- paste(ifelse(as.numeric(control@rage) < dims(object@stock.n)$min, "    Catchability independent of size for all ages\n\n",
        paste("    Catchability independent of size for ages >  ",control@rage,"\n\n",sep=" ")))

    set5 <- paste(ifelse(as.numeric(control@qage) < dims(object@stock.n)$min, "    Catchability independent of age for all ages\n\n",
        paste("    Catchability independent of age for ages >  ",control@qage,"\n\n",sep=" ")))

    set6 <- "Terminal population estimation :\n\n"
    set7 <- paste(ifelse(control@shk.f, paste("    Survivor estimates shrunk towards the mean F\n",
        "   of the final  ",control@shk.yrs,"years or the ",control@shk.ages,"oldest ages.\n\n",
        "   S.E. of the mean to which the estimates are shrunk =  ", control@fse,"\n",sep=" "),
        "    Final estimates not shrunk towards mean F\n"))
    set8 <- ifelse(as.numeric(control@min.nse)==0, "\n", paste("\n", "   Minimum standard error for population\n",
        "   estimates derived from each fleet = ",control@min.nse,"\n\n", sep=" "))
    set9 <- "   prior weighting not applied\n\n"

# cat(set1, set2, set3, set4, set5, set6, set7, set8, set9)

    regwtstitle <- "Regression weights\n"
    ### Calculation of time series weighting
    yr.range <- max(dims(object@harvest)$minyear,(dims(object@harvest)$maxyear-9)):dims(object@harvest)$maxyear
    regWt <- FLQuant(dimnames=list(age = 'all', year = yr.range))
    for(y in yr.range) 
      regWt[,as.character(y)] <- (1-((max(yr.range)-y)/control@tsrange)^control@tspower)^control@tspower
    regwts <- matrix(round(c(regWt),3),dims(regWt)$age,dimnames=list(age="all",year=yr.range))
# cat(regwtstitle)
# print(regwts)

    FMtitle <- "\n\n Fishing mortalities\n"
    FM      <- matrix(round(c(trim(object@harvest,year=yr.range)),3), dims(object@harvest)$age,
        dimnames=list(age=dims(object@harvest)$min:dims(object@harvest)$max, year=yr.range))
# cat(FMtitle)
# print(FM)

    PNtitle <- "\n\n XSA population number (Thousand)\n"
    PN <- (t(matrix(round(c(trim(object@stock.n,year=yr.range)),0), dims(object@stock.n)$age,
        dimnames=list(age=dims(object@stock.n)$min:dims(object@stock.n)$max, year=yr.range))))
# cat(PNtitle)
# print(PN)

    nextyear  <- dims(object@survivors)$maxyear
    survtitle <- paste("\n\n Estimated population abundance at 1st Jan ",nextyear,"\n")
    survivors <- t(matrix(round(c(object@survivors[,as.character(nextyear)])),
        dimnames=list(age=dims(object@survivors)$min:dims(object@survivors)$max, year=nextyear)))
# cat(survtitle)
# print(survivors)

    ## tuning info
    fleetname <- list()
    logQs     <- list()
    mlqbtitle <- list()
    mlqb      <- list()

    for (f in 1:length(object@index)) {
        fleetname[[f]] <- paste("\n\n Fleet: ",indices[[f]]@name,"\n\n","Log catchability residuals.\n\n")
        logQs[[f]] <- matrix(round(c(object@index.res[[f]]),3), nrow=dims(object@index.res[[f]])$age,
            dimnames=list(age=dimnames(object@index.res[[f]])$age, year=dimnames(object@index.res[[f]])$year))

#       print(fleetname[[f]])
#       print(logQs[[f]])

        if (control@rage < dims(object@index[[f]])$max){
          mlqbtitle[[f]] <- paste("\n\n Mean log catchability and standard error of ages with catchability \n",
              "independent of year class strength and constant w.r.t. time \n\n")
    q.tab <- rbind(
      Mean_Logq=round(c(log(object@q.hat[[f]])),4),
      S.E_Logq=round(c(apply(object@index.res[[f]],1,sd,na.rm=TRUE)), 4))

    colnames(q.tab) <- dimnames(object@q.hat[[f]])$age

          if (dims(object@index[[f]])$min <= control@rage ) {
              mlqb[[f]] <- q.tab[,as.character((control@rage+1):max(as.numeric(dimnames(object@index[[f]])$age)))]
          } else {mlqb[[f]] <- q.tab}
        }
        # print reg stats powermodel, note that maximum printed age is min(rage, max(age in tun series))
        if (dims(object@index[[f]])$min <= control@rage ) {
            mlqbtitle[[f]] <- paste("\n Regression statistics \n", "Ages with q dependent on year class strength \n")
            mlqb[[f]] <- paste(cbind((matrix(object@q2.hat[[f]][as.character(dims(object@index[[f]])$min:(min(control@rage,dims(object@index[[f]])$max)))],dimnames=list(age=paste("Age ",dims(object@index[[f]])$min:(min(control@rage,dims(object@index[[f]])$max)),sep=""),"slope"))),
            (matrix(object@q.hat[[f]][as.character(dims(object@index[[f]])$min:(min(control@rage,dims(object@index[[f]])$max)))],dimnames=list(age=dims(object@index[[f]])$min:(min(control@rage,dims(object@index[[f]])$max)),"intercept")))))
        }
#        cat(mlqbtitle[[f]])
#        print(mlqb[[f]])
    }

    tysurvtitle <- "\n\n Terminal year survivor and F summaries: \n "
    header      <- list()
    tysurvtext  <- list()
# cat(tysurvtitle)
    agevec <- sort(unique(object@diagnostics$age))
    for ( age in 1:length(agevec)){
#    cat("age: ", age, "\n")
        header[[age]] <- paste("\n ,Age ",agevec[age], " Year class =",  max(object@diagnostics$year) - agevec[age] ," \n\n","source \n", sep="")
        weights <- object@diagnostics[(object@diagnostics$age==agevec[age]) & (object@diagnostics$year== max(object@diagnostics$year)),]
        # calc surivors and scaled wts
        weights$survivors <- round(c(exp(weights$nhat)))
        weights$scaledWts <- round(c(weights$w / sum(weights$w)) ,3)
        row.names(weights) <- weights$source
        tysurvtext[[age]] <- weights[ ,c("scaledWts","survivors","yrcls") ]

# cat(header[[age]])
# print(tysurvtext[[age]])
    }

## send the text to the screen

if(sections[1]){
  cat(titledat)
  print(as.data.frame(idx.info))
}
if(sections[2]){
  cat(set1, set2, set3, set4, set5, set6, set7, set8, set9)
}
if(sections[3]){
  cat(regwtstitle)
  print(regwts)
}
if(sections[4]){
  cat(FMtitle)
  print(FM)
}
if(sections[5]){
  cat(PNtitle)
  print(PN)
}
if(sections[6]){
  cat(survtitle)
  print(survivors)
}
if(sections[7]){
  for(f in 1:length(object@index)) {
    cat(fleetname[[f]])
    print(logQs[[f]])
    cat(mlqbtitle[[f]])
    print(mlqb[[f]])
  }
}
if(sections[8]){
  cat(tysurvtitle)
  for ( age in 1:length(agevec)){
    cat(header[[age]])
    print(tysurvtext[[age]])
  }
}
    invisible()
}
) # }}}
