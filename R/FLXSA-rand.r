#' rand
#' @description 
#' Uses the internal standard errors to conduct a Monte Carlo simulation of the terminal 
#' N-at-age in the last year
#'           
#' @aliases 
#' rand
#' 
#' @param n number of iters to generate
#' @param mean an \code{FLXSA} object 
#' @return sd  an \code{FLStock} object
#' 
#' @export
#' @docType methods
#' @rdname rand
#' 
#' @examples
#' \dontrun{
#' data(ple4)
#' data(ple4.indices)
#' xsa =FLXSA(ple4,ple4.indices)
#' ple4=rand(100,ple4,xsa)}
setGeneric('rand', function(n, mean, sd) standardGeneric('rand'))

setMethod("rand", signature(n='numeric', mean="FLStock", sd="FLXSA"),
    function(n,mean,sd){
      
        mean=mean+sd
        yrs =ac(range(mean)["maxyear"])
        
        n. =stock.n(  sd)[-dims(sd)$max,yrs]
        ctc =catch.n(mean)[-dims(sd)$max,yrs]
        m   =      m(mean)[-dims(sd)$max,yrs]
        surv=sd@survivors[-1,ac(range(sd)["maxyear"]+1)]
        dimnames(surv)=dimnames(n.)
        
        surv=rlnorm(n,log(surv),sd@se.ext)
        n.  =(surv+ctc)*exp(m/2)
        f   =log(n./surv)-m
        
        n.=setPlusGroup(n.,range(mean)["plusgroup"])
        f =setPlusGroup(f, range(mean)["plusgroup"])
        
        stock.n(mean)=propagate(stock.n(mean),n)
        harvest(mean)=propagate(harvest(mean),n)
        
        pg=ac(range(mean)["plusgroup"])
        
        n.[pg]=catch.n(mean)[pg,yrs]%*%
                   (m(mean)[pg,yrs]%+%f[pg])/(f[pg]*(1-exp(-m(mean)[pg,yrs]%-%f[pg])))
        
        stock.n(mean)[,yrs]=n.
        harvest(mean)[,yrs]=f
       
        mean})