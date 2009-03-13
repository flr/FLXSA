# FLXSA.R - Comparison of FLXSA results for 3 stocks

# Author: Iago Mosqueira, CEFAS
# Maintainer: Laurie Kell, CEFAS
# Additions:
# Last Change: 21 May 2008 11:43
# $Id: FLXSA.R,v 1.28 2008/05/21 09:46:03 imosqueira Exp $

# Reference:
# Notes:

# TODO lun 26 feb 2007 16:15:37 CET iagoazti:

library(FLXSA)

# start test
setCon()
zz <- startTest("FLXSATest.txt")
tagTest("FLXSA testing ...")

# test function
foo <- function(stock, indices, f.ref, n.ref, fse, qage, shk.yrs, shk.ages) {
	for(i in 1:16) {
		control <- FLXSA.control(tol=1e-09, maxit=30, min.nse=0.3, fse=fse[i],
			rage=0, qage=qage[i], shk.n=1, shk.f=1, shk.yrs=shk.yrs[i], shk.ages=shk.ages[i],
			window=100, tsrange=100, tspower=100)
		xsa <- FLXSA(stock, indices, control)
		checkEqual(iter(f.ref, i), xsa@harvest, check.attributes=FALSE,
			tolerance = 0.05, msg=paste(name(stock), 'f', i))
		checkEqual(iter(n.ref, i), xsa@stock.n, check.attributes=FALSE,
			tolerance = 0.05, msg=paste(names(stock), 'n', i))
	}
}

# cod4
data(cod4)

foo(cod4.stock, cod4.indices, cod4.f.ref, cod4.n.ref, fse=rep(c(1.5, 2.5), 8),
	qage=rep(rep(c(6,7),4),2), shk.yrs=rep(rep(c(5,10),2),4), shk.ages=rep(c(4,6),8))

# her4
data(her4)

foo(her4.stock, her4.indices, her4.f.ref, her4.n.ref, fse=rep(c(1.5, 2.5), 8),
	qage=rep(rep(c(6,7),4),2), shk.yrs=rep(rep(c(5,10),2),4), shk.ages=rep(c(4,6),8))

# ple7a
data(ple7a)

foo(ple7a.stock, ple7a.indices, ple7a.f.ref, ple7a.n.ref, fse=rep(c(1.5, 2.5), 8),
	qage=rep(rep(c(6,7),4),2), shk.yrs=rep(rep(c(5,10),2),4), shk.ages=rep(c(4,6),8))

# vim: ft=r
