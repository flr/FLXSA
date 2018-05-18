# test.R - DESC
# /test.R

# Copyright European Union, 2017
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


context("Comparison 32-64 bit FLXSA")

data(ple4)
data(ple4.indices)

# NEW run
nre <- FLXSA(ple4, ple4.indices)

# OLD run
load('ple4FLXSA.RData')

test_that("ple4xsa 32 and 64 are equal", {
  expect_equivalent(ple4xsa, nre)
})

test_that("ple4 after XSA for 32 and 64 are equal", {

  nstk <- ple4 + nre
  ostk <- ple4 + ple4xsa
  
  expect_equal(ostk, nstk)
})
