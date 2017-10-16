# test-cod4.R - DESC
# /test-cod4.R

# Copyright European Union, 2017
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


library(FLXSA)

data(cod4)
cod4res2 <- FLXSA(cod4.stock, cod4.indices)

data(her4)
her4res2 <- FLXSA(her4.stock, her4.indices)

save(cod4res, her4res, file='res.RData', compress='xz')

load('res.RData')

all.equal(cod4res, cod4res2)

cl <- lapply(slotNames(cod4res), function(x) slot(cod4res, x))
cl2 <- lapply(slotNames(cod4res2), function(x) slot(cod4res2, x))

for(i in seq(length(cl)))
  print(identical(cl[[i]], cl2[[i]]))
