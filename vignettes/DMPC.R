## ------------------------------------------------------------------------
library(DMPC)
data("restroom")
DMPC.OOB.error(restroom, "Environment")

## ------------------------------------------------------------------------
per_variable_obj= DMPC.split.by.variable(restroom, "Environment")
per_variable_obj

## ------------------------------------------------------------------------
prevalences=DMPC.prevalence(per_variable_obj)
prevalences

## ------------------------------------------------------------------------
DMPC.best.prevalence(prevalences, "Environment")

## ------------------------------------------------------------------------
prevalence.60 = prevalences$`60`
prevalence.60

## ------------------------------------------------------------------------
DMPC.FDR(restroom, "Environment", bootstrap = 100, method.dist = "bray")

## ------------------------------------------------------------------------
DMPC.FDR(restroom, "Environment", bootstrap = 100, method.dist = "bray")

