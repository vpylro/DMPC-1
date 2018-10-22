#' Random Forest with full dataset
#'
#' This fucntion runs randonForests on full dataset. If OOB error >0 then proceed with folowing analysis.
#'
#'@param physeq full object phyloseq
#'@param variable "variable to run the classification"
#'@keywords OOB
#'@examples DMPC.OOB.error(restroom, "Environment")
#'@export
DMPC.OOB.error= function(physeq, variable){
  set.seed(2125)
  if ((phyloseq::taxa_are_rows(physeq))==TRUE) {
    train=t(phyloseq::otu_table((physeq)))
  } else {(train=phyloseq::otu_table(physeq))}
  #train=otu_table(physeq)
  # Make one column for our outcome/response variable
  met_data <- as(object = phyloseq::sample_data(physeq), Class = "data.frame")
  response <- factor(met_data[, variable])
  training.set <- data.frame(response, train)
  train.model = randomForest::randomForest(response ~ ., data = training.set, importance = TRUE)
  res=train.model$err.rate[500,1]
  if (res ==0) {print("OOB error rate is zero. Your dataset present large differences. No prevalence cutoff in necessary")
  } else {return(res)}
}
