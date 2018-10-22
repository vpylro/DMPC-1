#'FDR
#'
#'Testing the likelihood of introducing bias by calculations of different prevalence cutoffs.
#'Randomizes (bootstrap) variable distribuition, performs prevalences cutoff and perMANOVA.
#'
#'@param physeq Original phyloseq object used for determining the best prevalence cutoff.
#'@param pseudo.var Variable to run the model
#'@param bootstrap  Number of randomizations
#'@param method.dist Distance/dissimilarity. Default is Bray-Curtis. Suported mesurements from vegan::vegdist()
#'@keywords FDR permanova prevalence
#'@examples
#' phylist=DMPC.split.by.variable(restroom, "Environment")
#' prev=DMPC.prevalence(phylist)
#' DMPC.best.prevalence(prev, "Environment")
#' DMPC.FDR(restroom, "Environment", bootstrap = 100, method.dist = "bray")
#'
#' @importFrom phyloseq "otu_table"
#' @importFrom phyloseq "sample_data"
#'
#'@export

DMPC.FDR = function (physeq, pseudo.var, bootstrap, method.dist = "bray"){
  split0 = DMPC::DMPC.split.by.variable(physeq, pseudo.var)
  cc = sapply(split0, microbiome::prevalence)
  vv = sort(sapply(cc, max))
  max_cut = as.integer((vv[1] * 100) - 0.001)
  cutoff = as.factor(seq(5, max_cut, by = 5))
  bootstrap = bootstrap
  FDR=NULL
  physeq2 = physeq
  for (c in cutoff) {
    pvalue=NULL
    tab=data.frame(replicate(bootstrap,sample(phyloseq::sample_data(physeq)[[pseudo.var]],replace=TRUE)))
    rownames(tab) = phyloseq::sample_names(physeq)
    phyloseq::sample_data(physeq2) = tab
    df90.2 = as(phyloseq::sample_data(physeq2), "data.frame")
    for (g in seq_len(ncol(phyloseq::sample_data(physeq2)))) {
      split = list()
      split = DMPC::DMPC.split.by.variable(physeq, pseudo.var)
      varphy=lapply(split, function(y) microbiome::core(y,detection=0, prevalence=as.numeric(c)/100))
      prev.cut = do.call(phyloseq::merge_phyloseq, varphy)
      if (any(phyloseq::sample_sums(prev.cut) ==0) ==FALSE) {
        merged = prev.cut
      }
      if ((phyloseq::taxa_are_rows(merged) == TRUE) == TRUE) {
        otu = t(otu_table(merged))
      } else {
        otu = otu_table(merged) }
      Variable = as.factor(df90.2[,g])
      d = vegan::vegdist(otu, method = method.dist)
      adonis = vegan::adonis(d ~ Variable, data = df90.2)
      pvalue[[length(pvalue) + 1]]=adonis[[1]][[6]][[1]]
    }
    FDR[c]=rbind(length(which(pvalue<=0.05)))
    FDRv=((FDR*100)/bootstrap)/100
    print("FDR value", quote = FALSE)
    print(FDRv)
  }
  finalFDR=as.data.frame(FDRv)
  rownames(finalFDR)=paste("Prevalence",rownames(finalFDR), "%", sep = "")
  return(finalFDR)
}
