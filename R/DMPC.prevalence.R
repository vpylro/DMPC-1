#'Prevalences cutoffs
#'
#'First calculates which is the highest prevalence possible (no empty OTU table), calculates prevalence for each cutoff by increments of 5%.
#' After prevalence calculation, merging of phyloseq objects for each variable and cutoff is done when no sample is lost.
#' When any sample_sums==0 FALSE.
#'
#'@param phylist list of phyloseq objects (one for each variable)
#'@keywords prevalence
#'@examples phylist=DMPC.split.by.variable(restroom, "Environment")
#' prev=DMPC.prevalence(phylist)
#'@export
DMPC.prevalence = function(phylist) {
  cc=lapply(phylist, microbiome::prevalence)
  vv=sort(sapply(cc, max))
  max_cut=as.integer((vv[1]*100)-0.001)
  merged=list()
  var.phy=list()
  cutoff= as.factor(seq(5,max_cut, by=5))
  for (i in cutoff){
    for (j in names(phylist)){
      var.phy[[j]]=microbiome::core(phylist[[j]], detection = 0, prevalence= as.numeric(i)/100)
      prev.cut=do.call(phyloseq::merge_phyloseq, var.phy)
      if (any(phyloseq::sample_sums(prev.cut)==0)==FALSE){
        merged[[i]]=prev.cut
      }
    }
  }
  return(merged)
}
