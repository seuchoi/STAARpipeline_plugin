STAARpipeline_LOVO<-function(genofile,genes,category,Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff,rv_num_cutoff,silent=TRUE){

  ## leave one variant out analysis start
  results <- c()
  variant_info<-variantInfo(genofile)
  varids<-paste(variant_info$chr,variant_info$pos,variant_info$ref,variant_info$alt,sep=":")

  for (lvar in 0:length(varids)){
  pvalues <- NA

  if(lvar==0){

  try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
  rmvar<-"None"
  }else{
  cat("  remove number",lvar,"variant out of",length(varids),"\n")
  Geno2<-Geno[,-c(lvar)]
  Anno.Int.PHRED.sub2<-Anno.Int.PHRED.sub[-c(lvar),]
  rmvar<-varids[lvar]
  try(pvalues <- STAAR(Geno2,obj_nullmodel,Anno.Int.PHRED.sub2,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)

  }
	if(class(pvalues)=="list")
	{
		results_temp <- as.vector(genes[1,])
		results_temp[3] <- category
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[1,1])
    results_temp[4] <- rmvar
    results_temp[5] <- pvalues$num_variant


		results_temp <- c(results_temp,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
		pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
		pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)

		results <- rbind(results,results_temp)
	}
  }

	if(!is.null(results))
	{
		colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
		colnames(results)[1:5] <- c("Gene name","Chr","Category","rm_variant","#SNV")
		colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
	}
  return(results)
}
