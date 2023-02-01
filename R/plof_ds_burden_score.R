plof_ds_burden_score<-function (chr, gene_name, genofile, obj_nullmodel, genes, rare_maf_cutoff = 0.01, 
    rv_num_cutoff = 2, QC_label = "annotation/filter", variant_type = c("SNV",
        "Indel", "variant"), geno_missing_imputation = c("mean",
        "minor"), Annotation_dir = "annotation/info/FunctionalAnnotation",
    Annotation_name_catalog, Use_annotation_weights = c(TRUE,
        FALSE), Annotation_name = NULL, silent = FALSE)
{
    variant_type <- match.arg(variant_type)
    geno_missing_imputation <- match.arg(geno_missing_imputation)
    phenotype.id <- as.character(obj_nullmodel$id_include)
    filter <- seqGetData(genofile, QC_label)
    if (variant_type == "variant") {
        SNVlist <- filter == "PASS"
    }
    if (variant_type == "SNV") {
        SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    if (variant_type == "Indel") {
        SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    position <- as.numeric(seqGetData(genofile, "position"))
    variant.id <- seqGetData(genofile, "variant.id")
    rm(filter)
    gc()
    kk <- which(genes[, 1] == gene_name)
    sub_start_loc <- genes[kk, 3]
    sub_end_loc <- genes[kk, 4]
    is.in <- (SNVlist) & (position >= sub_start_loc) & (position <=
        sub_end_loc)
    variant.id.gene <- variant.id[is.in]
    seqSetFilter(genofile, variant.id = variant.id.gene, sample.id = phenotype.id)
    GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,
        Annotation_name_catalog$dir[which(Annotation_name_catalog$name ==
            "GENCODE.EXONIC.Category")]))
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,
        Annotation_name_catalog$dir[which(Annotation_name_catalog$name ==
            "GENCODE.Category")]))
    MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,
        Annotation_name_catalog$dir[which(Annotation_name_catalog$name ==
            "MetaSVM")]))
    variant.id.gene <- seqGetData(genofile, "variant.id")
    lof.in.plof <- (GENCODE.EXONIC.Category == "stopgain") |
        (GENCODE.EXONIC.Category == "stoploss") | (GENCODE.Category ==
        "splicing") | (GENCODE.Category == "exonic;splicing") |
        (GENCODE.Category == "ncRNA_splicing") | (GENCODE.Category ==
        "ncRNA_exonic;splicing") | ((GENCODE.EXONIC.Category ==
        "nonsynonymous SNV") & (MetaSVM_pred == "D"))
    variant.id.gene <- variant.id.gene[lof.in.plof]
    seqSetFilter(genofile, variant.id = variant.id.gene, sample.id = phenotype.id)
    id.genotype <- seqGetData(genofile, "sample.id")
    id.genotype.merge <- data.frame(id.genotype, index = seq(1,
        length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,
        id.genotype.merge, by = c(phenotype.id = "id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match, , drop = FALSE]
    if (!is.null(dim(Geno))) {
        if (dim(Geno)[2] > 0) {
            if (geno_missing_imputation == "mean") {
                Geno <- matrix_flip_mean(Geno)$Geno
            }
            if (geno_missing_imputation == "minor") {
                Geno <- matrix_flip_minor(Geno)$Geno
            }
        }
    }
    Anno.Int.PHRED.sub <- NULL
    Anno.Int.PHRED.sub.name <- NULL
    if (variant_type == "SNV") {
        if (Use_annotation_weights) {
            for (k in 1:length(Annotation_name)) {
                if (Annotation_name[k] %in% Annotation_name_catalog$name) {
                  Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,
                    Annotation_name[k])
                  Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,
                    Annotation_name_catalog$dir[which(Annotation_name_catalog$name ==
                      Annotation_name[k])]))
                  if (Annotation_name[k] == "CADD") {
                    Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
                  }
                  if (Annotation_name[k] == "aPC.LocalDiversity") {
                    Annotation.PHRED.2 <- -10 * log10(1 - 10^(-Annotation.PHRED/10))
                    Annotation.PHRED <- cbind(Annotation.PHRED,
                      Annotation.PHRED.2)
                    Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,
                      paste0(Annotation_name[k], "(-)"))
                  }
                  Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,
                    Annotation.PHRED)
                }
            }
            Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
            colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
        }
    }
    pvalues <- 0
    try(pvalues <- STAAR_burden_score(Geno, obj_nullmodel, Anno.Int.PHRED.sub,
        rare_maf_cutoff = rare_maf_cutoff, rv_num_cutoff = rv_num_cutoff),
        silent = silent)
    if(!is.null(pvalues)){
    pvalues[["Category"]]<-"plof_ds"
    pvalues[["Chr"]]<-chr
    pvalues[["Gene Name"]]<-as.character(genes[kk, 1])
    }

    seqResetFilter(genofile)
    return(pvalues)

}
