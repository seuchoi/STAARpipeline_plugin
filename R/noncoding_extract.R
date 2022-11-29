
noncoding_extract<-function(chr,gene_name,genofile,obj_nullmodel,genes,
                   rare_maf_cutoff=0.01,rv_num_cutoff=2,
                   QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                   Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                   Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,silent=FALSE){

    genes=genes_info[genes_info[,2]==chr,]

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- as.character(obj_nullmodel$id_include)

	#####################################
	#   Gene Info
	## get SNV id
	filter <- seqGetData(genofile, QC_label)
	if(variant_type=="variant")
	{
		SNVlist <- filter == "PASS"
	}

	if(variant_type=="SNV")
	{
		SNVlist <- (filter == "PASS") & isSNV(genofile)
	}

	if(variant_type=="Indel")
	{
		SNVlist <- (filter == "PASS") & (!isSNV(genofile))
	}

	variant.id <- seqGetData(genofile, "variant.id")

	rm(filter)
	gc()

	########################################
	#   Downstream

	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
	is.in <- (GENCODE.Category=="downstream")&(SNVlist)
	variant.id.downstream <- variant.id[is.in]

	seqSetFilter(genofile,variant.id=variant.id.downstream,sample.id=phenotype.id)

	rm(variant.id.downstream)
	gc()

	GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
	GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
	variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))

	variant.id.SNV <- seqGetData(genofile, "variant.id")
	variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)

	rm(GENCODE.Info)
	gc()

	rm(variant_gene_num)
	gc()

	Gene <- as.character(unlist(GENCODE.Info.split))

	rm(GENCODE.Info.split)
	gc()

	seqResetFilter(genofile)

	### Gene
	is.in <- which(Gene==gene_name)
	variant.is.in <- variant.id.SNV[is.in]

	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]

	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}

  
    ## SHC Add
    id.variant <- seqGetData(genofile,"variant.id")
      
    carrier0<-NULL
    for (vv in 1:length(id.variant)){
    idnums<-which(Geno[,vv]>=1)
    if(length(idnums)>0){
    id.carrier<-phenotype.id.merge[idnums,"phenotype.id"]
    var.carrier<-id.variant[vv]
    carrier1<-data.frame(sample.id=id.carrier,variant.id=var.carrier)        
    carrier0<-rbind(carrier0,carrier1)    
    }
    }
    if(!is.null(carrier0)){
    carrier0$type<-"downstream"
    downstream_carriers<-carrier0
    }else{
    downstream_carriers<-NULL
    }

  
	seqResetFilter(genofile)

	########################################
	#   Upstream

	is.in <- (GENCODE.Category=="upstream")&(SNVlist)
	variant.id.upstream <- variant.id[is.in]

	seqSetFilter(genofile,variant.id=variant.id.upstream,sample.id=phenotype.id)

	rm(variant.id.upstream)
	gc()

	GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
	GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
	variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))

	variant.id.SNV <- seqGetData(genofile, "variant.id")
	variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)

	rm(GENCODE.Info)
	gc()

	rm(variant_gene_num)
	gc()

	Gene <- as.character(unlist(GENCODE.Info.split))

	rm(GENCODE.Info.split)
	gc()

	seqResetFilter(genofile)

	### Gene
	is.in <- which(Gene==gene_name)
	variant.is.in <- variant.id.SNV[is.in]

	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]

	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
    
    }
	}
  
    ## SHC Add
    id.variant <- seqGetData(genofile,"variant.id")
      
    carrier0<-NULL
    for (vv in 1:length(id.variant)){
    idnums<-which(Geno[,vv]>=1)
    if(length(idnums)>0){
    id.carrier<-phenotype.id.merge[idnums,"phenotype.id"]
    var.carrier<-id.variant[vv]
    carrier1<-data.frame(sample.id=id.carrier,variant.id=var.carrier)        
    carrier0<-rbind(carrier0,carrier1)    
    }
    }
    if(!is.null(carrier0)){
    carrier0$type<-"upstream"
    upstream_carriers<-carrier0
    }else{
    upstream_carriers<-NULL
    }

	seqResetFilter(genofile)

	########################################################
	#                UTR

	is.in <- ((GENCODE.Category=="UTR3")|(GENCODE.Category=="UTR5")|(GENCODE.Category=="UTR5;UTR3"))&(SNVlist)
	variant.id.UTR <- variant.id[is.in]

	rm(GENCODE.Category)
	gc()

	seqSetFilter(genofile,variant.id=variant.id.UTR,sample.id=phenotype.id)

	rm(variant.id.UTR)
	gc()

	GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
	GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[(]")

	rm(GENCODE.Info)
	gc()

	Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[1]))

	rm(GENCODE.Info.split)
	gc()

	variant.id.SNV <- seqGetData(genofile, "variant.id")

	seqResetFilter(genofile)

	### Gene
	is.in <- which(Gene==gene_name)
	variant.is.in <- variant.id.SNV[is.in]

	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]

	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}

  
    ## SHC Add
    id.variant <- seqGetData(genofile,"variant.id")
      
    carrier0<-NULL
    for (vv in 1:length(id.variant)){
    idnums<-which(Geno[,vv]>=1)
    if(length(idnums)>0){
    id.carrier<-phenotype.id.merge[idnums,"phenotype.id"]
    var.carrier<-id.variant[vv]
    carrier1<-data.frame(sample.id=id.carrier,variant.id=var.carrier)        
    carrier0<-rbind(carrier0,carrier1)    
    }
    }
    if(!is.null(carrier0)){
    carrier0$type<-"UTR"
    UTR_carriers<-carrier0
    }else{
    UTR_carriers<-NULL
    }

    
	seqResetFilter(genofile)

	#############################################
	#   Promoter-CAGE

	## Promoter
	varid <- seqGetData(genofile, "variant.id")
	txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
	promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)

	# Subsetting Promoters that within +/-3kb of TSS and have CAGE signals
	CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
	CAGEBvt <- CAGEAnno!=""
	CAGEidx <- which(CAGEBvt,useNames=TRUE)
	seqSetFilter(genofile,variant.id=varid[CAGEidx])
	seqSetFilter(genofile,promGobj,intersect=TRUE)
	CAGEpromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
	CAGEGene <- unlist(lapply(strsplit(CAGEpromgene,"\\(|\\,|;|-"),`[[`,1))
	##obtain variants info
	CAGEvchr <- as.numeric(seqGetData(genofile,"chromosome"))
	CAGEvpos <- as.numeric(seqGetData(genofile,"position"))
	CAGEvref <- as.character(seqGetData(genofile,"$ref"))
	CAGEvalt <- as.character(seqGetData(genofile,"$alt"))
	dfPromCAGEVarGene <- data.frame(CAGEvchr,CAGEvpos,CAGEvref,CAGEvalt,CAGEGene)

	## get SNV id
	filter <- seqGetData(genofile, QC_label)
	if(variant_type=="variant")
	{
		SNVlist <- filter == "PASS"
	}

	if(variant_type=="SNV")
	{
		SNVlist <- (filter == "PASS") & isSNV(genofile)
	}

	if(variant_type=="Indel")
	{
		SNVlist <- (filter == "PASS") & (!isSNV(genofile))
	}

	variant.id <- seqGetData(genofile, "variant.id")
	variant.id.SNV <- variant.id[SNVlist]

	dfPromCAGEVarGene.SNV <- dfPromCAGEVarGene[SNVlist,]
	dfPromCAGEVarGene.SNV$CAGEvpos <- as.character(dfPromCAGEVarGene.SNV$CAGEvpos)
	dfPromCAGEVarGene.SNV$CAGEvref <- as.character(dfPromCAGEVarGene.SNV$CAGEvref)
	dfPromCAGEVarGene.SNV$CAGEvalt <- as.character(dfPromCAGEVarGene.SNV$CAGEvalt)

	seqResetFilter(genofile)

	rm(dfPromCAGEVarGene)
	gc()

	### Gene
	is.in <- which(dfPromCAGEVarGene.SNV[,5]==gene_name)
	variant.is.in <- variant.id.SNV[is.in]

	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]

	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}

  
    ## SHC Add
    id.variant <- seqGetData(genofile,"variant.id")
      
    carrier0<-NULL
    for (vv in 1:length(id.variant)){
    idnums<-which(Geno[,vv]>=1)
    if(length(idnums)>0){
    id.carrier<-phenotype.id.merge[idnums,"phenotype.id"]
    var.carrier<-id.variant[vv]
    carrier1<-data.frame(sample.id=id.carrier,variant.id=var.carrier)        
    carrier0<-rbind(carrier0,carrier1)    
    }
    }
    if(!is.null(carrier0)){
    carrier0$type<-"promoter_CAGE"
    promoter_CAGE_carriers<-carrier0
    }else{
    promoter_CAGE_carriers<-NULL
    }

	seqResetFilter(genofile)

	##################################################
	#       Promoter-DHS

	# Subsetting Promoters that within +/-3kb of TSS and have rOCRs signals
	rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
	rOCRsBvt <- rOCRsAnno!=""
	rOCRsidx <- which(rOCRsBvt,useNames=TRUE)
	seqSetFilter(genofile,variant.id=varid[rOCRsidx])

	seqSetFilter(genofile,promGobj,intersect=TRUE)
	rOCRspromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
	rOCRsGene <- unlist(lapply(strsplit(rOCRspromgene,"\\(|\\,|;|-"),`[[`,1))
	## obtain variants info
	rOCRsvchr <- as.numeric(seqGetData(genofile,"chromosome"))
	rOCRsvpos <- as.numeric(seqGetData(genofile,"position"))
	rOCRsvref <- as.character(seqGetData(genofile,"$ref"))
	rOCRsvalt <- as.character(seqGetData(genofile,"$alt"))
	dfPromrOCRsVarGene <- data.frame(rOCRsvchr,rOCRsvpos,rOCRsvref,rOCRsvalt,rOCRsGene)

	## get SNV id
	filter <- seqGetData(genofile, QC_label)
	if(variant_type=="variant")
	{
		SNVlist <- filter == "PASS"
	}

	if(variant_type=="SNV")
	{
		SNVlist <- (filter == "PASS") & isSNV(genofile)
	}

	if(variant_type=="Indel")
	{
		SNVlist <- (filter == "PASS") & (!isSNV(genofile))
	}

	variant.id <- seqGetData(genofile, "variant.id")
	variant.id.SNV <- variant.id[SNVlist]

	dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist,]
	dfPromrOCRsVarGene.SNV$rOCRsvpos <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvpos)
	dfPromrOCRsVarGene.SNV$rOCRsvref <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvref)
	dfPromrOCRsVarGene.SNV$rOCRsvalt <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvalt)

	seqResetFilter(genofile)

	rm(dfPromrOCRsVarGene)
	gc()

	### Gene
	is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
	variant.is.in <- variant.id.SNV[is.in]

	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]

	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}

  
    ## SHC Add
    id.variant <- seqGetData(genofile,"variant.id")
      
    carrier0<-NULL
    for (vv in 1:length(id.variant)){
    idnums<-which(Geno[,vv]>=1)
    if(length(idnums)>0){
    id.carrier<-phenotype.id.merge[idnums,"phenotype.id"]
    var.carrier<-id.variant[vv]
    carrier1<-data.frame(sample.id=id.carrier,variant.id=var.carrier)        
    carrier0<-rbind(carrier0,carrier1)    
    }
    }
    if(!is.null(carrier0)){
    carrier0$type<-"promotor_DHS"
    promoter_DHS_carriers<-carrier0
    }else{
    promoter_DHS_carriers<-NULL
    }

	seqResetFilter(genofile)

	###########################################
	#        Enhancer-CAGE

	#Now extract the GeneHancer with CAGE Signal Overlay
	genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
	genehancer <- genehancerAnno!=""

	CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
	CAGE <- CAGEAnno!=""
	CAGEGeneHancervt <- CAGEAnno!=""&genehancerAnno!=""
	CAGEGeneHanceridx <- which(CAGEGeneHancervt,useNames=TRUE)
	seqSetFilter(genofile,variant.id=varid[CAGEGeneHanceridx])

	# variants that covered by whole GeneHancer without CAGE overlap.
	genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
	enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
	enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
	enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
	enhancervpos <- as.numeric(seqGetData(genofile,"position"))
	enhancervref <- as.character(seqGetData(genofile,"$ref"))
	enhancervalt <- as.character(seqGetData(genofile,"$alt"))
	dfHancerVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

	## get SNV id
	filter <- seqGetData(genofile, QC_label)
	if(variant_type=="variant")
	{
		SNVlist <- filter == "PASS"
	}

	if(variant_type=="SNV")
	{
		SNVlist <- (filter == "PASS") & isSNV(genofile)
	}

	if(variant_type=="Indel")
	{
		SNVlist <- (filter == "PASS") & (!isSNV(genofile))
	}

	variant.id <- seqGetData(genofile, "variant.id")
	variant.id.SNV <- variant.id[SNVlist]

	dfHancerVarGene.SNV <- dfHancerVarGene[SNVlist,]
	dfHancerVarGene.SNV$enhancervpos <- as.character(dfHancerVarGene.SNV$enhancervpos)
	dfHancerVarGene.SNV$enhancervref <- as.character(dfHancerVarGene.SNV$enhancervref)
	dfHancerVarGene.SNV$enhancervalt <- as.character(dfHancerVarGene.SNV$enhancervalt)

	seqResetFilter(genofile)

	rm(dfHancerVarGene)
	gc()

	### Gene
	is.in <- which(dfHancerVarGene.SNV[,5]==gene_name)
	variant.is.in <- variant.id.SNV[is.in]

	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]

	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}

  
    ## SHC Add
    id.variant <- seqGetData(genofile,"variant.id")
      
    carrier0<-NULL
    for (vv in 1:length(id.variant)){
    idnums<-which(Geno[,vv]>=1)
    if(length(idnums)>0){
    id.carrier<-phenotype.id.merge[idnums,"phenotype.id"]
    var.carrier<-id.variant[vv]
    carrier1<-data.frame(sample.id=id.carrier,variant.id=var.carrier)        
    carrier0<-rbind(carrier0,carrier1)    
    }
    }
    if(!is.null(carrier0)){
    carrier0$type<-"enhancer_CAGE"
    enhancer_CAGE_carriers<-carrier0
    }else{
    enhancer_CAGE_carriers<-NULL
    }

	seqResetFilter(genofile)

	##################################################
	#       Enhancer-DHS

	rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
	rOCRs <- rOCRsAnno!=""
	rOCRsGeneHancervt <- rOCRsAnno!=""&genehancerAnno!=""
	rOCRsGeneHanceridx <- which(rOCRsGeneHancervt,useNames=TRUE)
	seqSetFilter(genofile,variant.id=varid[rOCRsGeneHanceridx])
	# variants that covered by whole GeneHancer without rOCRs overlap.

	genehancerSet <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
	enhancerGene <- unlist(lapply(strsplit(genehancerSet,"="),`[[`,4))
	enhancer2GENE <- unlist(lapply(strsplit(enhancerGene,";"),`[[`,1))
	enhancervchr <- as.numeric(seqGetData(genofile,"chromosome"))
	enhancervpos <- as.numeric(seqGetData(genofile,"position"))
	enhancervref <- as.character(seqGetData(genofile,"$ref"))
	enhancervalt <- as.character(seqGetData(genofile,"$alt"))
	dfHancerVarGene <- data.frame(enhancervchr,enhancervpos,enhancervref,enhancervalt,enhancer2GENE)

	rm(varid)
	gc()

	## get SNV id
	filter <- seqGetData(genofile, QC_label)
	if(variant_type=="variant")
	{
		SNVlist <- filter == "PASS"
	}

	if(variant_type=="SNV")
	{
		SNVlist <- (filter == "PASS") & isSNV(genofile)
	}

	if(variant_type=="Indel")
	{
		SNVlist <- (filter == "PASS") & (!isSNV(genofile))
	}

	variant.id <- seqGetData(genofile, "variant.id")
	variant.id.SNV <- variant.id[SNVlist]

	dfHancerVarGene.SNV <- dfHancerVarGene[SNVlist,]
	dfHancerVarGene.SNV$enhancervpos <- as.character(dfHancerVarGene.SNV$enhancervpos)
	dfHancerVarGene.SNV$enhancervref <- as.character(dfHancerVarGene.SNV$enhancervref)
	dfHancerVarGene.SNV$enhancervalt <- as.character(dfHancerVarGene.SNV$enhancervalt)

	seqResetFilter(genofile)

	rm(dfHancerVarGene)
	gc()

	### Gene
	is.in <- which(dfHancerVarGene.SNV[,5]==gene_name)
	variant.is.in <- variant.id.SNV[is.in]

	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]

	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}
  
    ## SHC Add
    id.variant <- seqGetData(genofile,"variant.id")
      
    carrier0<-NULL
    for (vv in 1:length(id.variant)){
    idnums<-which(Geno[,vv]>=1)
    if(length(idnums)>0){
    id.carrier<-phenotype.id.merge[idnums,"phenotype.id"]
    var.carrier<-id.variant[vv]
    carrier1<-data.frame(sample.id=id.carrier,variant.id=var.carrier)        
    carrier0<-rbind(carrier0,carrier1)    
    }
    }
    if(!is.null(carrier0)){
    carrier0$type<-"enhancer_DHS"
    enhancer_DHS_carriers<-carrier0
    }else{
    enhancer_DHS_carriers<-NULL
    }

	seqResetFilter(genofile)


  
  	############################################
  	#           results

  comb_carriers<-rbind(UTR_carriers,upstream_carriers,downstream_carriers,enhancer_CAGE_carriers,enhancer_DHS_carriers,promoter_CAGE_carriers,promoter_DHS_carriers)
  comb_carriers$index<-c(1:nrow(comb_carriers))
  all.variant.id<-comb_carriers$variant.id
  seqSetFilter(genofile,variant.id=all.variant.id,sample.id=phenotype.id)

  all.id.variant <- seqGetData(genofile,"variant.id")
  freq0<-seqAlleleFreq(genofile,1L)
  freq1<-data.frame(variant.id=all.id.variant,alt.freq=freq0)
  varinfo0<-variantInfo(genofile)
  varinfo0<-merge(varinfo0,freq1,by="variant.id")


  comb_carriers2<-merge(comb_carriers,varinfo0,by="variant.id")
  comb_carriers2<-comb_carriers2[order(comb_carriers2$index),]
  comb_carriers3<-subset(comb_carriers2,alt.freq<rare_maf_cutoff)
  comb_carriers4<-subset(comb_carriers3,select=-c(index))
  comb_carriers4<-data.frame(comb_carriers4,row.names = NULL)

seqResetFilter(genofile)

return(comb_carriers4)
}
