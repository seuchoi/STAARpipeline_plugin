

coding_extract<-function(chr,gene_name,genofile,obj_nullmodel,genes,
                   rare_maf_cutoff=0.01,rv_num_cutoff=2,
                   QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                   Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                   Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,silent=FALSE){


    genes=genes_info[genes_info[,2]==chr,]

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- as.character(obj_nullmodel$id_include)

	## get SNV id, position, REF, ALT (whole genome)
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

	position <- as.numeric(seqGetData(genofile, "position"))
	variant.id <- seqGetData(genofile, "variant.id")

	rm(filter)
	gc()

	### Gene
	kk <- which(genes[,1]==gene_name)

	sub_start_loc <- genes[kk,3]
	sub_end_loc <- genes[kk,4]

	is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
	variant.id.gene <- variant.id[is.in]

	rm(position)
	gc()

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Gencode
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
	## Meta.SVM.Pred
	MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

	################################################
	#           Coding
	################################################
	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.coding <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.EXONIC.Category=="nonsynonymous SNV")|(GENCODE.EXONIC.Category=="synonymous SNV")
	variant.id.gene <- variant.id.gene[lof.in.coding]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Gencode
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
	## Meta.SVM.Pred
	MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

	## Annotation
	Anno.Int.PHRED.sub <- NULL
	Anno.Int.PHRED.sub.name <- NULL

	if(variant_type=="SNV")
	{
		if(Use_annotation_weights)
		{
			for(k in 1:length(Annotation_name))
			{
				if(Annotation_name[k]%in%Annotation_name_catalog$name)
				{
					Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
					Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

					if(Annotation_name[k]=="CADD")
					{
						Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
					}

					if(Annotation_name[k]=="aPC.LocalDiversity")
					{
						Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
						Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
						Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
					}
					Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
				}
			}

			Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
			colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
		}
	}

	################################################
	#                  1) plof_ds
	################################################
	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
	variant.id.gene.category <- variant.id.gene[lof.in.plof]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	#id.genotype.match <- rep(0,length(id.genotype))

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
    carrier0$type<-"plof_ds"
    plof_ds_carriers<-carrier0
    }else{
    plof_ds_carriers<-NULL    
    }

	#####################################################
	#                      plof
	#####################################################
	# variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
	variant.id.gene.category <- variant.id.gene[lof.in.plof]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

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
    carrier0$type<-"plof"
    plof_carriers<-carrier0
    }else{
    plof_carriers<-NULL    
    }


	#############################################
	#             synonymous
	############################################# #### ERROR SHC
	lof.in.synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
	variant.id.gene.category <- variant.id.gene[lof.in.synonymous]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

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
    carrier0$type<-"synonymous"
    synonymous_carriers<-carrier0
    }else{
    synonymous_carriers<-NULL    
    }

	#################################################
	#        missense
	#################################################
	lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
	variant.id.gene.category <- variant.id.gene[lof.in.missense]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

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
    carrier0$type<-"missense"
    missense_carriers<-carrier0
    }else{
    missense_carriers<-NULL    
    }

	#################################################
	#         disruptive missense
	#################################################
	lof.in.dmissense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D")
	variant.id.gene.category <- variant.id.gene[lof.in.dmissense]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

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
    carrier0$type<-"disruptive_missense"
    disruptive_missense_carriers<-carrier0
    }else{
    disruptive_missense_carriers<-NULL    
    }


    comb_carriers<-rbind(plof_carriers,plof_ds_carriers,disruptive_missense_carriers,missense_carriers,synonymous_carriers)
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
