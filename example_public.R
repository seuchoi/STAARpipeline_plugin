## A example script to perform conditional anlaysis
## git clone https://github.com/seuchoi/STAARpipeline_plugin.git

### load your packages
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(data.table)
library(readr)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(R.utils)
sourceDirectory("/your/github/directory/STAARpipeline_plugin/","*.R")

### load your files
genofile=seqOpen("your gdsfile")
obj_nullmodel=get(load("your null model file.RData"))
obj_nullmodel$Sigma_i <- as(obj_nullmodel$Sigma_i, 'dsCMatrix') # optional
Annotation_name_catalog<-get(load("/your annotation catalog file/Annotation_name_catalog.Rdata"))
Annotation_name<-c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription","aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
method_cond = "optimal"

## add your parameters
chr=xx;gene_name="yyy";gene_name_cond="zzz";category="yyy-cat";cond_category="zzz-cat";genofile=genofile;obj_nullmodel=obj_nullmodel;rare_maf_cutoff=0.01;rv_num_cutoff=2;QC_label="annotation/filter";variant_type="SNV";geno_missing_imputation=c("minor");Annotation_dir="annotation/info/FunctionalAnnotation";Annotation_name_catalog;Use_annotation_weights=TRUE;Annotation_name=Annotation_name

## extract variant list in gene zzz
res<-coding_extract(chr=chr,gene_name=gene_name_cond,genofile=genofile,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=silent)


## performing conditional analysis
chr=xx;gene_name="yyy";gene_name_cond="zzz";category="yyy-cat";cond_category="zzz-cat";genofile=genofile;obj_nullmodel=obj_nullmodel;rare_maf_cutoff=0.01;rv_num_cutoff=2;QC_label="annotation/filter";variant_type="SNV";geno_missing_imputation=c("minor");Annotation_dir="annotation/info/FunctionalAnnotation";Annotation_name_catalog;Use_annotation_weights=TRUE;Annotation_name=Annotation_name

## perfroming original anlaysis
yyy_ori<-Gene_Centric_Noncoding(chr=14,gene_name="yyy",category="yyy-cat",genofile=genofile,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,QC_label="annotation/filter",variant_type="SNV",geno_missing_imputation=c("minor"),Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,Use_annotation_weights=TRUE,Annotation_name=Annotation_name)

## extract variant list
res<-coding_extract(chr=chr,gene_name=gene_name_cond,genofile=genofile,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=silent)

## extract variant for your category
varlist0<-subset(res,type==cond_category)
cond_variants<-varlist0[,c("chr","pos","ref","alt")]
known_loci<-cond_variants

## performing conditional analysis
results <- enhancer_DHS_cond_anyvar(chr, gene_name, genofile,
            obj_nullmodel, known_loci, rare_maf_cutoff = rare_maf_cutoff,
            rv_num_cutoff = rv_num_cutoff, method_cond = method_cond,
            QC_label = QC_label, variant_type = variant_type,
            geno_missing_imputation = geno_missing_imputation,
            Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog,
            Use_annotation_weights = Use_annotation_weights,
            Annotation_name = Annotation_name)
yyy_cond<-results


## create figure for original and conditional analyses
png("yyy_conditional_figure.png",width = 2500, height = 1000,res=100)
par(mar = c(20, 5, 2, 2))

STAARpipeline_pvalue_fig(result=yyy_ori,cond_result=yyy_cond)
dev.off()
