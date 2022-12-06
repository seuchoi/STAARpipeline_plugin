cd /restricted/projectnb/adsp-charge/seuchoi/rarevariant/script/STAAR/SNV/Coding/

module load R/4.2.1
module load staar/0.9.6.1
module load xsv/0.13.0

git clone https://github.com/seuchoi/STAARpipeline_plugin.git


library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(data.table)
library(readr)
library(SeqArray)
library(SeqVarTools)
#library(STAAR)
library(STAAR,lib.loc="/restricted/projectnb/adsp-charge/seuchoi/software/")
library(STAARpipeline,lib.loc="/restricted/projectnb/adsp-charge/seuchoi/software/")
#library(STAARpipeline)
library(R.utils)
sourceDirectory("/restricted/projectnb/adsp-charge/seuchoi/software/STAARpipeline_plugin/R/","*.R")


genofile=seqOpen("/restricted/projectnb/adsp-charge/seuchoi/rarevariant/data/gds/biallelic_chr14.gds")
obj_nullmodel=get(load("/restricted/projectnb/adsp-charge/seuchoi/rarevariant/data/17K_Null_wGRM_STAAR.RData"))
obj_nullmodel$Sigma_i <- as(obj_nullmodel$Sigma_i, 'dsCMatrix')
Annotation_name_catalog<-get(load("/restricted/projectnb/adsp-charge/seuchoi/rarevariant/result/STAAR/Annotation_name_catalog.Rdata"))
Annotation_name<-c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription","aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
method_cond = "optimal"

## performing conditional anlaysis
chr=14;gene_name="EIF2B2";gene_name_cond="PSEN1";category="enhancer_DHS";cond_category="plof_ds";genofile=genofile;obj_nullmodel=obj_nullmodel;rare_maf_cutoff=0.01;rv_num_cutoff=2;QC_label="annotation/filter";variant_type="SNV";geno_missing_imputation=c("minor");Annotation_dir="annotation/info/FunctionalAnnotation";Annotation_name_catalog;Use_annotation_weights=TRUE;Annotation_name=Annotation_name

## extract variant list
res<-coding_extract(chr=chr,gene_name=gene_name_cond,genofile=genofile,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=silent)




source("/restricted/projectnb/adsp-charge/seuchoi/software/STAARpipieline_plugin/R/coding_extract.R")
source("/restricted/projectnb/adsp-charge/seuchoi/software/STAARpipieline_plugin/R/noncoding_extract.R")
source("/restricted/projectnb/adsp-charge/seuchoi/software/STAARpipieline_plugin/R/STAARpipeline_cond.R")




## performing conditional analysis
chr=14;gene_name="EIF2B2";gene_name_cond="PSEN1";category="enhancer_DHS";cond_category="plof_ds";genofile=genofile;obj_nullmodel=obj_nullmodel;rare_maf_cutoff=0.01;rv_num_cutoff=2;QC_label="annotation/filter";variant_type="SNV";geno_missing_imputation=c("minor");Annotation_dir="annotation/info/FunctionalAnnotation";Annotation_name_catalog;Use_annotation_weights=TRUE;Annotation_name=Annotation_name

## extract variant list
res<-coding_extract(chr=chr,gene_name=gene_name_cond,genofile=genofile,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=silent)

######
varlist0<-subset(res,type==cond_category)
cond_variants<-varlist0[,c("chr","pos","ref","alt")]


known_loci<-cond_variants


results <- enhancer_DHS_cond_anyvar(chr, gene_name, genofile,
            obj_nullmodel, known_loci, rare_maf_cutoff = rare_maf_cutoff,
            rv_num_cutoff = rv_num_cutoff, method_cond = method_cond,
            QC_label = QC_label, variant_type = variant_type,
            geno_missing_imputation = geno_missing_imputation,
            Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog,
            Use_annotation_weights = Use_annotation_weights,
            Annotation_name = Annotation_name)

EIF2B2_cond<-results

EIF2B2_ori<-Gene_Centric_Noncoding(chr=14,gene_name="EIF2B2",category="enhancer_DHS",genofile=genofile,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,QC_label="annotation/filter",variant_type="SNV",geno_missing_imputation=c("minor"),Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,Use_annotation_weights=TRUE,Annotation_name=Annotation_name)

setwd("/restricted/projectnb/adsp-charge/seuchoi/rarevariant/result/STAAR/SNV/Noncoding/maf0.01/figure")
png("EIF2B2_conditional_figure.png",width = 2500, height = 1000,res=100)
par(mar = c(20, 5, 2, 2))

STAARpipeline_pvalue_fig(result=EIF2B2_ori,cond_result=EIF2B2_cond)
dev.off()
