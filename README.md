# STAARpipeline_plugin
Additional features of STAARpipeline

## software versions
```
R version 4.2.1
GENESIS_2.26.0
STAAR_0.9.6.1
STAARpipeline_0.9.6
```
## variant list extraction
this function extracts a list of variants used in the STAAR anlaysis.
```
coding_extract(chr=chr,gene_name=gene_name_cond,
            genofile=genofile,obj_nullmodel=obj_nullmodel,
            rare_maf_cutoff=rare_maf_cutoff,
            rv_num_cutoff=rv_num_cutoff,
            QC_label=QC_label,variant_type=variant_type,
            geno_missing_imputation=geno_missing_imputation,A
            nnotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
            Use_annotation_weights=Use_annotation_weights,
            Annotation_name=Annotation_name,silent=silent)
```

## conditional analaysis without regional limitaion
this function perfroms a conditional anlaysis without regional resctrictions. (currently, only implemented for gene centric noncoding analysis)
```
enhancer_DHS_cond_anyvar(chr, gene_name, genofile,
            obj_nullmodel, known_loci, rare_maf_cutoff = rare_maf_cutoff,
            rv_num_cutoff = rv_num_cutoff, method_cond = method_cond,
            QC_label = QC_label, variant_type = variant_type,
            geno_missing_imputation = geno_missing_imputation,
            Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog,
            Use_annotation_weights = Use_annotation_weights,
            Annotation_name = Annotation_name)
```
### leave one variant out analaysis
this function performs a leave one variant out anlaysis to understand the contibutaion of the each variant. (currently, this function implemtned for disruptive_missense and promoter_DHS)
```
promoter_DHS_LOVO(chr=chr,gene_name=gene_name,
          category=category, genofile=genofile,
          obj_nullmodel=obj_nullmodel, genes=genes,
          rare_maf_cutoff=0.01, rv_num_cutoff=2,
          QC_label="annotation/filter",
          variant_type="SNV",
          geno_missing_imputation=c("minor"),
          Annotation_dir="annotation/info/FunctionalAnnotation",
          Annotation_name_catalog=Annotation_name_catalog,
          Use_annotation_weights=TRUE,
          Annotation_name=Annotation_name)
