STAAR_burden_score<-function (genotype, obj_nullmodel, annotation_phred = NULL, rare_maf_cutoff = 0.01, rv_num_cutoff = 2){
    if (class(genotype)[1] != "matrix" && !(!is.null(attr(class(genotype),
        "package")) && attr(class(genotype), "package") == "Matrix")) {
        stop("genotype is not a matrix!")
    }
    if (dim(genotype)[2] == 1) {
        stop(paste0("Number of rare variant in the set is less than 2!"))
    }
    annotation_phred <- as.data.frame(annotation_phred)
    if (dim(annotation_phred)[1] != 0 & dim(genotype)[2] != dim(annotation_phred)[1]) {
        stop(paste0("Dimensions don't match for genotype and annotation!"))
    }
    if (!is.null(attr(class(genotype), "package")) && attr(class(genotype),
        "package") == "Matrix") {
        genotype <- as.matrix(genotype)
    }
    genotype <- matrix_flip(genotype)
    MAF <- genotype$MAF
    RV_label <- as.vector((MAF < rare_maf_cutoff) & (MAF > 0))
    Geno_rare <- genotype$Geno[, RV_label]
    rm(genotype)
    gc()
    annotation_phred <- annotation_phred[RV_label, , drop = FALSE]

    if (sum(RV_label) >= rv_num_cutoff) {
        G <- as(Geno_rare, "dgCMatrix")
        MAF <- MAF[RV_label]
        rm(Geno_rare)
        gc()
        annotation_rank <- 1 - 10^(-annotation_phred/10)
        w_1 <- dbeta(MAF, 1, 25)
        w_2 <- dbeta(MAF, 1, 1)
        if (dim(annotation_phred)[2] == 0) {
            w_B <- w_S <- as.matrix(cbind(w_1, w_2))
            w_A <- as.matrix(cbind(w_1^2/dbeta(MAF, 0.5, 0.5)^2,
                w_2^2/dbeta(MAF, 0.5, 0.5)^2))
        }else{
            w_B_1 <- annotation_rank * w_1
            w_B_1 <- cbind(w_1, w_B_1)
            w_B_2 <- annotation_rank * w_2
            w_B_2 <- cbind(w_2, w_B_2)
            w_B <- cbind(w_B_1, w_B_2)
            w_B <- as.matrix(w_B)
            w_S_1 <- sqrt(annotation_rank) * w_1
            w_S_1 <- cbind(w_1, w_S_1)
            w_S_2 <- sqrt(annotation_rank) * w_2
            w_S_2 <- cbind(w_2, w_S_2)
            w_S <- cbind(w_S_1, w_S_2)
            w_S <- as.matrix(w_S)
            w_A_1 <- annotation_rank * w_1^2/dbeta(MAF, 0.5,
                0.5)^2
            w_A_1 <- cbind(w_1^2/dbeta(MAF, 0.5, 0.5)^2, w_A_1)
            w_A_2 <- annotation_rank * w_2^2/dbeta(MAF, 0.5,
                0.5)^2
            w_A_2 <- cbind(w_2^2/dbeta(MAF, 0.5, 0.5)^2, w_A_2)
            w_A <- cbind(w_A_1, w_A_2)
            w_A <- as.matrix(w_A)
        }

        if (obj_nullmodel$relatedness) {
            if (!obj_nullmodel$sparse_kins) {
                P <- obj_nullmodel$P
                P_scalar <- sqrt(dim(P)[1])
                P <- P * P_scalar
                residuals.phenotype <- obj_nullmodel$scaled.residuals
                residuals.phenotype <- residuals.phenotype *
                  sqrt(P_scalar)
                pvalues <- STAAR_O_SMMAT(G, P, residuals.phenotype,
                  weights_B = w_B, weights_S = w_S, weights_A = w_A,
                  mac = as.integer(round(MAF * 2 * dim(G)[1])))
            }else{
                Sigma_i <- obj_nullmodel$Sigma_i
                Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
                cov <- obj_nullmodel$cov
                residuals.phenotype <- obj_nullmodel$scaled.residuals
                pvalues <- Burden_score_SMMAT_sparse(G, Sigma_i, Sigma_iX,
                  cov, residuals.phenotype, weights_B = w_B,
                  weights_S = w_S, weights_A = w_A, mac = as.integer(round(MAF * 2 * dim(G)[1])))
                #pvalues <- STAAR_O_SMMAT_sparse(G, Sigma_i, Sigma_iX,
                #  cov, residuals.phenotype, weights_B = w_B,
                #  weights_S = w_S, weights_A = w_A, mac = as.integer(round(MAF *
                #    2 * dim(G)[1])))
            }
        }
       else {
            X <- model.matrix(obj_nullmodel)
          working <- obj_nullmodel$weights
            sigma <- sqrt(summary(obj_nullmodel)$dispersion)
            if (obj_nullmodel$family[1] == "binomial") {
                fam <- 1
            } else if (obj_nullmodel$family[1] == "gaussian") {
                fam <- 0
            }
            residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values
            pvalues <- STAAR_O(G, X, working, sigma, fam, residuals.phenotype,
                weights_B = w_B, weights_S = w_S, weights_A = w_A,
                mac = as.integer(round(MAF * 2 * dim(G)[1])))
        }
        num_variant <- sum(RV_label)
        cMAC <- sum(G)
        num_annotation <- dim(annotation_phred)[2] + 1

        mat<-matrix(pvalues,ncol=3)
        pvalues_STAAR_B_1_25<-mat[c(1:(num_annotation)),]
        pvalues_STAAR_B_1_1<-mat[c((num_annotation+1):(2*num_annotation)),]
        pvalues_STAAR_B_1_25<-data.frame(pvalues_STAAR_B_1_25)
        pvalues_STAAR_B_1_25$weight<-c("Burden_1_25",names(annotation_phred))
        pvalues_STAAR_B_1_1<-data.frame(pvalues_STAAR_B_1_1)
        pvalues_STAAR_B_1_1$weight<-c("Burden_1_1",names(annotation_phred))
        cols<-c("Score","SE","Pvalue","weight")
        names(pvalues_STAAR_B_1_25)<-names(pvalues_STAAR_B_1_1)<-cols

        res<-list(num_variant = num_variant, cMAC = cMAC, RV_label = RV_label,
            results_STAAR_B_1_25 = pvalues_STAAR_B_1_25, results_STAAR_B_1_1 = pvalues_STAAR_B_1_1)



        }else{
        res<-NULL
        }

        return(res)
        }
