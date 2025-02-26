# MethylDetectR_run.R
# ✅ **檢查是否透過命令行傳遞參數**
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 4) {
    cat("✅ 檢測到 4 個命令行參數:\n")
    cat("  1. beta_matrix_path:", args[1], "\n")
    cat("  2. sample_sheet_path:", args[2], "\n")
    cat("  3. output_path:", args[3], "\n")
    cat("  4. project_root:", args[4], "\n")
    # ✅ 如果有傳遞 3 個參數，則使用命令行參數
    beta_matrix_path <- args[1]
    sample_sheet_path <- args[2]
    output_path <- args[3]
    project_root <- args[4]
    # ✅ gse從 `pd_file_path` 推測
    gse <- strsplit(basename(sample_sheet_path), "_")[[1]][1]

    MethylDetectR_dependencies <- file.path(project_root, "scripts", "MethylDetectR_dependencies")
    #測試用把MethylDetectR_dependencies直接寫死在這
    # MethylDetectR_dependencies <- "D:/SideProject/EpiAging/Python/LUCY_pipeline/scripts/MethylDetectR_dependencies/"
} else {
    stop("usage: Rscript MethylDetectR_run.R beta_matrix_path sample_sheet_path output_path")
}

cat("🚀 開始處理 GSE:", gse, "\n")
cat("📄 Beta Matrix:", beta_matrix_path, "\n")
cat("📄 Sample Sheet:", sample_sheet_path, "\n")
cat("📂 輸出檔案:", output_path, "\n")

run_MethylDetectR <- function(beta_matrix_path, sample_sheet_path, output_path) {
  # 1. Loading data
  message("1. Loading data")
  data <- read.csv(beta_matrix_path, header = TRUE, row.names = 1)
  sexageinfo <- read.csv(sample_sheet_path, header = TRUE)
  colnames(sexageinfo)[1] <- "ID"
  colnames(sexageinfo)[colnames(sexageinfo) == "Gender"] <- "Sex"
  sexageinfo <- sexageinfo[sexageinfo$ID %in% colnames(data), ]

  message("1.1 Loading Methylation data - rows to be CpGs and columns to be individuals")
  cpgs <- read.csv(file.path(MethylDetectR_dependencies, "Predictors_Shiny_by_Groups.csv"), header = T)
  
  coefficients <- read.delim(file.path(MethylDetectR_dependencies, "cage_coefficients_linear.tsv"), sep = "\t")
  coefficients_log <- read.delim(file.path(MethylDetectR_dependencies, "cage_coefficients_log.tsv"), sep = "\t")
  row.names(coefficients) <- coefficients[,"CpG_Site"]
  row.names(coefficients_log) <- coefficients_log[,"CpG_Site"]
  
  # Prepare CpGs for Bernabeu age predictor
  # Coefficients for log models
  coef_log_2 <- coefficients_log[rownames(coefficients_log)[grep('_2', rownames(coefficients_log))],,drop=FALSE]
  coef_log_2_simp <- gsub('_2', '', rownames(coef_log_2))
  coef_log <- coefficients_log[which(!(rownames(coefficients_log) %in% rownames(coef_log_2))),,drop=FALSE]
  
  # Coefficients for non-log models
  coef_2 <- coefficients[rownames(coefficients)[grep('_2', rownames(coefficients))],,drop=FALSE]
  coef_2_simp <- gsub('_2', '', rownames(coef_2))
  coef1 <- coefficients[rownames(coefficients)[which(!(rownames(coefficients) %in% rownames(coef_2)))],,drop=FALSE]
  
  # Beta means - for Bernabeu predictor
  means <- read.delim(file.path(MethylDetectR_dependencies, "cpg_meanbeta_gs20k.tsv"))
  row.names(means) <- means[,"cpg"]
  
  # Total CpGs
  cpgs_linear <- union(rownames(coef1), rownames(coef_log))
  cpgs_squared <- union(coef_2_simp, coef_log_2_simp)
  all_cpgs <- union(cpgs_linear, cpgs_squared)
  
  message("2. Quality Control and data Preparation")
  message("2.1 Checking if Row Names are CpG Sites")
  if(ncol(data) > nrow(data)){
    message("It seems that individuals are rows - data will be transposed!")
    data <- t(data)
  }
  
  message("2.2 Subsetting CpG sites to those required for Predictor Calculation")
  coef <- data[intersect(rownames(data), cpgs$CpG_Site),]
  coef.cage <- data[intersect(rownames(data), all_cpgs),]
  
  message("2.3 Checking if Beta Values are Present")
  m_to_beta <- function(val) {
    beta <- 2^val/(2^val + 1)
    return(beta)
  }
  
  coef <- if((range(coef,na.rm=T)> 1)[[2]] == "TRUE") {
    message("Suspect that M Values are present. Converting to Beta Values")
    m_to_beta(coef)
  } else {
    message("Suspect that Beta Values are present")
    coef
  }
  
  coef.cage <- if((range(coef.cage,na.rm=T)> 1)[[2]] == "TRUE") {
    message("Suspect that M Values are present. Converting to Beta Values")
    m_to_beta(coef.cage)
  } else {
    message("Suspect that Beta Values are present")
    coef.cage
  }
  
  message("2.4 Convert NA Values to Mean for each Probe")
  na_to_mean <- function(methyl){
    methyl[is.na(methyl)] <- mean(methyl,na.rm=T)
    return(methyl)
  }
  
  coef <- t(apply(coef,1,function(x) na_to_mean(x)))
  coef.cage <- t(apply(coef.cage,1,function(x) na_to_mean(x)))
  
  message("2.5 Find CpGs not present in uploaded file, add these with mean Beta Value for CpG site from Training Sample")
  coef <- if(nrow(coef) == length(unique(cpgs$CpG_Site))) {
    message("All sites present")
    coef
  } else if(nrow(coef)==0) {
    message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
  } else {
    missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef)),c("CpG_Site","Mean_Beta_Value")]
    message(paste(length(unique(missing_cpgs$CpG_Site)), "unique sites are missing - add to dataset with mean Beta Value from Training Sample", sep = " "))
    mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)),ncol = ncol(coef))
    row.names(mat) <- unique(missing_cpgs$CpG_Site)
    colnames(mat) <- colnames(coef)
    mat[is.na(mat)] <- 1
    missing_cpgs1 <- if(length(which(duplicated(missing_cpgs$CpG_Site))) > 1) {
      missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
    } else {
      missing_cpgs
    }
    ids = unique(row.names(mat))
    missing_cpgs1 = missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
    mat = mat*missing_cpgs1$Mean_Beta_Value
    coef = rbind(coef,mat)
  }
  
  coef.cage <- if(nrow(coef.cage)==length(all_cpgs)) {
    message("All sites present")
    coef.cage
  } else if(nrow(coef.cage)==0) {
    message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
  } else {
    missing_cpgs <- all_cpgs[which(!(all_cpgs %in% rownames(coef.cage)))]
    message(paste(length(missing_cpgs), "unique sites are missing - adding to dataset with mean Beta Value from GS training sample (N = 18,413)", sep = " "))
    mat = matrix(nrow=length(missing_cpgs), ncol = ncol(coef.cage))
    row.names(mat) <- missing_cpgs
    colnames(mat) <- colnames(coef.cage)
    mat[is.na(mat)] <- 1
    ids <- unique(row.names(mat))
    missing_cpg_means <- means[ids,]
    mat <- mat*missing_cpg_means[,"mean"]
    coef.cage <- rbind(coef.cage,mat)
  }
  
  message("3. Calculating the Predictors")
  loop <- unique(cpgs$Predictor)
  out <- data.frame()
  for(i in loop) {
    tmp <- coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% i,"CpG_Site"]),]
    tmp_coef = cpgs[cpgs$Predictor %in% i, ]
    if(nrow(tmp_coef) > 1) {
      tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
      out[colnames(coef),i] = colSums(tmp_coef$Coefficient*tmp)
    } else {
      tmp2 = as.matrix(tmp)*tmp_coef$Coefficient
      out[colnames(coef),i] = colSums(tmp2)
    }
  }
  
  out$'Epigenetic Age (Zhang)' <- out$'Epigenetic Age (Zhang)' + 65.79295
  out$ID <- row.names(out)
  out <- out[,c(ncol(out),1:(ncol(out)-1))]
  
  message("4. Obtaining cAge predictions - Bernabeu et al.")
  message("4.1. Preparing data for prediction")
  
  scores <- coef.cage[rownames(coef1),]
  scores_quadratic <- coef.cage[coef_2_simp,]**2
  rownames(scores_quadratic) <- paste0(rownames(scores_quadratic), "_2")
  scores_linear <- rbind(scores, scores_quadratic)
  
  scores_log <- coef.cage[rownames(coef_log),]
  scores_quadratic_log <- coef.cage[coef_log_2_simp,]**2
  rownames(scores_quadratic_log) <- paste0(rownames(scores_quadratic_log), "_2")
  scores_log <- rbind(scores_log, scores_quadratic_log)
  
  coefficients <- coefficients[rownames(scores_linear),]
  pred_linear <- scores_linear * coefficients[,"Coefficient"]
  pred_linear_pp <- colSums(pred_linear,na.rm=T)
  pred_linear_pp <- pred_linear_pp + 30.7873138968084
  
  over20s <- names(pred_linear_pp[pred_linear_pp > 20])
  pred_linear_pp1 <- pred_linear_pp[over20s]
  under20s <- names(pred_linear_pp[pred_linear_pp < 20])
  
  coefficients_log <- coefficients_log[rownames(scores_log),]
  pred_log <- scores_log * coefficients_log[,"Coefficient"]
  pred_log_pp <- colSums(pred_log,na.rm=T)
  pred_log_pp <- pred_log_pp + 3.87249881475691
  pred_log_pp <- exp(pred_log_pp[under20s])
  predictions <- c(pred_log_pp, pred_linear_pp1)
  results <- data.frame("ID" = names(predictions), "Epigenetic Age (Bernabeu)" = predictions,check.names=F)
  
  out <- merge(out,results,by="ID")
  out <- out[,c(1,2,ncol(out),3:(ncol(out)-1))]
  
  if(!exists("sexageinfo")){
    out$'True Age' <- NA
    out$Sex <- NA
  } else {
    ids = out$ID
    sexageinfo = sexageinfo[match(ids, sexageinfo$ID),]
    out$'True Age' <- sexageinfo$Age
    out$Sex <- sexageinfo$Sex
    message("4. Sex and Age Info Added")
  }
  
  out <- out[,c(1,ncol(out),(ncol(out)-1), 2:c(ncol(out)-2))]
  
  if(!exists("covariates")){
    NULL
  } else {
    ids = out$ID
    covariates = covariates[match(ids, covariates$ID),]
    out <- merge(out, covariates, by = "ID")
    for(i in names(out)[5:10]){
      out[,i] <- resid(lm(as.formula(paste(paste(bt(i),"~"),paste(bt(paste(c(names(covariates)[-which(names(covariates) %in% "ID")]))), collapse = "+"))), na.action = na.exclude, data = out))
    }
    out <- out[,c(1:10)]
    message("5. Covariates")
  }
  
  write.csv(out, output_path, row.names = FALSE)
  message("Analysis Finished! Thank you for using our application. Output File is called \"out\"")
  
  return(0)
}


run_MethylDetectR(beta_matrix_path, sample_sheet_path, output_path)