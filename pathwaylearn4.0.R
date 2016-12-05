################################################################################################################
##      ____        __  __                           __
##     / __ \____ _/ /_/ /_ _      ______ ___  __   / /   ___  ____ __________
##    / /_/ / __ `/ __/ __ \ | /| / / __ `/ / / /  / /   / _ \/ __ `/ ___/ __ \
##   / ____/ /_/ / /_/ / / / |/ |/ / /_/ / /_/ /  / /___/  __/ /_/ / /  / / / /
##  /  /   /____,_/\__/_/ /_/|__/|__/\__,_/\__/  /__/_____/\___/\_/_/  /_/ /_/ KEGG
## /__/
##             Author: David Venuto
###############################################################################################################

load("/home/yfwang/dvenuto/learned_Lungs.RData")

predict.pathways <-
  function(train_df,
           cat_1_end_train,
           alph = 1,
           cond1_name = "condition_1",
           cond2_name = "condition_2",
           qnorm = FALSE,
           voom_norm = FALSE,
           folds = 5) {
    
    ### LOAD PACKAGES #############################################################################################
   
    library(ROCR)
    library(glmnet)
    library(edgeR)
    
    ### QUANTILE NORMALIZATION SCRIPT #############################################################################
        
    quantile_normalisation <- function(df) {
      df_rank <- apply(df, 2, rank, ties.method = "min")
      df_sorted <- data.frame(apply(df, 2, sort))
      df_mean <- apply(df_sorted, 1, mean)
      
      index_to_mean <- function(my_index, my_mean) {
        return(my_mean[my_index])
      }
      df_final <-
        apply(df_rank, 2, index_to_mean, my_mean = df_mean)
      rownames(df_final) <- rownames(df)
      return(df_final)
    }
    
    ### PARAMETER DECLERATION ######################################################################################
    
    if (qnorm == TRUE) {
      train_df <- (quantile_normalisation(train_df))
    }
    
    if (voom_norm == TRUE) {
      train_df <- DGEList(counts = train_df)
      train_df <- calcNormFactors(train_df)
      train_df <- data.frame(cpm(train_df))
    }
    
    ### VARIABLE DECLERATION ########################################################################################
    
    cat1_list_train <- list()
    cat2_list_train <- list()
    up_down <- list()
    PermuteList <- list()
    ROC_permute_cor <- list()
    MIM_cond1 <- list()
    ROC_permute <- list()
    MIM_cond2 <- list()
    Kegg <- Kegg
    coe_list <- list()
    MIM_change <- list()
    avg_change <- NULL
    genes_list <- genes_list
    names(genes_list) <- Kegg$V2
    rowlist <- 1
    means <- list()
    merged_cat1_train <- data.frame(train_df[, 1:cat_1_end_train])
    merged_cat2_train <-
      data.frame(train_df[, (cat_1_end_train + 1):ncol(train_df)])
    train_list <- list()
    pred_list <- list()
    
    if(cat_1_end_train/folds < 10){
      print("AUC curve cannot be calculated with less than 10 category 1 observations per fold")
    }
    else if(((ncol(train_df)-(cat_1_end_train))/folds) < 10){
      print("AUC curve cannot be calculated with less than 10 category 2 observations per fold")
    }
    
    ### BUILD PATHWAY TRAINING STRUCTURE #############################################################################
    
    for (i in 1:nrow(Kegg)) {
      pathway_genes_cat1 <- data.frame(genes_list[i])
      colnames(pathway_genes_cat1) <- i
      row.names(pathway_genes_cat1) <- pathway_genes_cat1[, 1]
      pathway_genes_cat1 <-
        merge(pathway_genes_cat1, merged_cat1_train, by = "row.names")
      row.names(pathway_genes_cat1) <- pathway_genes_cat1$Row.names
      pathway_genes_cat1$Row.names <- NULL
      
      pathway_genes_cat2 <- data.frame(genes_list[i])
      colnames(pathway_genes_cat2) <- i
      row.names(pathway_genes_cat2) <- pathway_genes_cat2[, 1]
      pathway_genes_cat2 <-
        merge(pathway_genes_cat2, merged_cat2_train, by = "row.names")
      row.names(pathway_genes_cat2) <- pathway_genes_cat2$Row.names
      pathway_genes_cat2$Row.names <- NULL
      
      cat1_list_train[[i]] <- pathway_genes_cat1
      cat2_list_train[[i]] <- pathway_genes_cat2
    }
            
    ### LEARN AND PREDICT PATHWAYS ##############################################################################
    
    ROCList <- data.frame(AUCval=0)
    
    for (i in 1:nrow(Kegg)) {
      pathway_train_df <- data.frame(genes_list[i])
      colnames(pathway_train_df) <- i
      row.names(pathway_train_df) <- pathway_train_df[, 1]
      pathway_train_df <-
        merge(pathway_train_df, train_df, by = "row.names")
      row.names(pathway_train_df) <- pathway_train_df$Row.names
      pathway_train_df$Row.names <- NULL
      pathway_train_df[, 1] <- NULL 
      #Holdout Data
      HoldoutSize <- round((length(colnames(pathway_train_df))*0.2)/2, digits = 0)
      Cat1_Cols <- sample(1:(cat_1_end_train), HoldoutSize, replace=F)
      Cat2_Cols <- sample((cat_1_end_train+1):(ncol(pathway_train_df)),HoldoutSize, replace=F)
      pathway_test_df <- pathway_train_df[,c(Cat1_Cols,Cat2_Cols)]
      pathway_train_df <- pathway_train_df[,-c(Cat1_Cols,Cat2_Cols)]         
      if (rowMeans(pathway_test_df) != 0 &&
          nrow(pathway_test_df) > 1 &&
          rowMeans(pathway_train_df) != 0 &&
          nrow(pathway_train_df) > 1) {
        y <-
          (c(
            rep(cond1_name, each = cat_1_end_train-HoldoutSize),
            rep(cond2_name, each = ncol(pathway_train_df) - cat_1_end_train + HoldoutSize)
          ))
	try({
        fit_RR <-
          glmnet(t(pathway_train_df),
                 y,
                 family = "binomial",
                 alpha = alph)
        cvfit <- cv.glmnet(
          t(pathway_train_df),
          y,
          family = "binomial",
          type.measure = "class",
          alpha = alph,
          nfolds = folds
        )
	AUC <- cv.glmnet(
          t(pathway_train_df),
          y,
          type.measure = "auc",
          alpha = alph,
          family = "binomial",
          nfolds = folds
        )
        gene_mat <- coef(cvfit, s = "lambda.min")
        coe_list[[i]] <- (gene_mat)
        pathway_test_df <- data.matrix((pathway_test_df))
	pred <- 0 
        pred <-
          predict(
            AUC,
	    type="response",
            newx = t(pathway_test_df),
            s = AUC$lambda.min
          )
	labels <- c(rep(0,HoldoutSize),rep(1,HoldoutSize))
	AUC_pred <- prediction((pred),labels)
	pref <- performance(AUC_pred,"auc")
	pred <-
          predict(
            AUC,
            type="class",
            newx = t(pathway_test_df),
            s = AUC$lambda.min
          )	
	pred_list[[i]] <- pred
        train_list[[i]] <- pathway_train_df
        ROCObj <- data.matrix(unlist(pref@y.values))
        row.names(ROCObj) <- colnames(ROCObj)
        colnames(ROCObj) <- c("AUCval")
        ROCObj <- data.frame(ROCObj)
        ROCList <- rbind(ROCList, ROCObj)
        row.names(ROCList)[nrow(ROCList)] <-
          as.character(Kegg$V2[i])
        names(train_list)[i] <- as.character(Kegg$V2[i])
        names(pred_list)[i] <- as.character(Kegg$V2[i]) 
	names(coe_list)[i] <- as.character(Kegg$V2[i])
	})
      }	
      else{
        print("Not enough observations")
      }
    }

## FIND SUCCESS RATE #####################################################################################
    
    suc = 0
    fail = 0
    for (i in 1:nrow(Kegg)) {
      tab = data.frame(pred_list[i])
      if (nrow(tab) != 0) {
        for (n in 1:HoldoutSize) {
          if (tab[n, 1] == cond1_name) {
            suc = suc + 1
          }
          else{
            fail = fail + 1
          }
        }
        for (n in (HoldoutSize + 1):nrow(tab)) {
          if (tab[n, 1] == cond2_name) {
            suc = suc + 1
          }
          else{
            fail = fail + 1
          }
        }
      }
    }
    
    suc <- suc / (suc + fail)
    
    ### BUILD OUTPUT ###########################################################################################
    
   PathwayLearnData = function(suc,
                                genes_list,
                                pred_list,
                                train_list,
                                ROClist,
                                coef_list
                                ) {
      SuccessRate = suc
      GeneList = genes_list
      Predictions = pred_list
      TrainingData = train_list
      TrainData = function()
        TrainingData
      PredicitonData = function()
        Predictions
      GetSuccessRate = function()
        SuccessRate
      GetGeneList = function()
        GeneList
      ROCList = function()
        ROClist
      Coef_list = function()
        coef_list
      return(
        list(
          TrainData = TrainData,
          PredicitonData = PredicitonData,
          GetSuccessRate = GetSuccessRate,
          GetGeneList = GetGeneList,
          ROClist = ROClist,
          Coef_list = Coef_list
        )
      )
    }

    output <-
      PathwayLearnData(
        suc,
        genes_list,
        pred_list,
        train_list,
        ROCList,
        coe_list
      )
    return(output)

  }


train <- CDS_Colon_Cleaned

  learned_Colon <- predict.pathways(
    train_df = train,
    cat_1_end_train = 227,
    alph = 0,
    cond1_name = "case",
    cond2_name = "control",
    qnorm = FALSE,
    voom_norm = TRUE,
    folds = 5
  )
save.image("/home/yfwang/dvenuto/Leaned_AUC_Colon_RIDGE.RData")
