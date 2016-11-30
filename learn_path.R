################################################################################################################
##      ____        __  __                           __
##     / __ \____ _/ /_/ /_ _      ______ ___  __   / /   ___  ____ __________
##    / /_/ / __ `/ __/ __ \ | /| / / __ `/ / / /  / /   / _ \/ __ `/ ___/ __ \
##   / ____/ /_/ / /_/ / / / |/ |/ / /_/ / /_/ /  / /___/  __/ /_/ / /  / / / /
##  /  /   /____,_/\__/_/ /_/|__/|__/\__,_/\__/  /__/_____/\___/\__,_/_/  /_/_/
## /__/
##             Author: David Venuto
###############################################################################################################

predict.pathways <-
  function(train_df,
           test_df,
           cat_1_end_train,
           cat_1_end_test,
           alph,
           cond1_name,
           cond2_name,
           qnorm) {

    ### LOAD PACKAGES #############################################################################################

    library(glmnet)
    library(pheatmap)

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

    ### VARIABLE DECLERATION ########################################################################################

    cat1_list_train <- list()
    cat2_list_train <- list()
    cat1_list_test <- list()
    cat2_list_test <- list()
    avg_change <- NULL
    genes_list <- genes_list
    means <- list()
    merged_cat1_train <- data.frame(train_df[, 1:cat_1_end_train])
    merged_cat2_train <-
      data.frame(train_df[, (cat_1_end_train + 1):ncol(train_df)])
    merged_cat1_test <- data.frame(train_df[, 1:cat_1_end_test])
    merged_cat2_test <-
      data.frame(train_df[, (cat_1_end_test + 1):ncol(train_df)])
    train_list <- list()
    pred_list <- list()

    ### BUILD PATHWAY TRAINING STRUCTURE #############################################################################

    pheatmap(cor(train_df))

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

    ### BUILD PATHWAY TESTING STRUCTURE ##############################################################################

    for (i in 1:nrow(Kegg)) {
      pathway_genes_cat1 <- data.frame(genes_list[i])
      colnames(pathway_genes_cat1) <- i
      row.names(pathway_genes_cat1) <- pathway_genes_cat1[, 1]
      pathway_genes_cat1 <-
        merge(pathway_genes_cat1, merged_cat1_test, by = "row.names")
      row.names(pathway_genes_cat1) <- pathway_genes_cat1$Row.names
      pathway_genes_cat1$Row.names <- NULL

      pathway_genes_cat2 <- data.frame(genes_list[i])
      colnames(pathway_genes_cat2) <- i
      row.names(pathway_genes_cat2) <- pathway_genes_cat2[, 1]
      pathway_genes_cat2 <-
        merge(pathway_genes_cat2, merged_cat2_test, by = "row.names")
      row.names(pathway_genes_cat2) <- pathway_genes_cat2$Row.names
      pathway_genes_cat2$Row.names <- NULL

      cat1_list_test[[i]] <- pathway_genes_cat1
      cat2_list_test[[i]] <- pathway_genes_cat2
    }

    ### LEARN AND PREDICT PATHWAYS ##############################################################################

    for (i in 1:nrow(Kegg)) {
      pathway_train_df <- data.frame(genes_list[i])
      colnames(pathway_train_df) <- i
      row.names(pathway_train_df) <- pathway_train_df[, 1]
      pathway_train_df <-
        merge(pathway_train_df, train_df, by = "row.names")
      row.names(pathway_train_df) <- pathway_train_df$Row.names
      pathway_train_df$Row.names <- NULL
      pathway_train_df[, 1] <- NULL

      pathway_test_df <- data.frame(genes_list[i])
      colnames(pathway_test_df) <- i
      row.names(pathway_test_df) <- pathway_test_df[, 1]
      pathway_test_df <-
        merge(pathway_test_df, test_df, by = "row.names")
      row.names(pathway_test_df) <- pathway_test_df$Row.names
      pathway_test_df$Row.names <- NULL
      pathway_test_df[, 1] <- NULL

      if (rowMeans(pathway_test_df) != 0 &&
          nrow(pathway_test_df) != 1 &&
          rowMeans(pathway_train_df) != 0 &&
          nrow(pathway_train_df) != 1) {
        y <-
          (c(
            rep(cond1_name, each = cat_1_end_train),
            rep(cond2_name, each = ncol(train_df) - cat_1_end_train)
          ))
        fit_RR <-
          glmnet(t(pathway_train_df),
                 y,
                 family = "binomial",
                 alpha = alph)
        cvfit = cv.glmnet(
          t(pathway_train_df),
          y,
          family = "binomial",
          type.measure = "class",
          alpha = alph
        )
        gene_mat <- coef(cvfit, s = "lambda.min")
        pred_list[[i]] <-
          predict(
            fit_RR,
            newx = t(pathway_test_df),
            type = "class",
            s = cvfit$lambda.min
          )
        train_list[[i]] <- pathway_train_df
      }
      else{

      }
    }

    ### FIND SUCCESS RATE #################################################################################

    suc = 0
    fail = 0
    for (i in 1:nrow(Kegg)) {
      tab = data.frame(pred_list[i])
      if (nrow(tab) != 0) {
        if (tab[1, 1] == cond1_name && tab[2, 1] == cond2_name) {
          suc = suc + 1
        }
        else{
          fail = fail + 1
        }
      }
    }

    suc <- suc / (suc + fail)

    ### BUILD OUTPUT ##############################################################################################

    PathwayLearnData = function(suc, genes_list, pred_list, train_list) {
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
      return(
        list(
          TrainData = TrainData,
          PredicitonData = PredicitonData,
          GetSuccessRate = GetSuccessRate,
          GetGeneList = GetGeneList
        )
      )
    }

    output <- PathwayLearnData(suc, genes_list, pred_list, train_list)
    return(output)

  }
