Internal.boxplot.stats <- function(x, # borrowed from grDevices::boxplot.stats()
                                   coef = 1.5, do.conf = TRUE, do.out = TRUE) {
    if (coef < 0)
        stop("'coef' must not be negative")
    nna <- !is.na(x)
    n <- sum(nna)
    stats <- stats::fivenum(x, na.rm = TRUE)
    iqr <- diff(stats[c(2, 4)])
    names(iqr) <- NULL
    if (coef == 0)
        do.out <- FALSE
    else {
        lower_limit <- (stats[2L] - coef * iqr)
        upper_limit <- (stats[4L] + coef * iqr)
        names(lower_limit) <- NULL
        names(upper_limit) <- NULL
        out <- if (!is.na(iqr)) {
            x < lower_limit | x > upper_limit
        }
        else !is.finite(x)
        if (any(out[nna], na.rm = TRUE))
            stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)
    }
    conf <- if (do.conf)
        stats[3L] + c(-1.58, 1.58) * iqr/sqrt(n)
    list(stats = stats, n = n, conf = conf,
         out = if (do.out) x[out & nna] else numeric(),
         iqr = iqr, lower_limit = lower_limit, upper_limit = upper_limit)
}

Internal.remove_NA <- function(input_data_num, data_label = NULL,
                               replace_with_zero = FALSE) {
    data_na_sample_idx <- apply(input_data_num, 1, function(x) any(is.na(x)))
    data_na_sample_sum <- sum(data_na_sample_idx)

    if (data_na_sample_sum > 0) {
        if (!is.null(data_label)) warning(paste0("    ", data_na_sample_sum, " sample(s) in ", data_label, " contain(s) NA."))

        if (replace_with_zero) {
            input_data_num[is.na(input_data_num)] <- 0
        } else {
            input_data_num <- input_data_num[!data_na_sample_idx,]
            if (nrow(input_data_num) == 0) stop(paste0("    Variable selection failed: ", data_label, " contain too many NA!"))
        }
    }

    input_data_num
}

# Internal.impute_infinite <- function(input_data) {
#     data_check <- sapply(input_data, function(one_col) {
#         idx_inf_positive <- one_col == Inf
#         idx_inf_negative <- one_col == -Inf
#
#         if (any(idx_inf_positive, na.rm = TRUE)) one_col[idx_inf_positive] <- max(one_col[is.finite(one_col)], na.rm = TRUE) * 4
#         if (any(idx_inf_negative, na.rm = TRUE)) one_col[idx_inf_negative] <- min(one_col[is.finite(one_col)], na.rm = TRUE) / 4
#         one_col
#     })
#     data.frame(data_check)
# }

Internal.compute_cor <- function(train_num, test_num = NULL,
                                 selectVar_corType   = c("pcor", "cor"),
                                 selectVar_corMethod = c("spearman", "pearson")) {

    if (selectVar_corType == "pcor") {
        train_num_noNA <- Internal.remove_NA(train_num, data_label = "training data")
        if (!is.null(test_num)) test_num_noNA <- Internal.remove_NA(test_num, data_label = "test data")

        if (ncol(train_num_noNA) > 500) message("      Your data have more than 500 variables. It may take some time to process large datasets.")

        train_cor <- data.frame(ppcor::pcor(train_num_noNA, method = selectVar_corMethod)$estimate)
        names(train_cor) <- names(train_num_noNA)
        row.names(train_cor) <- names(train_num_noNA)

        if (!is.null(test_num)) {
            test_cor <- data.frame(ppcor::pcor(test_num_noNA, method = selectVar_corMethod)$estimate)
            names(test_cor) <- names(test_num_noNA)
            row.names(test_cor) <- names(test_num_noNA)
        } else test_cor <- NULL

    } else {

        train_cor <- data.frame(cor(train_num, method = selectVar_corMethod, use = "complete.obs"))

        if (!is.null(test_num)) {
            test_cor <- data.frame(cor(test_num, method = selectVar_corMethod, use = "complete.obs"))
        } else test_cor <- NULL
    }

    cor_info <- list(variable_name = names(train_cor),
                     train_cor = train_cor, test_cor = test_cor)
}

Internal.select_variable <- function(cor_info, selectVar_minNum = NULL,
                                     selectVar_maxNum = NULL) {

    variable_name <- cor_info$variable_name

    train_cor <- cor_info$train_cor
    test_cor  <- cor_info$test_cor

    train_cor[is.na(train_cor)] <- 0
    if (!is.null(test_cor)) test_cor[is.na(test_cor)] <- 0

    pb <- pbapply::timerProgressBar(min = 0, max = length(variable_name),
                                    initial = 0, style = 3, width = 70,
                                    min_time = 20)
    selected_var <- lapply(1:length(variable_name), function(var_idx,
                                                             variable_name,
                                                             train_cor, test_cor,
                                                             selectVar_minNum, selectVar_maxNum) {
        pbapply::setTimerProgressBar(pb, var_idx)

        input_one_variable_name <- variable_name[[var_idx]]
        correlated_train <- abs(train_cor[input_one_variable_name])
        correlated_train <- correlated_train[!row.names(correlated_train) %in% input_one_variable_name,,drop = FALSE]
        correlated_train <- correlated_train[order(correlated_train[[1]], decreasing = TRUE),,drop = FALSE]

        correlated_train_name <- row.names(correlated_train)
        candidate_train_name  <- correlated_train_name[correlated_train[[1]] > 0.5]
        candidate_var_name    <- candidate_train_name

        if (!is.null(test_cor)) {
            correlated_test <- abs(test_cor[input_one_variable_name])
            correlated_test <- correlated_test[!row.names(correlated_test) %in% input_one_variable_name,,drop = FALSE]
            correlated_test <- correlated_test[order(correlated_test[[1]], decreasing = TRUE),,drop = FALSE]

            correlated_test_name <- row.names(correlated_test)
            candidate_test_name  <- correlated_test_name[correlated_test[[1]] > 0.5]
            candidate_var_name   <- intersect(candidate_var_name, candidate_test_name)
        }

        if (length(candidate_var_name) < selectVar_minNum) {

            if (is.null(test_cor)) {
                selected_var_name <- correlated_train_name[1:selectVar_minNum]
            } else {
                current_upper_limit <- min(length(candidate_train_name), length(candidate_test_name))
                current_upper_limit <- max(current_upper_limit, selectVar_minNum)

                while (1) {
                    candidate_test_name_tmp  <- correlated_test_name[1:current_upper_limit]
                    candidate_train_name_tmp <- correlated_train_name[1:current_upper_limit]
                    candidate_var_name_tmp <- intersect(candidate_train_name_tmp, candidate_test_name_tmp)
                    if (length(candidate_var_name_tmp) < selectVar_minNum) {
                        current_upper_limit <- current_upper_limit + 1
                    } else {
                        selected_var_name <- candidate_var_name_tmp
                        break
                    }
                }
            }

        } else if (length(candidate_var_name) > selectVar_maxNum) {

            if (is.null(test_cor)) {
                selected_var_name <- correlated_train_name[1:selectVar_maxNum]
            } else {
                upper_limit <- selectVar_maxNum

                while (1) {
                    candidate_test_name_tmp  <- correlated_test_name[1:upper_limit]
                    candidate_train_name_tmp <- correlated_train_name[1:upper_limit]
                    candidate_var_name_tmp <- intersect(candidate_train_name_tmp, candidate_test_name_tmp)

                    if (length(candidate_var_name_tmp) < selectVar_maxNum) {
                        upper_limit <- upper_limit + 1
                        selected_var_name <- candidate_var_name_tmp
                    } else if (length(candidate_var_name_tmp) == selectVar_maxNum) {
                        selected_var_name <- candidate_var_name_tmp
                        break
                    } else {
                        break
                    }
                }
            }

        } else {
            selected_var_name <- candidate_var_name
        }
        selected_var_name
    },
    variable_name = variable_name,
    train_cor = train_cor, test_cor = test_cor,
    selectVar_minNum = selectVar_minNum, selectVar_maxNum = selectVar_maxNum)
    pbapply::closepb(pb)
    names(selected_var) <- variable_name
    selected_var
}

Internal.compute_targetVal <- function(input_col, targetVal_method = c("median", "mean"),
                                       targetVal_removeOutlier = TRUE) {

    input_col <- input_col[is.finite(input_col)]
    if (targetVal_removeOutlier) outlier_res <- Internal.boxplot.stats(input_col, coef = 1.5)$out

    if (targetVal_method == "mean") {
        if (targetVal_removeOutlier) input_col[input_col %in% outlier_res] <- mean(input_col[!input_col %in% outlier_res], na.rm = TRUE)
        res <- mean(input_col, na.rm = TRUE)
    } else {
        if (targetVal_removeOutlier) input_col[input_col %in% outlier_res] <- median(input_col[!input_col %in% outlier_res], na.rm = TRUE)
        res <- median(input_col, na.rm = TRUE)
    }
    res
}

Internal.compute_errorRatio <- function(rawVal, sampleType,
                                        targetVal) {

    out <- sapply(1:length(rawVal), function(val_idx, rawVal, sampleType, targetVal) {

        current_targetVal <- targetVal[row.names(targetVal) == sampleType[val_idx],] [[1]]
        current_rawVal <- rawVal[val_idx]
        errorRatio <- ifelse(current_targetVal == 0, 0, (current_rawVal - current_targetVal) / current_targetVal)
        out <- c(y = errorRatio, y_target = current_targetVal, y_raw = current_rawVal)
        out
    }, rawVal = rawVal, sampleType = sampleType, targetVal = targetVal)

    out_df <- data.frame(t(out))

    if (all(out_df$y == 0)) {
        out_df$y <- rnorm(length(out_df$y), sd = 0.001, mean = 0.0005)
    }
    out_df
}

Internal.run_ensemble <- function(trainSet, testSet,
                                  mtry_percent = seq(0.2, 0.8, 0.2),
                                  nodesize_percent = seq(0.2, 0.8, 0.2),
                                  ..., return_base_res = FALSE) {

    if (!is.null(mtry_percent)) mtry <- round(mtry_percent * (ncol(trainSet) - 3))
    if (!is.null(nodesize_percent)) nodesize <- round(nodesize_percent * nrow(trainSet))

    rf_hyperparams <- expand.grid(mtry = unique(mtry),
                                  nodesize = unique(nodesize),
                                  ... = ...)

    pred_ensemble <- lapply(1:nrow(rf_hyperparams), function(idx) {
        current_hyperparams <- as.list(rf_hyperparams[idx,])

        folds_train <- caret::createFolds(1:length(trainSet$y), k = 5, returnTrain = TRUE)

        res_folds <- lapply(folds_train, function(train_idx, rf_params) {

            train_fold     <- trainSet[train_idx,]
            validate_fold  <- trainSet[-train_idx,]

            fold_formula <- c(formula = as.formula(y ~ .),
                              data = list(train_fold[!names(train_fold) %in% c("y_target", "y_raw")]),
                              current_hyperparams)
            RF_fold_mod  <- do.call(randomForest::randomForest, fold_formula)

            pred_fold <- predict(RF_fold_mod, validate_fold)

            pred_convert <- validate_fold$y_raw / (pred_fold + 1)
            pred_loss <- abs(pred_convert - validate_fold$y_target) / validate_fold$y_target
            pred_loss
        })

        mean_loss <- mean(unlist(res_folds), na.rm = TRUE)
        mod_weight <- 1/exp(mean_loss)

        base_formula <- c(formula = as.formula(y ~ .), data = list(trainSet[!names(trainSet) %in% c("y_target", "y_raw")]),
                          current_hyperparams)
        RF_base_mod <- do.call(randomForest::randomForest, base_formula)

        pred_test   <- predict(RF_base_mod, testSet)
        pred_test_convert <- testSet$y_raw / (pred_test + 1)

        if (return_base_res) {
            out <- list(mod_weight = mod_weight, pred_test_convert = pred_test_convert, RF_base_mod = RF_base_mod)
        } else {
            out <- list(mod_weight = mod_weight, pred_test_convert = pred_test_convert)
        }
        out
    })
    mod_weights <- sapply(pred_ensemble, function(x) x$mod_weight)
    mod_weights_norm <- mod_weights / sum(mod_weights, na.rm = TRUE)

    pred_test <- sapply(pred_ensemble, function(x) x$pred_test_convert)
    pred_norm <- apply(pred_test, 1, function(x) sum(x * mod_weights_norm, na.rm = TRUE))

    if (return_base_res) {
        output <- list(pred_norm = pred_norm, base_res = pred_ensemble, base_weights = mod_weights_norm, rf_hyperparams = rf_hyperparams)
    } else {
        output <- pred_norm
    }
    output
}

# save(FF4_qc, file = "FF4_qc.RData", compress = "xz", version = 2)
# devtools::document(roclets = c('rd', 'collate', 'namespace'))
# devtools::build_manual()
# usethis::use_cran_comments()
