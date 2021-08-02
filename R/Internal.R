compute_RSD <- function(input_data) {
    val_RSD <- sd(input_data, na.rm = T) / mean(input_data, na.rm = T)
    val_RSD
}

Internal.boxplot.stats <- function(x, coef = 1.5, do.conf = TRUE, do.out = TRUE) {
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

Internal.remove_NA <- function(input_data_num, data_label = NULL) {
    data_na_sample_idx <- apply(input_data_num, 1, function(x) any(is.na(x)))
    data_na_sample_sum <- sum(data_na_sample_idx)
    if (data_na_sample_sum == 1) {
        warning(paste0("  - One sample in ", data_label, " contains NA and has been removed."))
        input_data_num <- input_data_num[!data_na_sample_idx,]
    } else if (data_na_sample_sum > 1) {
        warning(paste0("  - ", data_na_sample_sum, " samples in ", data_label, " contain NA and have been removed."))
        input_data_num <- input_data_num[!data_na_sample_idx,]
    }
    if (nrow(input_data_num) == 0) stop(paste0("  - Variable selection failed: ", data_label, " contain too many NA!"))
    input_data_num
}

Internal.compute_cor <- function(train_num, test_num = NULL,
                                 correlation_type = c("pcor", "cor"),
                                 correlation_method = c("spearman", "pearson")) {

    message("  - Checking missing values...")
    train_num_noNA <- Internal.remove_NA(train_num, data_label = "training data")
    if (!is.null(test_num)) test_num_noNA <- Internal.remove_NA(test_num, data_label = "test data")

    message("  - Computing correlation coefficients...")
    if (ncol(train_num_noNA) > 500) message("  - Your data have more than 500 variables. It may take some time to process large datasets.")

    if (correlation_type == "pcor") {
        train_cor <- data.frame(ppcor::pcor(train_num_noNA, method = correlation_method)$estimate)
        names(train_cor) <- names(train_num_noNA)
        row.names(train_cor) <- names(train_num_noNA)
        if (!is.null(test_num)) {
            test_cor <- data.frame(ppcor::pcor(test_num_noNA, method = correlation_method)$estimate)
            names(test_cor) <- names(test_num_noNA)
            row.names(test_cor) <- names(test_num_noNA)
        } else test_cor <- NULL

    } else {
        train_cor <- data.frame(cor(train_num_noNA, method = correlation_method, use = "complete.obs"))

        if (!is.null(test_num)) {
            test_cor <- data.frame(cor(test_num_noNA, method = correlation_method, use = "complete.obs"))
        } else test_cor <- NULL
    }

    cor_info <- list(variable_name = names(train_cor),
                     train_cor = train_cor, test_cor = test_cor)
}

Internal.select_variable <- function(cor_info, min_var_num = NULL,
                                     max_var_num = NULL) {

    variable_name <- cor_info$variable_name
    pb <- pbapply::timerProgressBar(min = 0, max = length(variable_name),
                                    initial = 0, style = 3, width = 70,
                                    min_time = 30)
    selected_var <- lapply(1:length(variable_name), function(var_idx,
                                                             variable_name,
                                                             train_cor, test_cor,
                                                             min_var_num, max_var_num) {
        pbapply::setTimerProgressBar(pb, var_idx)
        input_one_variable_name <- variable_name[[var_idx]]
        correlated_train <- abs(train_cor[input_one_variable_name])
        correlated_train <- correlated_train[order(correlated_train[[1]], decreasing = TRUE),,drop = FALSE]

        correlated_train_name <- row.names(correlated_train)
        candidate_train_name  <- correlated_train_name[correlated_train[[1]] > 0.5]
        candidate_var_name    <- candidate_train_name

        if (!is.null(test_cor)) {
            correlated_test <- abs(test_cor[input_one_variable_name])
            correlated_test <- correlated_test[order(correlated_test[[1]], decreasing = TRUE),,drop = FALSE]

            correlated_test_name <- row.names(correlated_test)
            candidate_test_name  <- correlated_test_name[correlated_test[[1]] > 0.5]
            candidate_var_name   <- intersect(candidate_var_name, candidate_test_name)
        }

        if (length(candidate_var_name) < min_var_num) {

            if (is.null(test_cor)) {
                selected_var_name <- correlated_train_name[1:min_var_num]
            } else {
                current_upper_limit <- min(length(candidate_train_name), length(candidate_test_name))
                current_upper_limit <- max(current_upper_limit, min_var_num)

                while (1) {
                    candidate_test_name_tmp  <- correlated_test_name[1:current_upper_limit]
                    candidate_train_name_tmp <- correlated_train_name[1:current_upper_limit]
                    candidate_var_name_tmp <- intersect(candidate_train_name_tmp, candidate_test_name_tmp)
                    if (length(candidate_var_name_tmp) < min_var_num) {
                        current_upper_limit <- current_upper_limit + 1
                        # if (current_upper_limit > length(correlated_test_name)) {
                        #     message(input_one_variable_name, " - min candidate: limit greater than var length.")
                        #     break
                        # }
                    } else {
                        selected_var_name <- candidate_var_name_tmp
                        break
                    }
                }

            }

        } else if (length(candidate_var_name) > max_var_num) {

            if (is.null(test_cor)) {
                selected_var_name <- correlated_train_name[1:max_var_num]
            } else {
                upper_limit <- max_var_num

                while (1) {
                    candidate_test_name_tmp  <- correlated_test_name[1:upper_limit]
                    candidate_train_name_tmp <- correlated_train_name[1:upper_limit]
                    candidate_var_name_tmp <- intersect(candidate_train_name_tmp, candidate_test_name_tmp)

                    if (length(candidate_var_name_tmp) < max_var_num) {
                        upper_limit <- upper_limit + 1
                        # if (upper_limit > length(correlated_test_name)) {
                        #     message(input_one_variable_name, " - max candidate: limit greater than var length.")
                        #     break
                        # }
                    } else {
                        selected_var_name <- candidate_var_name_tmp
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
    train_cor = cor_info$train_cor, test_cor = cor_info$test_cor,
    min_var_num = min_var_num, max_var_num = max_var_num)
    pbapply::closepb(pb)
    names(selected_var) <- variable_name
    selected_var
}

select_variable <- function(train_num, test_num = NULL,
                            correlation_type = c("pcor", "cor"),
                            correlation_method = c("spearman", "pearson"),
                            min_var_num = NULL, max_var_num = NULL,
                            coerce_numeric = FALSE) {

    correlation_type   <- match.arg(correlation_type)
    correlation_method <- match.arg(correlation_method)

    # test_num <- test_num[!sapply(test_num, anyNA)]
    # test_num <- test_num[sapply(test_num, function(x) all(is.finite(x)))]
    # if (nrow(test_num) < 3) stop("  - Data of test samples have too many NA or infinite values. Maybe you would like to impute your data!")
    #
    # train_num <- train_num[!sapply(train_num, anyNA)]
    # train_num <- train_num[sapply(train_num, function(x) all(is.finite(x)))]
    # if (nrow(train_num) < 3) stop("  - Data of train samples have too many NA or infinite values. Maybe you would like to impute your data!")

    min_var_num <- ifelse(is.null(min_var_num), 1, min_var_num)
    max_var_num <- ifelse(is.null(max_var_num), ncol(train_num), max_var_num)

    min_var_num <- as.integer(min_var_num)
    max_var_num <- as.integer(max_var_num)

    if (min_var_num < 1) stop("  - min_var_num must be a positive integer!")
    if (max_var_num > ncol(train_num)) stop("  - max_var_num cannot be greater than variable number!")

    if(coerce_numeric) {
        train_num <- as.data.frame(sapply(train_num, as.numeric))
        idx_NA <- sapply(train_num, function(x) {
            all(is.na(x))
        })
        train_num <- train_num[,!idx_NA]

        if (!is.null(test_num)) {
            test_num <- as.data.frame(sapply(test_num, as.numeric))
            idx_NA <- sapply(test_num, function(x) {
                all(is.na(x))
            })
            test_num <- test_num[,!idx_NA]
        }
    } else {
        if (!all(sapply(train_num, is.numeric))) stop("  - The values of train samples should be numeric!")
        if (!is.null(test_num) & !all(sapply(test_num, is.numeric))) stop("  - The values of test samples should be numeric!")
    }

    if (!is.null(test_num)) {
        if (any(names(train_num) != names(test_num))) stop("  - Variables in training and test data cannot match!")
    }

    # stop_cl <- FALSE
    # if (is.null(cl)) {
    #     cl <- parallel::makeCluster(parallel.cores)
    #     stop_cl <- TRUE
    # }

    cor_info <- Internal.compute_cor(train_num = train_num, test_num = test_num,
                                     correlation_type = correlation_type,
                                     correlation_method = correlation_method)
    message("  - Selecting variables...")
    selected_var <- Internal.select_variable(cor_info = cor_info,
                                             min_var_num = min_var_num,
                                             max_var_num = max_var_num)

    # if (stop_cl) {
    #     parallel::stopCluster(cl)
    # }
    selected_var
}

Internal.compute_targetVal <- function(input_col, targetVal_method = c("median", "mean"),
                                       targetVal_removeOutlier = TRUE) {
    if (targetVal_removeOutlier) outlier_res <- Internal.boxplot.stats(input_col, coef = 1.5)$out

    if (targetVal_method == "mean") {
        if (targetVal_removeOutlier) input_col[input_col %in% outlier_res] <- mean(input_col[!input_col %in% outlier_res], na.rm = T)
        res <- mean(input_col, na.rm = T)
    } else {
        if (targetVal_removeOutlier) input_col[input_col %in% outlier_res] <- median(input_col[!input_col %in% outlier_res], na.rm = T)
        res <- median(input_col, na.rm = T)
    }
    res
}

compute_targetVal <- function(QC_data, col_sampleType, col_batchID,
                              targetVal_method = c("median", "mean"),
                              targetVal_batch = FALSE,
                              targetVal_removeOutlier = TRUE,
                              coerce_numeric = FALSE) {

    targetVal_method   <- match.arg(targetVal_method)

    sampleID <- as.factor(QC_data[[col_sampleType]])
    batchID  <- as.factor(QC_data[[col_batchID]])
    data_numeric <- QC_data[!names(QC_data) %in% c(col_sampleType, col_batchID)]

    if(coerce_numeric) {
        data_numeric <- as.data.frame(sapply(data_numeric, as.numeric))
        idx_NA <- sapply(data_numeric, function(x) {
            all(is.na(x))
        })
        data_numeric <- data_numeric[,!idx_NA]
    } else {
        if (!all(sapply(data_numeric, is.numeric))) stop("The values of the input dataset (except sampleType and batchID) should be numeric!")
    }

    if (targetVal_batch) {
        target_values <- aggregate(data_numeric, by = list(batch = batchID, sample = sampleID),
                                   FUN = Internal.compute_targetVal, targetVal_method = targetVal_method,
                                   targetVal_removeOutlier = targetVal_removeOutlier)
        batchID <- target_values$batch
        target_values$batch <- NULL
        target_values_list <- split(target_values, f = batchID)
    } else {
        target_values_list <- list(wholeDataset = aggregate(data_numeric,
                                                            by = list(sample = sampleID),
                                                            FUN =  Internal.compute_targetVal,
                                                            targetVal_method = targetVal_method,
                                                            targetVal_removeOutlier = targetVal_removeOutlier))
    }

    target_values <- lapply(target_values_list, function(x) {
        row.names(x) <- x$sample
        x$sample <- NULL
        x
    })
    target_values
}

Internal.compute_errorRatio <- function(train_samples, col_sampleType,
                                        targetVal_df, current_var) {
    out <- sapply(1:nrow(train_samples), function(row_idx, train_samples,
                                                  col_sampleType, targetVal_df,
                                                  current_var) {

        targetVal <- targetVal_df[row.names(targetVal_df) == train_samples[[col_sampleType]] [row_idx],] [[current_var]]
        rawVal <- train_samples[[current_var]][row_idx]
        errorRatio <- ifelse(targetVal == 0, 0, (rawVal - targetVal) / targetVal)
        out <- c(errorRatio = errorRatio, targetVal = targetVal, rawVal = rawVal)

    }, train_samples = train_samples, col_sampleType = col_sampleType,
    targetVal_df = targetVal_df, current_var = current_var)

    out_df <- data.frame(t(out))
}

Internal.run_ensemble <- function(trainSet, testSet,
                                  mtry_ratio = seq(0.2, 0.8, 0.2),
                                  nodesize_ratio = seq(0.2, 0.8, 0.2),
                                  ..., return_base_res = FALSE) {

    if (!is.null(mtry_ratio)) mtry <- round(mtry_ratio * (ncol(trainSet) - 3))
    if (!is.null(nodesize_ratio)) nodesize <- round(nodesize_ratio * nrow(trainSet))

    rf_hyperparams <- expand.grid(mtry = unique(mtry),
                                  nodesize = unique(nodesize),
                                  ... = ...)

    pred_ensemble <- lapply(1:nrow(rf_hyperparams), function(idx) {
        current_hyperparams <- as.list(rf_hyperparams[idx,])

        folds_train <- caret::createFolds(1:length(trainSet$y), k = 5, returnTrain = T)

        res_folds <- lapply(folds_train, function(train_idx, rf_params) {

            train_fold     <- trainSet[train_idx,]
            validate_fold  <- trainSet[-train_idx,]
            # train_fold[c("y_target", "y_raw")] <- NULL

            fold_formula <- c(formula = as.formula(y ~ .),
                              data = list(train_fold[!names(train_fold) %in% c("y_target", "y_raw")]),
                              current_hyperparams)
            RF_fold_mod  <- do.call(randomForest::randomForest, fold_formula)

            pred_fold <- predict(RF_fold_mod, validate_fold)

            pred_convert <- validate_fold$y_raw / (pred_fold + 1)
            pred_loss <- abs(pred_convert - validate_fold$y_target) / validate_fold$y_target
            pred_loss
        })

        mean_loss <- mean(unlist(res_folds), na.rm = T)
        mod_weight <- 1/exp(mean_loss)

        base_formula <- c(formula = as.formula(y ~ .), data = list(trainSet[!names(trainSet) %in% c("y_target", "y_raw")]),
                          current_hyperparams)
        RF_base_mod <- do.call(randomForest::randomForest, base_formula)

        pred_test   <- predict(RF_base_mod, testSet)
        pred_test_convert <- testSet$y_raw / (pred_test + 1)

        out <- list(mod_weight = mod_weight, pred_test_convert = pred_test_convert)
        # if (return_base_res) out <- c(out, RF_base_mod = list(RF_base_mod))
        # if (return_base_res) {
        #     out <- list(mod_weight = mod_weight, pred_test_convert = pred_test_convert, RF_base_mod = RF_base_mod)
        #
        # } else {
        #     out <- list(mod_weight = mod_weight, pred_test_convert = pred_test_convert)
        # }
        out
    })
    mod_weights <- sapply(pred_ensemble, function(x) x$mod_weight)
    mod_weights_norm <- mod_weights / sum(mod_weights, na.rm = TRUE)

    pred_test <- sapply(pred_ensemble, function(x) x$pred_test_convert)
    pred_norm <- apply(pred_test, 1, function(x) sum(x * mod_weights_norm, na.rm = TRUE))

    if (return_base_res) {
        output <- list(pred_norm = pred_norm, base_res = pred_ensemble, base_weights = mod_weights_norm)
    } else {
        output <- pred_norm
    }
    output
}

run_TIGER <- function(test_samples, train_samples,
                      col_sampleID, col_sampleType, col_batchID,
                      col_order = NULL, col_position = NULL,
                      targetVal_external = NULL, targetVal_method = c("mean", "median"),
                      targetVal_batch = FALSE, targetVal_removeOutlier = TRUE,
                      correlation_type = c("cor", "pcor"),
                      correlation_method = c("pearson", "spearman"),
                      min_var_num = 10, max_var_num = 30,
                      mtry_ratio = seq(0.2, 0.8, 0.2),
                      nodesize_ratio = seq(0.2, 0.8, 0.2),
                      ..., parallel.cores = 2) {

    message("+ Initialising...   ", Sys.time())

    targetVal_method   <- match.arg(targetVal_method)
    correlation_type   <- match.arg(correlation_type)
    correlation_method <- match.arg(correlation_method)

    # zero values?
    # Inf values?
    # NA values?
    # order the output

    for (col_idx in c(col_sampleID, col_sampleType, col_batchID)) {
        test_samples[[col_idx]] <- as.character(test_samples[[col_idx]])
        train_samples[[col_idx]] <- as.character(train_samples[[col_idx]])
    }

    if (!is.null(col_order)) {
        if (anyNA(test_samples[[col_order]])  | any(!is.finite(test_samples[[col_order]]))  ) stop("  - test samples: col_order should be numeric only!")
        if (anyNA(train_samples[[col_order]]) | any(!is.finite(train_samples[[col_order]])) ) stop("  - train samples: col_order should be numeric only!")
    }

    if (!is.null(col_position)) {
        if (anyNA(test_samples[[col_position]])  | any(!is.finite(test_samples[[col_position]]))  ) stop("  - test samples: col_position should be numeric only!")
        if (anyNA(train_samples[[col_position]]) | any(!is.finite(train_samples[[col_position]])) ) stop("  - train samples: col_position should be numeric only!")
    }

    if (!all(sapply(train_samples[!names(train_samples) %in% c(col_sampleID, col_sampleType, col_batchID)], is.numeric))) {
        stop("  - The values of train samples (except sampleType and batchID) should be numeric!")
    }

    if (!all(sapply(test_samples[!names(test_samples) %in%c(col_sampleID, col_sampleType, col_batchID)], is.numeric))) {
        stop("  - The values of test samples (except sampleType and batchID) should be numeric!")
    }

    batchID_train <- unique(train_samples[[col_batchID]])
    batchID_test <- unique(test_samples[[col_batchID]])

    if(!all(batchID_test %in% batchID_train)) stop("  - The batchID in train and test samples cannot match!")
    batchID <- batchID_test

    # Target value computation
    if (is.null(targetVal_external)) {
        message("+ Computing target values...   ", Sys.time())
        targetVal_list <- compute_targetVal(QC_data = train_samples[!names(train_samples) %in% c(col_sampleID, col_order, col_position)],
                                            col_sampleType = col_sampleType,
                                            col_batchID = col_batchID,
                                            targetVal_method = targetVal_method,
                                            targetVal_batch = targetVal_batch,
                                            targetVal_removeOutlier = targetVal_removeOutlier,
                                            coerce_numeric = FALSE)
    } else {
        message("+ External target values loaded.   ", Sys.time())
        targetVal_list <- targetVal_external
    }

    message("+ Checking variable names...   ", Sys.time())
    var_names <- names(targetVal_list[[1]])
    train_num <- train_samples[!names(train_samples) %in% c(col_sampleID, col_sampleType, col_batchID, col_order, col_position)]
    test_num <- test_samples[!names(test_samples) %in% c(col_sampleID, col_sampleType, col_batchID, col_order, col_position)]

    if (!all(var_names == names(train_num))) stop("  - Varibale names in the train and test samples cannot match!")
    if (!all(var_names == names(test_num)))  stop("  - Varibale names in the train and test samples cannot match!")

    # Variable selection
    message("+ Creating clusters...   ", Sys.time())
    cl <- parallel::makeCluster(parallel.cores)
    message("+ Selecting highly-correlated variables...   ", Sys.time())
    var_selected <- select_variable(train_num = train_num,
                                    test_num = test_num,
                                    correlation_type = correlation_type,
                                    correlation_method = correlation_method,
                                    min_var_num = min_var_num, max_var_num = max_var_num,
                                    coerce_numeric = TRUE)
    pbapply::pboptions(type = "timer", style = 3, char = "=", txt.width = 70)
    message("+ Data correction started.   ", Sys.time())

    # Original sample order backup
    test_samples <- cbind(original_idx = 1:nrow(test_samples), test_samples)
    parallel::clusterExport(cl = cl, varlist = c("Internal.compute_errorRatio", "Internal.run_ensemble"))

    res_var <- pbapply::pblapply(var_names, function(current_var, var_selected, targetVal_batch,
                                                     train_samples, test_samples, col_sampleID, col_sampleType,
                                                     col_batchID, col_order, col_position, batchID, mtry_ratio,
                                                     nodesize_ratio, ...) {
        # message("  - Current variable: ", current_var, "   ", Sys.time())

        train_X_selected_var <- train_samples[c(col_sampleID, col_sampleType, col_batchID,
                                                col_order, col_position, var_selected[[current_var]]) ]

        train_X_selected_var <- train_X_selected_var[!apply(train_X_selected_var, 1, anyNA),]
        if (any(table(train_X_selected_var[[col_batchID]]) < 5)) stop("  - Your train data contain too many NA. Maybe you would like to impute your data!")

        test_data <- cbind(y_raw = test_samples[[current_var]], test_samples)

        if (!targetVal_batch) {
            train_y_all <- Internal.compute_errorRatio(train_samples = train_samples[!names(train_samples) %in% c(col_sampleID, col_batchID, col_order, col_position)],
                                                       col_sampleType = col_sampleType,
                                                       targetVal_df = targetVal_list$wholeDataset,
                                                       current_var = current_var)

            train_data_all <- cbind(y_target = train_y_all$targetVal, y_raw = train_y_all$rawVal,
                                    y = train_y_all$errorRatio, train_X_selected_var)
        }

        res_batch_list <- lapply(batchID, function(current_batch) {
            # message("    - Current batch: ", current_batch, "   ", Sys.time())

            if (targetVal_batch) {
                train_X_batch <- train_X_selected_var[train_X_selected_var[[col_batchID]] == current_batch,]

                train_y_batch <- Internal.compute_errorRatio(train_samples = train_X_batch[!names(train_X_batch) %in% c(col_sampleID, col_batchID, col_order, col_position)],
                                                             col_sampleType = col_sampleType,
                                                             targetVal_df = targetVal_list[[current_batch]],
                                                             current_var = current_var)

                train_data <- cbind(y_target = train_y_batch$targetVal, y_raw = train_y_batch$rawVal,
                                    y = train_y_batch$errorRatio, train_X_batch)
            } else {
                train_data <- train_data_all[train_data_all[[col_batchID]] == current_batch,]
            }

            trainSet <- train_data[!names(train_data) %in% c(col_sampleID, col_sampleType, col_batchID)]
            testSet  <- test_data[test_data[[col_batchID]] == current_batch,]

            var_pred <- Internal.run_ensemble(trainSet = trainSet, testSet = testSet,
                                              mtry_ratio = seq(0.2, 0.8, 0.2),
                                              nodesize_ratio = seq(0.2, 0.8, 0.2),
                                              ... = ..., return_base_res = FALSE)

            if (targetVal_batch) {
                test_targetVal_all   <- do.call(targetVal_method, list(test_data$y_raw, na.rm = TRUE))
                test_targetVal_batch <- do.call(targetVal_method, list(testSet$y_raw,   na.rm = TRUE))
                var_pred <- var_pred * test_targetVal_all / test_targetVal_batch
                # var_pred <- var_pred * test_targetVal_batch / test_targetVal_all
            }

            names(var_pred) <- testSet$original_idx
            var_pred
        })
        res_batch <- do.call("c", res_batch_list)
        res_batch_order <- res_batch[order(as.numeric(names(res_batch)))]
        res_batch_df <- data.frame(res_batch_order)
        names(res_batch_df) <- current_var
        res_batch_df
    }, var_selected = var_selected, targetVal_batch = targetVal_batch,
    train_samples = train_samples, test_samples = test_samples, col_sampleID = col_sampleID,
    col_sampleType = col_sampleType, col_batchID = col_batchID, col_order = col_order,
    col_position = col_position, batchID = batchID, mtry_ratio = mtry_ratio,
    nodesize_ratio = nodesize_ratio, ... = ..., cl = cl)

    parallel::stopCluster(cl)

    message("+ Merging results...   ", Sys.time())
    check_order <- sapply(res_var[-1], function(x) {
        any(row.names(x) != row.names(res_var[[1]]))
    })
    if (any(check_order)) stop("Error occurs when merging data!")

    res_var_df <- do.call("cbind", res_var)

    test_samples[names(res_var_df)] <- res_var_df
    test_samples$original_idx <- NULL
    message("+ Completed.   ", Sys.time())
    test_samples
}
