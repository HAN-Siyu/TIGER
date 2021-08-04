#' Compute RSD (relative standard deviation)
#' @description This function compute the RSD (relative standard deviation) of the values in input_data. Missing values are removed before the computation automatically.
#' @param input_data a numeric vector
#' @details The RSD in this function is compuated with the formula:
#' \code{sd(input_data, na.rm = T) / mean(input_data, na.rm = T)}.
#' @examples
#' compute_RSD(c(1:10))
#' @export

compute_RSD <- function(input_data) {
    val_RSD <- sd(input_data, na.rm = T) / mean(input_data, na.rm = T)
    val_RSD
}

#' Compute the target values for ensemble learning architecture
#' @description
#' @param QC_num a numeric data.frame including the matabolite values of quality control (QC) samples. Row: sample. Column: metabolite variable. see Examples.
#' @param sampleType a vector corresponding to \code{QC_num} to specify the type of each sample. see Examples.
#' @param batchID a vector corresponding to \code{QC_num} to specify the batch of each sample. see Examples.
#' @param targetVal_method a character string specifying how the target values are computed. Can be \code{"mean"} (default) or \code{"median"}.
#' @param targetVal_batchWise logical. If \code{TRUE}, the target values will be computed based on each batch, otherwise, based on the whole dataset. Setting \code{TRUE} might be useful if your dataset has very obvious batch effects, but this may also make the algorithm less robust. Default: \code{FALSE}.
#' @param targetVal_removeOutlier logical. If \code{TRUE}, outliers will be removed before the computation. We recommend turning this off when the target values are computed based on batches. Default: \code{!targetVal_batchWise}.
#' @param coerce_numeric logical. If \code{TRUE}, values in \code{QC_num} will be coerced to numeric before the computation. Default: \code{FALSE}.
#' @importFrom stats fivenum
#' @details
#' @examples
#' @export

compute_targetVal <- function(QC_num, sampleType, batchID,
                              targetVal_method = c("mean", "median"),
                              targetVal_batchWise = FALSE,
                              targetVal_removeOutlier = !targetVal_batchWise,
                              coerce_numeric = FALSE) {

    targetVal_method   <- match.arg(targetVal_method)

    sampleType <- as.factor(sampleType)
    batchID  <- as.factor(batchID)

    if(coerce_numeric) {
        QC_num <- as.data.frame(sapply(QC_num, as.numeric))
        idx_NA <- sapply(QC_num, function(x) {
            all(is.na(x))
        })
        QC_num <- QC_num[,!idx_NA]
    } else {
        if (!all(sapply(QC_num, is.numeric))) stop("The values of the input dataset (QC_num) should be numeric!")
    }

    if (targetVal_batchWise) {
        target_values <- aggregate(QC_num, by = list(batch = batchID, sample = sampleType),
                                   FUN = Internal.compute_targetVal, targetVal_method = targetVal_method,
                                   targetVal_removeOutlier = targetVal_removeOutlier)
        batchID <- target_values$batch
        target_values$batch <- NULL
        target_values_list <- split(target_values, f = batchID)
    } else {
        target_values_list <- list(wholeDataset = aggregate(QC_num,
                                                            by = list(sample = sampleType),
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


#' Select variables for ensemble learning architecture
#' @description
#' @param train_num a numeric data.frame including the matabolite values of training samples (can be quality control samples). Row: sample. Column: metabolite variable. see Examples.
#' @param test_num a optional numeric data.frame including the matabolite values of test samples (can be subject samples). Row: sample. Column: metabolite variable. see Examples. If \code{NULL}, the variables will be selected based on \code{train_num} only.
#' @param train_batchID ignored if \code{correlation_batchWise = FALSE}. \code{NULL} or a vector corresponding to \code{train_num} to specify the batch of each sample. see Examples.
#' @param test_batchID ignored if \code{correlation_batchWise = FALSE}. \code{NULL} or a vector corresponding to \code{test_num} to specify the batch of each sample. see Examples.
#' @param correlation_batchWise logical. Specify whether the variable selection should be performed for each batch. Setting \code{TRUE} might be useful if your dataset has very obvious batch effects, but this may also make the algorithm less robust. Default: \code{FALSE}.
#' @param correlation_type a character string indicating correlation (\code{"cor"}, default) or partial correlation (\code{"pcor"}) is to be used. Can be abbreviated. Note: computing partial correlations of a large dataset can be very time-consuming.
#' @param correlation_method a character string indicating which correlation coefficient is to be computed. One of \code{"spearman"} (default) or \code{"pearson"}. Can be abbreviated.
#' @param min_var_num an integer specifying the minimum number of the selected variables. If \code{NULL}, no limited, but 1 at least. Default: 5.
#' @param max_var_num an integer specifying the maximum number of the selected variables. If \code{NULL}, no limited, but \code{ncol(train_num)} at most. Default: 10.
#' @param coerce_numeric logical. If \code{TRUE}, values in \code{train_num} and  \code{test_num} will be coerced to numeric before the computation. Default: \code{FALSE}.
#' @importFrom pbapply timerProgressBar
#' @importFrom pbapply setTimerProgressBar
#' @importFrom pbapply closepb
#' @importFrom ppcor pcor
#' @details
#' @examples
#' @export

select_variable <- function(train_num, test_num = NULL,
                            train_batchID = NULL, test_batchID = NULL,
                            correlation_batchWise = FALSE,
                            correlation_type = c("cor", "pcor"),
                            correlation_method = c("spearman", "pearson"),
                            min_var_num = 5, max_var_num = 10,
                            coerce_numeric = FALSE) {

    correlation_type   <- match.arg(correlation_type)

    min_var_num <- ifelse(is.null(min_var_num), 1, min_var_num)
    max_var_num <- ifelse(is.null(max_var_num), ncol(train_num), max_var_num)

    # min_var_num <- as.integer(min_var_num) + 1 # add 1 as the current variable itself is included here,
    # max_var_num <- as.integer(max_var_num) + 1 # but the current variable will be removed when training the model.
    min_var_num <- as.integer(min_var_num)
    max_var_num <- as.integer(max_var_num)

    if (min_var_num < 1) stop("    min_var_num must be a positive integer!")
    if (max_var_num > ncol(train_num)) stop("    max_var_num cannot be greater than variable number!")

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
        if (!all(sapply(train_num, is.numeric))) stop("    The values of train samples should be numeric!")
        if (!is.null(test_num) & !all(sapply(test_num, is.numeric))) stop("    The values of test samples should be numeric!")
    }

    if (!is.null(test_num)) {
        if (any(names(train_num) != names(test_num))) stop("    Variables in training and test data cannot match!")
    }

    if (correlation_batchWise) {
        train_num_list <- split(train_num, f = train_batchID)
        test_num_list  <- split(test_num,  f = test_batchID)

        batch_names_train <- sort(names(train_num_list))
        batch_names_test  <- sort(names(test_num_list))
        batch_names_len   <- length(batch_names_train)
        if (!all(batch_names_train == batch_names_test) | batch_names_len != length(batch_names_test)) stop("    Batch names of train and test samples cannot match!")

        selected_var_list <- lapply(1:batch_names_len, function(batch_name_idx) {

            one_batch_name <- batch_names_train[[batch_name_idx]]
            message("  - Current batch: ", one_batch_name, "  (", batch_name_idx, "/", batch_names_len, ")")

            message("    Computing correlation coefficients...")
            cor_info <- Internal.compute_cor(train_num = train_num_list[[one_batch_name]],
                                             test_num = test_num_list[[one_batch_name]],
                                             correlation_type = correlation_type,
                                             correlation_method = correlation_method)

            message("    Selecting variables...")
            selected_var <- Internal.select_variable(cor_info = cor_info,
                                                     min_var_num = min_var_num,
                                                     max_var_num = max_var_num)
        })
        names(selected_var_list) <- batch_names_train
    } else {
        message("  - Computing correlation coefficients...")
        cor_info <- Internal.compute_cor(train_num = train_num, test_num = test_num,
                                         correlation_type = correlation_type,
                                         correlation_method = correlation_method)
        message("  - Selecting variables...")
        selected_var <- Internal.select_variable(cor_info = cor_info,
                                                 min_var_num = min_var_num,
                                                 max_var_num = max_var_num)
        selected_var_list <- list(wholeDataset = selected_var)
    }

    selected_var_list
}

#'
#' @description
#'
#' @param test_samples (required)
#' @param train_samples (required)
#' @param col_sampleID  (required)
#' @param col_sampleType (required)
#' @param col_batchID (required)
#' @param col_order (optional)
#' @param col_position (optional)
#' @param targetVal_external (optional)
#' @param targetVal_method Default: \code{"mean"}.
#' @param targetVal_batchWise Default: \code{FALSE}.
#' @param targetVal_removeOutlier Default: \code{!targetVal_batchWise}.
#' @param correlation_type Default: \code{"cor"}.
#' @param correlation_method Default: \code{"pearson"}.
#' @param correlation_batchWise Default: \code{FALSE}.
#' @param min_var_num Default: \code{5}.
#' @param max_var_num Default: \code{10}.
#' @param mtry_percent (advanced) a numeric vector indicating the percentages of selected variables randomly sampled as candidates at each split when training the random forest models. Default: \code{seq(0.2, 0.8, 0.2)}. Providing more values will train more models, which will increase the processing time.
#' @param nodesize_percent (advanced) a numeric vector indicating the percentages of sample size as the minimum sizes of the terminal nodes in random forest models. Default: \code{seq(0.2, 0.8, 0.2)}. Providing more values will train more models, which will increase the processing time.
#' @param ... (advanced) optional arguments (except \code{mtry} and \code{nodesize}) to be passed to \code{\link[randomForest]{randomForest}} for model training. Arguments \code{mtry} and \code{nodesize} are determined by \code{mtry_percent} and \code{nodesize_percent}. Default values are used for other arguments in \code{\link[randomForest]{randomForest}}. Providing more arguments will train more models, which will increase the processing time.
#' @param parallel.cores Default: \code{2}.
#'
#' @details
#'
#' @section References:
#'
#' @importFrom caret createFolds
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom pbapply pboptions
#' @importFrom pbapply pblapply
#' @importFrom pbapply timerProgressBar
#' @importFrom pbapply setTimerProgressBar
#' @importFrom pbapply closepb
#' @importFrom ppcor pcor
#' @importFrom randomForest randomForest
#' @importFrom stats fivenum
#' @details
#' @examples
#' @export
#'
run_TIGER <- function(test_samples, train_samples,
                      col_sampleID, col_sampleType, col_batchID,
                      col_order = NULL, col_position = NULL,
                      targetVal_external = NULL, targetVal_method = c("mean", "median"),
                      targetVal_batchWise = FALSE, targetVal_removeOutlier = !targetVal_batchWise,
                      correlation_type = c("cor", "pcor"),
                      correlation_method = c("pearson", "spearman"),
                      correlation_batchWise = targetVal_batchWise,
                      min_var_num = 5, max_var_num = 10,
                      mtry_percent = seq(0.2, 0.8, 0.2),
                      nodesize_percent = seq(0.2, 0.8, 0.2),
                      ..., parallel.cores = 2) {

    message("+ Initialising...   ", Sys.time())

    targetVal_method   <- match.arg(targetVal_method)
    correlation_type   <- match.arg(correlation_type)
    correlation_method <- match.arg(correlation_method)

    for (col_idx in c(col_sampleID, col_sampleType, col_batchID)) {
        test_samples[[col_idx]]  <- as.character(test_samples[[col_idx]])
        train_samples[[col_idx]] <- as.character(train_samples[[col_idx]])
    }

    if (!is.null(col_order)) {
        if (anyNA(test_samples[[col_order]])  | any(!is.finite(test_samples[[col_order]]))  ) stop("    test samples: col_order should be numeric only!")
        if (anyNA(train_samples[[col_order]]) | any(!is.finite(train_samples[[col_order]])) ) stop("    train samples: col_order should be numeric only!")
    }

    if (!is.null(col_position)) {
        if (anyNA(test_samples[[col_position]])  | any(!is.finite(test_samples[[col_position]]))  ) stop("    test samples: col_position should be numeric only!")
        if (anyNA(train_samples[[col_position]]) | any(!is.finite(train_samples[[col_position]])) ) stop("    train samples: col_position should be numeric only!")
    }

    if (!all(sapply(train_samples[!names(train_samples) %in% c(col_sampleID, col_sampleType, col_batchID)], is.numeric))) {
        stop("    The values of train samples (except sampleType and batchID) should be numeric!")
    }

    if (!all(sapply(test_samples[!names(test_samples) %in%c(col_sampleID, col_sampleType, col_batchID)], is.numeric))) {
        stop("    The values of test samples (except sampleType and batchID) should be numeric!")
    }

    test_samples_bak <- test_samples

    batchID_train <- unique(train_samples[[col_batchID]])
    batchID_test  <- unique(test_samples[[col_batchID]])

    if(!all(batchID_test %in% batchID_train)) stop("    The batchID in train and test samples cannot match!")
    batchID <- batchID_test

    train_num <- train_samples[!names(train_samples) %in% c(col_sampleID, col_sampleType, col_batchID, col_order, col_position)]
    test_num  <- test_samples[!names(test_samples) %in% c(col_sampleID, col_sampleType, col_batchID, col_order, col_position)]

    if ((length(train_num) != length(test_num)) | any(!names(train_num) %in% names(test_num))) stop("    Varibale names in the train and test samples cannot match!")

    # Target value computation
    if (is.null(targetVal_external)) {
        message("+ Computing target values...   ", Sys.time())
        targetVal_list <- compute_targetVal(QC_num = train_num,
                                            sampleType = train_samples[[col_sampleType]],
                                            batchID    = train_samples[[col_batchID]],
                                            targetVal_method = targetVal_method,
                                            targetVal_batchWise = targetVal_batchWise,
                                            targetVal_removeOutlier = targetVal_removeOutlier,
                                            coerce_numeric = FALSE)
    } else {
        message("+ External target values loaded.   ", Sys.time())
        targetVal_list <- targetVal_external
    }

    # Variable selection
    message("+ Selecting highly-correlated variables...   ", Sys.time())

    # To check infinite values. - Deprecated. Infinite values are not allowed.
    # train_samples_list  <- split(train_samples, f = train_samples[[col_sampleType]])
    # train_samples_check <- lapply(train_samples_list, Internal.impute_infinite)
    # train_samples <- do.call("rbind", train_samples_check)
    # test_samples  <- Internal.impute_infinite(test_samples)

    # message("  - Checking variable names...")
    # var_names <- names(targetVal_list[[1]])
    #
    #
    # if (!all(var_names == names(train_num))) stop("    Varibale names in the train and test samples cannot match!")
    # if (!all(var_names == names(test_num)))  stop("    Varibale names in the train and test samples cannot match!")

    var_selected_list <- select_variable(train_num = train_num, test_num = test_num,
                                         train_batchID = train_samples[[col_batchID]],
                                         test_batchID  = test_samples[[col_batchID]],
                                         correlation_batchWise = correlation_batchWise,
                                         correlation_type   = correlation_type,
                                         correlation_method = correlation_method,
                                         min_var_num = min_var_num,
                                         max_var_num = max_var_num,
                                         coerce_numeric = TRUE)
    idx_test_na <- is.na(test_samples)
    idx_train_na <- is.na(train_samples)
    if (any(idx_test_na))  test_samples[idx_test_na]   <- 0
    if (any(idx_train_na)) train_samples[idx_train_na] <- 0

    message("+ Data correction started.   ", Sys.time())
    message("  - Creating clusters...")
    cl <- parallel::makeCluster(parallel.cores, outfile = "log")
    parallel::clusterExport(cl = cl, varlist = c("Internal.compute_errorRatio", "Internal.run_ensemble"))
    pbapply::pboptions(type = "timer", style = 3, char = "=", txt.width = 70)

    # Original sample order backup
    test_samples <- cbind(original_idx = 1:nrow(test_samples), test_samples)

    message("  - Correcting data...")
    res_var <- pbapply::pblapply(var_names, function(current_var, var_selected_list, targetVal_list,
                                                     targetVal_batchWise, correlation_batchWise,
                                                     train_samples, test_samples, col_sampleID, col_sampleType,
                                                     col_batchID, col_order, col_position, batchID, mtry_percent,
                                                     targetVal_method, nodesize_percent, ...) {
        message("into loop of var")
        # if (!correlation_batchWise) {
        #     train_X_selected_var <- train_samples[c(col_sampleID, col_sampleType, col_batchID,
        #                                             col_order, col_position,
        #                                             var_selected_list$wholeDataset[[current_var]]) ]
        #
        # }

        if (!targetVal_batchWise) {
            train_y_all <- Internal.compute_errorRatio(input_samples = train_samples[!names(train_samples) %in% c(col_sampleID, col_batchID, col_order, col_position)],
                                                       col_sampleType = col_sampleType,
                                                       targetVal_df = targetVal_list$wholeDataset,
                                                       current_var = current_var)

            train_data_all <- cbind(y_target = train_y_all$targetVal, y_raw = train_y_all$rawVal,
                                    y = train_y_all$errorRatio,
                                    train_samples[var_selected_list[[current_batch]][[current_var]]])
        }
        test_data <- cbind(y_raw = test_samples[[current_var]], test_samples)

        res_batch_list <- lapply(batchID, function(current_batch) {
            message("into loop of batch")
            # if (correlation_batchWise) {
            #     train_X_selected_var <- train_samples[c(col_sampleID, col_sampleType, col_batchID,
            #                                             col_order, col_position,
            #                                             var_selected_list[[current_batch]][[current_var]]) ]
            #
            # }

            if (targetVal_batchWise) {
                # train_X_batch <- train_X_selected_var[train_X_selected_var[[col_batchID]] == current_batch,]
                train_X_batch <- train_samples[train_samples[[col_batchID]] == current_batch,]

                train_y_batch <- Internal.compute_errorRatio(input_samples = train_X_batch[!names(train_X_batch) %in% c(col_sampleID, col_batchID, col_order, col_position)],
                                                             col_sampleType = col_sampleType,
                                                             targetVal_df = targetVal_list[[current_batch]],
                                                             current_var = current_var)

                trainSet <- cbind(y_target = train_y_batch$targetVal, y_raw = train_y_batch$rawVal,
                                  y = train_y_batch$errorRatio,
                                  train_X_batch[var_selected_list[[current_batch]][[current_var]]])
            } else {
                trainSet <- train_data_all[train_data_all[[col_batchID]] == current_batch,]
            }

            # trainSet <- train_data[!names(train_data) %in% c(col_sampleID, col_sampleType, col_batchID, current_var)]
            testSet  <- test_data[test_data[[col_batchID]] == current_batch,]

            message("into ensemble")
            cat("current var:", current_var, "current batch:", current_batch,
                "selected var:", names(trainSet))
            var_pred <- Internal.run_ensemble(trainSet = trainSet, testSet = testSet,
                                              mtry_percent = seq(0.2, 0.8, 0.2),
                                              nodesize_percent = seq(0.2, 0.8, 0.2),
                                              ... = ..., return_base_res = FALSE)

            if (targetVal_batchWise) {
                message("out ensemble, targetVal_batchWise - convert back")
                test_targetVal_all   <- do.call(targetVal_method, list(test_data$y_raw, na.rm = TRUE))
                test_targetVal_batch <- do.call(targetVal_method, list(testSet$y_raw,   na.rm = TRUE))
                var_pred <- var_pred * test_targetVal_all / test_targetVal_batch
            }
            message("assign names")
            names(var_pred) <- testSet$original_idx
            var_pred
        })
        message("out loop of batch")
        res_batch       <- do.call("c", res_batch_list)
        res_batch_order <- res_batch[order(as.numeric(names(res_batch)))]
        res_batch_df    <- data.frame(res_batch_order)

        names(res_batch_df) <- current_var
        res_batch_df

    }, var_selected_list = var_selected_list, targetVal_list = targetVal_list,
    targetVal_batchWise = targetVal_batchWise, correlation_batchWise = correlation_batchWise,
    train_samples = train_samples, test_samples = test_samples, col_sampleID = col_sampleID,
    col_sampleType = col_sampleType, col_batchID = col_batchID, col_order = col_order,
    col_position = col_position, batchID = batchID, mtry_percent = mtry_percent,
    nodesize_percent = nodesize_percent, targetVal_method = targetVal_method, ... = ..., cl = cl)

    parallel::stopCluster(cl)

    message("  - Merging results...")
    check_order <- sapply(res_var[-1], function(x) {
        any(row.names(x) != row.names(res_var[[1]]))
    })
    if (any(check_order)) stop("Error occurs when merging data!")

    res_var_df <- do.call("cbind", res_var)

    test_samples[names(res_var_df)] <- res_var_df
    test_samples$original_idx <- NULL

    test_samples[is.na(test_samples)] <- as.numeric(test_samples_bak[is.na(test_samples)])
    test_samples[test_samples < 0]    <- as.numeric(test_samples_bak[test_samples < 0])
    test_samples[test_samples_bak == 0] <- 0
    if (any(idx_test_na)) test_samples[idx_test_na] <- NA

    message("+ Completed.   ", Sys.time())
    test_samples
}
