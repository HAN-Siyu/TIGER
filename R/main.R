#' Compute RSD (relative standard deviation)
#' @description This function compute the RSD (relative standard deviation) of the values in input_data. Missing values are removed before the computation automatically.
#' @param input_data a numeric vector
#' @details The RSD in this function is compuated with the formula:
#'
#' \code{sd(input_data, na.rm = T) / mean(input_data, na.rm = T)}.
#' @examples
#' compute_RSD(c(1:10))
#' @export

compute_RSD <- function(input_data) {
    val_RSD <- sd(input_data, na.rm = T) / mean(input_data, na.rm = T)
    val_RSD
}

#' Compute target values for ensemble learning architecture
#' @description
#' @param QC_num a numeric data.frame including the matabolite values of quality control (QC) samples. Row: sample. Column: metabolite variable. See Examples.
#' @param sampleType a vector corresponding to \code{QC_num} to specify the type of each sample. See Examples.
#' @param batchID a vector corresponding to \code{QC_num} to specify the batch of each sample. See Examples.
#' @param targetVal_method a character string specifying how the target values are computed. Can be \code{"mean"} (default) or \code{"median"}.
#' @param targetVal_batchWise logical. If \code{TRUE}, the target values will be computed based on each batch, otherwise, based on the whole dataset. Setting \code{TRUE} might be useful if your dataset has very obvious batch effects, but this may also make the algorithm less robust. Default: \code{FALSE}.
#' @param targetVal_removeOutlier logical. If \code{TRUE}, outliers will be removed before the computation. Outliers are determined with 1.5 * IQR (interquartile range) rule. We recommend turning this off when the target values are computed based on batches. Default: \code{!targetVal_batchWise}.
#' @param coerce_numeric logical. If \code{TRUE}, values in \code{QC_num} will be coerced to numeric before the computation. Default: \code{FALSE}.
#' @importFrom stats fivenum
#' @details Code for checking outliers is adapted from \code{\link[grDevices]{boxplot.stats}}.
#' @examples
#' @export

compute_targetVal <- function(QC_num, sampleType, batchID,
                              targetVal_method = c("mean", "median"),
                              targetVal_batchWise = FALSE,
                              targetVal_removeOutlier = !targetVal_batchWise,
                              coerce_numeric = FALSE) {

    message("+ Computing target values...   ", Sys.time())

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
        if (!all(sapply(QC_num, is.numeric))) stop("  The values of the input dataset (QC_num) should be numeric!")
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
#' @param train_num a numeric data.frame including the matabolite values of training samples (can be quality control samples). Row: sample. Column: metabolite variable. See Examples.
#' @param test_num an optional numeric data.frame including the matabolite values of test samples (can be subject samples). Row: sample. Column: metabolite variable. See Examples. If \code{NULL}, the variables will be selected based on \code{train_num} only.
#' @param train_batchID \code{NULL} or a vector corresponding to \code{train_num} to specify the batch of each sample. Ignored if \code{correlation_batchWise = FALSE}. See Examples.
#' @param test_batchID \code{NULL} or a vector corresponding to \code{test_num} to specify the batch of each sample. Ignored if \code{correlation_batchWise = FALSE}. See Examples.
#' @param correlation_batchWise logical. Specify whether the variable selection should be performed for each batch. Setting \code{TRUE} might be useful if your dataset has very obvious batch effects, but this may also make the algorithm less robust. Default: \code{FALSE}.
#' @param correlation_type a character string indicating correlation (\code{"cor"}, default) or partial correlation (\code{"pcor"}) is to be used. Can be abbreviated. \strong{Note}: computing partial correlations of a large dataset can be very time-consuming.
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

    message("+ Selecting highly-correlated variables...   ", Sys.time())

    correlation_type   <- match.arg(correlation_type)

    min_var_num <- ifelse(is.null(min_var_num), 1, min_var_num)
    max_var_num <- ifelse(is.null(max_var_num), ncol(train_num), max_var_num)

    # min_var_num <- as.integer(min_var_num) + 1 # add 1 as the current variable itself is included here,
    # max_var_num <- as.integer(max_var_num) + 1 # but the current variable will be removed when training the model.
    min_var_num <- as.integer(min_var_num)
    max_var_num <- as.integer(max_var_num)

    if (min_var_num < 1) stop("  min_var_num must be a positive integer!")
    if (max_var_num > ncol(train_num)) stop("  max_var_num cannot be greater than variable number!")

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
        if (!all(sapply(train_num, is.numeric))) stop("  The values of train samples should be numeric!")
        if (!is.null(test_num) & !all(sapply(test_num, is.numeric))) stop("  The values of test samples should be numeric!")
    }

    if (!is.null(test_num)) {
        if (!all(names(test_num) %in% names(train_num))) stop("  Variables in training and test data cannot match!")
        train_num <- train_num[names(test_num)]
    }

    if (correlation_batchWise) {
        train_num_list <- split(train_num, f = train_batchID)
        test_num_list  <- split(test_num,  f = test_batchID)

        batch_names_train <- sort(names(train_num_list))
        batch_names_test  <- sort(names(test_num_list))
        batch_names_len   <- length(batch_names_train)
        if (any(batch_names_train != batch_names_test) | batch_names_len != length(batch_names_test)) stop("    Batch names of train and test samples cannot match!")

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

#' Run TIGER to eliminate technical variation
#' @description
#'
#' @param test_samples (required) a data.frame containing the test samples (subject samples). This data.frame should contain columns of
#' \itemize{
#' \item sample ID (required): label for each sample,
#' \item sample type (required): indicating the type of each sample,
#' \item batch ID (required): the batch of each sample,
#' \item order information (optional): injection order or temporal information of each sample,
#' \item position information (optional): well position of each sample)columns of metabolite values,
#' \item metabolite values (required): values to be normalised.
#' }
#' Row: sample. Column: variable. See Examples.
#' @param train_samples (required) a data.frame containing the training samples (quality control samples). The columns in this data.frame should correspond to the columns in \code{test_samples}. \code{test_samples} and \code{train_samples} should have the identical column names.
#' @param col_sampleID  (required) a character string indicating the name of the column that specifies the sample ID of each sample. The values in this column will not affect the data correction process but can act as labels for different samples. See Examples.
#' @param col_sampleType (required) a character string indicating the name of the column that specifies the type (such as QC1, QC2, subject) of each sample. This column in the \code{train_samples} can be used to indicate different kinds of QC samples. \strong{QC samples of the same type should have the same type name.} See Examples.
#' @param col_batchID (required) a character string indicating the name of the column that specifies the batch ID of each sample. See Examples.
#' @param col_order (optional) \code{NULL} or a character string indicating the name of the column that contains the injection order or temporal information. This can explicitly ask the algorithm capture the technical variation introduced by injection order, which might be useful when the data have very obvious temporal drifts. If \code{NULL} (default), the input \code{train_samples} and \code{test_samples} should have \strong{no} column contains injection order information.
#' @param col_position (optional) \code{NULL} or a character string indicating the name of the column that contains the well position information. This can explicitly ask the algorithm capture the technical variation introduced by well position, which might be useful when the well position has a great impact during data acquisition. If \code{NULL} (default), the input \code{train_samples} and \code{test_samples} should have \strong{no} column contains well position information.
#' @param targetVal_external (optional) a list generated by function \code{\link{compute_targetVal}}. See Details.
#' @param targetVal_method a character string specifying how the target values are computed. Can be \code{"mean"} (default) or \code{"median"}.
#' @param targetVal_batchWise logical. If \code{TRUE}, the target values will be computed based on each batch, otherwise, based on the whole dataset. Setting \code{TRUE} might be useful if your dataset has very obvious batch effects, but this may also make the algorithm less robust. Default: \code{FALSE}.
#' @param targetVal_removeOutlier logical. If \code{TRUE}, outliers will be removed before the computation. Outliers are determined with 1.5 * IQR (interquartile range) rule. We recommend turning this off when the target values are computed based on batches. Default: \code{!targetVal_batchWise}.
#' @param correlation_type a character string indicating correlation (\code{"cor"}, default) or partial correlation (\code{"pcor"}) is to be used. Can be abbreviated. \strong{Note}: computing partial correlations of a large dataset can be very time-consuming.
#' @param correlation_method a character string indicating which correlation coefficient is to be computed. One of \code{"spearman"} (default) or \code{"pearson"}. Can be abbreviated.
#' @param correlation_batchWise logical. Specify whether the variable selection should be performed for each batch. Setting \code{TRUE} might be useful if your dataset has very obvious batch effects, but this may also make the algorithm less robust. Default: \code{FALSE}.
#' @param min_var_num an integer specifying the minimum number of the selected variables. If \code{NULL}, no limited, but 1 at least. Default: \code{5}.
#' @param max_var_num an integer specifying the maximum number of the selected variables. If \code{NULL}, no limited, but \code{ncol(train_num)} at most. Default: \code{10}.
#' @param mtry_percent (advanced) a numeric vector indicating the percentages of selected variables randomly sampled as candidates at each split when training the random forest models. Providing more values will train more models, which will increase the processing time. Default: \code{seq(0.2, 0.8, 0.2)}.
#' @param nodesize_percent (advanced) a numeric vector indicating the percentages of sample size used as the minimum sizes of the terminal nodes in random forest models. Providing more values will train more models, which will increase the processing time. Default: \code{seq(0.2, 0.8, 0.2)}.
#' @param ... (advanced) optional arguments (except \code{mtry} and \code{nodesize}) to be passed to \code{\link[randomForest]{randomForest}} for model training. Arguments \code{mtry} and \code{nodesize} are determined by \code{mtry_percent} and \code{nodesize_percent}. Default values are used for other arguments in \code{\link[randomForest]{randomForest}}. Providing more arguments will train more models, which will increase the processing time.
#' @param parallel.cores an integer (== -1 or >= 1) specifying the number of cores for parallel computation. Setting \code{-1} to run with all cores. Default: \code{2}.
#'
#' @details
#'
#' @section References:
#' Han S. \emph{et al}. TIGER: technical variation elimination for metabolomics data using ensemble learning architecture. (\emph{Submitted})
#'
#' @importFrom caret createFolds
#' @importFrom parallel detectCores
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
                      ..., parallel.cores = 2, logName) {

    message("+ Initialising...   ", Sys.time())

    targetVal_method   <- match.arg(targetVal_method)
    correlation_type   <- match.arg(correlation_type)
    correlation_method <- match.arg(correlation_method)

    for (col_idx in c(col_sampleID, col_sampleType, col_batchID)) {
        test_samples[[col_idx]]  <- as.character(test_samples[[col_idx]])
        train_samples[[col_idx]] <- as.character(train_samples[[col_idx]])
    }

    if (!is.null(col_order)) {
        if (anyNA(test_samples[[col_order]])  | any(!is.finite(test_samples[[col_order]]))  ) stop("  test samples: col_order should be numeric only!")
        if (anyNA(train_samples[[col_order]]) | any(!is.finite(train_samples[[col_order]])) ) stop("  train samples: col_order should be numeric only!")
    }

    if (!is.null(col_position)) {
        if (anyNA(test_samples[[col_position]])  | any(!is.finite(test_samples[[col_position]]))  ) stop("  test samples: col_position should be numeric only!")
        if (anyNA(train_samples[[col_position]]) | any(!is.finite(train_samples[[col_position]])) ) stop("  train samples: col_position should be numeric only!")
    }

    if (!all(sapply(train_samples[!names(train_samples) %in% c(col_sampleID, col_sampleType, col_batchID)], is.numeric))) {
        stop("  The values of train samples (except sampleType and batchID) should be numeric!")
    }

    if (!all(sapply(test_samples[!names(test_samples) %in%c(col_sampleID, col_sampleType, col_batchID)], is.numeric))) {
        stop("  The values of test samples (except sampleType and batchID) should be numeric!")
    }

    # To check infinite values. - Deprecated. Infinite values are not allowed.
    # train_samples_list  <- split(train_samples, f = train_samples[[col_sampleType]])
    # train_samples_check <- lapply(train_samples_list, Internal.impute_infinite)
    # train_samples <- do.call("rbind", train_samples_check)
    # test_samples  <- Internal.impute_infinite(test_samples)

    test_samples_bak <- test_samples

    batchID_train <- unique(train_samples[[col_batchID]])
    batchID_test  <- unique(test_samples[[col_batchID]])

    if(!all(batchID_test %in% batchID_train)) stop("  The batchID in train and test samples cannot match!")
    batchID <- batchID_test

    train_num <- train_samples[!names(train_samples) %in% c(col_sampleID, col_sampleType, col_batchID, col_order, col_position)]
    test_num  <- test_samples[!names(test_samples)   %in% c(col_sampleID, col_sampleType, col_batchID, col_order, col_position)]

    var_names <- names(test_num)

    # Target value computation
    if (is.null(targetVal_external)) {
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
    var_selected_list <- select_variable(train_num = train_num, test_num = test_num,
                                         train_batchID = train_samples[[col_batchID]],
                                         test_batchID  = test_samples[[col_batchID]],
                                         correlation_batchWise = correlation_batchWise,
                                         correlation_type   = correlation_type,
                                         correlation_method = correlation_method,
                                         min_var_num = min_var_num,
                                         max_var_num = max_var_num,
                                         coerce_numeric = TRUE)
    idx_test_na  <- is.na(test_samples)
    idx_train_na <- is.na(train_samples)
    if (any(idx_test_na))  test_samples[idx_test_na]   <- 0
    if (any(idx_train_na)) train_samples[idx_train_na] <- 0

    message("+ Data correction started.   ", Sys.time())
    message("  - Creating clusters...")
    parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
    cl <- parallel::makeCluster(parallel.cores, outfile = logName)
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
        message("debug: into loop of var")
        if (!targetVal_batchWise) {
            train_y_all <- Internal.compute_errorRatio(rawVal     = train_samples[[current_var]],
                                                       sampleType = train_samples[[col_sampleType]],
                                                       targetVal  = targetVal_list$wholeDataset[current_var])

        }

        test_data <- cbind(y_raw = test_samples[[current_var]], test_samples)

        res_batch_list <- lapply(batchID, function(current_batch) {
            message("debug: into loop of batch")

            train_X_batch <- train_samples[train_samples[[col_batchID]] == current_batch,]

            if (targetVal_batchWise) {
                train_y_batch <- Internal.compute_errorRatio(rawVal     = train_X_batch[[current_var]],
                                                             sampleType = train_X_batch[[col_sampleType]],
                                                             targetVal  = targetVal_list[[current_batch]][current_var])


            } else {
                train_y_batch <- train_y_all[train_samples[[col_batchID]] == current_batch,]
            }

            if (correlation_batchWise) {
                train_X_selected <- train_X_batch[ c(col_order, col_position,
                                                     var_selected_list[[current_batch]][[current_var]]) ]
            } else {
                train_X_selected <- train_X_batch[ c(col_order, col_position,
                                                     var_selected_list$wholeDataset[[current_var]]) ]
            }

            trainSet <- cbind(train_y_batch, train_X_selected)
            testSet  <- test_data[test_data[[col_batchID]] == current_batch,]

            message("debug: into ensemble")
            cat("debug: current var:", current_var, "current batch:", current_batch,
                "shape:", nrow(trainSet), "*", ncol(trainSet),
                "selected var:", names(trainSet), "\n")
            var_pred <- Internal.run_ensemble(trainSet = trainSet, testSet = testSet,
                                              mtry_percent = seq(0.2, 0.8, 0.2),
                                              nodesize_percent = seq(0.2, 0.8, 0.2),
                                              ... = ..., return_base_res = FALSE)

            if (targetVal_batchWise) {
                message("debug: out ensemble, targetVal_batchWise - convert back")
                test_targetVal_all   <- do.call(targetVal_method, list(test_data$y_raw, na.rm = TRUE))
                test_targetVal_batch <- do.call(targetVal_method, list(testSet$y_raw,   na.rm = TRUE))
                var_pred <- var_pred * test_targetVal_all / test_targetVal_batch
            }
            message("debug: assign names")
            names(var_pred) <- testSet$original_idx
            var_pred
        })
        message("debug: out loop of batch")
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
    if (any(check_order)) stop("    Error occurs when merging data!")

    res_var_df <- do.call("cbind", res_var)

    test_samples[names(res_var_df)] <- res_var_df
    test_samples$original_idx <- NULL

    idx_norm_na <- is.na(test_samples)
    if (any(idx_norm_na)) test_samples[idx_norm_na] <- as.numeric(test_samples_bak[idx_norm_na])

    idx_norm_zero <- test_samples < 0
    if (any(idx_norm_zero)) test_samples[idx_norm_zero] <- as.numeric(test_samples_bak[idx_norm_zero])

    test_samples[test_samples_bak == 0] <- 0
    if (any(idx_test_na)) test_samples[idx_test_na] <- NA

    message("+ Completed.   ", Sys.time())
    test_samples
}
