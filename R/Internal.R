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
        warning(paste0("One sample in ", data_label, " contains NA and has been removed."))
        input_data_num <- input_data_num[!data_na_sample_idx,]
    } else if (data_na_sample_sum > 1) {
        warning(paste0(data_na_sample_sum, " samples in ", data_label, " contain NA and have been removed."))
        input_data_num <- input_data_num[!data_na_sample_idx,]
    }
    if (nrow(input_data_num) == 0) stop(paste0("Variable selection failed: ", data_label, " contain too many NA!"))
    input_data_num
}

Internal.compute_cor <- function(train_num, test_num = NULL,
                                 cor_type = c("pcor", "cor"),
                                 cor_method = c("pearson", "spearman")) {

    train_num_noNA <- Internal.remove_NA(train_num, data_label = "training data")
    if (!is.null(test_num)) test_num_noNA <- Internal.remove_NA(test_num, data_label = "test data")

    if (cor_type == "pcor") {
        train_cor <- data.frame(ppcor::pcor(train_num_noNA, method = cor_method)$estimate)

        if (!is.null(test_num)) {
            test_cor <- data.frame(ppcor::pcor(test_num_noNA, method = cor_method)$estimate)
        } else test_cor <- NULL

    } else {
        train_cor <- data.frame(cor(train_num_noNA, method = cor_method, use = "complete.obs"))

        if (!is.null(test_num)) {
            test_cor <- data.frame(cor(test_num_noNA, method = cor_method, use = "complete.obs"))
        } else test_cor <- NULL
    }

    cor_info <- list(variable_name = names(train_cor),
                     train_cor = train_cor, test_cor = test_cor)
}

Internal.select_variable <- function(cor_info, min_var_num = NULL,
                                     max_var_num = NULL, cl) {

    selected_var <- parallel::parLapply(cl, cor_info$variable_name, function(input_one_variable_name,
                                                                             train_cor, test_cor,
                                                                             min_var_num, max_var_num) {

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
    train_cor = cor_info$train_cor, test_cor = cor_info$test_cor,
    min_var_num = min_var_num, max_var_num = max_var_num)

    names(selected_var) <- cor_info$variable_name
    selected_var
}

select_variable <- function(train_num, test_num = NULL,
                            cor_type = c("pcor", "cor"),
                            cor_method = c("pearson", "spearman"),
                            min_var_num = NULL, max_var_num = NULL,
                            parallel.cores = 2, cl = NULL, output_log = NULL) {

    cor_type   <- match.arg(cor_type)
    cor_method <- match.arg(cor_method)

    if (!is.null(test_num)) {
        if (any(names(train_num) != names(test_num))) stop("Variables in training and test data cannot match!")
    }

    min_var_num <- ifelse(is.null(min_var_num), 1, min_var_num)
    max_var_num <- ifelse(is.null(max_var_num), ncol(train_num), max_var_num)

    min_var_num <- as.integer(min_var_num)
    max_var_num <- as.integer(max_var_num)

    if (min_var_num < 1) stop("min_var_num must be a positive integer!")
    if (max_var_num > ncol(train_num)) stop("max_var_num cannot be greater than variable number!")

    stop_cl <- FALSE
    if (is.null(cl)) {
        cl <- parallel::makeCluster(parallel.cores, outfile = output_log)
        stop_cl <- TRUE
    }

    cor_info <- Internal.compute_cor(train_num = train_num, test_num = test_num,
                                     cor_type = cor_type, cor_method = cor_method)

    selected_var <- Internal.select_variable(cor_info = cor_info, min_var_num = min_var_num,
                                             max_var_num = max_var_num, cl = cl)

    if (stop_cl) {
        parallel::stopCluster(cl)
    }
    selected_var
}

Internal.compute_targetVal <- function(input_col, compute_method = c("mean", "median"),
                                       remove_outlier = TRUE) {
    if (remove_outlier) outlier_res <- Internal.boxplot.stats(input_col, coef = 1.5)$out

    if (compute_method == "mean") {
        if (remove_outlier) input_col[input_col %in% outlier_res] <- mean(input_col[!input_col %in% outlier_res], na.rm = T)
        res <- mean(input_col, na.rm = T)
    } else {
        if (remove_outlier) input_col[input_col %in% outlier_res] <- median(input_col[!input_col %in% outlier_res], na.rm = T)
        res <- median(input_col, na.rm = T)
    }
    res
}

compute_targetVal <- function(QC_data, col_sampleID, col_batchID,
                              compute_method = c("mean", "median"),
                              batch_based = FALSE, remove_outlier = TRUE) {
    compute_method   <- match.arg(compute_method)

    sampleID <- as.factor(QC_data[[col_sampleID]])
    batchID  <- as.factor(QC_data[[col_batchID]])
    data_numeric <- QC_data[-c(col_sampleID, col_batchID)]
    if (!all(sapply(data_numeric, is.numeric))) stop("The values of the input dataset (except sampleID and batchID) should be numeric!")

    if (batch_based) {
        target_values <- aggregate(data_numeric, by = list(batch = batchID, sample = sampleID),
                                   FUN = Internal.compute_targetVal, compute_method = compute_method,
                                   remove_outlier = remove_outlier)
        batchID <- target_values$batch
        target_values$batch <- NULL
        target_values_list <- split(target_values, f = batchID)
    } else {
        target_values_list <- list(wholeDataset = aggregate(data_numeric,
                                                            by = list(sample = sampleID),
                                                            FUN =  Internal.compute_targetVal,
                                                            compute_method = compute_method,
                                                            remove_outlier = remove_outlier))
    }

    target_values <- lapply(target_values_list, function(x) {
        row.names(x) <- x$sample
        x$sample <- NULL
        x
    })
    target_values
}

tar_mean_batch <- compute_targetVal(QC_data = input_dataset_bak,
                                    col_sampleID = 1, col_batchID = 2,
                                    compute_method = "mean",
                                    batch_based = F, remove_outlier = T)

