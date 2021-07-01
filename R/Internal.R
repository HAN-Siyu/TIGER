compute_RSD <- function(input_data) {
    cv <- sd(input_data, na.rm = T) / mean(input_data, na.rm = T)
    cv
}

.boxplot.stats <- function(x, coef = 1.5, do.conf = TRUE, do.out = TRUE) {
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

compute_average <- function(input_qc_data, qc_idx = 1, var_start_idx = 5) {
    set.seed(1)
    input_qc_data[[qc_idx]] <- as.character(input_qc_data[[qc_idx]])

    qc_data <- input_qc_data[,c(qc_idx, var_start_idx:ncol(input_qc_data))]

    names(qc_data)[1] <- "QC_ID"

    var_start_idx <- 2

    for(i in var_start_idx:ncol(qc_data)) {
        qc_data[!is.finite(qc_data[,i]), i] <- NA
    }

    ### start - calculate pseudo-true values
    qc_id <- qc_data$QC_ID
    qc_types <- as.character(unique(qc_id))
    qc_num <- qc_data[,var_start_idx:ncol(qc_data)]

    qc_meanVal <- sapply(qc_types, function(qc_type, qc_num, qc_id) {
        qc_type_num <- qc_num[qc_id %in% qc_type,]
        qc_remove_outlier <- sapply(qc_type_num, function(input_col) {
            outlier_res <- .boxplot.stats(input_col, coef = 1.5)$out
            input_col[input_col %in% outlier_res] <- mean(input_col[!input_col %in% outlier_res], na.rm = T)
            input_col
        })

        qc_value <- colMeans(qc_remove_outlier, na.rm = T)
    }, qc_num = qc_num, qc_id = qc_id)

    qc_reference <- data.frame(t(qc_meanVal))
}

remove_NA <- function(input_data_num, data_label = NULL) {
    data_na_sample_idx <- (apply(input_data_num, 1, function(x) any(is.na(x))))
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

compute_cor <- function(train_num, test_num = NULL,
                        cor_type = c("pcor", "cor"),
                        cor_method = c("pearson", "spearman")) {
    cor_type   <- match.arg(cor_type)
    cor_method <- match.arg(cor_method)

    if (!is.null(test_num)) {
        if (any(names(train_num) != names(test_num))) stop("Error: Variables in training and test data cannot match!")
    }

    train_num_noNA <- remove_NA(train_num, data_label = "training data")
    if (!is.null(test_num)) test_num_noNA <- remove_NA(test_num, data_label = "test data")

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

select_variable <- function(cor_info, min_var_num = NULL,
                            max_var_num = NULL, cl) {

    min_var_num <- ifelse(is.null(min_var_num), 1, min_var_num)
    max_var_num <- ifelse(is.null(max_var_num), length(cor_info$variable_name), max_var_num)

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

scale_train_test <- function(train_num, test_num = NULL) {
    train_mean <- colMeans(train_num, na.rm = T)
    train_sd <- apply(train_num, 2, sd, na.rm = T)

    train_scaled <- apply(train_num, 1, function(x) (x - train_mean) / train_sd)
    train_scaled <- data.frame(t(train_scaled))
    output <- list(train_scaled = train_scaled)

    if (!is.null(test_set)) {
        test_scaled <- apply(test_num, 1, function(x) (x - train_mean) / train_sd)
        test_scaled <- data.frame(t(test_scaled))
        output <- c(output, test_scaled = list(test_scaled))
    }
    output
}
