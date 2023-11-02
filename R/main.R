#' Compute RSD (relative standard deviation)
#'
#' @description This function computes the RSD (relative standard deviation) of the values in \code{input_data}. Missing values are removed before the computation automatically.
#'
#' @param input_data a numeric vector
#'
#' @importFrom stats sd
#'
#' @details The RSD in this function is computed by:
#'
#' \code{sd(input_data, na.rm = TRUE) / mean(input_data, na.rm = TRUE)}.
#'
#' @return The RSD of the values in \code{input_data} is computed, as a numeric of length one.
#'
#' @examples
#' RSD_1 <- compute_RSD(c(1:10))
#'
#' data(FF4_qc) # load demo dataset
#'
#' # RSD of QC:
#' RSD_2 <- sapply(FF4_qc[FF4_qc$sampleType == "QC", -c(1:5)], compute_RSD)
#' quantile(RSD_2)
#'
#' # RSD of different types of QC samples:
#' # (each metabolote has its own RSD)
#' RSD_3 <- aggregate(FF4_qc[-c(1:5)], by = list(Type = FF4_qc$sampleType),
#'                    FUN = compute_RSD)
#' @export

compute_RSD <- function(input_data) {
    val_RSD <- sd(input_data, na.rm = TRUE) / mean(input_data, na.rm = TRUE)
    val_RSD
}

#' Compute target values for ensemble learning architecture
#'
#' @description This function provides an advanced option to calculate the target values of one reference dataset (i.e. \code{QC_num}, numeric values of quality control samples). The generated target values (a list) can be further passed to argument \code{targetVal_external} in function \code{\link{run_TIGER}} such that TIGER can align the \code{test_samples} with the reference dataset. This is useful for longitudinal datasets correction and cross-kit adjustment. See case study section of our original paper for detailed explanation.
#'
#' @param QC_num a numeric data.frame including the metabolite values of quality control (QC) samples. Missing values and infinite values will not be taken into account. Row: sample. Column: metabolite variable. See Examples.
#' @param sampleType a vector corresponding to \code{QC_num} to specify the type of each QC sample. QC samples of the \strong{same type} should have the \strong{same type name}. See Examples.
#' @param batchID a vector corresponding to \code{QC_num} to specify the batch of each sample. Ignored if \code{targetVal_batchWise = FALSE}. See Examples.
#' @param targetVal_method a character string specifying how the target values are computed. Can be \code{"mean"} (default) or \code{"median"}. See Details.
#' @param targetVal_batchWise logical. If \code{TRUE}, the target values will be computed based on each batch, otherwise, based on the whole dataset. Setting \code{TRUE} might be useful if your dataset has very obvious batch effects, but this may also make the algorithm less robust. See Details. Default: \code{FALSE}.
#' @param targetVal_removeOutlier logical. If \code{TRUE}, outliers will be removed before the computation. Outliers are determined with 1.5 * IQR (interquartile range) rule. We recommend turning this off when the target values are computed based on batches. See Details. Default: \code{!targetVal_batchWise}.
#' @param coerce_numeric logical. If \code{TRUE}, values in \code{QC_num} will be coerced to numeric before the computation. The columns cannot be coerced will be removed (with warnings). See Examples. Default: \code{FALSE}.
#'
#' @importFrom stats fivenum
#' @importFrom stats aggregate
#' @importFrom stats median
#'
#' @details See \code{\link{run_TIGER}}.
#'
#' @return
#' If \code{targetVal_batchWise = FALSE}, the function returns a list of length one containing the target values computed on the whole dataset.
#'
#' If \code{targetVal_batchWise = TRUE}, a list containing the target values computed on different batches is returned. The length of the returned list equals the number of batch specified by \code{batchID}.
#'
#' @examples
#' data(FF4_qc) # load demo dataset
#' QC_num <- FF4_qc[-c(1:5)] # only contain numeric metabolite values.
#'
#' # target values computed on the whole dataset:
#' tarVal_1 <- compute_targetVal(QC_num = QC_num,
#'                               sampleType = FF4_qc$sampleType,
#'                               batchID = FF4_qc$plateID,
#'                               targetVal_method = "mean",
#'                               targetVal_batchWise = FALSE,
#'                               targetVal_removeOutlier = TRUE)
#'
#' # target values computed on batches:
#' tarVal_2 <- compute_targetVal(QC_num = QC_num,
#'                               sampleType = FF4_qc$sampleType,
#'                               batchID = FF4_qc$plateID,
#'                               targetVal_method = "mean",
#'                               targetVal_batchWise = TRUE,
#'                               targetVal_removeOutlier = FALSE)
#'
#' # If coerce_numeric = TRUE,
#' # columns cannot be coerced to numeric will be removed (with warnings):
#' tarVal_3 <- compute_targetVal(QC_num = FF4_qc[-c(4:5)],
#'                               sampleType = FF4_qc$sampleType,
#'                               batchID = FF4_qc$plateID,
#'                               targetVal_method = "mean",
#'                               targetVal_batchWise = TRUE,
#'                               targetVal_removeOutlier = FALSE,
#'                               coerce_numeric = TRUE)
#' identical(tarVal_2, tarVal_3)  # identical to tarVal_2
#'
#' \dontrun{
#'
#' # will throw errors if input data have non-numeric columns
#' # and coerce_numeric = FALSE:
#'
#' tarVal_4 <- compute_targetVal(QC_num = FF4_qc,
#'                               sampleType = FF4_qc$sampleType,
#'                               batchID = FF4_qc$plateID,
#'                               targetVal_method = "mean",
#'                               targetVal_batchWise = TRUE,
#'                               targetVal_removeOutlier = FALSE,
#'                               coerce_numeric = FALSE)
#' }
#' @export

compute_targetVal <- function(QC_num, sampleType, batchID = NULL,
                              targetVal_method = c("mean", "median"),
                              targetVal_batchWise = FALSE,
                              targetVal_removeOutlier = !targetVal_batchWise,
                              coerce_numeric = FALSE) {

    message("+ Computing target values...   ", Sys.time())

    targetVal_method <- match.arg(targetVal_method)

    sampleType <- as.factor(sampleType)
    batchID    <- as.factor(batchID)

    if(coerce_numeric) {
        QC_num <- as.data.frame(sapply(QC_num, as.numeric))
        idx_NA <- sapply(QC_num, function(x) {
            all(is.na(x))
        })

        if (sum(idx_NA) > 0) {
            QC_num <- QC_num[,!idx_NA]
            warning("  ", sum(idx_NA), " column(s) in QC_num removed due to non-numeric values." )
        }
    } else {
        if (!all(sapply(QC_num, is.numeric))) stop("  The values in QC_num should be numeric!")
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
#'
#' @description This function provides an advanced option to select metabolite variables from external dataset(s). The selected variables (as a list) can be further passed to argument \code{selectVar_external} in function \code{\link{run_TIGER}} for a customised data correction.
#'
#' @param train_num a numeric data.frame \strong{only} including the metabolite values of training samples (can be quality control samples). Information such as injection order or well position need to be excluded. Row: sample. Column: metabolite variable. See Examples.
#' @param test_num an optional numeric data.frame including the metabolite values of test samples (can be subject samples). If provided, the column names of \code{test_num} should correspond to the column names of \code{train_num}. Row: sample. Column: metabolite variable. If \code{NULL}, the variables will be selected based on \code{train_num} only. See Examples.
#' @param train_batchID \code{NULL} or a vector corresponding to \code{train_num} to specify the batch of each sample. Ignored if \code{selectVar_batchWise = FALSE}. See Examples.
#' @param test_batchID \code{NULL} or a vector corresponding to \code{test_num} to specify the batch of each sample. Ignored if \code{selectVar_batchWise = FALSE}. See Examples.
#' @param selectVar_corType a character string indicating correlation (\code{"cor"}, default) or partial correlation (\code{"pcor"}) is to be used. Can be abbreviated. See Details. \strong{Note}: computing partial correlations of a large dataset can be very time-consuming.
#' @param selectVar_corMethod a character string indicating which correlation coefficient is to be computed. One of \code{"spearman"} (default) or \code{"pearson"}. Can be abbreviated. See Details.
#' @param selectVar_corUse Indicate how to deal with missing values when computing correlation (\code{selectVar_corType = "cor"}). The value will be passed to the parameter \code{use} from \code{\link[stats]{cor}}. Support values: \code{"complete.obs"} (default), \code{"everything"}, \code{"all.obs"}, \code{"na.or.complete"}, or \code{"pairwise.complete.obs"}. See \code{\link[stats]{cor}} for details.
#' @param selectVar_minNum an integer specifying the minimum number of the selected variables. If \code{NULL}, no limited, but 1 at least. See Details. Default: 5.
#' @param selectVar_maxNum an integer specifying the maximum number of the selected variables. If \code{NULL}, no limited, but \code{ncol(train_num) - 1} at most. See Details. Default: 10.
#' @param selectVar_batchWise (advanced) logical. Specify whether the variable selection should be performed based on each batch. Default: \code{FALSE}. \strong{Note}: if \code{TRUE}, batch ID of each sample are required. The support of batch-wise variable selection is provided for data requiring special processing (for example, data with strong batch effects). But in most case, batch-wise variable selection is not necessary. Setting \code{TRUE} might make the algorithm less robust. See Details.
#' @param coerce_numeric logical. If \code{TRUE}, values in \code{train_num} and  \code{test_num} will be coerced to numeric before the computation. The columns cannot be coerced will be removed (with warnings). See Examples. Default: \code{FALSE}.
#'
#' @importFrom pbapply timerProgressBar
#' @importFrom pbapply setTimerProgressBar
#' @importFrom pbapply closepb
#' @importFrom ppcor pcor
#' @importFrom stats cor
#'
#' @details See \code{\link{run_TIGER}}.
#'
#' @return
#' If \code{selectVar_batchWise = FALSE}, the function returns a list of length one containing the selected variables computed on the whole dataset.
#'
#' If \code{selectVar_batchWise = TRUE}, a list containing the selected variables computed on different batches is returned. The length of the returned list equals the number of batch specified by \code{test_batchID} and/or \code{train_batchID}.
#'
#' @examples
#'
#' data(FF4_qc) # load demo dataset
#'
#' # QC as training samples; QC1, QC2 and QC3 as test samples:
#' train_samples <- FF4_qc[FF4_qc$sampleType == "QC",]
#' test_samples  <- FF4_qc[FF4_qc$sampleType != "QC",]
#'
#' # Only numeric data of metabolite variables are allowed:
#' train_num = train_samples[-c(1:5)]
#' test_num  = test_samples[-c(1:5)]
#'
#' # If the selection is performed on the whole dataset:
#' # based on training samples only:
#' selected_var_1 <- select_variable(train_num = train_num,
#'                                   test_num  = NULL,
#'                                   selectVar_batchWise = FALSE)
#'
#' # also consider test samples:
#' selected_var_2 <- select_variable(train_num = train_num,
#'                                   test_num  = test_num,
#'                                   selectVar_batchWise = FALSE)
#'
#' # If the selection is based on different batches:
#' # (In selectVar_batchWise, batch ID is required.)
#' selected_var_3 <- select_variable(train_num = train_num,
#'                                   test_num  = NULL,
#'                                   train_batchID = train_samples$plateID,
#'                                   test_batchID  = NULL,
#'                                   selectVar_batchWise = TRUE)
#'
#' # If coerce_numeric = TRUE,
#' # columns cannot be coerced to numeric will be removed (with warnings):
#' # (In this example, columns of injection order and well position are excluded.
#' # Because we don't want to calculate the correlations between metabolites and
#' # injection order/well position.)
#' selected_var_4 <- select_variable(train_num = train_samples[-c(4,5)],
#'                                   train_batchID = train_samples$plateID,
#'                                   selectVar_batchWise = TRUE,
#'                                   coerce_numeric = TRUE)
#' identical(selected_var_3, selected_var_4)  # identical to selected_var_3
#'
#' \dontrun{
#'
#' # will throw errors if input data have non-numeric columns
#' # and coerce_numeric = FALSE:
#'
#' selected_var_5 <- select_variable(train_num = train_samples[-c(4,5)],
#'                                   coerce_numeric = FALSE)
#' }
#' @export

select_variable <- function(train_num, test_num = NULL,
                            train_batchID = NULL, test_batchID = NULL,
                            selectVar_corType   = c("cor", "pcor"),
                            selectVar_corMethod = c("spearman", "pearson"),
                            selectVar_corUse = "complete.obs",
                            selectVar_minNum = 5, selectVar_maxNum = 10,
                            selectVar_batchWise = FALSE,
                            coerce_numeric = FALSE) {

    message("+ Selecting highly-correlated variables...   ", Sys.time())

    selectVar_corType   <- match.arg(selectVar_corType)
    selectVar_corMethod <- match.arg(selectVar_corMethod)
    selectVar_corUse <- match.arg(selectVar_corUse, choices = c("complete.obs", "everything", "all.obs",
                                                                "na.or.complete", "pairwise.complete.obs"))

    if(coerce_numeric) {
        train_num <- as.data.frame(sapply(train_num, as.numeric))
        idx_NA <- sapply(train_num, function(x) {
            all(is.na(x))
        })

        if (sum(idx_NA) > 0) {
            train_num <- train_num[,!idx_NA]
            warning("  ", sum(idx_NA), " column(s) in train_num removed due to non-numeric values." )
        }

        if (!is.null(test_num)) {
            test_num <- as.data.frame(sapply(test_num, as.numeric))
            idx_NA <- sapply(test_num, function(x) {
                all(is.na(x))
            })
            if (sum(idx_NA) > 0) {
                test_num <- test_num[,!idx_NA]
                warning("  ", sum(idx_NA), " column(s) in test_num removed due to non-numeric values." )
            }
        }
    } else {
        if (!all(sapply(train_num, is.numeric))) stop("  The values in train_num should be numeric!")
        if (!is.null(test_num) & !all(sapply(test_num, is.numeric))) stop("  The values in test_num should be numeric!")
    }

    if (!is.null(test_num)) {
        if (!all(names(test_num) %in% names(train_num))) stop("  Variables in training and test data cannot match!")
        train_num <- train_num[names(test_num)]
    }

    selectVar_minNum <- ifelse(is.null(selectVar_minNum), 1, selectVar_minNum)
    selectVar_maxNum <- ifelse(is.null(selectVar_maxNum), (ncol(train_num) - 1), selectVar_maxNum)

    selectVar_minNum <- max(as.integer(selectVar_minNum), 1)
    selectVar_maxNum <- max(as.integer(selectVar_maxNum), 1)

    selectVar_maxNum <- min(selectVar_maxNum, ncol(train_num) - 1)
    selectVar_minNum <- min(selectVar_minNum, selectVar_maxNum)

    if (selectVar_batchWise) {
        train_num_list    <- split(train_num, f = train_batchID)
        batch_names_train <- sort(names(train_num_list))
        batch_names_len   <- length(batch_names_train)

        if (!is.null(test_num)) {
            test_num_list     <- split(test_num,  f = test_batchID)
            batch_names_test  <- sort(names(test_num_list))

            if (any(batch_names_train != batch_names_test) | batch_names_len != length(batch_names_test)) stop("  Batch names of train and test samples cannot match: their column names should be identical!")
        } else test_num_list <- NULL

        selected_var_list <- lapply(1:batch_names_len, function(batch_name_idx) {

            one_batch_name <- batch_names_train[[batch_name_idx]]
            message("  - Current batch: ", one_batch_name, "  (", batch_name_idx, "/", batch_names_len, ")")

            message("      Computing correlation coefficients...")
            cor_info <- Internal.compute_cor(train_num = train_num_list[[one_batch_name]],
                                             test_num = test_num_list[[one_batch_name]],
                                             selectVar_corType = selectVar_corType,
                                             selectVar_corMethod = selectVar_corMethod,
                                             selectVar_corUse = selectVar_corUse)

            message("      Selecting variables...")
            selected_var <- Internal.select_variable(cor_info = cor_info,
                                                     selectVar_minNum = selectVar_minNum,
                                                     selectVar_maxNum = selectVar_maxNum)
        })
        names(selected_var_list) <- batch_names_train
    } else {
        message("  - Computing correlation coefficients...")
        cor_info <- Internal.compute_cor(train_num = train_num, test_num = test_num,
                                         selectVar_corType = selectVar_corType,
                                         selectVar_corMethod = selectVar_corMethod)
        message("  - Selecting variables...")
        selected_var <- Internal.select_variable(cor_info = cor_info,
                                                 selectVar_minNum = selectVar_minNum,
                                                 selectVar_maxNum = selectVar_maxNum)
        selected_var_list <- list(wholeDataset = selected_var)
    }

    selected_var_list
}

#' Run TIGER to eliminate technical variation
#' @description Use TIGER algorithm to eliminate the technical variation in metabolomics data. TIGER supports targeted and untargeted metabolomics data and is competent to perform both intra- and inter-batch technical variation removal.
#'
#' @param test_samples (required) a data.frame containing the samples to be corrected (for example, subject samples). This data.frame should contain columns of
#' \itemize{
#' \item sample ID (required): name or label for each sample,
#' \item sample type (required): indicating the type of each sample,
#' \item batch ID (required): the batch of each sample,
#' \item order information (optional): injection order or temporal information of each sample,
#' \item position information (optional): well position of each sample,
#' \item metabolite values (required): values to be normalised. Infinite values are not allowed.
#' }
#' Row: sample. Column: variable. See Examples.
#' @param train_samples (required) a data.frame containing the quality control (QC) samples used for model training. The columns in this data.frame should correspond to the columns in \code{test_samples}. And \code{test_samples} and \code{train_samples} should have the identical column names.
#' @param col_sampleID  (required) a character string indicating the name of the column that specifies the sample ID of each sample. The values in this column will not affect the data correction process but can act as labels for different samples. See Examples.
#' @param col_sampleType (required) a character string indicating the name of the column that specifies the type (such as QC1, QC2, subject) of each sample. This column can be used to indicate different kinds of QC samples in \code{train_samples}. QC samples of the \strong{same type} should have the \strong{same type name}. See Examples.
#' @param col_batchID (required) a character string indicating the name of the column that specifies the batch ID of each sample. See Examples.
#' @param col_order (optional) \code{NULL} or a character string indicating the name of the column that contains the injection order or temporal information (numeric values). This can explicitly ask the algorithm capture the technical variation introduced by injection order, which might be useful when your data have very obvious temporal drifts. If \code{NULL} (default), \code{train_samples} and \code{test_samples} should have \strong{No} column contains injection order information.
#' @param col_position (optional) \code{NULL} or a character string indicating the name of the column that contains the well position information (numeric values). This can explicitly ask the algorithm capture the technical variation introduced by well position, which might be useful when the well position has a great impact during data acquisition. If \code{NULL} (default), \code{train_samples} and \code{test_samples} should have \strong{No} column contains well position information.
#' @param targetVal_external (optional) a list generated by function \code{\link{compute_targetVal}}. See Details.
#' @param targetVal_method a character string specifying how target values are to be computed. Can be \code{"mean"} (default) or \code{"median"}. Ignored if a list of external target values has been assigned to  \code{targetVal_external}.
#' @param targetVal_batchWise logical. If \code{TRUE}, the target values will be computed based on each batch, otherwise, based on the whole dataset. Setting \code{TRUE} might be useful if your dataset has very obvious batch effects, but this may also make the algorithm less robust. Default: \code{FALSE}. Ignored if a list of external target values has been assigned to  \code{targetVal_external}.
#' @param targetVal_removeOutlier logical. If \code{TRUE}, outliers will be removed before the computation. Outliers are determined with 1.5 * IQR (interquartile range) rule. We recommend turning this off when the target values are computed based on batches. Default: \code{!targetVal_batchWise}. Ignored if a list of external target values has been assigned to  \code{targetVal_external}.
#' @param selectVar_external (optional) a list generated by function \code{\link{select_variable}}. See Details.
#' @param selectVar_corType a character string indicating correlation (\code{"cor"}, default) or partial correlation (\code{"pcor"}) is to be used. Can be abbreviated. Ignored if a list of selected variables has been assigned to \code{selectVar_external}. \strong{Note}: computing partial correlations of a large dataset can be very time-consuming.
#' @param selectVar_corMethod a character string indicating which correlation coefficient is to be computed. One of \code{"spearman"} (default) or \code{"pearson"}. Can be abbreviated. Ignored if a list of selected variables has been assigned to \code{selectVar_external}.
#' @param selectVar_corUse Indicate how to deal with missing values when computing correlation (\code{selectVar_corType = "cor"}). The value will be passed to the parameter \code{use} from \code{\link[stats]{cor}}. Support values: \code{"complete.obs"} (default), \code{"everything"}, \code{"all.obs"}, \code{"na.or.complete"}, or \code{"pairwise.complete.obs"}. See \code{\link[stats]{cor}} for details.
#' @param selectVar_minNum an integer specifying the minimum number of the selected metabolite variables (injection order and well position are not regarded as metabolite variables). If \code{NULL}, no limited, but 1 at least. Default: \code{5}. Ignored if a list of selected variables has been assigned to \code{selectVar_external}.
#' @param selectVar_maxNum an integer specifying the maximum number of the selected metabolite variables (injection order and well position are not regarded as metabolite variables). If \code{NULL}, no limited, but no more than the number of all available metabolite variables. Default: \code{10}. Ignored if a list of selected variables has been assigned to \code{selectVar_external}.
#' @param selectVar_batchWise (advanced) logical. Specify whether the variable selection should be performed based on each batch. Default: \code{FALSE}. Ignored if a list of selected variables has been assigned to \code{selectVar_external}. \strong{Note}: the support of batch-wise variable selection is provided for data requiring special processing (for example, data with strong batch effects). But in most case, batch-wise variable selection is not necessary. Setting \code{TRUE} can make the algorithm less robust.
#' @param mtry_percent (advanced) a numeric vector indicating the percentages of selected variables randomly sampled as candidates at each split when training random forest models (base learners). \strong{Note}: providing more arguments will include more base learners into the ensemble model, which will increase the processing time. Default: \code{seq(0.2, 0.8, 0.2)}.
#' @param nodesize_percent (advanced) a numeric vector indicating the percentages of sample size used as the minimum sizes of the terminal nodes in random forest models (base learners). \strong{Note}: providing more arguments will include more base learners into the ensemble model, which will increase the processing time. Default: \code{seq(0.2, 0.8, 0.2)}.
#' @param ... (advanced) optional arguments (except \code{mtry} and \code{nodesize}) to be passed to \code{\link[randomForest]{randomForest}} for model training. Arguments \code{mtry} and \code{nodesize} are determined by \code{mtry_percent} and \code{nodesize_percent}. See \code{\link[randomForest]{randomForest}} and Examples. \strong{Note}: providing more arguments will include more base learners into the ensemble model, which will increase the processing time.
#' @param parallel.cores an integer (== -1 or >= 1) specifying the number of cores for parallel computation. Setting \code{-1} to run with all cores. Default: \code{2}.
#'
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
#' @importFrom stats aggregate
#' @importFrom stats as.formula
#' @importFrom stats cor
#' @importFrom stats fivenum
#' @importFrom stats median
#' @importFrom stats predict
#' @importFrom stats quantile
#' @importFrom stats rnorm
#' @importFrom stats sd
#'
#' @details
#' TIGER can effectively process the datasets with its default setup. The following hyperparameters are provided to customise the algorithm and achieve the best possible performance. These hyperparameters are also practical for some special purposes (such as cross-kit adjustment, longitudinal dataset correction) or datasets requiring special processing (for example, data with very strong temporal drifts or batch effects). We recommend users to examine the normalised result with different metrics, such as RSD (relative standard deviation), MAPE (mean absolute percentage error) and PCA (principal component analysis), especially when using advanced options of TIGER.
#'
#' \strong{Hyperparameters for target value computation}
#'
#' \itemize{
#' \item \code{targetVal_external}
#'
#' TIGER by default captures and eliminate the technical variation within the input dataset, and the target values are automatically computed from \code{train_samples}. The target values can also be calculated from a reference dataset using function \code{\link{compute_targetVal}} and then passed to this function as an argument. This will enable TIGER to align \code{test_samples} with the reference dataset. In this case, \code{train_samples} is still the accompanying QC samples of \code{test_samples}. And argument \code{targetVal_external} accepts external target values (a list). If the list of external target values is provided, values in \code{targetVal_method}, \code{targetVal_batchWise} and \code{targetVal_removeOutlier} will be ignored.
#'
#' \item \code{targetVal_method}
#'
#' The target values can be the mean or median values of metabolite values. The target values of different kinds of QC samples are computed separately. \code{"mean"} is recommended here, but the optimal selection can differ for different datasets.
#'
#' \item \code{targetVal_batchWise}
#'
#' The target values can be computed from the whole dataset or from different batches. By default, the target values are computed based on the whole dataset. Computing based on batches (\code{targetVal_batchWise = TRUE}) is only recommended when the samples has very strong batch effects. For example, we set this as \code{TRUE} when normalising WaveICA's Amide dataset in our original paper.
#'
#' \item  \code{targetVal_removeOutlier}
#'
#' If computing is based on the whole dataset (\code{targetVal_batchWise = TRUE}), users can remove the outliers in each metabolite by setting \code{targetVal_removeOutlier} as \code{TRUE}. This can weaken the impact of extreme values. If \code{targetVal_batchWise = FALSE}, it is generally not recommended to remove outliers, as we assume the input data have strong batch effects and contain extreme values—we hope TIGER can take these into account. Code for checking outliers is adapted from \code{\link[grDevices]{boxplot.stats}}.
#' }
#'
#' \strong{Hyperparameters for variable selection}
#'
#' \itemize{
#' \item \code{selectVar_external}:
#'
#' This argument accepts a list of selected variables generated by \code{\link{select_variable}}. This is helpful when you want to use the same selected variables to correct several datasets. You can also pass a self-defined list to this argument, as long as the self-defined list has similar data structure as the one generated by \code{\link{select_variable}}.
#'
#' \item \code{selectVar_corType} and \code{selectVar_corMethod}:
#'
#' TIGER supports Pearson product-moment correlation (\code{"pearson"}) and Spearman's rank correlation (\code{"spearman"}) to compute correlation coefficients (\code{"cor"}) or partial correlation coefficients (\code{"por"}) for variable selection. See \code{\link[stats]{cor}} and \code{\link[ppcor]{pcor}} for further details.
#'
#' \item \code{selectVar_minNum} and \code{selectVar_maxNum}:
#'
#' For an objective metabolite to be corrected, the intersection of its top \emph{t} highly-correlated metabolites calculated from training and test samples are selected to train the ensemble model. The highly-correlated metabolites are the ones with correlation coefficients greater than 0.5 (the objective metabolite itself will not be regarded as its highly-correlated metabolite). Arguments \code{selectVar_minNum} and \code{selectVar_maxNum} are used to avoid selecting too many or too few metabolites. Selecting too many metabolites can lower the process, sometimes even lower the accuracy.
#'
#' \item \code{selectVar_batchWise}:
#'
#' Advanced option designed for special cases. Setting it \code{TRUE} might be useful when your data have very obvious batch effects.
#' }
#'
#' \strong{Hyperparameters for model construction}
#'
#' \itemize{
#' \item \code{mtry_percent}, \code{nodesize_percent} and \code{...}:
#'
#' Advanced options to specify \code{mtry}, \code{nodesize} and other related arguments in \code{\link[randomForest]{randomForest}} for a customised ensemble learning architecture. See Examples.
#' }
#'
#' @return This function returns a data.frame with the same data structure as the input \code{test_samples}, but the metabolite values are the normalised/corrected ones. \code{NA} and zeros in the original \code{test_samples} will not be changed or normalised.
#'
#' @section Reference:
#' Han S. \emph{et al}. TIGER: technical variation elimination for metabolomics data using ensemble learning architecture. \emph{Briefings in Bioinformatics} (2022) bbab535. \doi{10.1093/bib/bbab535}.
#'
#' @examples
#' \donttest{
#' data(FF4_qc) # load demo dataset
#'
#' # QC as training samples; QC1, QC2 and QC3 as test samples:
#' train_samples <- FF4_qc[FF4_qc$sampleType == "QC",]
#' test_samples  <- FF4_qc[FF4_qc$sampleType != "QC",]
#'
#' # col_sampleID includes labels. You can assign names for different samples:
#' train_samples$sampleID <- "train"
#' test_samples$sampleID  <- "test"
#'
#' # Use default setting and
#' # include injection order and well position into feature set:
#' test_norm_1 <- run_TIGER(test_samples = test_samples,
#'                          train_samples = train_samples,
#'                          col_sampleID  = "sampleID",     # input column name
#'                          col_sampleType = "sampleType",  # input column name
#'                          col_batchID = "plateID",        # input column name
#'                          col_order = "injectionOrder",   # input column name
#'                          col_position = "wellPosition",  # input column name
#'                          parallel.cores = 2)
#'
#' # If the information of injection order and well position is not available,
#' # or you don't want to use them:
#' train_data <- train_samples[-c(4:5)]  # remove the two columns
#' test_data  <- test_samples[-c(4:5)]   # remove the two columns
#'
#' test_norm_2 <- run_TIGER(test_samples = test_data,
#'                          train_samples = train_data,
#'                          col_sampleID  = "sampleID",
#'                          col_sampleType = "sampleType",
#'                          col_batchID = "plateID",
#'                          col_order = NULL,                # set NULL
#'                          col_position = NULL,             # set NULL
#'                          parallel.cores = 2)
#'
#' # If use external target values and selected variables with
#' # customised settings:
#' target_val <- compute_targetVal(QC_num = train_samples[-c(1:5)],
#'                                 sampleType = train_samples$sampleType,
#'                                 batchID = train_samples$plateID,
#'                                 targetVal_method = "median",
#'                                 targetVal_batchWise = TRUE)
#'
#' select_var <- select_variable(train_num = train_samples[-c(1:5)],
#'                               test_num = test_samples[-c(1:5)],
#'                               train_batchID = train_samples$plateID,
#'                               test_batchID = test_samples$plateID,
#'                               selectVar_corType = "pcor",
#'                               selectVar_corMethod = "spearman",
#'                               selectVar_minNum = 10,
#'                               selectVar_maxNum = 30,
#'                               selectVar_batchWise = TRUE)
#'
#' test_norm_3 <- run_TIGER(test_samples = test_samples,
#'                          train_samples = train_samples,
#'                          col_sampleID  = "sampleID",
#'                          col_sampleType = "sampleType",
#'                          col_batchID = "plateID",
#'                          col_order = "injectionOrder",
#'                          col_position = "wellPosition",
#'                          targetVal_external = target_val,
#'                          selectVar_external = select_var,
#'                          parallel.cores = 2)
#'
#' # The definitions of other hyperparameters correspond to
#' # randomForest::randomForest().
#' # If want to include more hyperparameters into model training,
#' # put hyperparameter values like this:
#' mtry_percent <- c(0.4, 0.8)
#' nodesize_percent <- c(0.4, 0.8)
#' replace <- c(TRUE, FALSE)
#' ntree <- c(100, 200, 300)
#'
#' test_norm_4 <- run_TIGER(test_samples = test_data,
#'                          train_samples = train_data,
#'                          col_sampleID  = "sampleID",
#'                          col_sampleType = "sampleType",
#'                          col_batchID = "plateID",
#'                          mtry_percent = mtry_percent,
#'                          nodesize_percent = nodesize_percent,
#'                          replace = replace,
#'                          ntree = ntree,
#'                          parallel.cores = 2)
#'
#' # test_norm_4 is corrected by the ensemble model consisted of base learners
#' # trained with (around) 24 different hyperparameter combinations:
#' expand.grid(mtry_percent, nodesize_percent, replace, ntree)
#'
#' # Note: mtry and nodesize are calculated by mtry_percent and nodesize_percent,
#' #       duplicated hyperparameter combinations, if any, will be removed.
#' #       Thus, the total number of hyperparameter combinations can be less than 24.
#' #       This is determined by the shape of your input datasets.
#' }
#' @export
#'
run_TIGER <- function(test_samples, train_samples,
                      col_sampleID, col_sampleType, col_batchID,
                      col_order = NULL, col_position = NULL,

                      targetVal_external = NULL, targetVal_method = c("mean", "median"),
                      targetVal_batchWise = FALSE, targetVal_removeOutlier = !targetVal_batchWise,

                      selectVar_external = NULL, selectVar_corType = c("cor", "pcor"),
                      selectVar_corMethod = c("pearson", "spearman"),
                      selectVar_corUse = "complete.obs",
                      selectVar_minNum = 5, selectVar_maxNum = 10,
                      selectVar_batchWise = FALSE,

                      mtry_percent = seq(0.2, 0.8, 0.2),
                      nodesize_percent = seq(0.2, 0.8, 0.2),
                      ..., parallel.cores = 2) {

    message("+ Initialising...   ", Sys.time())

    selectVar_corType   <- match.arg(selectVar_corType)
    selectVar_corMethod <- match.arg(selectVar_corMethod)
    targetVal_method    <- match.arg(targetVal_method)

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
        if (length(targetVal_list) == 1 & names(targetVal_list[1]) == "wholeDataset") {
            targetVal_batchWise <- FALSE
        } else {
            if (all(names(targetVal_list) %in% batchID_train)) {
                targetVal_batchWise <- TRUE
            } else stop("  Batch ID of targetVal_external cannot match your training data!")
        }
    }

    # Variable selection
    if (is.null(selectVar_external)) {
        var_selected_list <- select_variable(train_num = train_num, test_num = test_num,
                                             train_batchID = train_samples[[col_batchID]],
                                             test_batchID  = test_samples[[col_batchID]],
                                             selectVar_batchWise = selectVar_batchWise,
                                             selectVar_corType   = selectVar_corType,
                                             selectVar_corMethod = selectVar_corMethod,
                                             selectVar_corUse = selectVar_corUse,
                                             selectVar_minNum = selectVar_minNum,
                                             selectVar_maxNum = selectVar_maxNum,
                                             coerce_numeric = FALSE)
    } else {
        message("+ External selected variables loaded.   ", Sys.time())
        var_selected_list <- selectVar_external
        if (length(var_selected_list) == 1 & names(var_selected_list[1]) == "wholeDataset") {
            selectVar_batchWise <- FALSE
        } else {
            if (all(names(var_selected_list) %in% batchID_test)) {
                selectVar_batchWise <- TRUE
            } else stop("  Batch ID of selectVar_external cannot match your test data!")
        }
    }

    idx_test_na  <- is.na(test_samples)
    idx_train_na <- is.na(train_samples)
    if (any(idx_test_na))  test_samples[idx_test_na]   <- 0
    if (any(idx_train_na)) train_samples[idx_train_na] <- 0

    message("+ Data correction started.   ", Sys.time())
    parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
    cl <- parallel::makeCluster(parallel.cores)
    parallel::clusterExport(cl = cl, varlist = c("Internal.compute_errorRatio", "Internal.run_ensemble"), envir = environment())
    pbapply::pboptions(type = "timer", style = 3, char = "=", txt.width = 70)

    # Original sample order backup
    test_samples <- cbind(original_idx = 1:nrow(test_samples), test_samples)

    res_var <- pbapply::pblapply(var_names, function(current_var, var_selected_list, targetVal_list,
                                                     targetVal_batchWise, selectVar_batchWise,
                                                     train_samples, test_samples, col_sampleID, col_sampleType,
                                                     col_batchID, col_order, col_position, batchID, mtry_percent,
                                                     targetVal_method, nodesize_percent, ...) {
        if (!targetVal_batchWise) {
            train_y_all <- Internal.compute_errorRatio(rawVal     = train_samples[[current_var]],
                                                       sampleType = train_samples[[col_sampleType]],
                                                       targetVal  = targetVal_list$wholeDataset[current_var])

        }

        test_data <- cbind(y_raw = test_samples[[current_var]], test_samples)

        res_batch_list <- lapply(batchID, function(current_batch) {

            train_X_batch <- train_samples[train_samples[[col_batchID]] == current_batch,]

            if (targetVal_batchWise) {
                train_y_batch <- Internal.compute_errorRatio(rawVal     = train_X_batch[[current_var]],
                                                             sampleType = train_X_batch[[col_sampleType]],
                                                             targetVal  = targetVal_list[[current_batch]][current_var])


            } else {
                train_y_batch <- train_y_all[train_samples[[col_batchID]] == current_batch,]
            }

            if (selectVar_batchWise) {
                train_X_selected <- train_X_batch[ c(col_order, col_position,
                                                     var_selected_list[[current_batch]][[current_var]]) ]
            } else {
                train_X_selected <- train_X_batch[ c(col_order, col_position,
                                                     var_selected_list$wholeDataset[[current_var]]) ]
            }

            trainSet <- cbind(train_y_batch, train_X_selected)
            testSet  <- test_data[test_data[[col_batchID]] == current_batch,]

            var_pred <- Internal.run_ensemble(trainSet = trainSet, testSet = testSet,
                                              mtry_percent = mtry_percent,
                                              nodesize_percent = nodesize_percent,
                                              ... = ..., return_base_res = FALSE)

            if (targetVal_batchWise) {
                test_targetVal_all   <- do.call(targetVal_method, list(test_data$y_raw, na.rm = TRUE))
                test_targetVal_batch <- do.call(targetVal_method, list(testSet$y_raw,   na.rm = TRUE))
                var_pred <- var_pred * test_targetVal_all / test_targetVal_batch
            }

            names(var_pred) <- testSet$original_idx
            var_pred
        })

        res_batch       <- do.call("c", res_batch_list)
        res_batch_order <- res_batch[order(as.numeric(names(res_batch)))]
        res_batch_df    <- data.frame(res_batch_order)

        names(res_batch_df) <- current_var
        res_batch_df

    }, var_selected_list = var_selected_list, targetVal_list = targetVal_list,
    targetVal_batchWise = targetVal_batchWise, selectVar_batchWise = selectVar_batchWise,
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
