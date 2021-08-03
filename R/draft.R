# library(ranger)
# library(randomForest)

# refer
# load("/Volumes/Work/Projects/Helmholtz/Helmholtz_normaliation/Data/norm_RF/qc_class_refer_forEvalCV.RData")

# strat
# load("/Volumes/Work/Projects/Helmholtz/Helmholtz_normaliation/Data/norm_RF/qc_class_strat_forEvalCV.RData")

# qc_Refer
# load("/Volumes/Work/Projects/Helmholtz/Helmholtz_normaliation/Data/norm_RF/qc_class_qcRefer_forEvalCV.RData")

# SERRF demo dataset
# load("/Volumes/Work/Projects/Helmholtz/Helmholtz_normaliation/Data/norm_RF/SERRF_test/SERRF_demo_class.RData")

# Amide dataset
# load("/Volumes/Work/Projects/Helmholtz/Helmholtz_normaliation/Data/norm_Amide/Amide_class.RData")

source_file <- '/Volumes/Work/Projects/Helmholtz/Helmholtz_normaliation/scripts/metabolomics/norm_RF_Internal.R'
if(file.exists(source_file)) {
    source(source_file)
} else {
    # source("~/R_tmp/HGMU/norm_RF_Internal.R")
    source("norm_RF_Internal.R")
}

train_test_kNNreg <- function(input_qc_class = NULL, train_set = NULL, test_set = NULL,
                              reference_data = NULL, parallel.cores = 14, batch_wise = T,
                              min_candidate_len = 10, max_candidate_len = 30,
                              scale_mode = c("none", "together", "separate", "train"),
                              feature_selection_mode = c("all", "plate"),
                              cor_type = c("pcor", "cor"),
                              cor_method = c("pearson", "spearman"),
                              plate_based_refer = F, cv_mode = F,
                              # using_train_std = T,  scale_all = F,
                              checkOutlier_batchVal = F, checkOutlier_recoverNorm = F,
                              well_position = F, injection_order = T, algorithm = "kd_tree",
                              k_range = c(1:10, 15, 20, 30, 40, 50), k_batch = 10, sep_char = "_",
                              predict_test = T, output_pred_test = T, output_file_label = NULL,
                              rf_stack_mode = c("weight_recover_y", "mean", "weight_target_y", "linear"),
                              use_methods = c("TIGER", "kNN", "kNN_stack", "rf", "ranger", "xgboost"),
                              xgboost_params = NULL, rf_params = NULL, log_file = "log") {
    message("+ Processing...   ", Sys.time())
    library(parallel)
    library(pbapply)
    # browser()
    scale_all <- F
    using_train_std <- T
    std <- T
    scale_mode <- match.arg(scale_mode)
    if("none" %in% scale_mode) std <- F
    if("together" %in% scale_mode) scale_all <- T
    if("separate" %in% scale_mode) using_train_std <- F
    if("train" %in% scale_mode) using_train_std <- T

    select_all <- F
    feature_selection_mode <- match.arg(feature_selection_mode)
    if("all" %in% feature_selection_mode) {
        select_all <- T
    }
    message("  - Parameter:")
    message("    scale_mode: ", scale_mode)
    message("    feature_selection_mode: ", feature_selection_mode)
    message("    plate_based_refer: ", plate_based_refer)

    if(!is.null(input_qc_class)) {
        train_set <- do.call("rbind", input_qc_class$qc_train)
        test_set <- do.call("rbind", input_qc_class$qc_test)
        reference_data <- input_qc_class$qc_reference
    }

    train_set$QC_ID <- as.character(train_set$QC_ID)
    train_set$Plate_ID <- as.character(train_set$Plate_ID)

    test_set$QC_ID <- as.character(test_set$QC_ID)
    test_set$Plate_ID <- as.character(test_set$Plate_ID)

    # train_set <- train_set[order(train_set$Plate_ID),]
    # test_set <- test_set[order(test_set$Plate_ID),]

    if(!well_position) {
        train_set <- train_set[!names(train_set) %in% "Well_Position"]
        test_set <- test_set[!names(test_set) %in% "Well_Position"]
    }

    if(!injection_order) {
        train_set <- train_set[!names(train_set) %in% "Injection_Order"]
        test_set <- test_set[!names(test_set) %in% "Injection_Order"]
    }

    if(is.null(log_file)) {
        cl <- parallel::makeCluster(parallel.cores)
    } else {
        cl <- parallel::makeCluster(parallel.cores, outfile = log_file)
    }

    # browser()
    message("  - Calculating correlation...   ", Sys.time())
    train_test_data <- generate_train_test_proc(train_set = train_set,
                                                test_set = test_set,
                                                reference_data = reference_data,
                                                cor_type = cor_type,
                                                cor_method = cor_method,
                                                scale_all = F, select_feature = select_all,
                                                min_candidate_len = min_candidate_len,
                                                max_candidate_len = max_candidate_len,
                                                std = F, using_train_std = F)

    var_names <- intersect(names(train_set), names(test_set))
    var_names <- intersect(names(reference_data), var_names)

    pbapply::pboptions(type = "timer", style = 1, char = "=")
    message("  - Normalising...   ", Sys.time())
    if (batch_wise) {
        var_res <- pbapply::pblapply(var_names, function(var_name, train_test_data, k_batch, k_range,
                                                         algorithm, reference_data, using_train_std,
                                                         checkOutlier_batchVal, use_methods, scale_all,
                                                         xgboost_params, rf_params, rf_stack_mode, std,
                                                         cor_type, cor_method) {
            # var_name <- "Var18"; var_name <- var_names[3];
            source_file <- '/Volumes/Work/Projects/Helmholtz/Helmholtz_normaliation/scripts/metabolomics/norm_RF_Internal.R'
            if(file.exists(source_file)) {
                source(source_file)
            } else {
                # source("~/R_tmp/HGMU/norm_RF_Internal.R")
                source("norm_RF_Internal.R")
            }

            input_test <- train_test_data$test_raw
            input_train <- train_test_data$train_raw
            selected_var_names <- train_test_data$correlation$selected_var
            if (all(input_test$Plate_ID %in% input_train$Plate_ID)) {
                plate_names <- unique(input_test$Plate_ID)
            } else stop("plate ID cannot match!")



            res_plate <- lapply(plate_names, function(plate_name) {
                # plate_name <- plate_names[[2]];
                # message("  - Current variable: ", var_name, ", ",
                # "Current plate: ", plate_name, ", ", Sys.time())
                input_train_plate <- input_train[input_train$Plate_ID == plate_name,]
                input_test_plate  <- input_test[input_test$Plate_ID == plate_name,]

                if (checkOutlier_batchVal) {
                    input_train_plate_num <- input_train_plate[!names(input_train_plate) %in% c("QC_ID", "Plate_ID", "Well_Position", "Injection_Order")]
                    input_train_plate_checked <- as.data.frame(sapply(input_train_plate_num, function(x) {
                        outlier_res <- .boxplot.stats(x)$out

                        replace_val <- median(x[!(x %in% outlier_res)], na.rm = T)

                        x[x %in% outlier_res] <- replace_val
                        x
                    }))
                    input_train_plate <- cbind(input_train_plate[names(input_train_plate) %in% c("QC_ID", "Plate_ID", "Well_Position", "Injection_Order")],
                                               input_train_plate_checked)
                }

                if(is.null(selected_var_names)) select_plate <- T else select_plate <- F

                train_test_plate_data <- generate_train_test_proc(train_set = input_train_plate,
                                                                  test_set = input_test_plate,
                                                                  reference_data = reference_data,
                                                                  scale_all = scale_all,
                                                                  cor_type = cor_type,
                                                                  cor_method = cor_method,
                                                                  min_candidate_len = min_candidate_len,
                                                                  max_candidate_len = max_candidate_len,
                                                                  select_feature = select_plate,
                                                                  std = std, using_train_std = using_train_std)

                if(std) {
                    input_train_plate_scale <- train_test_plate_data$train_set_std
                    input_test_plate_scale <- train_test_plate_data$test_set_std
                } else {
                    input_train_plate_scale <- train_test_plate_data$train_raw
                    input_test_plate_scale <- train_test_plate_data$test_raw
                }

                if (select_plate) selected_var_names <- train_test_plate_data$correlation$selected_var

                input_train_plate_scale <- input_train_plate_scale[ selected_var_names[[var_name]] ]
                train_X <- input_train_plate_scale[!names(input_train_plate_scale) %in% c(var_name, "QC_ID", "Plate_ID")]
                # train_X <- input_train_plate_scale[!names(input_train_plate_scale) %in% c("QC_ID", "Plate_ID")]
                test_X <- input_test_plate_scale[names(input_test_plate_scale) %in% colnames(train_X)]

                train_nonNA_idx <- sapply(train_X, function(x) any(is.finite(x)))
                test_nonNA_idx <- sapply(test_X, function(x) any(is.finite(x)))
                nonNA_idx <- train_nonNA_idx & test_nonNA_idx

                train_X <- train_X[nonNA_idx]
                test_X <- test_X[nonNA_idx]

                train_X <- train_X[names(train_X) %in% names(test_X)]

                # message("    - feature num: ", ncol(train_X) - 1, " ; sample num: ", nrow(train_X))
                # message("    - current feature: ", var_name, ": ", paste(names(train_X), collapse = " "))

                test_y_raw <- train_test_plate_data$test_raw[[var_name]]
                # test_y_true <- train_test_plate_data$test_true[[var_name]]

                train_y_raw <- train_test_plate_data$train_raw[[var_name]]
                train_y_true <- train_test_plate_data$train_true[[var_name]]



                ######
                if(plate_based_refer) {
                    if(!all(train_y_raw == 0)) {

                        # outlier_res_batch_refer <- .boxplot.stats(train_y_raw, coef = 1.5)$out
                        # train_y_raw_batch_refer <- train_y_raw
                        # train_y_raw_batch_refer[train_y_raw_batch_refer %in% outlier_res_batch_refer] <-
                        #     mean(train_y_raw_batch_refer[!train_y_raw_batch_refer %in% outlier_res_batch_refer], na.rm = T)
                        # train_y_raw_batch_refer <- mean(train_y_raw_batch_refer, na.rm = T)
                        # if(train_y_raw_batch_refer != 0) {
                        #     train_y <- (train_y_raw - train_y_raw_batch_refer) / train_y_raw_batch_refer
                        # } else {
                        #     train_y <- scale(train_y_raw, scale = F)
                        #     if(attr(train_y, "scaled:center") != 0) {
                        #         train_y <- train_y / attr(train_y, "scaled:center")
                        #     } else train_y <- 0; train_y_raw_batch_refer <- 0
                        # }


                        train_y_raw_batch_refer <- mean(train_y_raw, na.rm = T)
                        train_y <- train_y_raw - train_y_raw_batch_refer
                        if(train_y_raw_batch_refer != 0) {
                            train_y <- train_y / train_y_raw_batch_refer
                        } else train_y <- 0



                        # train_y <- train_y_raw - median(train_y_raw, na.rm = T)
                        #
                        # if(median(train_y_raw, na.rm = T) != 0) {
                        #     train_y <- train_y / median(train_y_raw, na.rm = T)
                        # } else train_y <- 0

                    } else {
                        train_y_raw_batch_refer <- 0
                        train_y <- 0
                    }
                } else {
                    train_y <- train_test_plate_data$train_target[[var_name]]
                }


                ######



                if ("TIGER" %in% use_methods) {
                    # message("rf_stack")
                    input_train_data <- cbind(y = train_y, train_X)
                    mtry <- round(seq(0.2, 0.8, 0.2) * ncol(train_X))
                    nodesize <- round(seq(0.2, 0.8, 0.2) * nrow(train_X))
                    # samplesize <- round(seq(0.6, 1, 0.2) * nrow(train_X))
                    rf_params <- expand.grid(mtry = unique(mtry),
                                             nodesize = unique(nodesize))
                    # samplesize <- round(seq(0.3, 0.8, 0.2) * nrow(train_X))
                    # rf_params_tmp <- expand.grid(mtry = mtry,
                    #                               samplesize = samplesize)
                    # rf_params_tmp_2 <- lapply(seq(0.3, 0.8, 0.2) , function(nodesize_ratio, rf_params_tmp) {
                    #     rf_params_tmp$nodesize <- round(rf_params_tmp$samplesize * nodesize_ratio)
                    #     rf_params_tmp
                    # }, rf_params_tmp = rf_params_tmp)
                    # rf_params <- do.call("rbind", rf_params_tmp_2)
                    # # rf_params <- rbind(rf_params_tmp, rf_params_tmp_2)
                    # rf_params <- unique(rf_params)

                    ensemble_y_pred_rf <- lapply(1:nrow(rf_params), function(idx, input_train_data,
                                                                             rf_params, test_X) {

                        if (cv_mode) {
                            folds_train <- caret::createFolds(1:length(input_train_data$y), k = 5, returnTrain = T)

                            # if reference == 0, put raw values back.

                            res_folds <- lapply(folds_train, function(train_idx, rf_params) {
                                input_test_fold_raw <- train_y_raw[-train_idx]
                                input_test_fold_true <- train_y_true[-train_idx]

                                input_train_fold <- input_train_data[train_idx,]
                                input_test_fold  <- input_train_data[-train_idx,]
                                RF_fold_mod <- randomForest::randomForest(y ~., data = input_train_fold,
                                                                          ntree = 500, replace = T,
                                                                          mtry = rf_params$mtry[idx],
                                                                          # samplesize = rf_params$samplesize[idx],
                                                                          nodesize = rf_params$nodesize[idx])

                                pred_rf_fold <- predict(RF_fold_mod, input_test_fold)
                                pred_rf_fold_df <- data.frame(pred = pred_rf_fold,
                                                              target = input_test_fold$y,
                                                              raw = input_test_fold_raw,
                                                              true = input_test_fold_true)
                            }, rf_params = rf_params)

                            train_y_pred_rf <- do.call("rbind", res_folds)

                            RF_mod <- randomForest::randomForest(y ~., data = input_train_data,
                                                                 ntree = 500, replace = T,
                                                                 mtry = rf_params$mtry[idx],
                                                                 # samplesize = rf_params$samplesize[idx],
                                                                 nodesize = rf_params$nodesize[idx])
                            test_y_pred_rf <- predict(RF_mod, test_X)
                            out_pred <- list(pred_train = train_y_pred_rf,
                                             pred_test = test_y_pred_rf)

                        } else {
                            RF_mod <- randomForest::randomForest(y ~., data = input_train_data,
                                                                 ntree = 500, replace = T,
                                                                 mtry = rf_params$mtry[idx],
                                                                 # samplesize = rf_params$samplesize[idx],
                                                                 nodesize = rf_params$nodesize[idx])
                            train_y_pred_rf <- predict(RF_mod, input_train_data)
                            train_y_pred_rf <- data.frame(pred = train_y_pred_rf,
                                                          target = input_train_data$y,
                                                          raw = train_y_raw,
                                                          true = train_y_true)

                            test_y_pred_rf <- predict(RF_mod, test_X)
                            out_pred <- list(pred_train = train_y_pred_rf, pred_test = test_y_pred_rf)
                        }
                        out_pred

                    }, input_train_data = input_train_data, rf_params = rf_params, test_X)

                    ensemble_train_y_pred <- sapply(ensemble_y_pred_rf, function(x, input_train_data, train_y_true, train_y_raw, rf_stack_mode) {
                        train_y_pred <- x$pred_train
                        train_normFactor <- train_y_pred$pred + 1
                        train_recover <- train_y_pred$raw / train_normFactor
                        train_error_ratio <- (abs(train_recover - train_y_pred$true)/train_y_pred$true)
                        out <- train_error_ratio
                    }, input_train_data = input_train_data, train_y_true = train_y_true, train_y_raw = train_y_raw, rf_stack_mode = rf_stack_mode)
                    # ensemble_train_y_pred[!is.finite(ensemble_train_y_pred)] <- NA

                    # mod_weight <- 1/colMeans(ensemble_train_y_pred, na.rm = T)
                    mod_weight <- 1/exp(colMeans(ensemble_train_y_pred, na.rm = T))
                    mod_weight_normlised <- mod_weight/sum(mod_weight)

                    ensemble_test_y_pred <- lapply(ensemble_y_pred_rf, function(x) x$pred_test)
                    ensemble_test_y_pred <- as.data.frame(do.call("rbind", ensemble_test_y_pred))
                    if (nrow(ensemble_test_y_pred) != length(mod_weight_normlised)) stop("ERROR: mod weights cannot match!")
                    ensemble_test_y_pred_weighted <- sapply(ensemble_test_y_pred, function(x, mod_weight_normlised) {
                        sum(x * mod_weight_normlised)
                    }, mod_weight_normlised = mod_weight_normlised)
                    test_y_pred <- ensemble_test_y_pred_weighted

                    test_y_pred <- test_y_pred
                } else stop("ERROR: Wrong method name!")

                if(plate_based_refer) {

                    test_raw_refer <- input_test[[var_name]]
                    # outlier_res_test_raw_refer <- .boxplot.stats(test_raw_refer, coef = 1.5)$out
                    # test_raw_refer[test_raw_refer %in% outlier_res_test_raw_refer] <-
                    #     mean(test_raw_refer[!test_raw_refer %in% outlier_res_test_raw_refer], na.rm = T)
                    test_raw_refer <- mean(test_raw_refer, na.rm = T)

                    test_y_raw_batch <- test_y_raw
                    # outlier_res_test_y_raw_batch <- .boxplot.stats(test_y_raw_batch, coef = 1.5)$out
                    # test_y_raw_batch[test_y_raw_batch %in% outlier_res_test_y_raw_batch] <-
                    #     mean(test_y_raw_batch[!test_y_raw_batch %in% outlier_res_test_y_raw_batch], na.rm = T)
                    test_y_raw_batch <- mean(test_y_raw_batch, na.rm = T)

                    if(test_y_raw_batch != 0) {
                        validate_recover <- test_y_raw / ( (1 + test_y_pred) * test_y_raw_batch ) * test_raw_refer
                    } else validate_recover <- test_y_raw

                    # validate_recover <- test_y_raw / ( test_y_pred + mean(test_y_raw, na.rm = T) ) * median(input_test[[var_name]], na.rm = T)
                } else {
                    normFactor <- test_y_pred + 1
                    validate_recover <- test_y_raw / normFactor
                }


                out <- data.frame(Plate_ID = input_test_plate_scale$Plate_ID, QC_ID = input_test_plate_scale$QC_ID, validate_recover)
                names(out)[3] <- var_name
                out
            })
            res_plate_df <- do.call("rbind", res_plate)
        },
        train_test_data = train_test_data, k_batch = k_batch, k_range = k_range,
        algorithm = algorithm, reference_data = reference_data,
        scale_all = scale_all, std = std, cor_type = cor_type, cor_method = cor_method,
        using_train_std = using_train_std, checkOutlier_batchVal = checkOutlier_batchVal,
        use_methods = use_methods, xgboost_params = xgboost_params, rf_params = rf_params,
        rf_stack_mode = rf_stack_mode, cl = cl)

        formatted_pred_test <- var_res[[1]]
        formatted_pred_test <- cbind(label = paste(formatted_pred_test$QC_ID, formatted_pred_test$Plate_ID, sep = sep_char),
                                     formatted_pred_test)

        for(idx in 2:length(var_res)) {
            if(!all(formatted_pred_test$Plate_ID == var_res[[idx]]$Plate_ID)) stop("ERROR: Plate_ID!")
            if(!all(formatted_pred_test$QC_ID == var_res[[idx]]$QC_ID)) stop("ERROR: QC_ID!")
            formatted_pred_test <- cbind(formatted_pred_test, var_res[[idx]][-c(1,2)])
        }
        row.names(formatted_pred_test) <- NULL
        formatted_pred_test <- formatted_pred_test[!names(formatted_pred_test) %in% c("Plate_ID", "QC_ID")]

        output <- formatted_pred_test
        if (output_pred_test) {
            tmp <- data.frame(t(output))
            # if na or nan or negative values, put original values back.
            na_idx <- is.na(tmp)
            tmp[na_idx] <- 0
            write.table(tmp, file = paste0("normalized by - ", use_methods, output_file_label, ".csv"),
                        row.names = T, col.names = F, quote = F, sep = ",")
        }

    } else {
        var_res <- parallel::parLapply(cl, var_names, function(var_name, train_test_data, k_range, reference_data, using_train_std, algorithm)
        {
            # print(var_name)
            source_file <- '/Volumes/Work/Projects/Helmholtz/Helmholtz_normaliation/scripts/metabolomics/norm_RF_Internal.R'
            if(file.exists(source_file)) {
                source(source_file)
            } else {
                # source("~/R_tmp/HGMU/norm_RF_Internal.R")
                source("norm_RF_Internal.R")
            }

            input_test <- train_test_data$test_set_std
            input_train <- train_test_data$train_raw
            selected_var_names <- train_test_data$correlation$selected_var

            dummy_plate_id <- data.frame(model.matrix(~Plate_ID,
                                                      data = rbind(input_train["Plate_ID"],
                                                                   input_test["Plate_ID"])))[-1]
            dummy_plate_train <- dummy_plate_id[1:nrow(input_train),]
            dummy_plate_test  <- dummy_plate_id[(nrow(input_train) + 1):nrow(dummy_plate_id),]

            set.seed(1)
            folds <- caret::createFolds(y = input_train$Plate_ID, k = 5, returnTrain = T)

            res_tune <- lapply(k_range, function(k) {
                res_folds <- lapply(1:5, function(fold_idx) {
                    input_fold <- folds[[fold_idx]]

                    train_set <- input_train[input_fold,]
                    validation_set <- input_train[-input_fold,]

                    train_set_categorical <- dummy_plate_train[input_fold,]
                    validation_set_categorical <- dummy_plate_train[-input_fold,]

                    train_validate_data <- generate_train_test_proc(train_set = train_set,
                                                                    test_set = validation_set,
                                                                    reference_data = reference_data,
                                                                    cor_type = cor_type,
                                                                    cor_method = cor_method,
                                                                    min_candidate_len = 99999,
                                                                    max_candidate_len = 99999,
                                                                    select_feature = F,
                                                                    std = T, using_train_std = using_train_std)

                    input_train_fold <- train_validate_data$train_set_std
                    input_valid_fold <- train_validate_data$test_set_std

                    input_predictor_fold <- input_train_fold[ selected_var_names[[var_name]] ]
                    # input_predictor_fold <- input_predictor_fold[!names(input_predictor_fold) %in% var_name]
                    input_train_fold <- cbind(y = train_validate_data$train_target[[var_name]], input_predictor_fold)
                    # if(var_name %in% names(input_train_fold)) stop("ERROR: Target found in predictor!")

                    input_train_fold <- cbind(train_set_categorical, input_train_fold)
                    y <- input_train_fold[["y"]]
                    X <- input_train_fold[!names(input_train_fold) %in% c("Plate_ID", "QC_ID", "y")]
                    X <- as.matrix(X)

                    input_valid_fold <- cbind(validation_set_categorical, input_valid_fold)
                    validate_set_y_raw <- train_validate_data$test_raw[[var_name]]
                    validate_set_y_true <- train_validate_data$test_true[[var_name]]
                    validate_set_qc_id <- input_valid_fold$QC_ID
                    input_validate_X <- input_valid_fold[colnames(input_valid_fold) %in% colnames(X)]
                    input_validate_X <- as.matrix(input_validate_X)

                    kNNreg_mod <- FNN::knn.reg(train = X, test = input_validate_X,
                                               y = y, k = k, algorithm = algorithm)

                    validate_pred <- kNNreg_mod$pred
                    normFactor <- validate_pred + 1
                    validate_recover <- validate_set_y_raw / normFactor

                    perf_validate_recover <- summary_prediction(preds = validate_recover,
                                                                trues = validate_set_y_true,
                                                                raws = validate_set_y_raw,
                                                                qc_id = validate_set_qc_id,
                                                                label = paste0("fold_", fold_idx))[-c(10:19)]

                })
                res_folds_df <- do.call("rbind", res_folds)
                if(length(unique(res_folds_df$type)) > 1) {
                    res_folds_df <- res_folds_df[res_folds_df$type == "Overall",]
                }
                res_folds_mean <- data.frame(t(colMeans(res_folds_df[sapply(res_folds_df, class) == "numeric"], na.rm = T)))
                res_folds_mean <- cbind(var = var_name, k = k, res_folds_mean)
            })
            res_tune_df <- do.call("rbind", res_tune)
            perf_tune <- res_tune_df[which.min(res_tune_df$errDiffPer_absMean_recover),]
        }, train_test_data = train_test_data, k_range = k_range, reference_data = reference_data, using_train_std = using_train_std, algorithm = algorithm)

        perf_validate <- do.call("rbind", var_res)
        output <- list(perf_validation = perf_validate)

        if(predict_test) {
            best_k <- as.list(perf_validate$k)
            names(best_k) <- perf_validate$var

            test_pred <- parallel::parLapply(cl, var_names, function(var_name, train_test_data, best_k) {
                input_test <- train_test_data$test_set_std
                input_train <- train_test_data$train_set_std
                selected_var_names <- train_test_data$correlation$selected_var

                dummy_plate_id <- data.frame(model.matrix(~Plate_ID,
                                                          data = rbind(input_train["Plate_ID"],
                                                                       input_test["Plate_ID"])))[-1]
                dummy_plate_train <- dummy_plate_id[1:nrow(input_train),]
                dummy_plate_test  <- dummy_plate_id[(nrow(input_train) + 1):nrow(dummy_plate_id),]

                input_predictor <- input_train[ selected_var_names[[var_name]] ]
                # input_predictor <- input_predictor[!names(input_predictor) %in% var_name]
                input_train <- cbind(y = train_test_data$train_target[[var_name]], input_predictor)
                # if(var_name %in% names(input_train)) stop("ERROR: Target found in predictor!")

                input_train <- cbind(dummy_plate_train, input_train)
                y <- input_train[["y"]]
                X <- input_train[!names(input_train) %in% c("Plate_ID", "QC_ID", "y")]
                X <- as.matrix(X)

                input_test <- cbind(dummy_plate_test, input_test)
                test_y_raw <- train_test_data$test_raw[[var_name]]
                # validate_set_y_true <- train_test_data$test_true[[var_name]]
                # validate_set_qc_id <- input_test$QC_ID
                test_X <- input_test[colnames(input_test) %in% colnames(X)]
                test_X <- as.matrix(test_X)

                kNNreg_test_mod <- FNN::knn.reg(train = X, test = test_X,
                                                y = y, k = best_k[[var_name]], algorithm = "kd_tree")

                test_pred <- kNNreg_test_mod$pred
                normFactor <- test_pred + 1
                test_recover <- test_y_raw / normFactor
                as.data.frame(test_recover)
            },
            train_test_data = train_test_data, best_k = best_k)
            test_pred_df <- do.call("cbind", test_pred)
            test_pred_df <- data.frame(t(test_pred_df))
            names(test_pred_df) <- train_test_data$test_raw$QC_ID
            test_pred_df <- cbind(label = var_names, test_pred_df)
            row.names(test_pred_df) <- NULL
            if (output_pred_test) {
                write.table(test_pred_df, file = paste0("normalized by - kNN", output_file_label, ".csv"),
                            row.names = F, quote = F, sep = ",")
            }
            output <- c(output, predicted_test = list(test_pred_df))
        }
    }
    message("+ Completed.   ", Sys.time())
    parallel::stopCluster(cl)
    output
}


out <- train_test_kNNreg(input_qc_class = input_qc_class, parallel.cores = parallel::detectCores()-2,
                         batch_wise = T, min_candidate_len = 10, max_candidate_len = 30,
                         scale_mode = "none", feature_selection_mode = "all",
                         cor_type = "pcor", cor_method = "spearman",
                         plate_based_refer = F, cv_mode = T, well_position = F,
                         injection_order = F, sep_char = "_", use_methods = "TIGER",
                         rf_stack_mode = "weight_recover_y", log_file = "log",
                         output_file_label = "_pcor_spearman")

out <- train_test_kNNreg(input_qc_class = input_qc_class, parallel.cores = parallel::detectCores()-2,
                         batch_wise = T, min_candidate_len = 10, max_candidate_len = 10,
                         scale_mode = "none", feature_selection_mode = "all",
                         cor_type = "pcor", cor_method = "pearson",
                         plate_based_refer = F, cv_mode = T, well_position = F,
                         injection_order = T, sep_char = "_", use_methods = "TIGER",
                         rf_stack_mode = "weight_recover_y", log_file = "log",
                         output_file_label = "_pcor_pearson_fixed10")

out <- train_test_kNNreg(input_qc_class = input_qc_class, parallel.cores = parallel::detectCores()-2,
                         batch_wise = T, min_candidate_len = 10, max_candidate_len = 10,
                         scale_mode = "none", feature_selection_mode = "all",
                         cor_type = "cor", cor_method = "spearman",
                         plate_based_refer = F, cv_mode = T, well_position = F,
                         injection_order = T, sep_char = "_", use_methods = "TIGER",
                         rf_stack_mode = "weight_recover_y", log_file = "log",
                         output_file_label = "_cor_spearman_fixed10")

out <- train_test_kNNreg(input_qc_class = input_qc_class, parallel.cores = parallel::detectCores()-2,
                         batch_wise = T, min_candidate_len = 10, max_candidate_len = 10,
                         scale_mode = "none", feature_selection_mode = "all",
                         cor_type = "cor", cor_method = "pearson",
                         plate_based_refer = F, cv_mode = T, well_position = F,
                         injection_order = T, sep_char = "_", use_methods = "TIGER",
                         rf_stack_mode = "weight_recover_y", log_file = "log",
                         output_file_label = "_cor_pearson_fixed10")

# out <- train_test_kNNreg(input_qc_class = input_qc_class, parallel.cores = 30,
#                          batch_wise = T, min_candidate_len = 5, max_candidate_len = 10,
#                          scale_mode = "none", feature_selection_mode = "all",
#                          plate_based_refer = T, cv_mode = T, well_position = F,
#                          injection_order = T, sep_char = "_", use_methods = "TIGER",
#                          rf_stack_mode = "weight_recover_y", log_file = NULL,
#                          output_file_label = "_stackWeightCV_fea5.10_order")
#
# out <- train_test_kNNreg(input_qc_class = input_qc_class, parallel.cores = 30,
#                          batch_wise = T, min_candidate_len = 10, max_candidate_len = 30,
#                          scale_mode = "none", feature_selection_mode = "all",
#                          plate_based_refer = T, cv_mode = T, well_position = F,
#                          injection_order = F, sep_char = "_", use_methods = "TIGER",
#                          rf_stack_mode = "weight_recover_y", log_file = NULL,
#                          output_file_label = "_stackWeightCV_fea10.30")

###### fill in NA ######

# bak <- read.table("normalized by - normFactor.csv", sep = ",",
#                   row.names = 1)
# bak_bak <- bak
#
# table(which(is.na(bak)) == which(is.na(bak_bak)))
#
# TIGER_stackWeightCV_fea5.10_order <- data.frame(t(TIGER_stackWeightCV_fea5.10_order))
#
# bak[103,104] == TIGER_stackWeightCV_fea5.10_order[103,104]
# bak[(is.na(bak))] <- TIGER_stackWeightCV_fea5.10_order[(is.na(bak))]
# anyNA(bak)
# write.table(bak, row.names = T, col.names = F, quote = F,
#             file = "normalized by - TIGER_stackWeightCV_fea5.10.csv", sep = ",")

###
#
# bak <- read.table("normalized by - TIGER_bak_stackWeightTrain.csv", sep = ",",
#                   row.names = 1)
#
# tmp <- sapply(bak[1,], function(x) {
#     strsplit(x, ".", fixed = T)[[1]][[1]]
# })
#
# tmp <- data.frame(val = tmp, ori = unlist(bak[1,]))
# bak[1,] <- tmp$val
#
# write.table(bak, row.names = T, col.names = F, quote = F,
#                         file = "normalized by - TIGER_bak_stackWeightTrain.csv", sep = ",")
