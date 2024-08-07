#' @import dplyr
#' @import ggplot2

# Constructor for class 'sineFit'
#' @title sineFitData
#' @export
sineFitData = setClass("sineFitData", slots = c(Results = "data.frame",
                                                RSS_per_period = "data.frame",
                                                FStat_per_period = "data.frame",
                                                FStat_PVal_per_period = "data.frame",
                                                RSquared_per_period = "data.frame",
                                                Adj_RSquared_per_period = "data.frame",
                                                Input_data = "data.frame",
                                                Row_data = "data.frame"))

# Modify show function for sineFit objects
#' @export
setMethod("show", "sineFitData",
          function(object){
            num_samples = nrow(slot(object, "Input_data"))
            num_genes = ncol(slot(object, "Input_data")) - 1
            rowData_names = colnames(slot(object, "Row_data"))
            if(length(rowData_names) > 2){
              term_string = paste0("\nAdditional linear model coefficients: ",
                                   rowData_names[3:length(rowData_names)], collapse = ", ")
            } else {
              term_string = ""
            }
            rss_tb = slot(object, "RSS_per_period")
            min_per = min(rss_tb$Period)
            max_per = max(rss_tb$Period)
            if(min_per == max_per){
              period_string = paste0("Generated with fixed period: ", min_per)
            } else {
              period_string = paste0("Generated with variable period: Min. period = ", min_per, ", Max. period = ", max_per)
            }
            # Message:
            cat(paste0("A sineFitData object:\n", num_samples, " samples, ", num_genes, " genes\n",
                       period_string, term_string))
          }
)


# Accessor method for results
#' @name results
#' @title results
#' @export
#' @param object A sinFitData object
results = function(object) slot(object, "Results")
# setGeneric("results", function(object)
#   standardGeneric("results"))
# setMethod("results", signature(object = "sineFitData"), results.sineFit)
#setGeneric("results", def = results)
setMethod("results", signature(object = "sineFitData"), results)


# Define plotting function for RSS
#' @name plot_RSS
#' @title plot_RSS
#' @export
#' @param object A sineFitData object
#' @param gene Gene name
plot_RSS = function(object, gene){
  rss_tb = slot(object, "RSS_per_period")
  min_per = min(rss_tb$Period)
  max_per = max(rss_tb$Period)
  ggplot(rss_tb, aes(x = Period, y = .data[[gene]])) +
    geom_line() +
    ggtitle(paste0("RSS vs period: ", gene)) +
    ylab("RSS") +
    scale_x_continuous(breaks = seq(min_per, max_per, by = 1))
}
#setGeneric("plot_RSS", def = plot_RSS)
setMethod("plot_RSS", signature(object = "sineFitData"), plot_RSS)

# Define plotting function for F-test (F-statistic or p-value)
#' @name plot_FStat
#' @title plot_FStat
#' @export
#' @param object A sineFitData object
#' @param gene Gene name
#' @param plot_pvalues Logical to indicate if p-values should be plotted
plot_FStat = function(object, gene, plot_pvalues = FALSE){
  if(plot_pvalues){
    fstat_tb = slot(object, "FStat_PVal_per_period")
    plot_title = paste0("F-Statistic p-values vs period: ", gene)
    y_lab = "p-value"
  } else {
    fstat_tb = slot(object, "FStat_per_period")
    plot_title = paste0("F-Statistic vs period: ", gene)
    y_lab = "F-Statistic"
  }
  min_per = min(fstat_tb$Period)
  max_per = max(fstat_tb$Period)
  ggplot(fstat_tb, aes(x = Period, y = .data[[gene]])) +
    geom_line() +
    ggtitle(plot_title) +
    ylab(y_lab) +
    scale_x_continuous(breaks = seq(min_per, max_per, by = 1))
}
#setGeneric("plot_FStat", def = plot_FStat)
setMethod("plot_FStat", signature(object = "sineFitData"), plot_FStat)

# Define plotting function for R-Squared (or adjusted R-Squared)
#' @name plot_RSquared
#' @title plot_RSquared
#' @export
#' @param object A sineFitData object
#' @param gene Gene name
#' @param plot_adj_rsquared Logical to indicate if adjusted R-Squared should be plotted
plot_RSquared = function(object, gene, plot_adj_rsquared = FALSE){
  if(plot_adj_rsquared){
    rsq_tb = slot(object, "Adj_RSquared_per_period")
    plot_title = paste0("Adjusted R-Squared vs period: ", gene)
    y_lab = "Adjusted R-Squared"
  } else {
    rsq_tb = slot(object, "RSquared_per_period")
    plot_title = paste0("R-Squared vs period: ", gene)
    y_lab = "R-Squared"
  }
  min_per = min(rsq_tb$Period)
  max_per = max(rsq_tb$Period)
  ggplot(rsq_tb, aes(x = Period, y = .data[[gene]])) +
    geom_line() +
    ggtitle(plot_title) +
    ylab(y_lab) +
    scale_x_continuous(breaks = seq(min_per, max_per, by = 1))
}
#setGeneric("plot_RSquared", def = plot_RSquared)
setMethod("plot_RSquared", signature(object = "sineFitData"), plot_RSquared)


# Function to get a set of points on sine wave
#' @name get_sine_points
#' @title get_sine_points
#' @param x Co-ordinates on x-axis
#' @param vert Vertical offset
#' @param amp Amplitude
#' @param per Period
#' @param horiz Horizontal offset
get_sine_points = function(x, vert, amp, per, horiz){
  y = vert + amp * sin((2*pi*x/per) + horiz)
  return(y)
}


# Function to get peak phase (i.e. time at which sine wave is at max) for given amplitude, period and offsets
#' @name get_peak_phase
#' @title get_peak_phase
#' @param vert Vertical offset
#' @param amp Amplitude
#' @param per Period
#' @param horiz Horizontal offset
get_peak_phase = function(vert, amp, per, horiz){
  x = seq(0, 24, length.out = 288)
  sine_points = get_sine_points(x, vert, amp, per, horiz)
  peak_phase = x[which.max(sine_points)]
  return(peak_phase)
}


# Function to get RSS from QR decomposition
#' @name get_RSS
#' @title get_RSS
#' @param qr_obj QR decomposition object
#' @param y_mat Matrix to use as second argument to qr()
get_RSS = function(qr_obj, y_mat){
  colSums((qr.resid(qr_obj, y_mat))^2)
}


# Function to get sine-fits for genes using QR decomposition (fixed 24 hour period by default)
#' @name sineFit
#' @title sineFit
#' @export
#' @param in_data Data frame with sample IDs and expression values
#' @param rowData Data frame with sample IDs, time points and any other covariates
#' @param min_per Minimum period to generate fit
#' @param max_per Maximum period to generate fit
sineFit = function(in_data, rowData, min_per = 24, max_per = 24){

  # Sort in_data and rowData by sample ID
  in_data = in_data[order(in_data[[1]]), ]
  rowData = rowData[order(rowData[[1]]), ]

  # Check in_data and row data match
  stopifnot("Sample IDs in in_data and rowData do not match" = in_data[[1]] == rowData[[1]])

  # Generate periods to test
  period_vec = seq(min_per, max_per, length.out = 1 + (max_per - min_per) * 12) # every 5 mins

  # Get basis functions phi_1 and phi_2 for each period
  phi_1_df = sapply(period_vec, function(p) sin(2*pi*(rowData[[2]])/p)) %>% as.data.frame() %>% round(., 7)
  phi_2_df = sapply(period_vec, function(p) cos(2*pi*(rowData[[2]])/p)) %>% as.data.frame() %>% round(., 7)

  # Get a list of design matrices including additional covariates in rowData
  if(ncol(rowData) > 2){
    x_list = lapply(1:length(period_vec),
                    function(x) {bind_cols(
                      tibble(intercept = 1, phi_1 = phi_1_df[[x]], phi_2 = phi_2_df[[x]]),
                      rowData[, 3:ncol(rowData)])
                    }
    )
  } else {
    x_list = lapply(1:length(period_vec), function(x) {
      tibble(intercept = 1, phi_1 = phi_1_df[[x]], phi_2 = phi_2_df[[x]])
    }
    )
  }

  # Get a list of QR decompositions (one for each design matrix)
  qr_x_list = lapply(x_list, qr)

  # Get RSS per gene per period)
  rss_per_period = sapply(qr_x_list, get_RSS, y_mat = as.matrix(in_data[, -1])) %>% as.matrix()

  # Get TSS for all genes
  lm_mean_devs = sweep(in_data[, -1], MARGIN = 2, STATS = colMeans(in_data[, -1]), FUN = "-")
  lm_tss = colSums(lm_mean_devs^2)

  # Get F-statistic and p-value for each period for each gene
  df_full = nrow(in_data) - ncol(x_list[[1]]) # all design matrices in x_list have same number of columns, just using the first one
  df_reduced = ncol(x_list[[1]]) - 1
  # Get (TSS - RSS) * df_full per gene per period
  fstat_numerator = sweep(rss_per_period, MARGIN = 1, lm_tss, FUN = "-") * (-1 * df_full)
  # Get RSS * df_reduced per gene per period
  fstat_denominator = rss_per_period * df_reduced
  # Calculate F-statistics for each period
  fstat_per_period = fstat_numerator/fstat_denominator
  # Get F-test p-values
  fstat_pvals = sapply(1:ncol(fstat_per_period), function(p) pf(fstat_per_period[, p], df_reduced, df_full, lower.tail = FALSE)) %>%
    as.matrix()

  # Get R-squared (1 - RSS/TSS)
  rsquared_per_period = sweep(rss_per_period, MARGIN = 1, lm_tss, FUN = "/") * (-1) + 1
  # Get adjusted R-squared (using same formula as lm; see summary.lm)
  adj_rsquared_per_period = 1 - ((1 - rsquared_per_period) * ((nrow(in_data) - 1) / df_full ))

  # Get period which minimises RSS for each gene
  min_rss_period = apply(rss_per_period, 1, which.min)
  # Get coefficients for genes using the period which minimises RSS
  qr_coef = mapply(qr.coef, qr_x_list[min_rss_period], in_data[names(min_rss_period)])

  # Get F-statistic and p-value for the period which minimises RSS for each gene
  lm_fstat = sapply(1:nrow(fstat_per_period), function(x) fstat_per_period[x, min_rss_period[x]])
  lm_fstat_pval = sapply(1:nrow(fstat_pvals), function(x) fstat_pvals[x, min_rss_period[x]])
  # Get adjusted p-values
  lm_fstat_padj = p.adjust(lm_fstat_pval)

  # Get RSS for the period which minimises RSS for each gene
  lm_rss = sapply(1:nrow(rss_per_period), function(x) rss_per_period[x, min_rss_period[x]])

  # Get R-Squared and adjusted R-Squared for the period which minimises RSS for each gene
  lm_rsquared = sapply(1:nrow(rsquared_per_period), function(x) rsquared_per_period[x, min_rss_period[x]])
  lm_adj_rsquared = sapply(1:nrow(adj_rsquared_per_period), function(x) adj_rsquared_per_period[x, min_rss_period[x]])

  # Calculate vertical shift, amplitude, phase shift and base
  lm_vert = qr_coef[1, ]
  beta_1 = qr_coef[2, ]
  beta_2 = qr_coef[3, ]

  lm_amp = sqrt(beta_1^2 + beta_2^2)
  lm_horiz = atan2(beta_2, beta_1)

  # Get peak phase (i.e. time at which sine wave is at max)
  lm_phase = mapply(get_peak_phase, lm_vert, lm_amp, period_vec[min_rss_period], lm_horiz)

  # Make tibble of results
  result_tb = tibble(Gene = names(lm_rss),
                     lm_Ftest_Pval = lm_fstat_pval, lm_Ftest_Adj_Pval = lm_fstat_padj,
                     lm_RSquared = lm_rsquared, lm_Adj_RSquared = lm_adj_rsquared,
                     lm_RSS = as.vector(lm_rss),
                     Base = lm_vert, Amp = lm_amp, Period = period_vec[min_rss_period],
                     Phase_Shift = lm_horiz, Peak_Phase = lm_phase)

  # Rearrange RSS per period
  rss_per_period = t(rss_per_period) %>% as.data.frame() %>% as_tibble() %>%
    mutate(Period = period_vec, .before =  everything())

  # Rearrange fstat_per_period
  fstat_per_period = t(fstat_per_period) %>% as.data.frame() %>% as_tibble() %>%
    mutate(Period = period_vec, .before =  everything())

  # Rearrange fstat_pvals
  fstat_pvals = t(fstat_pvals) %>% as.data.frame() %>% as_tibble() %>%
    mutate(Period = period_vec, .before =  everything())

  # Rearrange rsquared_per_period
  rsquared_per_period = t(rsquared_per_period) %>% as.data.frame() %>% as_tibble() %>%
    mutate(Period = period_vec, .before =  everything())

  # Rearrange adj_rsquared_per_period
  adj_rsquared_per_period = t(adj_rsquared_per_period) %>% as.data.frame() %>% as_tibble() %>%
    mutate(Period = period_vec, .before =  everything())

  # Return S4 sineFitData object
  sf = sineFitData(Results = result_tb,
                   RSS_per_period = rss_per_period,
                   FStat_per_period = fstat_per_period,
                   FStat_PVal_per_period = fstat_pvals,
                   RSquared_per_period = rsquared_per_period,
                   Adj_RSquared_per_period = adj_rsquared_per_period,
                   Input_data = in_data,
                   Row_data = rowData)

  return(sf)
}
#setGeneric("sineFit", def = sineFit)
setMethod("sineFit", signature(in_data = "data.frame", rowData = "data.frame"), sineFit)


# Function to plot genes with best sine fits
#' @name plot_Fit
#' @title plot_Fit
#' @export
#' @param object A sineFitData object
#' @param genes A vector of genes to plot
#' @param y_lower Optional lower y-axis limit
#' @param y_upper Optional upper y-axis limit
#' @param plot_title Optional plot title
plot_Fit = function(object, genes, y_lower = NULL, y_upper = NULL, plot_title = NULL){
  # object = wk00.L.sf.24; genes = c("BMAL1", "PER1", "PER3", "TNMD"); y_lower = 2; y_upper = 9

  # Get gene expression data
  gene_expr = slot(object, "Input_data") %>% select(1, all_of(genes))
  time_data = slot(object, "Row_data") %>% pull(2)
  # These are in the same order (checked when creating the sinFit object)
  gene_expr$time_point = time_data
  gene_data = pivot_longer(gene_expr, cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
    mutate_at("gene", as.factor)

  # Get sine wave variables
  fit_vars = results(object) %>% filter(Gene %in% genes)
  # Get points for sine waves
  x = seq(0, 24, length.out = 288)
  sine_points = mapply(get_sine_points,
                       fit_vars$Base, fit_vars$Amp, fit_vars$Period, fit_vars$Phase_Shift,
                       MoreArgs = list(x = x)) %>% as.data.frame()
  colnames(sine_points) = fit_vars$Gene
  sine_points$x = x
  sine_data = pivot_longer(sine_points, cols = all_of(genes), names_to = "gene", values_to = "y") %>%
    mutate_at("gene", as.factor)

  # Make plot:
  # Set plot title if not provided
  if(is.null(plot_title)){
    plot_title = paste0(deparse(substitute(object)), " : ", paste(genes, collapse = ", "))
  }
  # Set y-axis limits if not provided
  y_range = range(c(gene_data$expression, sine_data$y))
  if(is.null(y_lower)) y_lower = y_range[1]
  if(is.null(y_upper)) y_upper = y_range[2]
  g = ggplot(gene_data, aes(x = time_point, y = expression, colour = gene)) +
    geom_point() +
    geom_line(data = sine_data, aes(x = x, y = y, colour = gene)) +
    ylim(y_lower, y_upper) +
    xlab("Time") +
    ggtitle(plot_title) +
    scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24))

  return(g)
}
#setGeneric("plot_Fit", def = plot_Fit)
setMethod("plot_Fit", signature(object = "sineFitData"), plot_Fit)


# Function to plot gene fits for given periods
# N.B. Works even if given period(s) are not in the interval used when generating the object
# Intended to allow flexibility and time saving when investigating how fits look with the data
#' @name plot_PeriodFits
#' @title plot_PeriodFits
#' @export
#' @param object A sineFitData object
#' @param gene Gene to plot periods for
#' @param periods Periods to plot best fits for
#' @param y_lower Optional lower y-axis limit
#' @param y_upper Optional upper y-axis limit
#' @param plot_title Optional plot title
plot_PeriodFits = function(object, gene, periods, y_lower = NULL, y_upper = NULL, plot_title = NULL){
  # object = wk00.L.sf.24; gene = "BMAL1"; periods = c(18, 24, 30); y_lower = NULL; y_upper = NULL; plot_title = NULL

  # Designed to work for a single gene, so check argument
  stopifnot("The gene argument of plot_Periods() should be a single gene" = length(gene) == 1)

  # Get gene data
  in_data = slot(object, "Input_data") %>% select(1, all_of(gene))
  rowData = slot(object, "Row_data")

  # Get gene values for plot
  time_data = rowData %>% pull(2)
  # These are in the same order (checked when creating the sinFit object)
  gene_data = mutate(in_data, time_point = time_data)
  colnames(gene_data)[2] = "expression"

  # Get list of sineFit objects for requested periods
  sf_list = mapply(sineFit, periods, periods, MoreArgs = list(in_data = in_data, rowData = rowData))
  # Get sineFit results
  fit_vars = lapply(sf_list, results.sineFit) %>% bind_rows()  ### should this be results or results.sineFit?
  # Get points for sine waves
  x = seq(0, 24, length.out = 288)
  sine_points = mapply(get_sine_points,
                       fit_vars$Base, fit_vars$Amp, fit_vars$Period, fit_vars$Phase_Shift,
                       MoreArgs = list(x = x)) %>% as.data.frame()
  colnames(sine_points) = paste0("period_", round(periods, 2))
  sine_points$x = x
  sine_data = pivot_longer(sine_points, cols = starts_with("period_"), names_to = "Period", values_to = "y") %>%
    mutate(Period = gsub("period_", "", Period)) %>% mutate_at("Period", as.factor)

  # Make plot:
  # Set plot title if not provided
  if(is.null(plot_title)){
    plot_title = paste0(deparse(substitute(object)), " : ", gene)
  }
  # Set y-axis limits if not provided
  y_range = range(c(gene_data$expression, sine_data$y))
  if(is.null(y_lower)) y_lower = y_range[1]
  if(is.null(y_upper)) y_upper = y_range[2]
  g = ggplot(gene_data, aes(x = time_point, y = expression)) +
    geom_point() +
    geom_line(data = sine_data, aes(x = x, y = y, colour = Period)) +
    ylim(y_lower, y_upper) +
    xlab("Time") +
    ggtitle(plot_title) +
    scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24))

  return(g)
}
#setGeneric("plot_PeriodFits", def = plot_PeriodFits)
setMethod("plot_PeriodFits", signature(object = "sineFitData"), plot_PeriodFits)

