compute_regression_slope = function(pt_list, plot_only) {
  f_st = pt_list$f_st
  t = pt_list$t
  lin_model = lm(f_st ~ t)
  a = lin_model$coefficients[[1]] # Intercept
  b = lin_model$coefficients[[2]] # Slope
  reg_RMSE = sqrt(mean(lin_model$residuals^2))

  naive_intercept = head(f_st, 1)
  naive_slope = (tail(f_st, 1) - head(f_st, 1))/tail(t, 1)
  f_predictions = naive_intercept + naive_slope*t
  naive_RMSE = RMSE(f_st, f_predictions)

  max_idx = head(which(f_st == max(f_st)), 1)
  min_idx = head(which(f_st == min(f_st)), 1)

  if (min_idx == max_idx) {
    min_max_slope = 0
    min_max_intercept = 0
  } else {
    if (min_idx < max_idx) {
      f_start = f_st[min_idx]
      f_delta = f_st[max_idx] - f_st[min_idx]
      t_end = t[max_idx]
      t_start = t[min_idx]
    } else {
      f_start = f_st[max_idx]
      f_delta = f_st[min_idx] - f_st[max_idx]
      t_start = t[max_idx]
      t_end = t[min_idx]
    }
    t_delta = t_end - t_start
    min_max_slope = f_delta/t_delta
    min_max_intercept = f_start - t_start*min_max_slope
    min_max_predictions = min_max_intercept + min_max_slope*t
    min_max_RMSE = RMSE(f_st, min_max_predictions)
  }


  if (plot_only) {
    plot(t, f_st)
    abline(lin_model)
    lines(t, f_predictions, 'col' = 'red')
    lines(t, min_max_predictions, 'col' = 'blue')
    title("Slopes compared")
    labels = c(
      paste('Slope naive, RMSE:', round(naive_RMSE, 1)),
      paste('Slope regression, RMSE:', round(reg_RMSE, 1)),
      paste('Slope min/max, RMSE:', round(min_max_RMSE, 1))
    )
    colors = c("red", "black", "blue")
    legend(min(t), max(f_st), legend=labels, col=colors, lwd = 2, cex = 0.8, lty = 1)
  } else {
    return(data.frame(
      name = pt_list$name,
      a = a,
      b = b,
      reg_RMSE = reg_RMSE,
      naive_intercept = naive_intercept,
      naive_slope = naive_slope,
      naive_RMSE = naive_RMSE,
      min_max_slope = min_max_slope,
      min_max_RMSE = min_max_RMSE
    ))
  }
}

compute_regression_slopes = function(pt_file_names, pt_path, plot_only = FALSE) {
  results = NULL
  for (file_name in pt_file_names) {
    pt_list = get_pt_list(pt_path, file_name)
    row = compute_regression_slope(pt_list, plot_only)
    if (!plot_only){
      results = rbind(results, row)
    }
  }

  if (!plot_only){
    return(results)
  }
}

compute_EAC_ICCs = function(pt_file_names, pt_path, plot_only = FALSE){
  require(dplyr)
  min_gap_bwn_ICCs = 5 # at least 5 points no pitch tracking
  min_ICCs_length = 5 # Each ICCs must contain at least 5 points

  results = NULL
  for (file_name in pt_file_names){
    pt_list = get_pt_list(pt_path, file_name)
    f_st = pt_list$f_st
    t = pt_list$t

    t_max_idx = length(t)
    durations = t[2:t_max_idx] - t[1:(t_max_idx - 1)]
    samp_dur = min(durations)
    ICCs_gap = 5*samp_dur
    ICCs_borders = which(durations >= ICCs_gap)
    # if (plot_only){
    #   plot(t, f_st)
    #   abline(v = t[ICCs_borders])
    #   title("ICC borders")
    # }
    if (length(ICCs_borders) > 0){
      ICCs = NULL
      start_idx = 1
      ICC_idx = 1
      for (border in ICCs_borders){
        end_idx = border
        if (end_idx < start_idx) {
          stop(paste("Start index", start_idx, "must always be >= end index", end_idx))
        }

        if (end_idx - start_idx + 1 >= min_ICCs_length){
          idxs = start_idx:end_idx
          ICCs = rbind(ICCs, data.frame(
            t = t[idxs],
            f = f_st[idxs],
            ICC_idx = ICC_idx
          ))
          ICC_idx = ICC_idx + 1
        } else{
          warning(paste("ICC size to small", end_idx - start_idx, "SKIP current ICC"))
        }
        start_idx = end_idx + 1
      }
      if (border != t_max_idx){
        idxs = start_idx:t_max_idx
        ICCs = rbind(ICCs, data.frame(
          t = t[idxs],
          f = f_st[idxs],
          ICC_idx = ICC_idx
        ))
      }
    } else {
      warning("No ICC detected")
      if (length(f_st) > 0){
        ICCs = data.frame(t = t, f = f_st, ICC_idx = 1)
      } else {
        stop(paste("The file", file_name, "contains no pitch!"))
      }
    }

    ICCs$ICC_idx = as.factor(ICCs$ICC_idx)

    if (plot_only){
      library(ggplot2)
      print(ggplot(ICCs) +
        geom_point(aes(x = t, y = f), color = 'grey', data = data.frame(t = t, f = f_st)) +
        geom_point(aes(x = t, y = f, color = ICC_idx)) +
        ggtitle(paste("ICCs for", pt_list$name)) +
        xlab("Time in seconds") +
        ylab("F0 in semitones") +
        theme_minimal())
    }

    perc_excluded_points = ((length(t)-nrow(ICCs))/length(t))*100
    ICC_idxs = unique(ICCs$ICC_idx)
    num_ICCs = length(ICC_idxs)

    slopes = c()
    delta_F0s = c() # fmax - fmin
    tcs = c() # durations
    delta_tcs = c() # absolute distance in time between max and min F0

    for (idx in ICC_idxs){
      ICC = dplyr::filter(ICCs, ICC_idx == idx)
      #slope = lm(ICC$f ~ ICC$t)$coefficients[[2]]
      slope = (tail(ICC$f, 1) - head(ICC$f, 1))/(tail(ICC$t, 1) - head(ICC$t, 1))
      slopes = c(slopes, slope)

      delta_F0s = c(delta_F0s, max(ICC$f) - min(ICC$f))
      tcs = c(tcs, max(ICC$t) - min(ICC$t))
      delta_tc = abs(ICC$t[ICC$f == max(ICC$f)] - ICC$t[ICC$f == min(ICC$f)])
      delta_tcs = c(tcs, delta_tc)
    }

    num_pos_ICCs = length(which(slopes > 0))
    num_neg_ICCs = length(which(slopes < 0))

    if (num_pos_ICCs > 0){
      EAC = 0
      for (idx in which(slopes > 0)){
        EAC = EAC + delta_F0s[idx] * delta_tcs[idx]
      }
      t1 = head(ICCs$t, 1)
      t2 = tail(ICCs$t, 1)
      EAC = EAC/(mean(ICCs$f)*(t2-t1))
    } else {
      EAC = NA
    }


    results = rbind(results, data.frame(
      name = pt_list$name,
      num_ICCs = num_ICCs,
      num_pos_ICCs = num_pos_ICCs,
      num_neg_ICCs = num_neg_ICCs,
      mean_delta_F0 = mean(delta_F0s),
      sd_delta_F0 = sd(delta_F0s),
      mean_tc = mean(tcs),
      sd_tc = sd(tcs),
      mean_delta_tc = mean(delta_tcs),
      sd_delta_tc = sd(delta_tcs),
      EAC = EAC
    )
    )
  }

  if (!plot_only){
  return(results)
  }
}
