superposition_by_word = function(filenames, pt_path, tier_path, grouping_list){
  # TODO
  library(stringr)
  require(dplyr)
  #require(factoextra)
  if (!(is.list(grouping_list) && length(grouping_list) > 0)){
    stop("grouping_list needs to be a non empty list")
  } else{
    list_lengths = c()
    for (g in grouping_list){
      list_lengths = c(list_lengths, length(g))
    }

    if (sd(list_lengths) != 0){
      stop("Each sentence must have the same number of words")
    }
  }

  duration_df = NULL
  pts = list()

  # Check all words contain F0
  for (filename in filenames){
    print(filename)
    pt = read_PitchTier(paste0(pt_path, filename))
    filebase = strsplit(filename, "\\.")[[1]][1]

    # tg = read_TextGrid(paste0(tg_path, filebase, ".TextGrid"))
    # tier = tg$words
    tier = read_TextGrid_Tier(paste0(tier_path, filebase, ".csv"))
    tg_word_idxs = which(tier$label != "")
    t1_ref = tier$t1[tg_word_idxs]
    t2_ref = tier$t2[tg_word_idxs]
    labels = tier$label[tg_word_idxs]

    s = as.numeric(str_extract(strsplit(filebase, "_")[[1]][3], "(\\d)+"))

    word_idxs = grouping_list[[s]]
    word_num = 0
    for (i in word_idxs){
      word_num = word_num + 1
      pt_idxs = which(pt$t >= t1_ref[i] & pt$t < t2_ref[i])
      if (length(pt_idxs) > 1){
        # Must contain at least two pitch points
        duration_df = rbind(duration_df, data.frame(filename = filebase, word_idx = i, word_num = word_num, label = labels[i], duration = max(pt$t[pt_idxs]) - min(pt$t[pt_idxs]), num_points = length(pt_idxs)))
      } else{
        warning(paste("The word", labels[i], i, "in", filebase, "contains no F0!\n"))
      }

      l = list(t = pt$t[pt_idxs], f = pt$f[pt_idxs], filename = filebase, label = labels[i])

      if (is.null(pts[word_num][[1]])){
        pts[[word_num]] = list(l)
      } else{
        pts[[word_num]] = append(pts[[word_num]], list(l))
      }
    }
  }


  results = NULL
  avg_dur = mean(duration_df$duration)
  duration_df$compression_rate = avg_dur/duration_df$duration
  num_points = median(duration_df$num_points)
  tiny_padding = 0.001 # 1 ms padding to avoid NA
  t = seq(0, (avg_dur - 0.001), length=num_points)

  for (w in 1:length(pts)){
    pt_list = pts[[w]]
    sel_df = dplyr::filter(duration_df, word_num == w)

    for (pt in pt_list){
      if(length(pt$f) < 2){
        results = rbind(results, c(rep(NA, num_points), compression_rate, w, pt$filename, pt$label))
      } else {
        compression_rate = sel_df[sel_df$filename == pt$filename, "compression_rate"]
        # Start at 0
        pt$t = pt$t - min(pt$t)

        # Scale it
        pt$t = pt$t * compression_rate

        # Resample it
        interpol = approxfun(pt$t, pt$f)
        f = interpol(t)

        # pitch to st
        f = f2st(f)

        if (any(is.na(f))){
          stop("This may not happen!")
        }
        result = data.frame(t(c(f, compression_rate)))
        names(result) =  c(paste0('p', 1:num_points), "compression_rate")
        result$word_num = w
        result$filename = pt$filename
        result$label = pt$label
        results = rbind(results, result)
      }
    }
  }
  results = as.data.frame(results)

  numeric_columns = names(results)[1:(length(names(results))-3)]
  for (col in numeric_columns){
    results[[col]] = as.numeric(results[[col]])
  }
  results$word_num = as.factor(results$word_num)
  return(results)

}

eigenshape_analysis = function(df, method, compression_col_name = "compression_rate"){
  df = add_meta_data(df)
  idxs = grep('p[0-9]', names(df))
  num_landmarks = length(idxs)
  if (method == "Eigenshape_ZR"){
    num_angles = num_landmarks - 2
  } else {
    num_angles = num_landmarks - 1
  }
  if (length(idxs) == 0){
    stop("DF needs to have at least one point")
  }
  if (!compression_col_name %in% names(df)){
    stop(paste("DF needs to have a column named", compression_col_name))
  }

  results = NULL
  for (i in 1:nrow(df)){
    speaker = df$speaker[i]
    if (any(speaker %in% c("DF","MG"))){
      pitch_floor = 70
    } else {
      pitch_floor = 120
    }

    measures_size = 0.75/pitch_floor
    xs = seq(0, (num_landmarks-1)*measures_size, measures_size)
    xs = xs * df[i, compression_col_name]
    ys = as.numeric(df[i, idxs])

    angles = xyToPhi(xs, ys, method)
    result = data.frame(t(angles))
    names(result) =  c(paste0('p', 1:(num_angles)))
    result$word_num = df$word_num[i]
    result$filename = df$filename[i]
    result$label = df$labeldf$word_num
    results = rbind(results, result)
  }

  return(results)
}

xyToPhi = function(xs, ys, method, move_right_to_left = FALSE){
  #' This generates the phi function from x, y data
  #' As currently written, phi goes from the second point (after 0) to -2 π
  #' This follows Zahn and Roskies 1972, but introduces redundant information
  #'
  #' @param xs: time, regular intervals
  #' @param ys: F0, in semitones
  #' @param method: either "angle_to_angle" or "reference_angle", angle to angle computes between two succesive angles, reference angle implements Zahn & Roskies, 1972

  if (length(xs) != length(ys)){
    stop("xs and ys must be of the same length")
  }

  if (!method %in% c("angle_to_angle", "reference_angle", "Eigenshape_ZR")){
    stop("Only the methods 'angle_to_angle', 'reference_angle' (as in Zahn & Roskies, 1972) and Eigenshape as described here https://www.palass.org/publications/newsletter/palaeomath-101/palaeomath-part-24-centre-cannot-hold-i-z-r-fourier-analysis")
  }

  num_points = length(xs)

  # Array for the computed angles
  angle_degrees = c()

  # Boolean indicating the method
  is_angle_2_angle = method == "angle_to_angle"
  is_angle_2_ref = method == "reference_angle"
  is_Eigenshape_ZR = method == "Eigenshape_ZR"

  total_iterations = num_points - 1
  # ZR uses a fixed reference point, in this case the first landmark
  if (is_angle_2_ref){
    # Variable to save rotation
    angle_rotation = 0 # actually not needed for open outlines

    # Get the reference
    baseline_angle = (atan2(ys[2] - ys[1], xs[2] - xs[1])*180)/pi
  } else if (is_Eigenshape_ZR){
    compute_c = function(i, xs, ys){
      dx_1 = xs[i+1] - xs[i]
      dy_1 = ys[i+1] - ys[i]

      dx_2 = xs[i+2] - xs[i+1]
      dy_2 = ys[i+2] - ys[i+1]

      c = ((dx_1*dx_2) + (dy_1*dy_2))/sqrt((dx_1^2 + dy_1^2) * (dx_2^2 + dy_2^2))
      return(c)
    }

    c1 = compute_c(1, xs, ys)
    total_iterations = num_points - 2
  }

  for (i in 1:total_iterations){
    if (is_angle_2_angle){
      # Recompute the angle and slope for each point
      angle = (atan2(ys[i] - ys[i+1], xs[i] - xs[i+1])*180)/pi
    } else if (is_Eigenshape_ZR){
      c = compute_c(i, xs, ys)

      if (c^2 >= 1){
        s = 0
      } else{
        s = sqrt(1-(c^2))
      }

      angle = (atan2(s, c1)*180)/pi

      dx_1 = xs[i+1] - xs[i]
      dy_1 = ys[i+1] - ys[i]

      dx_2 = xs[i+2] - xs[i+1]
      dy_2 = ys[i+2] - ys[i+1]

      if (((dx_1*dy_2) - (dx_2*dy_1)) < 0){
        angle = angle * -1
      }

    } else {
      if (is_angle_2_ref){
        # The reference is the first landmark
        # Update each iteration reference slope and angle
        angle = (atan2(ys[i+1] - ys[i], xs[i+1] - xs[i])*180)/pi


        # In the case of 180 degrees, you are probably at the other side of the shape
        # E.g. -180 + -90 degrees => -270, not 90 degrees
        angle_rounded = round(angle, 1)
        if (is.na(angle_rounded)){
          print(paste(ys[i+1], ys[i], xs[i+1], xs[i]))
        }
        if (abs(floor(angle_rounded)) == 180){
          # Check to flip direction
          if (xs[i+1] > xs[i]){
            # moving from left to right
            if (move_right_to_left){
              angle_rotation = angle_rotation + 180
            }
            move_right_to_left = FALSE
          } else if (xs[i+1] < xs[i]){
            # moving from right to left
            if (!move_right_to_left){
              angle_rotation = angle_rotation - 180
            }
            move_right_to_left = TRUE
          }
          # Avoid tha
          if (move_right_to_left){
            if (angle_rotation == 0){
              multiplyer = 0
            } else{
              multiplyer = round(abs(180/angle_rotation)) - 1
            }
            if (angle_rotation < 0){
              angle = multiplyer*-180 - abs(angle)
            } else {
              angle = multiplyer*180 + abs(angle)
            }
          }
        }
      } else{
        # Not 180 degrees, e.g. 78 degrees
        if (move_right_to_left){
          angle = angle_rotation - angle
        }
      }

      angle = angle - baseline_angle
    }
    angle = round(angle, 1)
    angle_degrees = c(angle_degrees, angle)
  }

  return(angle_degrees)
}

plot_shape_function = function(xs, ys, angles, method,
                               save = FALSE, plot3.y_start = 0,
                               plot3.y_lim = range(-360, 360),
                               plot3.vlines = c(),
                               plot3.use_radians = FALSE,
                               plot3.x_ticks_at = c (),
                               plot3.x_tick_labels = c(),
                               prefix_name = 'plot'){
  # Boolean indicating the method

  ref_color = "#f2c62c"
  comp_color = "#3695d2"
  move_right_to_left = FALSE
  is_angle_2_angle = method == "angle_to_angle"
  is_angle_2_ref = method == "reference_angle"

  num_landmarks = length(xs)

  if (is_angle_2_ref){
    # Get the reference for plotting
    m1 = (ys[2]-ys[1])/(xs[2]-xs[1])
  }

  for (i in 1:(num_landmarks - 1)){
    angle = angles[i]
    # Compute
    if (is_angle_2_angle){
      # Recompute the angle and slope for each point
      m1 = (ys[i+1]-ys[i])/(xs[i+1]-xs[i])
    } else if (is_angle_2_ref){
      # The reference is the first landmark
      # Update each iteration reference slope and angle
      m2 = (ys[i+1]-ys[i])/(xs[i+1]-xs[i])
    }

    if (save){
      pdf(paste0(prefix_name, "_", i,".pdf"))
    }
    layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))

    # Plot shape
    plot(xs, ys)
    if (is_angle_2_ref){
      # Reference
      points(xs[1:2], ys[1:2], col = ref_color, type='l', lwd=2)
      # Moving points
      points(xs[i:(i+1)], ys[i:(i+1)], col = comp_color, type='l', lwd=2)

      if (i <= 2){
        points(xs[1:(i+1)], ys[1:(i+1)], pch = 19)
      } else{
        points(xs[i:(i+1)], ys[i:(i+1)], pch = 19, col=comp_color)
        points(xs[1:2], ys[1:2], pch = 19, col=ref_color)
      }
    } else if (is_angle_2_angle){
      # Moving points
      points(xs[i], ys[i], pch = 19, col = ref_color)
      points(xs[i+1], ys[i+1], pch = 19, col = comp_color)
    }
    title = paste("Comparing point", i, "and", i+1)
    if (move_right_to_left){
      title = paste0(title, "\nMoving from right to left")
    } else {
      title = paste0(title, "\nMoving from left to right")
    }
    title(title)

    # Plot angle
    t = seq(-1, 1, 0.01)
    plot(t, (m1*t + 0.05), col=comp_color, type='l', ylim=range(-1, 1), xaxt='n', yaxt='n', ann=FALSE, lwd=2)
    abline(v=0, col=ref_color, lwd=2)
    title(paste("Angle between lines:", angle, "degrees"))

    #plot(NA,NA, xlim=c(-1,1), ylim=c(-1,1))


    degrees = (angles + plot3.y_start)
    if (plot3.use_radians){
      degrees = (degrees/180)*pi
      ylab = "Angles in radians"
    } else {
      ylab = "Angles in degrees"
    }

    # Plot angle development
    plot(
      1:i,
      degrees[1:i],
      xlim=range(1:num_landmarks),
      ylim=plot3.y_lim,
      ylab=ylab,
      xlab="Iterations",
      axes=FALSE,
      pch = 19
    )
    title(paste("Current iteration:", i))
    if (length(plot3.vlines) > 0){
      abline(v = plot3.vlines)
    }

    if (length(plot3.x_ticks_at) > 0){
      if (length(plot3.x_tick_labels) > 0){
        axis(side = 2, at = plot3.x_ticks_at, plot3.x_tick_labels)
      } else {
        axis(side = 2, at = plot3.x_ticks_at)
      }
    } else{
      axis(side = 2)
    }

    axis(side = 1)

    if (save){
      dev.off()
    }
  }
}

# Methods on combined measures
pca_analysis = function(results, plot = TRUE, title_prefix="", scale = TRUE, prefix = "", center = TRUE, compression_rate_col = "compression_rate", colors = c(), labels = c(), return_plots = FALSE){
  if (!compression_rate_col %in% names(results)){
    stop(paste("DF must contain the following column:", compression_rate_col))
  }

  # Add some meta data used for plotting
  results = add_meta_data(results)

  # Remove missing values from analysis
  results_wo_NA = na.omit(results)
  filenames = results_wo_NA$filename
  emotions = results_wo_NA$emotion

  # Add the points + compression to the analysis
  pca_col_idxs = 1:(max(grep("p[0-9]", names(results))))
  pca_col_idxs = c(pca_col_idxs, which(names(results) == compression_rate_col))

  pca_df = results_wo_NA[,pca_col_idxs]

  # Do PCA
  pca = prcomp(pca_df, center = center, scale. = scale)

  pc_df = data.frame(emotion = emotions, PC1 = pca$x[,2], PC2 = pca$x[,2])

  # Compute explained variance
  eig = (pca$sdev)^2
  variance = eig*100/sum(eig)
  cumvar = cumsum(variance)

  if (plot) {
    if (title_prefix != "") {
      title_prefix = paste(title_prefix, "")
    }

    # Compute total variance
    PC1_variance = cumvar[1]
    PC2_variance = round(cumvar[2] - PC1_variance, 1)
    PC1_variance = round(PC1_variance, 1)

    # PCA plot
    library(ggbiplot)
    library(ggplot2)
    emotions = results_wo_NA$emotion
    if (length(labels) == length(levels(emotions))) {
      levels(emotions) = labels
    }
    pca_plot = ggbiplot(pca, obs.scale = 1, var.scale = 1, groups = emotions,
                 ellipse = FALSE, circle = FALSE, varname.size=0, var.axes = F) +

      #geom_point(aes(colour=emotions), size = 1) +
      ggtitle(paste0(title_prefix, "PCA for all words together (num points: ", max(pca_col_idxs) - 1, ")")) +
      theme_minimal() +
      theme(legend.title = element_blank())
    if (length(colors) == length(levels(emotions))){
      pca_plot = pca_plot + scale_color_manual(values=colors)
    }
    if (plot & !return_plots) {
      print(pca_plot)
    }
    # print(
    #   ggplot(pc_df) +
    #                 geom_point(aes(x = PC1, y = PC2, colour=emotions), size = 1) +
    #                 xlab(paste("PC1 (explains", PC1_variance, "% of variance)")) +
    #                 ylab(paste("PC2 (explains", PC2_variance, "% of variance)")) +
    #                 ggtitle(paste0(title_prefix, "PCA for all words together (num points: ", max(pca_col_idxs) - 1, ")")) +
    #                 theme_minimal()
    # )



    # Display variable contributions to first PC
    var = factoextra::get_pca_var(pca)
    most_variance = abs(var$coord[, 1])
    threshold = .5
    most_variance = most_variance[most_variance > threshold]
    variance_df = data.frame(value=as.numeric(most_variance), name = names(most_variance))
    variance_df$name = factor(variance_df$name, levels = variance_df$name)
    variance_plot = ggplot(variance_df) +
      geom_bar(aes(x=name, y = value), fill="#02a3dd", stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90)) +
      ggtitle(paste0(title_prefix, "Contribution of each variable to PC1 (explains ", PC1_variance, " % of all variance)")) +
      xlab(paste("All variables contributing more than", threshold)) +
      ylab("Contributions of the variables")

    if (plot & !return_plots) {
      print(variance_plot)
    }
  }
  if (prefix != ""){
    prefix = paste0(prefix, "_")
  }
  if (length(unique(results_wo_NA$word_num)) > 1){
    PC1 = data.frame(PC1 = pca$x[,1], filename = filenames, word_num = paste0("PC1_", prefix, results_wo_NA$word_num))
    PC1 = tidyr::spread(PC1, word_num, PC1)

    PC2 = data.frame(PC2 = pca$x[,2], filename = filenames, word_num = paste0("PC2_", prefix, results_wo_NA$word_num))
    PC2 = tidyr::spread(PC2, word_num, PC2)

    features = merge(PC1, PC2, by = "filename")
  } else {
    features = data.frame(PC1 = pca$x[,1],PC2 = pca$x[,2], filename = filenames)
    names(features)[1:2] = c(paste0(prefix, "PC1"), paste0(prefix, "PC2"))
  }
  returned_list = list(cumvar = cumvar, pca = pca, features = features)
  if (return_plots) {
    returned_list$variance_plot = variance_plot
    returned_list$pca_plot = pca_plot
  }
  return(returned_list)
}

plot_morphospace = function(pca_list, results, row_idxs,
                            baseline_row_idx = NULL, baseline_label = "baseline",
                            number_of_lines = 10, relative_width = 0.2, relative_height = 0.2,
                            colors = NULL, labels = NULL, save_to = NULL,
                            title = "PCA morphospace"){
  if (length(row_idxs) == 0){
    stop("Need to specify row indexes")
  }

  if (!is.null(colors) && length(colors) != length(row_idxs)){
    stop("Color needs to be specified for each row")
  }

  if (!is.null(labels)){
    labels = as.character(labels)
    if(length(labels) != length(row_idxs)){
      stop("Color needs to be specified for each row")
    }
  }

  pca_scores = pca_list$pca$x

  pc1 = pca_scores[,1]
  pc2 = pca_scores[,2]
  xlim = range(pc1)
  ylim = range(pc2)

  max_width = diff(xlim)*relative_width
  max_height = diff(ylim)*relative_height
  scale_to = min(max_width, max_height)

  if (is.null(baseline_row_idx)){
    baseline_row_idx = head(which((floor(pc1) == 0) & (floor(pc2) == 0)), 1)
    x_ref = 0
    y_ref = 0
  } else{
    x_ref = pc1[baseline_row_idx]
    y_ref = pc2[baseline_row_idx]
  }

  # Compute arbitrary time
  col_idxs_points = grep("p[0-9]", names(results))

  # min_val = min(as.numeric(results[col_idxs_points, 1]))
  # max_val = max(as.numeric(results[col_idxs_points, 1]))
  #
  # for (col_idx in col_idxs_points){
  #   test_min = min(as.numeric(results[col_idxs_points, col_idx]))
  #   test_max = max(as.numeric(results[col_idxs_points, col_idx]))
  #
  #   if (test_max > max_val){
  #     max_val = test_max
  #   }
  #   if (test_min < min_val){
  #     min_val = test_min
  #   }
  # }

  t = col_idxs_points - 1 # e.g. 0, 1, 2, 3, 4
  compression_col_idx = max(col_idxs_points) + 1

  t_ref = normalize(results[baseline_row_idx,compression_col_idx]*t)
  #t_ref = normalize(t)
  f_ref = normalize(as.numeric(results[baseline_row_idx, col_idxs_points]))
  ref = as.matrix(data.frame(t_ref, f_ref))

  grid_list = compute_grid(ref, ref, number_of_lines, scale_to, x_ref, y_ref)

  # General plot in PCA space
  cumvar = pca_list$cumvar
  xlab = paste0("PC1 (", round(cumvar[1], 1), "%)")
  ylab = paste0("PC2 (", round(cumvar[2]-cumvar[1], 1), "%)")

  if (!is.null(save_to)){
    cairo_pdf(save_to)
  }

  plot(x = c(), xlim=xlim, ylim=ylim, ylab = ylab, xlab = xlab)
  abline(v=0, col="#eeeeee")
  abline(h=0, col="#eeeeee")
  title(title)

  plot_grid(grid_list)

  c = 0
  for (idx in row_idxs){
    c = c + 1
    if (idx != baseline_row_idx){
      x = pc1[idx]
      y = pc2[idx]
      t_comp = normalize(results[idx, compression_col_idx]*t)
      f_comp = normalize(as.numeric(results[idx, col_idxs_points]))
      comp = as.matrix(data.frame(t_comp, f_comp))
      grid_list = compute_grid(ref, comp, number_of_lines, scale_to, x, y)
      if (!is.null(colors)){
        color = colors[c]
      } else {
        color = "#424143"
      }
      plot_grid(grid_list, color = color)
    }
  }

  if (!(is.null(labels) && is.null(colors))){
    if (length(which(row_idxs == baseline_row_idx)) == 1){
      idx = which(row_idxs == baseline_row_idx)
      colors[idx] = "grey"
    } else{
      colors = c("grey", colors)
      labels = c(baseline_label, labels)
    }
    legend(min(pc1), max(pc2), legend=labels, col=colors, lwd=2, cex=0.8, lty=1)
  }

  if (!is.null(save_to)){
    dev.off()
  }
}


plot_point_to_morphospace = function(results, row_idxs, colors = NULL, save_to = NULL){
  if (length(row_idxs) != 2){
    stop("Need to specify 2 row indexes")
  }

  if(!is.null(colors) && length(row_idxs) != length(colors)){
    stop("Each row needs to have a color")
  }

  if (is.null(colors)){
    c1 = "#1aa5df"
    c2 = "#414042"
  } else{
    c1 = colors[1]
    c2 = colors[2]
  }

  # Compute indexes of points
  col_idxs_points = grep("p[0-9]", names(results))

  t = col_idxs_points - 1 # e.g. 0, 1, 2, 3, 4
  t_norm =  normalize(t)
  compression_col_idx = max(col_idxs_points) + 1

  #t_ref = normalize(results[row_idxs[1],compression_col_idx]*t)
  f_ref = normalize(as.numeric(results[row_idxs[1], col_idxs_points]))
  ref = as.matrix(data.frame(t_norm, f_ref))

  #t_comp = normalize(results[row_idxs[2], compression_col_idx]*t)
  f_comp = normalize(as.numeric(results[row_idxs[2], col_idxs_points]))
  comp = as.matrix(data.frame(t_norm, f_comp))

  if (!is.null(save_to)){
    cairo_pdf(save_to)
  }
  layout(matrix(c(1,3,2,4), 2, 2, byrow = TRUE))

  n_grid = compute_grid(ref, ref, 20, scale = FALSE, translate = FALSE)
  plot_grid(n_grid, new_plot = TRUE)
  lines(ref[col_idxs_points,],lty=3,lwd=4, col=c2)
  title('Reference')

  n_grid = compute_grid(comp, comp, 20, scale = FALSE, translate = FALSE)
  plot_grid(n_grid, new_plot = TRUE)
  points(comp, asp=1,pch=20,col=c1)
  lines(comp[col_idxs_points,],lwd=2, col=c1)
  title('Comparison')

  plot(comp, asp=1,pch=20,col=c1, ylim=range(-1.5, 1.5), xlim=range(-1.5, 1.5), ylab="", xlab="")
  lines(comp,lwd=2,col=c1)
  lines(ref,lty=3,lwd=2, col=c2)
  # Draw arrows
  for (i in 1:dim(comp)[1]){
    arrows(ref[i,1],ref[i,2],comp[i,1],comp[i,2],length=0.05, lwd=1.2, col="grey50")
  }
  title("Transforming reference to comparison")

  n_grid = compute_grid(ref, comp, 20, scale_to = 3, translate = FALSE)
  plot_grid(n_grid, new_plot = TRUE, color = c2)
  title("Deformation grid")

  if (!is.null(save_to)){
    dev.off()
  }
}

compute_grid<-function(matr, matt, n, scale_to, x, y, scale = TRUE, translate = TRUE){
  # Define the range of the graph to estimate the size of the grid.
  xm<-min(matt[,1])
  ym<-min(matt[,2])
  xM<-max(matt[,1])
  yM<-max(matt[,2])
  rX<-xM-xm
  rY<-yM-ym

  # Calculate the coordinates of the intersections between the lines of the grid.
  a<-seq(xm-1/5*rX, xM+1/5*rX, length=n)
  b<-seq(ym-1/5*rX, yM+1/5*rX, by=(xM-xm)*7/(5*(n-1)))
  m<-round(0.5+(n-1)*(2/5*rX+ yM-ym)/(2/5*rX+ xM-xm))
  # Baseline grid
  M<-as.matrix(expand.grid(a,b))
  # Grid
  n_grid = tps2d(M,matr,matt)
  if (scale){
    scalar = scale_to/diff(range(n_grid[,2]))
    n_grid = n_grid*scalar
  }
  if (translate){
    n_grid[,1] = n_grid[,1] + x
    n_grid[,2] = n_grid[,2] + y
  }
  return(list(grid = n_grid, n = n, m = m))
}

plot_grid = function(grid_list, color = "grey", new_plot = FALSE){
  n_grid = grid_list$grid
  n = grid_list$n
  m = grid_list$m
  if (new_plot){
    plot(n_grid, cex=0.2,asp=1, col=color, ylab="", xlab="")
  } else {
    points(n_grid, cex=0.2,asp=1, col=color)
  }
  for (i in 1:m){
    lines(n_grid[(1:n)+(i-1)*n,], col=color)
  }
  for (i in 1:n){
    lines(n_grid[(1:m)*n-i+1,], col=color)
  }
}


# producing 2D grids of deformation according to our needs
tps2d<-function(M, matr, matt){
  p<-dim(matr)[1]
  q<-dim(M)[1]
  n1<-p+3
  P<-matrix(NA, p, p)
  for (i in 1:p){
    for (j in 1:p){
      r2<-sum((matr[i,]-matr[j,])^2)
      P[i,j]<- r2*log(r2)
    }
  }
  P[which(is.na(P))]<-0
  Q<-cbind(1, matr)
  L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,3,3)))
  m2<-rbind(matt, matrix(0, 3, 2))
  coefx<-solve(L)%*%m2[,1]
  coefy<-solve(L)%*%m2[,2]
  fx<-function(matr, M, coef){Xn<-numeric(q)
  for (i in 1:q){
    Z<-apply((matr-matrix(M[i,], p, 2, byrow=T))^2, 1, sum)
    Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*
      M[i,2]+sum(coef[1:p]*(Z*log(Z)))
  }
  Xn
  }
  matg<-matrix(NA, q, 2)
  matg[,1]<-fx(matr, M, coefx)
  matg[,2]<-fx(matr, M, coefy)
  matg
}


normalize = function(X, from = -1, to = 1, max_val = NULL, min_val = NULL){
  #' Normalize into -1 to 1 space
  #' @param X, numeric array
  #' @param from starting at
  #' @param to ends at
  #'

  if (is.null(min_val)){
    min_val = min(X)
  }

  if (is.null(max_val)){
    max_val = max(X)
  }

  return(sum(abs(from), abs(to))*((X-min_val)/(max_val - min_val)) + from)
}

equaly_space = function(X, num){
  #' Get equally spaced landmarks
  #' @param x, numeric array
  #' @param num, number of landmarks
  return((X[seq(1,length(X),length=(num+1))])[-1])
}


compute_minimal_harmonics = function(
  df,
  precision_threshold = 0.01, steps = c(-10:-1, 1:40), step_size = 0.2, plot = FALSE
){
  #' Computes the minimal amount of harmonics to reconstruct the shape within a certain precision threshold
  #' df,

  if (!any(c("filename", "f", "t") %in% names(df))){
    stop("Dataframe must contain the column filename")
  }

  results = NULL

  for (file in as.character(unique(df$filename))){
    smallest_harmonic = NA

    sel_df = dplyr::filter(df, filename == file)
    trajectory = sel_df$f # get pitch
    ts = sel_df$t # Get time

    X.k_backup = fft(trajectory) # FFT

    level_min = min(trajectory)
    scale_to = max(trajectory) - level_min
    RMSEs = c()

    max_freq = floor(length(X.k_backup)/2)
    for (iter in 1:max_freq){
      N   <- length(ts)
      i   <- complex(real = 0, imaginary = 1)
      x.n <- rep(0,N)           # create vector to keep the trajectory

      X.k = X.k_backup
      X.k = X.k[1:iter]
      ks  <- 0:(length(X.k)-1)

      for(n in 0:(N-1)) {       # compute each time point x_n based on freqs X.k
        x.n[n+1] <- sum(X.k * exp(i*2*pi*ks*n/N)) / N
      }

      # Get the estimation
      y_reconstructed = Re(x.n)

      # Scale them so they exactly match
      y_diff = max(y_reconstructed) - min(y_reconstructed)
      if (y_diff != 0){
        scalar = scale_to/(max(y_reconstructed) - min(y_reconstructed))
        y_reconstructed = scalar*y_reconstructed
      }

      # Do a first estimation for the y translation
      botom_level = min(y_reconstructed)
      y_reconstructed = y_reconstructed + level_min - botom_level

      # Optimalize y translation
      RMSE_check = c()
      for (i in steps){
        RMSE_check = c(RMSE_check, RMSE(y_reconstructed + i*step_size, trajectory))
      }
      optimal_y_translation = steps[min(RMSE_check) == RMSE_check]*step_size
      y_reconstructed = y_reconstructed + optimal_y_translation

      # Save the RMSE
      RMSEs = c(RMSEs, RMSE(y_reconstructed, trajectory))

      if (plot){
        plot(trajectory, type='l', col='red')
        lines(y_reconstructed, type='l')
        title(paste("Num harmonics:", iter))
      }
    }

    good_RMSEs = which(RMSEs < scale_to*precision_threshold)

    if (length(good_RMSEs) > 0){
      smallest_harmonic = min(good_RMSEs)
    }

    results = rbind(results, data.frame(filename = file, smallest_harmonic = smallest_harmonic, max_freq = max_freq, relative_harmonic = smallest_harmonic/max_freq))
  }

  return(results)
}

compute_intsint_features = function(tg_path){
  filenames = list.files(tg_path)
  filenames = filenames[grepl('.TextGrid', filenames)]

  results = NULL

  for (filename in filenames){
    base_name = str_split(filename, "\\.")[[1]][1]
    tg = read_TextGrid(paste0(tg_path, filename))
    results = rbind(results, data.frame(filename = base_name, INTSINT_count = length(tg$Intsint$label)))
  }

  return(results)
}

frequency_analysis = function(f_list, t_list, padding, foi){
  if (length(f_list) !=  length(t_list)){
    stop("f_list and t_list must be equally long!")
  }

  # Internal helper functions
  polyremoval = function(x){
    # Compute peudoinverse
    invxcov = MASS::ginv(length(x))[[1]]
    # Multiply pseudoinverse with total sum and substract it from each item in array
    beta = sum(x)*invxcov
    return(x - beta)
  }

  hanning = function(n){
    # Taken from Matlab
    # compute taper
    N   = n+1
    tap = 0.5*(1-cos((2*pi*(1:n))/N))

    # make symmetric
    halfn = floor(n/2)
    tap[(n+1-halfn):n] = rev(tap[1:halfn])

    return(tap)
  }

  # Number of frequencies of interest
  nfoi = length(foi)

  # Number of trials
  ntrail = length(f_list)

  # Empty array for all trials
  powspctrm = matrix(ncol = length(foi), nrow = ntrail)

  for (itrial in 1:ntrail){
    # Get time and F0 for each trial
    f = f_list[[itrial]]
    t = t_list[[itrial]]

    # Remove polynomial fit from the data
    f = polyremoval(f)

    # Sample rate
    fsample = 1/mean(diff(t))

    # Total samples
    ndatsample = length(f)

    # Total time in seconds of input data
    dattime = ndatsample / fsample

    # Compute padding
    postpad    = ceiling((padding - dattime) * fsample)
    endnsample = round(padding * fsample) #  total number of samples of padded data
    endtime    = padding # total time in seconds of padded data

    freqboi = round(foi * endtime) + 1
    freqboi = unique(freqboi)
    freqoi  = (freqboi-1)/endtime; # boi - 1 because 0 Hz is included in fourier output
    nfreqboi = length(freqboi)
    nfreqoi  = length(freqoi)

    # Compute hanning
    tap = hanning(ndatsample)
    tap  = tap / norm(as.matrix(tap), 'F')

    # Set ntaper
    ntaper = rep(ndatsample, nfreqoi)

    # determine phase-shift so that for all frequencies angle(t=0) = 0
    timedelay = t[1]
    if (timedelay != 0){
      angletransform = as.complex(rep(0, nfreqoi))
      for (ifreqoi in 1:nfreqoi){
        missedsamples = round(timedelay * fsample);
        # determine angle of freqoi if oscillation started at 0
        # the angle of wavelet(cos,sin) = 0 at the first point of a cycle, with sine being in upgoing flank, which is the same convention as in mtmconvol
        anglein = missedsamples * ((2*pi/fsample) * freqoi[ifreqoi])
        angletransform[ifreqoi] = atan2(sin(anglein), cos(anglein))
      }
    }
    # complex to real
    angletransform = as.numeric(angletransform)

    # Empty spectrum
    spectrum = matrix(ncol = nfreqoi, nrow = ntaper)
    for (itap in 1:ndatsample){
      # Pad data and apply window
      padded_data = c(tap*f, rep(0, postpad))

      # FFT
      dum = fft(padded_data)

      # Select FOIs
      dum = dum[freqboi]

      # Remove time delay
      if (timedelay != 0){
        dum = dum*exp(-1i*angletransform)
      }

      dum = dum * sqrt(2/endnsample)
      spectrum[itap,] = dum
    }

    # Add values to power spectrum
    for (ifoi in 1:nfoi){
      powdum = abs(spectrum[,ifoi])^2;
      powspctrm[itrial, ifoi] = mean(powdum)
    }
  }
  # First frequency needs to be devided by 2
  powspctrm[,1] = powspctrm[,1]/2

  return(powspctrm)
}

