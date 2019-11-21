read_pac_file = function(pac_filename){
  #' read data from .PAC file from FujiParaEditor.exe
  #'
  #' @param pac_filename: string
  #' @return: nsamps, samp_dur, baseline, df
  #'
  #' @examples
  #' read_pac_file('test_file.pac')

  if (!file.exists(pac_filename)) {
    stop("File does not exist!")
  } else{
    pac_file = file(pac_filename)
    open(pac_file)
    lines = readLines(pac_file, n = 11, warn = FALSE)
    close(pac_file)

    nsamps = as.numeric(lines[7])
    num_phrase = as.numeric(lines[8])
    num_accent = as.numeric(lines[9])
    baseline = as.numeric(lines[10])
    samp_dur = as.numeric(lines[11])

    start_line = 21
    end_line = start_line + num_phrase + num_accent -1

    pac_file = file(pac_filename)
    open(pac_file)
    lines = readLines(pac_file, n = end_line, warn = FALSE)
    close(pac_file)

    # Read csv at line 21, 'T1', 'T2', 'A', 'param'
    T1 = c()
    T2 = c()
    A  = c()
    param = c()

    for (line in lines[start_line:end_line]){
      values = strsplit(line, " ")[[1]]
      values = values[values != ""]
      T1 = c(T1, as.numeric(values[1]))
      T2 = c(T2, as.numeric(values[2]))
      A = c(A, as.numeric(values[3]))
      param = c(param, as.numeric(values[4]))
    }

    df = data.frame(T1, T2, A,  param)

    if (num_phrase > 0){
      phrase_idxs = 1:num_phrase
    } else {
      phrase_idxs = c()
    }

    if (num_accent > 0){
      accent_idxs = (num_phrase+1):(num_phrase + num_accent)
    } else {
      accent_idxs = c()
    }

    return(list(
      nsamps = nsamps,
      samp_dur = samp_dur,
      baseline = baseline,
      df = df,
      phrase_idxs = phrase_idxs,
      accent_idxs = accent_idxs
    ))
  }
}

get_accent_contour = function (pac_filename, force_beta = NULL, pac_list = NULL){
  #' calculate accent curve from .PAC file, output from FujiParaEditor.exe
  #'
  #' @param pac_filename: string
  #' @return: phrase_contour
  #'
  if (is.null(pac_list)){
    pac_list = read_pac_file(pac_filename)
  }
  nsamps = pac_list$nsamps
  samp_dur = pac_list$samp_dur
  baseline = pac_list$baseline
  df = pac_list$df
  tmax = nsamps * samp_dur

  fs = 100
  accent_contour = rep(0, as.integer(tmax * fs))

  requireNamespace("dplyr")
  library(dplyr)

  row_idxs = pac_list$accent_idxs
  for (row_idx in row_idxs){
    if (is.numeric(force_beta)){
      beta = force_beta
    } else {
      beta = df[row_idx,'param']
    }
    T1 = df[row_idx,'T1']
    dur =  df[row_idx,'T2'] - T1
    A =  df[row_idx,'A']

    accent_contour = accent_contour + A * Ga2(beta, T1, dur, tmax=tmax, fs=fs)
  }
  return(accent_contour)
}

get_phrase_contour = function (pac_filename, force_alpha = NULL, pac_list = NULL){
  #' calculate phrase curve from .PAC file, output from FujiParaEditor.exe
  #'
  #' @param pac_filename: string
  #' @return: phrase_contour

  if (is.null(pac_list)){
    pac_list = read_pac_file(pac_filename)
  }
  nsamps = pac_list$nsamps
  samp_dur = pac_list$samp_dur
  df = pac_list$df
  tmax = nsamps * samp_dur

  fs = 100
  phrase_contour = rep(0, as.integer(tmax * fs))

  row_idxs = pac_list$phrase_idxs
  for (row_idx in row_idxs){
    if (is.numeric(force_alpha)){
      alpha = force_alpha
    } else {
      alpha =df[row_idx,'param'] # eigenvalue of the ith phrase control mechanism
    }
    T1 = df[row_idx,'T1']
    A =  df[row_idx,'A']

    phrase_contour = phrase_contour + A * Gp(alpha, T1, tmax=tmax, fs=fs)
  }

  return(phrase_contour)
}

get_baseline = function(pac_filename, pac_list = NULL){
  if (is.null(pac_list)){
    pac_list = read_pac_file(pac_filename)
  }
  return(pac_list$baseline)
}


fuji2pitch = function(baseline, accent, phrase){
  #' Get reconstruction of pitch contour from fujisaki components
  #'
  #' @param baseline: float. Use get_baseline(pac_filename)
  #' @param accent: vector. Use get_accent_contour(pac_filename)
  #' @param phrase: vector. Use get_phrase_contour(pac_filename)
  #'
  #' @return: array, contour reconstruction
  #'

  logpitch = accent + log(baseline) + phrase
  pitch = exp(logpitch)

  return (pitch)
}

# functions of the components of the Fujisaki model, see Torres, Gurlekian, Mixdorff, & Pfitzinger (2014)
# DOI: 10.1186/s13636-014-0028-3
Ga_func = function(beta, tt, gamma=Inf){
  val = 1 - (1 + beta * tt) * exp(-beta * tt)
  val[val > gamma] = gamma
  return(val)
}

Ga = function(beta, T1, tmax=2, fs=100){
  val = rep(0,as.integer(tmax * fs))
  trange = seq(0, tmax, 1/fs)[1:len(val)]
  func_domain = trange - T1 >= 0
  val[func_domain] = Ga_func(beta, trange[func_domain] - T1)
  return(val)
}

Ga2 = function(beta, T1, dur, tmax=3, fs=100){
  val = Ga(beta, T1, tmax=tmax, fs=fs) - Ga(beta, T1 + dur, tmax=tmax, fs=fs)
  return(val)
}

Gp_func = function(alpha, tt){
  val = (alpha ** 2) * exp(-alpha * tt) * tt
  return(val)
}

Gp = function(alpha, T1, tmax=2, fs=100){
  val = rep(0, as.integer(tmax * fs))
  trange = seq(0, tmax, 1/fs)[1:length(val)]
  func_domain = trange - T1 >= 0
  val[func_domain] = Gp_func(alpha, trange[func_domain] - T1)
  return(val)
}

estimate_fujisaki_parameters = function (files, f0_ascii_path, host, passwd = NULL) {
  # Currently only works over SSH

  if (any(grepl('\\.', files))) {
    stop('None of the files may contain a "."')
  }

  if (!is.null(passwd)) {
    session = ssh::ssh_connect(host, passwd = passwd)
  } else{
    session = ssh::ssh_connect(host)
  }

  zipname = paste0(tempdir(), "/temp.zip")
  zipdir = paste0(tempdir(), "/temp/")
  required_files = c('Accent1Hz.fir', 'interpolation.exe')
  if (!all(dir.exists(zipdir), required_files %in% list.files(zipdir))) {
    download.file('http://public.beuth-hochschule.de/~mixdorff/thesis/files/fujisaki_analysis.zip', zipname)
    #old_wd = getwd()
    #setwd(tempdir())
    unzip(zipname, exdir = zipdir)
    #setwd(old_wd)
  }

  # Upload template
  template_dir = paste0(.libPaths(),"/contouR/templates/")
  template = paste0(template_dir,"ascii_to_PAC.bat")
  if (!file.exists(template)) {
    stop("Template file does not exist")
  }
  ssh::scp_upload(session, template)


  file_dir = paste0(f0_ascii_path, 'f0_ascii.zip')
  max_batch_size = 40
  time_out = 10
  batches = split(files, ceiling(seq_along(files)/max_batch_size))

  for (batch in batches) {
    # Upload files zip
    if (file.exists(file_dir)) {
      file.remove(file_dir)
    }
    filename_with_ext = paste0(f0_ascii_path, batch, '.f0_ascii')
    zip(file_dir, filename_with_ext, flags = "-j")
    ssh::scp_upload(session, file_dir, verbose = FALSE)

    # Upload the required files (interpolate.exe)
    for (f in required_files) {
      ssh::scp_upload(session, paste0(zipdir, f), verbose = FALSE)
    }

    # Execute script
    options(show.error.messages = FALSE)
    try({R.utils::withTimeout({ssh::ssh_exec_wait(session, 'ascii_to_PAC.bat > output.txt') }, timeout = time_out)})
    options(show.error.messages = TRUE)
    Sys.sleep(time_out)

    # Download results
    for (f in batch) {
      download_files = paste0(f, c('.4.txt', '.PAC'))
      for (df in download_files) {
        ssh::scp_download(session, paste0('f0_ascii/', df), f0_ascii_path, verbose = FALSE)
      }
    }

    ssh::ssh_exec_wait(session, 'del \f f0_ascii.zip')
    ssh::ssh_exec_wait(session, 'rmdir /Q /S f0_ascii')
  }
  ssh::scp_download(session, 'output.txt', f0_ascii_path)

  rm_files = c("ascii_to_PAC.bat", "output.txt")

  # Remove the required files
  for (f in c(rm_files)) {
    ssh::ssh_exec_wait(session, paste0('del \f ', f) )
  }


  for (f in batch) {
    for (ext in c('.4.txt', '.PAC')) {
      if (!file.exists(paste0(f0_ascii_path, f, ext))) {
        warning(paste("File", f, "is missing!"))
      }
    }
  }

  # Disconnect from session
  ssh::ssh_disconnect(session)
}

plot_pitch = function(pac_filename, original_F0 = c(), time = c()){
  #' Plot reconstructed pitch contour from fujisaki components
  #'
  #' @param pac_filename: string. Name of Pac file
  #' @param original_F0: array. F0 values from PitchTier
  #' @param time: array. Time stamps for the F0 values
  #'
  #' @return: plot
  #'
  accent = get_accent_contour(pac_filename)
  phrase = get_phrase_contour(pac_filename)
  baseline = get_baseline(pac_filename)

  fuji_F0 = fuji2pitch(baseline, accent, phrase)

  samp_rate = 0.01 # every 10 ms
  samp_time = seq(samp_rate, length(fuji_F0) * samp_rate, samp_rate)
  requireNamespace("ggplot2")
  library(ggplot2)
  df = data.frame(t = samp_time, f=fuji_F0, method = "Fujisaki")
  title = "Fujisaki reconstructed contour"
  if (length(original_F0) != 0){
    interpol = approxfun(time, original_F0)
    original_F0 = interpol(samp_time)
    df = rbind(df, data.frame(t = samp_time, f=original_F0, method = "Original"))
    rmse = RMSE(fuji_F0, original_F0)
    title = paste("Original vs. Fujisaki reconstructed contour. RMSE:", rmse)
  }

  ggplot(df) + geom_point(aes(x = t, y = f, color = method)) + ggtitle(title) + ylab('F0 in Hz') + xlab('Time in ms')
}

write_f0_ascii = function(t, f, fname){
  #' Extract a pitch contour for a sentence and save it in a file format that is accepted by FujiParamEditor
  #' @param t: time
  #' @param f: F0 values
  #' @param fname: name of the f0_ascii file

  if (!grepl('.f0_ascii', fname)){
    err0("Filename (", fname, ") needs to have the .f0_ascii extension")
  }

  # f = f[!is.na(t)]
  # t = t[!is.na(t)]
  #
  # t = t[!is.na(f)]
  # f = f[!is.na(f)]
  #
  # if (any(length(t) <= 2, length(f))){
  #   err0("Filename (", fname, ") must contain more than 2 not NA values")
  # }

  # Important addition
  interpol = approxfun(t, f)
  sampling_rate = .01
  num_samples = ceiling(max(t)/sampling_rate)
  max_t = sampling_rate * num_samples
  t = seq(from=0, to=max_t, by=sampling_rate)
  f = interpol(t)

  f[is.na(f)] = 0
  voiced = as.double(f != 0)
  df = data.frame(pitch=f, voiced = voiced, na1 = 1.0, na2 = voiced)
  write.table(df, fname, sep=" ", row.names= FALSE, col.names = FALSE)
}

read_f0_ascii = function(f0_ascii_file){
  #' reads f0 file that is in the format taken by FujiParamEditor
  #' @param f0_ascii_file: filepath
  #' @return: pitch, tt_pitch, dt

  df = read.table(f0_ascii_file, sep=' ', header = FALSE)
  names(df) = c('pitch', 'voiced', 'na', 'na')
  pitch = df$pitch
  pitch[pitch == 0] = NA
  tt_pitch = (1:length(pitch)) / 100

  return(list(pitch=pitch, tt_pitch=tt_pitch, dt=.01))
}

compute_optimal_hyper_parameter = function(
  full_path, pt,
  a_start, a_step, a_end,
  b_start, b_step, b_end
){
  #' Computes the RMSE between raw contour (pt) and PAC file for alpha (a) and beta (b) combinations
  #' @param full_path: the full path to the PAC file
  #' @param pt: PitchTier object from rPraat
  #' @param a_start: first alpha number
  #' @param a_step: steps between alpha start and end
  #' @param a_end: last alpha number
  #' @param b_start: first beta number
  #' @param b_step: steps between beta start and end
  #' @param b_end: last beta number
  #'

  # Get time and pitch from PitchTier
  t = pt$t
  f = pt$f

  # Read settings from PAC
  pac_list = read_pac_file(full_path)

  # Get the baseline
  baseline = pac_list$baseline

  # Generate possible alpha and beta values
  beta_values = seq(b_start, b_end, b_step)
  alpha_values = seq(a_start, a_end, a_step)

  # DF for plotting; this contains all RMSE for all alpha vs. beta combinations
  plot_df = NULL
  for(beta in beta_values){
    accent = get_accent_contour(force_beta = beta, pac_list = pac_list)
    for (alpha in alpha_values){
      phrase = get_phrase_contour(force_alpha = alpha, pac_list = pac_list)
      f_est = fuji2pitch(baseline, accent, phrase)
      t_est = seq(pt$t[1], tail(pt$t, 1), length=length(f_est))
      # Compute RMSE between raw pitch points and estimated contour
      # Important: This means that the number of comparisons between reconstructed and raw pitch contour,
      # depends on the points in the raw contour!
      plot_df = rbind(plot_df, data.frame(alpha = alpha, beta = beta, RMSE = RMSE_points_only(t, f, t_est, f_est)))
    }
  }

  return(plot_df)
}

optimize_PAC = function(PAC, pt_path, Fujisaki_path, a_start, a_step, a_end, b_start, b_step, b_end, plot = FALSE){
  # Read PitchTier
  name = strsplit(PAC, "\\.")[[1]][1]
  pt = read_PitchTier(paste0(pt_path, name, ".PitchTier"))

  full_path = paste0(Fujisaki_path, name, '.PAC')
  plot_df = compute_optimal_hyper_parameter(full_path, pt, a_start, a_step, a_end, b_start, b_step, b_end)

  # Take the 5 best (i.e. with the lowest RMSE) beta, alpha combinations
  top_5 = plot_df[tail(rev(order(plot_df$RMSE)), 5), ]
  top_5$file = name
  top_5$placement = 1:5


  #pp("Finished", name)

  if (plot){
    # Create a heat map of the RMSEs for all beta vs. alpha combinations
    RMSE_min = min(plot_df$RMSE)
    RMSE_max = max(plot_df$RMSE)

    cols = c(colorRampPalette(c("#e7f0fa", "#c9e2f6", "#95cbee", "#0099dc", "#4ab04a", "#ffd73e"))(10), colorRampPalette(c("#eec73a", "#e29421", "#e29421", "#f05336","#ce472e"), bias=2)(90))

    requireNamespace("scales")
    requireNamespace("ggplot2")
    library(ggplot2)
    ggplot(plot_df, aes(y=alpha, x=beta, fill=RMSE)) +
      geom_tile(colour="white", width=.9, height=.9) +
      theme_minimal() +
      scale_fill_gradientn(colours=cols,
                           limits=c(RMSE_min, RMSE_max),
                           breaks=seq(RMSE_min, RMSE_max, length=2),
                           na.value=rgb(246, 246, 246, max=255),
                           labels=floor(seq(RMSE_min, RMSE_max, length=2)),
                           guide=guide_colourbar(ticks=T, nbin=50, barheight=.5, label=T, barwidth=10)
      ) +
      theme(legend.position=c(.5, -.13),
            legend.direction="horizontal",
            legend.text=element_text(colour="grey20"),
            plot.margin=grid::unit(c(.5,.5,2,.5), "cm"),
            panel.grid=element_blank(),
            axis.text.y=element_text(size=10, family="Helvetica"),
            axis.text.x=element_text(size=10),
            title=element_text(family="Helvetica"),
      ) +
      ggtitle(paste("RSME for different hyper parameters", name))

    # Due to the amount of plots, they will be saved
    ggsave(paste0("heatmap_", name, ".pdf"))
  }

  return(top_5)
}

explore_ft_space = function(
  Fujisaki_path, pt_path,
  plot = TRUE,
  a_start = 0.5, a_step = 0.2, a_end = 4,
  b_start = 15, b_step = 1, b_end = 35
){
  #' Explores hyper parameters alpha and beta
  #' @param Fujisaki_path: absolute path to directory containing all PAC files (generated with Mixdorff tool)
  #' @param pt_path: absolute path to directory containing all PitchTiers
  #' @param plot: boolean if should display plots; if TRUE, it will save in current directory
  #' @param a_start: first alpha number
  #' @param a_step: steps between alpha start and end
  #' @param a_end: last alpha number
  #' @param b_start: first beta number
  #' @param b_step: steps between beta start and end
  #' @param b_end: last beta number
  #'

  require(grid)
  require(stringr)
  require(dplyr)
  require(scales)

  # Get all PAC files
  files = list.files(Fujisaki_path)
  files = files[grepl(".PAC", files)]

  # Empty df for scores
  top_scores = NULL
  for (PAC in files){
    top_scores = rbind(top_scores, optimize_PAC(PAC, pt_path, Fujisaki_path, a_start, a_step, a_end, b_start, b_step, b_end, plot))
  }

  # Add meta data to top_scores
  top_scores = cbind(top_scores, stringr::str_split_fixed(top_scores$file, "_", 3))
  names(top_scores)[6:8] = c("speaker", "emotion", "sentence_ID")

  if (plot){
    # Take only the best estimations
    place_1 = dplyr::filter(top_scores, placement == 1)

    # Statistically compare RMSE across emotions, surprise significantly higher (Kruskal-Wallis)
    comparisons = list(c("SUR", "SAD"), c("SUR", "NEU"), c("SUR", "HAP"), c("SUR", "FER"), c("SUR", "DIS"), c("SUR", "ANG"))
    requireNamespace("ggpubr")
    library("ggpubr")
    ggboxplot(place_1, x = "emotion", y = "RMSE",
              color = "emotion")+
      stat_compare_means(comparisons = comparisons)+ # Add pairwise comparisons p-value
      stat_compare_means(label.y = 50)     # Add global p-value

    # Compare how much of the corpus is below a certain RMSE value
    distr_RMSE = NULL
    xs = 1:150
    emotions = c("SUR", "ANG", "DIS", "FER", "HAP", "SAD", "NEU", "AVERAGE", "AVG-SUR")

    for (emo in emotions){
      if (emo == "AVERAGE"){
        sel_df = place_1
      } else if (emo == "AVG-SUR"){
        sel_df = dplyr::filter(place_1, emotion != "SUR")
      }else{
        sel_df = dplyr::filter(place_1, emotion == emo)
      }
      for (x in xs){
        y = length(which(sel_df$RMSE < x))/nrow(sel_df)
        distr_RMSE = rbind(distr_RMSE, data.frame(x = x, y = y, emotion = emo))
      }
    }
    AVERAGE = dplyr::filter(distr_RMSE, emotion == "AVERAGE")
    AVG_SUR = dplyr::filter(distr_RMSE, emotion == "AVG-SUR")
    distr_RMSE = dplyr::filter(distr_RMSE, emotion %in% emotions[1:7])
    ggplot(distr_RMSE) +
      geom_line(aes(x=x, y=y, color = emotion)) +
      geom_line(aes(x=x, y=y), data = AVERAGE, color="darkgrey", lwd=1) +
      geom_line(aes(x=x, y=y), data = AVG_SUR, lwd=1) +
      geom_vline(xintercept = 20, color="grey") +
      theme_minimal() +
      theme(
        plot.title=element_text(family="IBM Plex Sans Medium", size=16),
        legend.position="bottom",
        legend.title = element_blank()
      ) +
      scale_y_continuous(labels = scales::percent) +
      ylab('Amount of corpus below RMSE') +
      xlab('RMSE thesholds in Hz') +
      ggtitle('Amount of corpus below RMSE theshold')

    ggsave('RMSE_threshold.pdf', device = cairo_pdf)
  }

  return(top_scores)
}

plot_optimal_estimation = function(
  name, pt_path, Fujisaki_path,
  save=FALSE,
  a_start = 0.5, a_step = 0.2, a_end = 4,
  b_start = 15, b_step = 1, b_end = 35
){
  #' Plots all possible
  #' @param name: base file name, e.g. "SL_DIS_VX1b"
  #' @param Fujisaki_path: absolute path to directory containing all PAC files (generated with Mixdorff tool)
  #' @param pt_path: absolute path to directory containing all PitchTiers
  #' @param save: boolean, indicating to save it
  #' @param a_start: first alpha number
  #' @param a_step: steps between alpha start and end
  #' @param a_end: last alpha number
  #' @param b_start: first beta number
  #' @param b_step: steps between beta start and end
  #' @param b_end: last beta number
  #'
  require(dplyr)
  pt = read_PitchTier(paste0(pt_path, name, ".PitchTier"))

  full_path = paste0(Fujisaki_path, name, '.PAC')
  param_df = compute_optimal_hyper_parameter(full_path, pt, a_start, a_step, a_end, b_start, b_step, b_end)
  pac_list = read_pac_file(paste0(Fujisaki_path, name, '.PAC'))
  baseline = pac_list$baseline

  plot_df = NULL

  for (i in 1:nrow(param_df)){
    alpha = param_df[[i, 'alpha']]
    beta = param_df[[i, 'beta']]
    RMSE = param_df[[i, 'RMSE']]

    accent = get_accent_contour(pac_list = pac_list, force_beta = beta)
    phrase = get_phrase_contour(pac_list = pac_list, force_alpha = alpha)

    fuji_F0 = fuji2pitch(baseline, accent, phrase)
    samp_rate = 0.01 # every 10 ms
    samp_time = seq(samp_rate, length(fuji_F0) * samp_rate, samp_rate)
    plot_df = rbind(plot_df, data.frame(
      alpha = alpha,
      beta = beta,
      RMSE = RMSE,
      f = fuji_F0,
      t = samp_time
    ))
  }

  # Convert F0 to semitones
  plot_df$f = f2st(plot_df$f)
  pt$f = f2st(pt$f)

  library(ggplot2)
  p = ggplot() +
    geom_line(aes(x = t, y = f, color = RMSE, group=interaction(alpha, beta)), data = plot_df) +
    scale_color_gradient2(high="#dddddd") +
    theme_minimal() +
    geom_line(aes(x = t, y = f), color='black', data = dplyr::filter(plot_df, RMSE==sort(unique(plot_df$RMSE))[2]))+
    geom_point(aes(x = t, y = f), data = pt) +
    ylab("F0 in semitones") +
    xlab("Time in seconds") +
    ggtitle(paste("Optimal contour estimation and raw contour for", name))

  if (save){
    ggsave(paste0('optimal_estimation_', name, '.pdf'), useDingbats = FALSE)
  }

  return(p)
}

compute_phrase_feature = function(alpha, beta, file, Fujisaki_path, plot = FALSE){
  # Compute the phrase
  pac_list = read_pac_file(paste0(Fujisaki_path, file, '.PAC'))

  # Feature 1
  num_phrases = length(pac_list$phrase_idxs)

  phrase_df = pac_list$df[pac_list$phrase_idxs, ]

  # Feature 2
  Ap = mean(phrase_df$A)

  if (num_phrases == 1){
    # Feature 3
    T0 = phrase_df[['T1']]
    phrase = get_phrase_contour(force_alpha = alpha, pac_list = pac_list)

    n_xs = length(phrase)
    xs = 1:n_xs

    # Feature 4: regression slope
    lin_model = lm(phrase ~ xs)
    reg_slope = lin_model$coefficients[[2]]
    reg_RMSE = sqrt(mean(lin_model$residuals^2))
    if (plot){
      plot(phrase)
      abline(lin_model)
      title(paste("RMSE linear regression slope:", reg_RMSE))
    }

    # Feature 5: naive slope
    intercept = phrase[1]
    naive_slope = (phrase[n_xs] - phrase[1])/(n_xs - 1)
    y_values = intercept + naive_slope*xs
    naive_RMSE = RMSE(y_values, phrase, remove_NA = FALSE)
    if (plot){
      plot(phrase)
      lines(y_values)
      title(paste("RMSE naive slope:", naive_RMSE))
    }

  } else {
    T0 = NA
    reg_slope = NA
    reg_RMSE = NA
    naive_slope = NA
    naive_RMSE = NA
  }

  return(data.frame(
    name = file,
    num_phrases = num_phrases,
    Ap = Ap,
    T0 = T0,
    reg_slope = reg_slope,
    reg_RMSE = reg_RMSE,
    naive_slope = naive_slope,
    naive_RMSE = naive_RMSE
  ))
}

compute_phrase_features = function(top_scores, Fujisaki_path, plot = FALSE){
  #' Compute features on phrase command
  #'
  #' @param top_scores: resulting DF of explore_ft_space(...)
  #' @param plot: boolean indicating whether to plot
  #' @return results: dataframe containing all features
  require(dplyr)
  # Select only the optimal settings & extract relevant columns
  place_1 = dplyr::filter(top_scores, placement == 1)
  place_1 = place_1[, c(1,2,4)]

  results = NULL
  for (r in 1:nrow(place_1)){
    row = place_1[r, ]
    results = rbind(results, compute_phrase_feature(row$alpha, row$beta, row$file, Fujisaki_path))
  }

  return(results)
}

compute_accent_features = function(top_scores, Fujisaki_path){
  #' Compute features on accent command
  #'
  #' @param top_scores: resulting DF of explore_ft_space(...)
  #' @return results: dataframe containing all features

  # Select only the optimal settings & extract relevant columns
  require(dplyr)
  place_1 = dplyr::filter(top_scores, placement == 1)
  place_1 = place_1[, c(1,2,4)]

  results = NULL
  for (row in 1:nrow(place_1)){
    # Compute the phrase
    name = as.character(place_1[row, 3])
    pac_list = read_pac_file(paste0(Fujisaki_path, name, '.PAC'))

    # Feature 1: Number of accents
    num_accents = length(pac_list$accent_idxs)

    if (num_accents == 0){
      Aa_mean = NA
      Aa_1 = NA
      Aa_last = NA
    } else {
      accent_df = pac_list$df[pac_list$accent_idxs, ]

      # Feature 2: Average accent amplitude
      Aa_mean = mean(accent_df$A)

      # Feature 3: First accent amplitude
      Aa_1 = head(accent_df$A, 1)

      # Feature 4: Last accent amplitude
      Aa_last = tail(accent_df$A, 1)
    }

    results = rbind(results, data.frame(name, num_accents, Aa_mean, Aa_1, Aa_last))
  }

  return(results)
}

