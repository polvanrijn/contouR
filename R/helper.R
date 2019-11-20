# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# This contains all kinds of helper functions

resample = function(X, num_landmarks, output = FALSE){
  if (output){
    pp("Downsampling from", length(X), "to", num_landmarks)
  }
  return((X[seq(1,length(X),length=num_landmarks)])[-1])
}

add_meta_data = function(df, ID_col_name = "filename"){
  require(stringr)
  if (ID_col_name %in% names(df)){
    if (!any(names(df) %in% c("speaker", "emotion", "sentence_ID"))){
      df = cbind(df, stringr::str_split_fixed(df[[ID_col_name]], "_", 3))
      names(df)[(ncol(df)-2):ncol(df)] = c("speaker", "emotion", "sentence_ID")
      df$sentence = as.numeric(stringr::str_extract(df$sentence_ID, "(\\d)+"))
    }
  } else{
    err0("Does not contain the column '", ID_col_name, "'")
  }
  return(df)
}

RMSE = function(m, o, output = FALSE, remove_NA = TRUE){
  if (length(m) != length(o)){
    stop("Both arrays must have the same size!")
  }

  NA_idx = rep(FALSE, length(m))
  if (length(which(is.na(m))) > 0){
    NA_idx = is.na(m) | NA_idx
  }
  if (length(which(is.na(o))) > 0){
    NA_idx = is.na(o) | NA_idx
  }

  if (length(which(NA_idx)) > 0 && remove_NA){
    o = o[!NA_idx]
    m = m[!NA_idx]
    if (output){
      pp(length(which(NA_idx)), "NAs removed!")
    }
  }

  rmse = sqrt(mean((m - o)^2))

  if (output){
    pp("RMSE:", rsme)
  }

  return(rmse)
}

RMSE_points_only = function(t, f, t_est, f_est, output = FALSE, remove_NA = TRUE){
  interpol = approxfun(t_est, f_est)
  f_est = interpol(t)
  return(RMSE(f, f_est, output, remove_NA))
}

combine_features = function(df1, df2, key = "filename"){
  if (!(key %in% names(df1) && key %in% names(df2))){
    stop(paste("Key", key, "must be present in all dataframes"))
  }
  # d2_ordered = df2[match(df1[[key]], df2[[key]]), ]
  # col_bool_idx = !names(d2_ordered) == key
  # return(cbind(df1, d2_ordered[, col_bool_idx]))
  return(merge(df1, df2, by=key, incomparables = NA, all = TRUE))
}

# Alias
len = function(x){
  return(length(x))
}

pp = function(..., sep = " ", collapse = NULL){
  return(print(paste(..., sep = sep, collapse = collapse)))
}

pp0 = function(..., collapse = NULL){
  return(print(paste0(..., collapse = collapse)))
}

err0 = function(..., collapse = NULL){
  stop(paste0(..., collapse = collapse))
}

warn0 = function(..., collapse = NULL){
  war = paste0(..., '\n', collapse = collapse)
  warning(war)
}

f2st = function(f){
  require(hqmisc)
  requireNamespace("hqmisc")
  return(hqmisc::f2st(f, base = 70))
}

get_pt_list = function(pt_path, file_name){
  name = stringr::str_replace(file_name, ".PitchTier", "")
  pt = read_PitchTier(paste0(pt_path, file_name))
  f = pt$f
  f_st = f2st(f)
  return(list(f_st = f_st, t = pt$t, name = name))
}

read_PitchTier = function(full_path, cache=TRUE){
  f_name = tail(strsplit(full_path, '\\/')[[1]], 1)
  splitted_fn = strsplit(f_name, '\\.')[[1]]
  splitted_path = strsplit(full_path, '\\.')[[1]]
  splitted_path = paste(splitted_path[1:(length(splitted_path) - 1)], collapse = ".")
  if (length(splitted_fn) != 2){
    err0("Filename (", full_path, ") may only contain a single dot (.)!!!")
  }
  csv_path = paste0(splitted_path[1], '.csv')
  if (all(cache, file.exists(csv_path))){
    return(read.csv(csv_path))
  } else{
    loadNamespace("rPraat")
    library(rPraat)
    pt = pt.read(full_path)
    pt_df = data.frame(t = pt$t, f = pt$f)
    write.csv(pt_df, csv_path, row.names = FALSE)
    return(pt_df)
  }
}

get_pt_filnames = function(path){
  files = list.files(path)
  return(files[grep('.PitchTier', files)])
}

cache_pts = function(path){
  #' This is a function that helps you cache all your PTs. This is usefull as there is some weird bug in the rPraat library
  pts = get_pt_filnames(path)

  for (pt_name in pts){
    success = FALSE
    while (!success){
      tryCatch({
        read_PitchTier(paste0(path, pt_name))
        success = TRUE
      })
    }
  }
}


read_TextGrid = function(full_path){
  f_name = tail(strsplit(full_path, '\\/')[[1]], 1)
  splitted_fn = strsplit(f_name, '\\.')[[1]]
  if (length(splitted_fn) != 2) {
    err0("Filename (", full_path, ") may only contain a single dot (.)!!!")
  }

  loadNamespace("rPraat")
  library(rPraat)
  return(tg.read(full_path, 'auto'))
}

significance_test = function(df, comp_col, ref_col = "emotion", significance_test_method = "by_group", global_sig_value = TRUE){
  library(ggpubr)
  if (!any(c("speaker", "emotion", "sentence") %in% names(df))){
    df = add_meta_data(df)
  }

  if (!significance_test_method %in% c("by_group", "sig_only")){
    stop("Method not supported")
  }

  p = ggboxplot(df, x = ref_col, y = comp_col, color = ref_col) + ggtitle(paste("Comparing", comp_col, "across", ref_col))

  if (significance_test_method == "sig_only"){
    sig_tab = compare_means(as.formula(paste(comp_col, "~", ref_col)), data=df, method = "t.test",p.adjust.method = "bonferroni")
    sig_tab = as.data.frame(sig_tab)
    sig_idxs = which(sig_tab$p.signif != "ns" & sig_tab$group2 != "NEU")
    # Only significant values
    my_comparisons <- mapply(c, sig_tab[sig_idxs,2], sig_tab[sig_idxs,3], SIMPLIFY=FALSE)

    p = p + stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
  } else if (significance_test_method == "by_group") {
    p = p + stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") # Add pairwise comparisons p-value
  }

  if (global_sig_value){
    p = p + stat_compare_means(label.y = max(df[[comp_col]]) + sd(df[[comp_col]]))     # Add global p-value
  }

  return(p)
}
