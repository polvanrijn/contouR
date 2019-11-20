compute_distribution_features = function(pt_file_names, pt_path){
  library(moments)
  require(EnvStats)
  require(DescTools)
  requireNamespace("EnvStats")
  requireNamespace("DescTools")
  results = NULL
  for (file_name in pt_file_names) {
    pt_list = get_pt_list(pt_path, file_name)
    f_st = pt_list$f_st

    sk = EnvStats::skewness(f_st)
    ku = EnvStats::kurtosis(f_st)
    #ent = DescTools::Entropy(f_st)
    ent = entropy(f_st)
    results = rbind(results, data.frame(name = pt_list$name, skewness = sk, kurtosis = ku, entropy = ent, sd = sd(f_st)))
  }

  return(results)
}


entropy = function(x, scale = F){
  k = length(unique(x))
  freqs = table(x)/length(x)
  entrop = -sum(freqs * log2(freqs))
  if (scale) {
    entrop = entrop/log2(k)
  }
  return(entrop)
}
