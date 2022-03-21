suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(readr))
suppressMessages(library(smoother))

cmd_is_available = function(programs) {
  missing_programs = c()
  for(prg in programs) {
    if (Sys.which(prg) == "") {
      missing_programs = c(missing_programs, prg)
    }
  }
  if(length(missing_programs)>0) {
    error_message = paste0("Some executables required to run the pipeline were not found:\n    ", paste(missing_programs, collapse="\n    "))
    stop(error_message)
  }
}

cmd_run = function(cmd, intern=F) {
  start_time = Sys.time()
  writeLines(paste(">", cmd))
  x = system(cmd, intern=intern)
  end_time = Sys.time()
  writeLines(paste0(round(as.numeric(difftime(end_time, start_time, units="secs"))), "s"))

  if(intern) return(x)
  return(NULL)
}

repliseq_read = function(path) {
  repliseq_df = as.data.frame(t(readr::read_tsv(path, col_names=F))) %>%
    dplyr::rename(repliseq_chrom="V1", repliseq_start="V2", repliseq_end="V3") %>%
    dplyr::mutate(repliseq_start=as.numeric(repliseq_start), repliseq_end=as.numeric(repliseq_end)) %>%
    reshape2::melt(id.vars=c("repliseq_chrom", "repliseq_start", "repliseq_end"), variable.name="repliseq_fraction", value.name="repliseq_value") %>%
    dplyr::mutate(repliseq_fraction=as.numeric(gsub("V", "", repliseq_fraction))-3, repliseq_value=as.numeric(repliseq_value))
  repliseq_df.keep = repliseq_df %>%
    dplyr::arrange(repliseq_start) %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end) %>%
    dplyr::summarize(repliseq_isna=mean(is.na(repliseq_value))>0.5) %>%
    dplyr::group_by(repliseq_chrom) %>%
    dplyr::mutate(first_value_position=min(which(!repliseq_isna)), last_value_position=max(which(!repliseq_isna)), keep_value=dplyr::between(1:n(), first_value_position[1], last_value_position[1])) %>%
    dplyr::filter(keep_value) %>%
    dplyr::ungroup() %>%
    dplyr::select(repliseq_chrom, repliseq_start, repliseq_end)
  repliseq_df.f = repliseq_df %>%
    dplyr::inner_join(repliseq_df.keep, by=c("repliseq_chrom", "repliseq_start", "repliseq_end"))

  repliseq_df.f
}

repliseq_summarize = function(repliseq_df, window=5) {
  th.repliseq_value_norm = 0.1
  repliseq_df = repliseq_df %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end) %>%
    dplyr::mutate(repliseq_value_norm=((repliseq_value)/max(repliseq_value))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end) %>%
    do((function(z) {
      zz<<-z

      i.max = z$repliseq_fraction[which.max(z$repliseq_value_norm)]
      lb = which(z$repliseq_value_norm < th.repliseq_value_norm & z$repliseq_fraction < i.max)
      lb = z$repliseq_fraction[ifelse(length(lb), lb[length(lb)]+1, 1)]
      ub = which(z$repliseq_value_norm < th.repliseq_value_norm & z$repliseq_fraction > i.max)
      ub = z$repliseq_fraction[ifelse(length(ub), ub[1]-1, nrow(z))]
      # dplyr::between(z$repliseq_fraction, lb, ub)

      z$repliseq_value_in_scope = T
      z
    })(.)) %>%
    dplyr::ungroup()

  repliseq_time_df = repliseq_df %>%
    dplyr::mutate(repliseqTime_chrom=repliseq_chrom, repliseqTime_start=repliseq_start, repliseqTime_end=repliseq_end) %>%
    dplyr::group_by(repliseqTime_chrom, repliseqTime_start, repliseqTime_end) %>%
    dplyr::mutate(repliseqTime_avg=weighted.mean(repliseq_fraction[repliseq_value_in_scope], repliseq_value_norm[repliseq_value_in_scope], na.rm=T)) %>%
    dplyr::mutate(lb=dplyr::between(repliseq_fraction, floor(repliseqTime_avg[1]-4), ceiling(repliseqTime_avg[1]-2)), repliseqTime_min=ifelse(any(lb), weighted.mean(repliseq_fraction[lb], repliseq_value_norm[lb]), NA_real_)) %>%
    dplyr::mutate(ub=dplyr::between(repliseq_fraction, floor(repliseqTime_avg[1]+2), ceiling(repliseqTime_avg[1]+4)), repliseqTime_max=ifelse(any(ub), weighted.mean(repliseq_fraction[ub], repliseq_value_norm[ub]), NA_real_)) %>%
    dplyr::summarize(repliseqTime_avg=repliseqTime_avg[1], repliseqTime_min=repliseqTime_min[1], repliseqTime_max=repliseqTime_max[1], repliseqTime_lb=min(which(repliseq_value_in_scope)), repliseqTime_ub=max(which(repliseq_value_in_scope))) %>%
    dplyr::group_by(repliseqTime_chrom) %>%
    dplyr::mutate(repliseqTime_avg=smoother::smth.gaussian(repliseqTime_avg, window=window), repliseqTime_min=smoother::smth.gaussian(repliseqTime_min, window=window), repliseqTime_max=smoother::smth.gaussian(repliseqTime_max, window=window)) %>%
    dplyr::mutate(repliseqTime_avg=zoo::na.fill(repliseqTime_avg, "extend"), repliseqTime_min=zoo::na.fill(repliseqTime_min, "extend"), repliseqTime_max=zoo::na.fill(repliseqTime_max, "extend")) %>%
    dplyr::ungroup()

  repliseq_time_df
}

asda = function() {
  library(readr)
  library(ggplot2)
  library(Biostrings)
  devtools::load_all('~/Workspace/breaktools/')
  rdc_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_pnas_mm10.tsv") %>%
    tidyr::crossing(rdc_offset=c(-1, 0, 1)) %>%
    dplyr::mutate(rdc_region_start=rdc_start+(rdc_end-rdc_start)*rdc_offset, rdc_region_end=rdc_end+(rdc_end-rdc_start)*rdc_offset)
    dplyr::mutate(rdc_region_start=(rdc_region_end+rdc_region_start)/2-5000, rdc_region_end=rdc_region_start+10000)
  rdc_ranges = rdc_df %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_region_start, end=rdc_region_end) %>% GenomicRanges::makeGRangesFromDataFrame()
  rdc_df$seq = get_seq("~/Workspace/genomes/mm10/mm10.fa", rdc_ranges)$sequence

  r1 = ".*(A{12,}(.{0,5})T{12,}).*"
  r2 = ".*(T{12,}(.{0,5})A{12,}).*"
  rdc_df$r1_motif = ifelse(grepl(r1, rdc_df$seq,  perl=T, ignore.case=T), gsub(r1, "\\1", rdc_df$seq,  perl=T, ignore.case=T), "")
  rdc_df$r1_gap = ifelse(nchar(rdc_df$r1_motif)>0, gsub(r1, "\\2", rdc_df$r1_motif,  perl=T, ignore.case=T), NA_character_)
  rdc_df$r2_motif = ifelse(grepl(r2, rdc_df$seq,  perl=T, ignore.case=T), gsub(r2, "\\1", rdc_df$seq,  perl=T, ignore.case=T), "")
  rdc_df$r2_gap = ifelse(nchar(rdc_df$r2_motif)>0, gsub(r2, "\\2", rdc_df$r2_motif,  perl=T, ignore.case=T), NA_character_)
  rdc_df$gap = pmax(nchar(rdc_df$r1_gap), nchar(rdc_df$r2_gap), na.rm=T)
  rdc_df$polA_count = sapply(gregexpr("(A{2,}(.{0,1})){10}A{2,}", rdc_df$seq), function(x) sum(x>0))
  rdc_df$polT_count = sapply(gregexpr("(T{2,}(.{0,1})){10}T{2,}", rdc_df$seq), function(x) sum(x>0))

  table(rdc_df$rdc_offset, rdc_df$polA_count>0)
  table(rdc_df$rdc_offset, rdc_df$polT_count>0)
  table(rdc_df$rdc_offset, !is.na(rdc_df$r1_gap))
  table(rdc_df$rdc_offset, !is.na(rdc_df$r2_gap))

  ggplot(rdc_df) +
    geom_boxplot(aes(x=1, fill=as.factor(rdc_offset), y=gap))

}