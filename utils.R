suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(readr))
suppressMessages(library(smoother))

calc_bamsizes = function(i)
{
  source("utils.R")
  suppressMessages(library(dplyr))

  path_bamsize = file.path(path_bam, stringr::str_glue("{bam}_size.txt", bam=basename(binned_arguments_df$SAMPLE_BAM[i])))
  if(!file.exists(path_bamsize)) {
    cmd_bamsize = stringr::str_glue("samtools view -c {bam} > {size}", bam=binned_arguments_df$SAMPLE_BAM[i], size=path_bamsize)
    cmd_run(cmd_bamsize)

    df_bamsize = data.frame(SAMPLE_BAM=binned_arguments_df$SAMPLE_BAM[i], SAMPLE_SIZE=as.numeric(readLines(path_bamsize)))
    readr::write_tsv(df_bamsize, path=path_bamsize)
  }
}

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
  repliseq_time_df = repliseq_df %>%
    dplyr::mutate(repliseqTime_chrom=repliseq_chrom, repliseqTime_start=repliseq_start, repliseqTime_end=repliseq_end) %>%
    dplyr::group_by(repliseqTime_chrom, repliseqTime_start, repliseqTime_end) %>%
    dplyr::mutate(
      repliseq_value_norm=repliseq_value,
      repliseq_value_norm=repliseq_value_norm-min(repliseq_value_norm, na.rm=T),
      repliseq_value_norm=repliseq_value_norm/max(repliseq_value_norm, na.rm=T),
      repliseq_value_norm=repliseq_value_norm^2,
    ) %>%
    dplyr::summarize(repliseqTime_avg=ifelse(any(!is.na(repliseq_value) & repliseq_value>0), weighted.mean(repliseq_fraction, repliseq_value_norm, na.rm=T), mean(repliseq_fraction))) %>%
    dplyr::group_by(repliseqTime_chrom) %>%
    dplyr::mutate(repliseqTime_avg=smoother::smth.gaussian(repliseqTime_avg, window=window)) %>%
    dplyr::mutate(repliseqTime_avg=zoo::na.fill(repliseqTime_avg, "extend")) %>%
    dplyr::ungroup()

  repliseq_time_df
}