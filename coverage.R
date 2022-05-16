#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
options(dplyr.summarise.inform=F)
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(spatstat))
suppressMessages(require(parallel))
source("utils.R")

make_bedgraph = function(i) {
  suppressMessages(library(dplyr))
  source("utils.R")

  r = paste0(i, "/", nrow(binned_arguments_df))
  writeLines(paste0("\n==========================================\n", r, ": ", binned_arguments_df$SAMPLE_NAME[i], " / bin", binned_arguments_df$SAMPLE_BINSIZE[i], "\n=========================================="))

  path_tiles = R.utils::getAbsolutePath(file.path(path_output, stringr::str_glue("bedgraph/genome_{format(binsize, scientific=F)}tiles.bed", binsize=binned_arguments_df$SAMPLE_BINSIZE[i])), expandTilde=T)
  path_tiles_tmp = R.utils::getAbsolutePath(file.path(path_output, stringr::str_glue("bedgraph/genome_{format(binsize, scientific=F)}tiles_tmp.bed", binsize=binned_arguments_df$SAMPLE_BINSIZE[i])), expandTilde=T)
  path_bedgraph = R.utils::getAbsolutePath(file.path(path_output, stringr::str_glue("bedgraph/{sample}_bin{format(binsize, scientific=F)}.bdg", sample=binned_arguments_df$SAMPLE_NAME[i], binsize=binned_arguments_df$SAMPLE_BINSIZE[i])), expandTilde=T)
  path_chromsizes = R.utils::getAbsolutePath(file.path(path_output, stringr::str_glue("bedgraph/{sample}_{format(binsize, scientific=F)}.sizes", sample=binned_arguments_df$SAMPLE_NAME[i], binsize=binned_arguments_df$SAMPLE_BINSIZE[i])), expandTilde=T)

  if(!file.exists(path_bedgraph)) {
    cmd_chromizes = stringr::str_glue("samtools idxstats '{bam}' | cut -f 1-2 | grep -vP '\\t0' > '{chromsizes}'", bam=binned_arguments_df$SAMPLE_BAM[i], chromsizes=path_chromsizes)
    cmd_run(cmd_chromizes)
    cmd_makewindows = stringr::str_glue("bedtools makewindows -w {format(binsize, scientific=F)} -s {format(binsize, scientific=F)} -g '{chromsizes}' > '{tiles}'; bedtools sort -faidx '{chromsizes}' -i '{tiles}' > {tiles_tmp}; mv '{tiles_tmp}' '{tiles}'", chromsizes=path_chromsizes, binsize=binned_arguments_df$SAMPLE_BINSIZE[i], tiles=path_tiles, tiles_tmp=path_tiles_tmp)
    cmd_run(cmd_makewindows)
    cmd_coverage = stringr::str_glue("bedtools intersect -abam '{tiles}' -b '{bam}' -c -bed -g '{chromsizes}' -sorted > '{bedgraph}'", sample=binned_arguments_df$SAMPLE_NAME[i], chromsizes=path_chromsizes, bam=binned_arguments_df$SAMPLE_BAM[i], tiles=path_tiles, bedgraph=path_bedgraph)
    cmd_run(cmd_coverage)
    file.remove(path_tiles)
    file.remove(path_chromsizes)
  } else {
    writeLines(paste0("Bedgraph file '", path_bedgraph, "' exists, skipping..."))
  }
}

write_results = function(i, bedgraph_df) {
  suppressMessages(library(dplyr))
  source("utils.R")

  z = bedgraph_df %>%
    dplyr::inner_join(bedgraph_output_df[i,,drop=F], by=c("SAMPLE_CONDITION", "SAMPLE_BINSIZE")) %>%
    dplyr::arrange(repliseq_chrom, repliseq_start, repliseq_end, dplyr::desc(repliseq_fraction))
  writeLines(paste0("\n==========================================\n", i, "/", nrow(bedgraph_output_df), ": ", z$SAMPLE_CONDITION[1], " / bin", z$SAMPLE_BINSIZE[1], "\n=========================================="))

  z_path_igv = file.path(path_output, stringr::str_glue("results/{cond}_repliseq{format(binsize, scientific=F)}.igv", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))
  z_path_mat = file.path(path_output, stringr::str_glue("results/{cond}_repliseq{format(binsize, scientific=F)}.mat", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))
  z_path_repliseqTime = file.path(path_output, stringr::str_glue("results/{cond}_repliseqTime{format(binsize, scientific=F)}.tsv", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))
  z_path_bedgraph_repliseqTime = file.path(path_output, stringr::str_glue("results/{cond}_repliseqTime{format(binsize, scientific=F)}.bedgraph", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))
  z_path_repliseq = file.path(path_output, stringr::str_glue("results/{cond}_repliseq{format(binsize, scientific=F)}.tsv", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))

  z_wide = z %>%
    dplyr::mutate(repliseq_dtype="repliseq") %>%
    dplyr::arrange(SAMPLE_CONDITION, repliseq_chrom, repliseq_start, dplyr::desc(repliseq_fraction)) %>%
    dplyr::mutate(SAMPLE_FRACTION_COL=gsub("([A-Z]+)_(\\d+)", "\\1\\2", paste0(SAMPLE_PHASE, "_", SAMPLE_FRACTION))) %>%
    dplyr::mutate(SAMPLE_FRACTION_COL=factor(SAMPLE_FRACTION_COL, unique(SAMPLE_FRACTION_COL))) %>%
    dplyr::mutate(repliseq_value=round(repliseq_value, 3)) %>%
    reshape2::dcast(repliseq_chrom+repliseq_start+repliseq_end+repliseq_dtype~SAMPLE_FRACTION_COL, fill=0, value.var="repliseq_value")

  repliseq_df = z %>% dplyr::select(dplyr::matches("^repliseq_"))
  repliseqTime_df = repliseq_summarize(repliseq_df)

  writeLines(paste0("Writing RepliSeq TSV file '", z_path_repliseq, "'..."))
  readr::write_tsv(repliseq_df, file=z_path_repliseq)

  writeLines(paste0("Writing RepliSeq matrix file '", z_path_mat, "'..."))
  write.table(t(z_wide %>% dplyr::select(-repliseq_dtype)), file=z_path_mat, row.names=F, col.names=F, quote=F, sep="\t", na="")

  writeLines(paste0("Writing RepliSeq IGV file '", z_path_igv, "'..."))
  readr::write_lines("#type=GENE_EXPRESSION", file=z_path_igv)
  readr::write_tsv(z_wide, file=z_path_igv, append=T, col_names=T)

  writeLines(paste0("Writing RepliSeq time TSV file '", z_path_repliseqTime, "'..."))
  readr::write_tsv(repliseqTime_df, file=z_path_repliseqTime)

  repliseqTime_df %>%
    dplyr::select(repliseqTime_chrom, repliseqTime_start, repliseqTime_end, repliseqTime_avg) %>%
    readr::write_tsv(z_path_bedgraph_repliseqTime, col_names=F)
}

pipeline_coverage = function(path_metadata, path_output, binsizes, threads=1)
{
  path_metadata = "~/Workspace/Datasets/B400_RS_002/B400_RS_002_repliseq.tsv"
  path_output = "~/Workspace/Datasets/B400_RS_002"
  binsizes = 50e3
  threads = 1

  path_metadata = "~/Workspace/Datasets/B400_RS_001/B400_RS_001_repliseq.tsv"
  path_output = "~/Workspace/Datasets/B400_RS_001"
  binsizes = 50e3
  threads = 1

  path_metadata = "~/Workspace/Datasets/B400_RS_ALL/B400_RS_ALL_repliseq.tsv"
  path_output = "~/Workspace/Datasets/B400_RS_ALL"
  binsizes = 50e3
  threads = 1

  #
  # Check external executables
  #
  cmd_is_available(c("bedtools", "samtools"))

  path_bam = file.path(path_output, "alignments")
  if(!dir.exists(path_bam)) stop(paste0("Path '", path_bam, "' doesn't exist"))

  allowed_phases = c("G1", "S", "G2")
  arguments_df = suppressMessages(readr::read_tsv(path_metadata)) %>%
    dplyr::distinct(SAMPLE_NAME, SAMPLE_CONDITION, SAMPLE_PHASE, SAMPLE_FRACTION) %>%
    dplyr::mutate(SAMPLE_BAM=R.utils::getAbsolutePath(file.path(path_bam, paste0(SAMPLE_NAME, ".bam")), expandTilde=T), SAMPLE_BAM_EXIST=file.exists(SAMPLE_BAM)) %>%
    dplyr::mutate(SAMPLE_PHASE_I=match(SAMPLE_PHASE, allowed_phases)) %>%
    dplyr::group_by(SAMPLE_CONDITION) %>%
    dplyr::arrange(SAMPLE_CONDITION, SAMPLE_PHASE_I, SAMPLE_FRACTION) %>%
    dplyr::mutate(repliseq_fraction=1:dplyr::n()) %>%
    data.frame()

  error_message = ""
  if(!all(arguments_df$SAMPLE_PHASE %in% allowed_phases)) {
    incorrect_phase_df = arguments_df %>% dplyr::filter(!(SAMPLE_PHASE %in% allowed_phases))
    error_missing = paste0("Incorrect phases in metadata file:\n    ", paste0(incorrect_phase_df$SAMPLE_NAME, " (", incorrect_phase_df$SAMPLE_PHASE, ")", collapse="\n    "))
    error_message = paste0(error_message, error_missing, "\n-------------------------------\n")
  }

  if(!all(arguments_df$SAMPLE_BAM_EXIST)) {
    missing_bam = arguments_df %>% dplyr::filter(!SAMPLE_BAM_EXIST)
    error_missing = paste0("Missing BAM files in '", path_bam, "':\n    ", paste(missing_bam$SAMPLE_BAM, collapse="\n    "))
    error_message = paste0(error_message, error_missing, "\n-------------------------------\n")
  }

  if(nchar(error_message)>0) {
    stop(error_message)
  }

  if(!dir.exists(path_output)) cmd_run(paste0("mkdir ", path_output))
  if(!dir.exists(file.path(path_output, "results"))) cmd_run(paste0("mkdir ", file.path(path_output, "results")))
  if(!dir.exists(file.path(path_output, "bedgraph"))) cmd_run(paste0("mkdir ", file.path(path_output, "bedgraph")))


  writeLines("Calculating coverage...")
  binned_arguments_df = arguments_df %>% tidyr::crossing(data.frame(SAMPLE_BINSIZE=binsizes)) %>% data.frame()

  if(threads > 0) {
    x = sapply(1:nrow(binned_arguments_df), FUN=make_bedgraph)
  } else {
    cl = parallel::makeCluster(threads, outfile="")
    parallel::clusterExport(cl, "binned_arguments_df", envir=environment())
    x = parallel::parSapply(cl , 1:nrow(binned_arguments_df), FUN=make_bedgraph)
    parallel::stopCluster(cl)
  }

  bedgraph_cols = readr::cols(
    repliseq_chrom=readr::col_character(),
    repliseq_start=readr::col_double(),
    repliseq_end=readr::col_double(),
    repliseq_value_abs=readr::col_double()
  )

  writeLines("Building Repliseq matrix...")
  bedgraph_raw_df = arguments_df %>%
    # dplyr::filter(SAMPLE_PHASE=="S") %>%
    tidyr::crossing(data.frame(SAMPLE_BINSIZE=binsizes)) %>%
    dplyr::filter(!is.na(SAMPLE_NAME)) %>%
    dplyr::group_by(SAMPLE_BINSIZE, SAMPLE_NAME, SAMPLE_CONDITION, SAMPLE_PHASE, SAMPLE_FRACTION, repliseq_fraction) %>%
    dplyr::do(readr::read_tsv(file.path(path_output, stringr::str_glue("bedgraph/{sample}_bin{format(binsize, scientific=F)}.bdg", sample=.$SAMPLE_NAME[1], binsize=.$SAMPLE_BINSIZE[1])), col_names=names(bedgraph_cols$cols), col_types=bedgraph_cols)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(grepl("^chr([0-9XY]+)$", repliseq_chrom))

  # libfactors_df = bedgraph_raw_df %>%
  #   dplyr::arrange(dplyr::desc(repliseq_value_abs)) %>%
  #   dplyr::distinct(SAMPLE_NAME, SAMPLE_CONDITION, repliseq_chrom, repliseq_start, repliseq_end, .keep_all=T) %>%
  #   dplyr::group_by(SAMPLE_NAME, SAMPLE_CONDITION, repliseq_chrom) %>%
  #   dplyr::summarize(SAMPLE_BASELINE=median(repliseq_value_abs)) %>%
  #   dplyr::group_by(SAMPLE_CONDITION, repliseq_chrom) %>%
  #   dplyr::mutate(library_factor=ifelse(SAMPLE_BASELINE==0, 0, min(SAMPLE_BASELINE)/SAMPLE_BASELINE)) %>%
  #   dplyr::ungroup()
  #
  # libfactors_df1 = bedgraph_raw_df %>%
  #   dplyr::filter(grepl("chr\\d+", repliseq_chrom)) %>%
  #   dplyr::group_by(SAMPLE_NAME, SAMPLE_CONDITION) %>%
  #   dplyr::summarize(SAMPLE_BASELINE=sum(repliseq_value_abs)) %>%
  #   dplyr::group_by(SAMPLE_CONDITION) %>%
  #   dplyr::mutate(library_factor=min(SAMPLE_BASELINE)/SAMPLE_BASELINE)  %>%
  #   data.frame()

  libfactors_df = bedgraph_raw_df %>%
    dplyr::filter(grepl("chr\\d+", repliseq_chrom)) %>%
    dplyr::arrange(repliseq_start) %>%
    dplyr::group_by(SAMPLE_CONDITION, SAMPLE_BINSIZE, SAMPLE_NAME, repliseq_chrom) %>%
    dplyr::summarize(SAMPLE_BASELINE=na.omit(unique(zoo::rollapply(repliseq_value_abs, FUN=max, width=1e8/50000, by=1e6/50000, na.rm=T)))) %>%
    dplyr::group_by(SAMPLE_CONDITION, SAMPLE_BINSIZE, SAMPLE_NAME) %>%
    dplyr::filter(quantile(SAMPLE_BASELINE, 0.3, na.rm=T)<=SAMPLE_BASELINE & SAMPLE_BASELINE <= quantile(SAMPLE_BASELINE, 0.7, na.rm=T)) %>%
    dplyr::summarize(SAMPLE_BASELINE=mean(SAMPLE_BASELINE, na.rm=T))  %>%
    dplyr::group_by(SAMPLE_CONDITION, SAMPLE_BINSIZE)%>%
    dplyr::mutate(library_factor=min(SAMPLE_BASELINE)/SAMPLE_BASELINE)  %>%
    data.frame()

  # ggplot(libfactors_df %>% dplyr::inner_join(libfactors_df1, by=c("SAMPLE_NAME")), aes(library_factor.y, library_factor.x)) +
  #   geom_point() +
  #   geom_abline(xintercept=0, slope=1) +
  #   labs(x="Simpla adjustment", y="Max value adjustment") +
  #   ggrepel::geom_text_repel(aes(label=SAMPLE_NAME))
  #
  #
  # ggplot2::ggplot(libfactors_df) +
  #   ggridges::geom_density_ridges(ggplot2::aes(x=repliseq_value_abs, y=SAMPLE_NAME))

  # bedgraph_raw_df %>% dplyr::filter(SAMPLE_NAME=="DMSO S1" & repliseq_chrom=="chr6")

  # libfactors_cor = libfactors_df %>%
  #   reshape2::dcast(repliseq_chrom~SAMPLE_NAME, value.var="SAMPLE_BASELINE") %>%
  #   dplyr::select(-repliseq_chrom) %>%
  #   cor()
  # corrplot::corrplot.mixed(corr=libfactors_cor, lower = 'shade', upper = 'ellipse', order = 'hclust')
  #
  # ggplot(bedgraph_df %>% dplyr::filter(repliseq_value_abs<500)) +
  #   ggridges::geom_density_ridges(aes(x=repliseq_value_abs, y=SAMPLE_NAME)) +
  #   facet_wrap(~SAMPLE_CONDITION, scales="free")

   # softmax <- function(z){
   #   expsum = sum(exp(z))
   #   exp(z) / expsum
   # }


  bedgraph_df = bedgraph_raw_df %>%
    dplyr::inner_join(libfactors_df, by=intersect(colnames(libfactors_df), colnames(bedgraph_raw_df))) %>%
    dplyr::mutate(repliseq_value=library_factor*repliseq_value_abs) %>%
    # dplyr::group_by(SAMPLE_BINSIZE, SAMPLE_CONDITION, repliseq_chrom, repliseq_start, repliseq_end) %>%
    # dplyr::mutate(repliseq_value=repliseq_value/quantile(repliseq_value[SAMPLE_PHASE=="S" & repliseq_value>0], 0.1, na.rm=T))) %>%
    # dplyr::mutate(repliseq_value=ifelse(is.infinite(repliseq_value) | is.na(repliseq_value) | repliseq_value<0, 0, repliseq_value)) %>%
    dplyr::group_by(SAMPLE_BINSIZE, SAMPLE_CONDITION, repliseq_chrom, repliseq_start, repliseq_end) %>%
    dplyr::mutate(repliseq_value=repliseq_value/sum(repliseq_value)) %>%
    dplyr::group_by(SAMPLE_BINSIZE, SAMPLE_CONDITION) %>%
    dplyr::mutate(repliseq_value=repliseq_value/quantile(repliseq_value, 0.99, na.rm=T), repliseq_value=tidyr::replace_na(repliseq_value, 0)) %>%
    dplyr::ungroup()
    # dplyr::arrange(SAMPLE_CONDITION, repliseq_chrom, repliseq_start, dplyr::desc(repliseq_fraction)) %>%
    # dplyr::group_by(SAMPLE_CONDITION, repliseq_chrom, repliseq_fraction) %>%
    # dplyr::mutate(repliseq_value=smoother::smth.gaussian(repliseq_value, window=3, tails=T)) %>%
    # dplyr::ungroup() %>%
    # dplyr::arrange(SAMPLE_CONDITION, repliseq_chrom, repliseq_start, dplyr::desc(repliseq_fraction)) %>%
    # dplyr::group_by(SAMPLE_CONDITION, repliseq_chrom, repliseq_start) %>%
    # dplyr::mutate(repliseq_value=smoother::smth.gaussian(repliseq_value, window=2, tails=T)) %>%
    # dplyr::ungroup()

  #
  # Write IGV and MAT files
  #

  writeLines("Writing results in matrix and IGV formats...")
  bedgraph_output_df = bedgraph_df %>% dplyr::distinct(SAMPLE_BINSIZE, SAMPLE_CONDITION)
  if(threads == 1) {
    x = sapply(1:nrow(bedgraph_output_df), FUN=write_results, bedgraph_df=bedgraph_df)
  } else {
    cl = parallel::makeCluster(threads, outfile="")
    parallel::clusterExport(cl, c("bedgraph_output_df"), envir=environment())
    x = parallel::parSapply(cl , 1:nrow(bedgraph_output_df), FUN=write_results, bedgraph_df=bedgraph_df)
    parallel::stopCluster(cl)
  }
}


pipeline_coverage_cli = function() {
  option_list = list(
    optparse::make_option(c("-m", "--metadata"), dest="metadata", type="character", default=NULL, help="File with metadata", metavar="character"),
    optparse::make_option(c("-o", "--output-dir"), dest="out", type="character", default=".", help="Output directory [default= %default]", metavar="character"),
    optparse::make_option("--binsizes", type="character", default="50000", help="Comma separated bin size [default= %default]", metavar="character"),
    optparse::make_option(c("-t", "--threads"), type="integer", default=1, help="Output directory [default= %default]", metavar="integer")
  )
  opt_parser = optparse::OptionParser(option_list=option_list, usage="./coverage.R [options]", description = "Use the aligned replieq to build final replication program matrix")
  opt = optparse::parse_args(opt_parser)

  #
  # Parse CLI arguments
  #
  path_metadata = opt$metadata
  path_output = opt$out
  threads = opt$threads
  binsizes = sapply(unlist(strsplit(opt$binsizes, ", *")), as.numeric)

  # binsizes = c(50e3, 100e3)
  # # # path_output = "~/Workspace/repliseq/zhao"
  # # # path_metadata = "~/Workspace/Datasets/zhao_bmc_repliseq_2020/zhao_metadata.tsv"
  # path_output = "~/Workspace/Datasets/Repliseq"
  # path_metadata = "~/Workspace/Datasets/Repliseq/raw/B400_RS_001_23911/23911_meta.tsv"

  #
  # Print arguments
  #
  writeLines(paste0("metadata=", path_metadata))
  writeLines(paste0("out=", path_output))
  writeLines(paste0("threads=", threads))
  writeLines(paste0("binsizes=", paste(binsizes, collapse=",")))

# path_output = "~/Workspace/Datasets/Repliseq"

  #
  # Create missing folders
  #
  if(!dir.exists(path_output)) cmd_run(paste0("mkdir ", path_output))

  pipeline_coverage(path_metadata=path_metadata, path_output=path_output, binsizes=binsizes, threads=threads)
}

pipeline_coverage_cli()