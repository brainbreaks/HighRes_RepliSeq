#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
options(dplyr.summarise.inform=F)
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(spatstat))
suppressMessages(require(parallel))
source("utils.R")


pipeline_coverage = function(path_metadata, path_output, binsizes, threads=1)
{
  #
  # Check external executables
  #
  cmd_is_available(c("bedtools", "samtools", "bowtie2"))

  path_bam = file.path(path_output, "alignments")
  if(!dir.exists(path_bam)) stop(paste0("Path '", path_bam, "' doesn't exist"))

  arguments_df = suppressMessages(readr::read_tsv(path_metadata)) %>%
    dplyr::distinct(SAMPLE_NAME, SAMPLE_CONDITION, repliseq_fraction=SAMPLE_FRACTION) %>%
    dplyr::mutate(SAMPLE_NUMBER=match(SAMPLE_NAME, unique(SAMPLE_NAME))) %>%
    dplyr::mutate(SAMPLE_BAM=file.path(path_bam, paste0(SAMPLE_NAME, ".bam"))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(SAMPLE_BAM_EXIST=file.exists(SAMPLE_BAM))

  error_message = ""
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
  make_bedgraph = function(i) {
    suppressMessages(library(dplyr))
    source("utils.R")

    r = paste0(i, "/", nrow(binned_arguments_df))
    writeLines(paste0("\n==========================================\n", r, ": ", binned_arguments_df$SAMPLE_NAME[i], " / bin", binned_arguments_df$SAMPLE_BINSIZE[i], "\n=========================================="))

    path_tiles = file.path(path_output, stringr::str_glue("bedgraph/genome_{format(binsize, scientific=F)}tiles.bed", binsize=binned_arguments_df$SAMPLE_BINSIZE[i]))
    path_tiles_tmp = file.path(path_output, stringr::str_glue("bedgraph/genome_{format(binsize, scientific=F)}tiles_tmp.bed", binsize=binned_arguments_df$SAMPLE_BINSIZE[i]))
    path_bedgraph = file.path(path_output, stringr::str_glue("bedgraph/{sample}_bin{format(binsize, scientific=F)}.bdg", sample=binned_arguments_df$SAMPLE_NAME[i], binsize=binned_arguments_df$SAMPLE_BINSIZE[i]))
    path_chromsizes = file.path(path_output, stringr::str_glue("bedgraph/{sample}_{format(binsize, scientific=F)}.sizes", sample=binned_arguments_df$SAMPLE_NAME[i], binsize=binned_arguments_df$SAMPLE_BINSIZE[i]))

    if(!file.exists(path_bedgraph)) {
      cmd_chromizes = stringr::str_glue("samtools idxstats {bam} | cut -f 1-2 | grep -vP '\\t0' > {chromsizes}", bam=binned_arguments_df$SAMPLE_BAM[i], chromsizes=path_chromsizes)
      cmd_run(cmd_chromizes)

      cmd_makewindows = stringr::str_glue("bedtools makewindows -w {format(binsize, scientific=F)} -s {format(binsize, scientific=F)} -g {chromsizes} > {tiles}; bedtools sort -faidx {chromsizes} -i {tiles} > {tiles_tmp}; mv {tiles_tmp} {tiles}", chromsizes=path_chromsizes, binsize=binned_arguments_df$SAMPLE_BINSIZE[i], tiles=path_tiles, tiles_tmp=path_tiles_tmp)
      cmd_run(cmd_makewindows)
      cmd_coverage = stringr::str_glue("bedtools intersect -abam {tiles} -b {bam} -c -bed -g {chromsizes} -sorted > {bedgraph}", sample=binned_arguments_df$SAMPLE_NAME[i], chromsizes=path_chromsizes, bam=binned_arguments_df$SAMPLE_BAM[i], tiles=path_tiles, bedgraph=path_bedgraph)
      cmd_run(cmd_coverage)
      file.remove(path_tiles)
      file.remove(path_tiles_tmp)
      file.remove(path_chromsizes)
      # cmd_run(paste("rm", path_tiles))
      # cmd_run(paste("rm", path_tiles_tmp))
      # cmd_run(paste("rm", path_chromsizes))
    } else {
      writeLines(paste0("Bedgraph file '", path_bedgraph, "' exists, skipping..."))
    }
  }

  # if(threads == 0) {
  x = sapply(1:nrow(binned_arguments_df), FUN=make_bedgraph)
  # } else {
  #   cl = parallel::makeCluster(threads, outfile="")
  #   parallel::clusterExport(cl, "binned_arguments_df", envir=environment())
  #   x = parallel::parSapply(cl , 1:nrow(binned_arguments_df), FUN=make_bedgraph)
  #   parallel::stopCluster(cl)
  # }

  bedgraph_cols = readr::cols(
    repliseq_chrom=readr::col_character(),
    repliseq_start=readr::col_double(),
    repliseq_end=readr::col_double(),
    repliseq_value_abs=readr::col_double()
  )

  writeLines("Building Repliseq matrix...")
  bedgraph_df = arguments_df %>%
    tidyr::crossing(data.frame(SAMPLE_BINSIZE=binsizes)) %>%
    dplyr::filter(!is.na(SAMPLE_NAME)) %>%
    dplyr::group_by(SAMPLE_BINSIZE, SAMPLE_NAME, SAMPLE_CONDITION, repliseq_fraction) %>%
    dplyr::do(readr::read_tsv(file.path(path_output, stringr::str_glue("bedgraph/{sample}_bin{format(binsize, scientific=F)}.bdg", sample=.$SAMPLE_NAME[1], binsize=.$SAMPLE_BINSIZE[1])), col_names=names(bedgraph_cols$cols), col_types=bedgraph_cols)) %>%
    dplyr::group_by(SAMPLE_BINSIZE, SAMPLE_NAME) %>%
    dplyr::mutate(library_size=dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(library_factor=min(library_size)/library_size) %>%
    dplyr::group_by(SAMPLE_BINSIZE, SAMPLE_CONDITION, repliseq_chrom, repliseq_start, repliseq_end) %>%
    dplyr::mutate(repliseq_start=repliseq_start) %>%
    dplyr::mutate(repliseq_value=library_factor*repliseq_value_abs, repliseq_value=repliseq_value/max(repliseq_value)) %>%
    dplyr::mutate(repliseq_value=tidyr::replace_na(repliseq_value, 0)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(repliseq_fraction)) %>%
    dplyr::filter(grepl("^chr([0-9XY]+)$", repliseq_chrom))

  writeLines("Smoothening repliseq matrix...")
  bedgraph_smooth_df = bedgraph_df %>%
    dplyr::group_by(SAMPLE_BINSIZE, SAMPLE_CONDITION) %>%
    dplyr::do((function(z){
      zz<<-z
      # adads()
      # z = bedgraph_df %>% dplyr::filter(SAMPLE_CONDITION=="DMSO" & SAMPLE_BINSIZE==100000)
      z_wide = z %>%
        dplyr::arrange(repliseq_chrom, repliseq_start) %>%
        reshape2::dcast(repliseq_chrom+repliseq_start+repliseq_end~repliseq_fraction, fill=0, value.var="repliseq_value")
      z_mat = z_wide %>% dplyr::select(-(repliseq_chrom:repliseq_end)) %>% as.matrix()
      z_fractions = colnames(z_mat)
      z_mat = as.matrix(spatstat.core::blur(spatstat.geom::as.im(z_mat), sigma=1))
      z_mat = t(apply(z_mat, 1, scales::rescale))
      z_wide[,z_fractions] = z_mat
      z_long = z_wide %>%
        reshape2::melt(measure.vars=z_fractions, value.name="repliseq_value") %>%
        dplyr::mutate(variable=as.numeric(as.character(variable)))
      z %>%
        dplyr::rename(repliseq_value_unsmoothed="repliseq_value") %>%
        dplyr::inner_join(z_long, by=c("repliseq_fraction"="variable", "repliseq_chrom", "repliseq_start", "repliseq_end"))
    })(.))

  #
  # Write IGV and MAT files
  #
  writeLines("Writing results in matrix and IGV formats...")
  bedgraph_smooth_output_df = bedgraph_smooth_df %>% dplyr::distinct(SAMPLE_BINSIZE, SAMPLE_CONDITION)

  write_results = function(i) {
    suppressMessages(library(dplyr))
    source("utils.R")

    z = bedgraph_smooth_df %>%
      dplyr::inner_join(bedgraph_smooth_output_df[i,,drop=F], by=c("SAMPLE_CONDITION", "SAMPLE_BINSIZE")) %>%
      dplyr::arrange(repliseq_chrom, repliseq_start, repliseq_end, dplyr::desc(repliseq_fraction))
    writeLines(paste0("\n==========================================\n", i, "/", nrow(bedgraph_smooth_output_df), ": ", z$SAMPLE_NAME[1], " / bin", z$SAMPLE_BINSIZE[1], "\n=========================================="))

    # z = bedgraph_smooth_df %>% dplyr::filter(SAMPLE_CONDITION=="DMSO")
    z_path_igv = file.path(path_output, stringr::str_glue("results/{cond}_repliseq{format(binsize, scientific=F)}.igv", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))
    z_unsmoothed_igv_path = file.path(path_output, stringr::str_glue("results/{cond}_smooth_repliseq{format(binsize, scientific=F)}.igv", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))

    z_path_mat = file.path(path_output, stringr::str_glue("results/{cond}_smooth_repliseq{format(binsize, scientific=F)}.mat", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))
    z_unsmoothed_path_mat = file.path(path_output, stringr::str_glue("results/{cond}_smooth_repliseq{format(binsize, scientific=F)}.mat", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))

    z_path_repliseqTime = file.path(path_output, stringr::str_glue("results/{cond}_repliseqTime{format(binsize, scientific=F)}.tsv", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))
    z_unsmoothed_path_repliseqTime = file.path(path_output, stringr::str_glue("results/{cond}_repliseqTime{format(binsize, scientific=F)}.tsv", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))

    z_path_repliseq = file.path(path_output, stringr::str_glue("results/{cond}_repliseq{format(binsize, scientific=F)}.tsv", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))

    z_wide = z %>%
      dplyr::mutate(repliseq_dtype="repliseq") %>%
      dplyr::arrange(repliseq_chrom, repliseq_start) %>%
      reshape2::dcast(repliseq_chrom+repliseq_start+repliseq_end+repliseq_dtype~repliseq_fraction, fill=0, value.var="repliseq_value")

    z_unsmoothed_wide = z %>%
      dplyr::mutate(repliseq_dtype="repliseq") %>%
      dplyr::arrange(repliseq_chrom, repliseq_start) %>%
      reshape2::dcast(repliseq_chrom+repliseq_start+repliseq_end+repliseq_dtype~repliseq_fraction, fill=0, value.var="repliseq_value_unsmoothed")

    repliseq_df = z %>% dplyr::select(dplyr::matches("^repliseq_"))
    repliseqTime_df = repliseq_summarize(repliseq_df)

    repliseq_unsmoothed_df = z %>% dplyr::select(dplyr::matches("^repliseq_")) %>% dplyr::mutate(repliseq_value=repliseq_value_unsmoothed)
    repliseqTime_unsmoothed_df = repliseq_summarize(repliseq_unsmoothed_df)

    writeLines(paste0("Writing RepliSeq TSV file '", z_path_repliseq, "'..."))
    readr::write_tsv(repliseq_df, file=z_path_repliseq)

    writeLines(paste0("Writing RepliSeq matrix file '", z_path_mat, "'..."))
    readr::write_tsv(as.data.frame(t(z_wide %>% dplyr::select(-repliseq_dtype))), file=z_path_mat, append=F, col_names=F)

    writeLines(paste0("Writing unsmoothed RepliSeq matrix file '", z_unsmoothed_path_mat, "'..."))
    readr::write_tsv(as.data.frame(t(z_unsmoothed_wide %>% dplyr::select(-repliseq_dtype))), file=z_unsmoothed_path_mat, append=F, col_names=F)

    writeLines(paste0("Writing RepliSeq IGV file '", z_unsmoothed_igv_path, "'..."))
    readr::write_lines("#type=GENE_EXPRESSION", file=z_unsmoothed_igv_path)
    readr::write_tsv(z_unsmoothed_wide, file=z_unsmoothed_igv_path, append=T, col_names=T)

    writeLines(paste0("Writing smoothed RepliSeq IGV file '", z_path_igv, "'..."))
    readr::write_lines("#type=GENE_EXPRESSION", file=z_path_igv)
    readr::write_tsv(z_wide, file=z_path_igv, append=T, col_names=T)

    writeLines(paste0("Writing RepliSeq time TSV file '", z_unsmoothed_path_repliseqTime, "'..."))
    readr::write_tsv(repliseqTime_unsmoothed_df, file=z_unsmoothed_path_repliseqTime)

    writeLines(paste0("Writing smoothed RepliSeq time TSV file '", z_path_repliseqTime, "'..."))
    readr::write_tsv(repliseqTime_df, file=z_path_repliseqTime)
  }

  if(threads == 1) {
    x = sapply(1:nrow(bedgraph_smooth_output_df), FUN=make_bedgraph)
  } else {
    cl = parallel::makeCluster(threads, outfile="")
    parallel::clusterExport(cl, c("bedgraph_smooth_output_df", "bedgraph_smooth_df"), envir=environment())
    x = parallel::parSapply(cl , 1:nrow(bedgraph_smooth_output_df), FUN=write_results)
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