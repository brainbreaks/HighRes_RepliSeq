setwd("~/Workspace/HighRes_RepliSeq")
library(dplyr)
library(openssl)
library(readr)
library(spatstat)
library(optparse)
source("utils.R")

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

cmd_run = function(cmd) {
  start_time = Sys.time()
  print(cmd)
  system(cmd)
  end_time = Sys.time()
  print(end_time - start_time)
}

pipeline = function() {
  option_list = list(
    optparse::make_option(c("-m", "--metadata"), dest="metadata", type="character", default=NULL, help="File with metadata", metavar="character"),
    optparse::make_option(c("-g", "--genome-dir"), dest="genome_dir", type="character", default=NULL, help="Path to genome annotation", metavar="character"),
    optparse::make_option(c("-o", "--output-dir"), dest="out", type="character", default=".", help="Output directory [default= %default]", metavar="character"),
    optparse::make_option("--binsizes", type="character", default="50000", help="Comma separated bin size [default= %default]", metavar="character"),
    optparse::make_option("--fastq-dir", dest="fastq", type="character", default=".", help="Path to folder with fastq files [default= %default]", metavar="character"),
    optparse::make_option("--md5-suffix", dest="md5", type="character", default=NULL, help="MD5 files suffix. Leave out this argument for not validating fastq files", metavar="character"),
    optparse::make_option(c("-t", "--threads"), type="integer", default=1, help="Output directory [default= %default]", metavar="integer")
  )
  opt_parser = optparse::OptionParser(option_list=option_list)
  opt = optparse::parse_args(opt_parser)

  #
  # Parse CLI arguments
  #
  path_metadata = opt$metadata
  path_database = opt$genome_dir
  path_output = opt$out
  path_fastq = opt$fastq
  threads = opt$threads
  binsizes = sapply(strsplit(opt$binsizes, ", *"), as.numeric)
  md5_suffix = opt$md5

  # threads = 32
  # binsizes = c(50e3, 100e3)
  # md5_suffix = NULL
  # path_database = "~/Workspace/genomes/mm10/mm10"
  # # path_output = "~/Workspace/repliseq/zhao"
  # # path_fastq = "~/Workspace/Datasets/zhao_bmc_repliseq_2020/zhao_raw"
  # # path_metadata = "~/Workspace/Datasets/zhao_bmc_repliseq_2020/zhao_metadata.tsv"
  # path_output = "~/Workspace/repliseq/wei"
  # path_fastq = "~/Workspace/Datasets/Repliseq/raw/B400_RS_001_23911"
  # path_metadata = "~/Workspace/Datasets/Repliseq/raw/B400_RS_001_23911/23911_meta.tsv"

  #
  # Print arguments
  #
  writeLines(paste0("metadata=", path_metadata))
  writeLines(paste0("genome=", path_database))
  writeLines(paste0("out=", path_output))
  writeLines(paste0("fastq=", path_fastq))
  writeLines(paste0("threads=", threads))
  writeLines(paste0("binsizes=", paste(binsizes, collapse=",")))
  writeLines(paste0("md5-suffix=", paste(md5_suffix, collapse=",")))



  #
  # Create missing folders
  #
  if(!dir.exists(path_output)) cmd_run(paste0("mkdir ", path_output))
  if(!dir.exists(file.path(path_output, "results"))) cmd_run(paste0("mkdir ", file.path(path_output, "results")))
  if(!dir.exists(file.path(path_output, "bedgraph"))) cmd_run(paste0("mkdir ", file.path(path_output, "bedgraph")))
  if(!dir.exists(file.path(path_output, "alignments"))) cmd_run(paste0("mkdir ", file.path(path_output, "alignments")))

  files_df = data.frame(FASTQ_PATH=list.files(path_fastq, recursive=T, full.names=T)) %>% dplyr::mutate(FASTQ_FILE=basename(FASTQ_PATH))
  meta_df = readr::read_tsv(path_metadata)

  #
  # Check external executables
  #
  cmd_is_available(c("bedtools", "samtools", "bowtie2"))

  #
  # Check for errors in arguments
  #
  error_message = ""
  if(!is.null(md5_suffix)) {
    files_df = files_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(md5=as.character(openssl::md5(file(FASTQ_PATH, open = "rb"))))

    files_df = files_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(md5=as.character(md5)) %>%
      dplyr::mutate(md5_expected=gsub(" .*", "", readr::read_tsv(paste0(FASTQ_PATH, md5_suffix), col_names=F)[[1]])) %>%
      dplyr::ungroup()

    error_md5 = paste0("HASH doesn't match for: ", paste(basename(files_df %>% dplyr::filter(md5!=md5_expected) %>% .$FASTQ_FILE), collapse=", "))
    error_message = paste0(error_message, error_md5, "\n-------------------------------\n")
  }

  duplicate_sample = meta_df %>%
    dplyr::group_by(SAMPLE_NAME) %>%
    dplyr::filter(length(unique(SAMPLE_CONDITION, SAMPLE_FRACTION))>1) %>%
    dplyr::ungroup() %>%
    .$SAMPLE_NAME
  if(length(duplicate_sample)>0) {
    error_dupsample = paste0("Duplicate sample names in '", path_metadata, "':\n    ", paste(duplicate_sample, collapse="\n    "))
    error_message = paste0(error_message, error_dupsample, "\n-------------------------------\n")
  }

  duplicate_fastq = meta_df %>%
    dplyr::group_by(FASTQ_FILE) %>%
    dplyr::filter(dplyr::n()>1) %>%
    dplyr::ungroup() %>%
    .$FASTQ_FILE
  if(length(duplicate_fastq)>0) {
    error_dupfile = paste0("Duplicate fasta files in '", path_fastq, "':\n    ", paste(duplicate_fastq, collapse="\n    "))
    error_message = paste0(error_message, error_dupfile, "\n-------------------------------\n")
  }

  arguments1_df = meta_df %>% dplyr::left_join(files_df, by="FASTQ_FILE")
  if(any(is.na(arguments1_df$FASTQ_PATH))) {
    missing_fasta = arguments1_df %>%
      dplyr::filter(is.na(FASTQ_PATH)) %>%
      dplyr::group_by(SAMPLE_NAME) %>%
      dplyr::summarize(path=paste0(SAMPLE_NAME[1], ": ", paste(FASTQ_FILE, collapse=", ")))
    error_missing = paste0("Missing fasta files in '", path_fastq, "':\n    ", paste(missing_fasta$path, collapse="\n    "))
    error_message = paste0(error_message, error_missing, "\n-------------------------------\n")
  }

  if(nchar(error_message)>0) {
    stop(error_message)
  }

  arguments_df = arguments1_df %>% reshape2::dcast(SAMPLE_NAME+SAMPLE_CONDITION+SAMPLE_FRACTION~FASTQ_PAIR, value.var="FASTQ_PATH")

  #
  # Run pipeline
  #
  writeLines("Aligning sequenced reads with reference...")
  for(i in 1:nrow(arguments_df)) {
    print(paste0(i, "/", nrow(arguments_df), ": ", arguments_df$SAMPLE_NAME[i]))

    path_bam = file.path(path_output, stringr::str_glue("alignments/{sample}.bam", sample=arguments_df$SAMPLE_NAME[i]))
    path_bai = paste0(path_bam, ".bai")
    path_tmpbam = file.path(path_output, stringr::str_glue("alignments/{sample}_tmp.bam", sample=arguments_df$SAMPLE_NAME[i]))

    if(!file.exists(path_bam)) {
      cmd_bowtie = stringr::str_glue("bowtie2 -q -p {threads} --phred33 --end-to-end --sensitive --no-mixed --no-discordant --reorder --no-unal -x {db} -1  {r1} -2  {r2} | samtools view -b -@ {threads} - > {bam}", bam=path_bam, db=path_database, r1=arguments_df$`1`[i], r2=arguments_df$`2`[i], sample=arguments_df$SAMPLE_NAME[i], threads=threads)
      cmd_run(cmd_bowtie)

      cmd_mapq = stringr::str_glue("samtools view -bSq 20 -@ {threads} {bam} > {tmp}; mv {tmp} {bam}", tmp=path_tmpbam, bam=path_bam, threads=threads)
      cmd_run(cmd_mapq)

      cmd_fixmate = stringr::str_glue("samtools fixmate -@ {threads} -m {bam} {tmp}; mv {tmp} {bam}", tmp=path_tmpbam, bam=path_bam, threads=threads)
      cmd_run(cmd_fixmate)

      cmd_sort = stringr::str_glue("samtools sort -@ {threads} -o {tmp} {bam}; mv {tmp} {bam}", tmp=path_tmpbam, bam=path_bam, threads=threads)
      cmd_run(cmd_sort)

      cmd_unique = stringr::str_glue("samtools markdup -@ {threads} -r {bam} {tmp}; mv {tmp} {bam}", tmp=path_tmpbam, bam=path_bam, threads=threads)
      cmd_run(cmd_unique)
    } else {
      writeLines(paste0("BAM file '", path_bam, "' exists, skipping..."))
    }

    if(!file.exists(path_bai)) {
      cmd_index = stringr::str_glue("samtools index -@ {threads} -b {bam}", bam=path_bam, threads=threads)
      cmd_run(cmd_index)
    } else {
      writeLines(paste0("BAM index file '", path_bai, "' exists, skipping..."))
    }
  }

  #
  # Calculate coverage
  #
  writeLines("Calculating coverage...")
  for(i in 1:nrow(arguments_df)) {
    print(paste0(i, "/", nrow(arguments_df), ": ", arguments_df$SAMPLE_NAME[i]))

    path_bam = file.path(path_output, stringr::str_glue("alignments/{sample}.bam", sample=arguments_df$SAMPLE_NAME[i]))
    path_chromsizes = file.path(path_output, stringr::str_glue("{sample}.sizes", sample=arguments_df$SAMPLE_NAME[i]))
    cmd_chromizes = stringr::str_glue("samtools idxstats {bam} | cut -f 1-2 | grep -vP '\\t0' > {chromsizes}", bam=path_bam, chromsizes=path_chromsizes)
    cmd_run(cmd_chromizes)

    for(binsize in binsizes) {
      path_tiles = file.path(path_output, stringr::str_glue("bedgraph/genome_{format(binsize, scientific=F)}tiles.bed", binsize=binsize))
      path_tiles_tmp = file.path(path_output, stringr::str_glue("bedgraph/genome_{format(binsize, scientific=F)}tiles_tmp.bed", binsize=binsize))
      path_bedgraph = file.path(path_output, stringr::str_glue("bedgraph/{sample}_bin{format(binsize, scientific=F)}.bdg", sample=arguments_df$SAMPLE_NAME[i], binsize=binsize))

      if(!file.exists(path_bedgraph)) {
      cmd_makewindows = stringr::str_glue("bedtools makewindows -w {format(binsize, scientific=F)} -s {format(binsize, scientific=F)} -g {chromsizes} > {tiles}; bedtools sort -faidx {chromsizes} -i {tiles} > {tiles_tmp}; mv {tiles_tmp} {tiles}", chromsizes=path_chromsizes, binsize=binsize, tiles=path_tiles, tiles_tmp=path_tiles_tmp)
      cmd_run(cmd_makewindows)
      cmd_coverage = stringr::str_glue("bedtools intersect -abam {tiles} -b {bam} -c -bed -g {chromsizes} -sorted > {bedgraph}", sample=arguments_df$SAMPLE_NAME[i], chromsizes=path_chromsizes, bam=path_bam, tiles=path_tiles, bedgraph=path_bedgraph)
      cmd_run(cmd_coverage)
      cmd_run(paste("rm", path_tiles))
      } else {
        writeLines(paste0("Bedgraph file '", path_bedgraph, "' exists, skipping..."))
      }
    }

    cmd_run(paste("rm", path_chromsizes))
  }

  writeLines("Building Repliseq matrix...")
  bedgraph_df = meta_df %>%
    dplyr::distinct(SAMPLE_NAME, SAMPLE_CONDITION, repliseq_fraction=SAMPLE_FRACTION) %>%
    tidyr::crossing(data.frame(SAMPLE_BINSIZE=binsizes)) %>%
    dplyr::filter(!is.na(SAMPLE_NAME)) %>%
    dplyr::group_by(SAMPLE_BINSIZE, SAMPLE_NAME, SAMPLE_CONDITION, repliseq_fraction) %>%
    dplyr::do(readr::read_tsv(file.path(path_output, stringr::str_glue("bedgraph/{sample}_bin{format(binsize, scientific=F)}.bdg", sample=.$SAMPLE_NAME[1], binsize=.$SAMPLE_BINSIZE[1])), col_names=c("repliseq_chrom", "repliseq_start", "repliseq_end", "repliseq_value_abs"))) %>%
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
  bedgraph_smooth_df %>%
    dplyr::group_by(SAMPLE_BINSIZE, SAMPLE_CONDITION) %>%
    dplyr::do((function(z){
      zz <<- z
      # z = bedgraph_smooth_df %>% dplyr::filter(SAMPLE_CONDITION=="DMSO")
      z_path_igv = file.path(path_output, stringr::str_glue("results/{cond}_repliseq{format(binsize, scientific=F)}.igv", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))
      z_path_mat = file.path(path_output, stringr::str_glue("results/{cond}_repliseq{format(binsize, scientific=F)}.mat", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))
      z_path_repliseq = file.path(path_output, stringr::str_glue("results/{cond}_repliseq{format(binsize, scientific=F)}.tsv", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))
      z_path_repliseqTime = file.path(path_output, stringr::str_glue("results/{cond}_repliseqTime{format(binsize, scientific=F)}.tsv", cond=z$SAMPLE_CONDITION[1], binsize=z$SAMPLE_BINSIZE[1]))

      z_wide = z %>%
        dplyr::mutate(repliseq_dtype="repliseq") %>%
        dplyr::arrange(repliseq_chrom, repliseq_start) %>%
        reshape2::dcast(repliseq_chrom+repliseq_start+repliseq_end+repliseq_dtype~repliseq_fraction, fill=0, value.var="repliseq_value")

      repliseq_df = z %>% dplyr::select(dplyr::matches("^repliseq_"))
      repliseqTime_df = repliseq_summarize(repliseq_df)

      writeLines(paste0("Writing RepliSeq matrix file '", z_path_mat, "'..."))
      readr::write_tsv(as.data.frame(t(z_wide %>% dplyr::select(-repliseq_dtype))), file=z_path_mat, append=F, col_names=F)

      writeLines(paste0("Writing RepliSeq IGV file '", z_path_igv, "'..."))
      readr::write_lines("#type=GENE_EXPRESSION", file=z_path_igv)
      readr::write_tsv(z_wide, file=z_path_igv, append=T, col_names=T)

      writeLines(paste0("Writing RepliSeq TSV file '", z_path_repliseq, "'..."))
      readr::write_tsv(repliseq_df, file=z_path_repliseq)

      writeLines(paste0("Writing RepliSeq time TSV file '", z_path_repliseqTime, "'..."))
      readr::write_tsv(repliseqTime_df, file=z_path_repliseqTime)

      data.frame()
    })(.))
}

pipeline()

