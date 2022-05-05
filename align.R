#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(readr))
source("utils.R")

pipeline_align = function(path_metadata, path_fastq, path_database, path_output, threads=1, md5_suffix=NULL)
{
  if(!dir.exists(file.path(path_output, "alignments"))) cmd_run(paste0("mkdir ", file.path(path_output, "alignments")))

  files_df = data.frame(FASTQ_PATH=list.files(path_fastq, recursive=T, full.names=T)) %>% dplyr::mutate(FASTQ_FILE=basename(FASTQ_PATH))
  meta_df = readr::read_tsv(path_metadata) %>%
    dplyr::mutate(SAMPLE_NUMBER=match(SAMPLE_NAME, unique(SAMPLE_NAME)))

  #
  # Check external executables
  #
  cmd_is_available(c("bedtools", "samtools", "bowtie2"))

  #
  # Check for errors in arguments
  #
  error_message = ""
  if(!is.null(md5_suffix)) {
    suppressMessages(library(openssl))
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

  arguments_df = arguments1_df %>%
    reshape2::dcast(SAMPLE_NAME+SAMPLE_NUMBER+SAMPLE_CONDITION+SAMPLE_FRACTION~FASTQ_PAIR, value.var="FASTQ_PATH") %>%
    dplyr::arrange(SAMPLE_NUMBER)

  #
  # Run pipeline
  #
  writeLines("Aligning sequenced reads with reference...")
  for(i in 1:nrow(arguments_df)) {
    writeLines(paste0(i, "/", nrow(arguments_df), ": ", arguments_df$SAMPLE_NAME[i]))

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
}

pipeline_align_cli = function()
{
    option_list = list(
    optparse::make_option(c("-m", "--metadata"), dest="metadata", type="character", default=NULL, help="File with metadata", metavar="character"),
    optparse::make_option(c("-g", "--genome-dir"), dest="genome_dir", type="character", default=NULL, help="Path to genome {align}", metavar="character"),
    optparse::make_option(c("-o", "--output-dir"), dest="out", type="character", default=".", help="Output directory [default= %default]", metavar="character"),
    optparse::make_option("--binsizes", type="character", default="50000", help="Comma separated bin size [default= %default]", metavar="character"),
    optparse::make_option("--fastq-dir", dest="fastq", type="character", default=".", help="Path to folder with fastq files [default= %default]", metavar="character"),
    optparse::make_option("--md5-suffix", dest="md5", type="character", default=NULL, help="MD5 files suffix. Leave out this argument for not validating fastq files", metavar="character"),
    optparse::make_option(c("-t", "--threads"), type="integer", default=1, help="Output directory [default= %default]", metavar="integer")
  )
  opt_parser = optparse::OptionParser(option_list=option_list, usage="./coverage.R [options]", description = "Align short reads from FASTQ files to reference genome")
  opt = optparse::parse_args(opt_parser)

  #
  # Parse CLI arguments
  #
  path_metadata = opt$metadata
  path_database = opt$genome_dir
  path_output = opt$out
  path_fastq = opt$fastq
  threads = opt$threads
  md5_suffix = opt$md5

  # threads = 32
  # path_database = "~/Workspace/genomes/mm10/mm10"
  # md5_suffix = NULL

  # path_output = "~/Workspace/repliseq/zhao"
  # path_fastq = "~/Workspace/Datasets/zhao_bmc_repliseq_2020/zhao_raw"
  # path_metadata = "~/Workspace/Datasets/zhao_bmc_repliseq_2020/zhao_metadata.tsv"

  # path_output = "~/Workspace/Datasets/Repliseq"
  # path_fastq = "~/Workspace/Datasets/Repliseq/raw/B400_RS_001_23911"
  # path_metadata = "~/Workspace/Datasets/Repliseq/raw/B400_RS_001_23911/23911_meta.tsv"

  # path_output = "~/Workspace/HighRes_RepliSeq/test_repliseq"
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
  writeLines(paste0("md5-suffix=", paste(md5_suffix, collapse=",")))



  #
  # Create missing folders
  #
  if(!dir.exists(path_output)) cmd_run(paste0("mkdir ", path_output))

  pipeline_align(path_metadata=path_metadata, path_fastq=path_fastq, path_database=path_database, path_output=path_output, threads=threads, md5_suffix=md5_suffix)
}

pipeline_align_cli()