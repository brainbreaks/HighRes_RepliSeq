#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(openssl))
suppressMessages(library(foreach))
source("utils.R")

pipeline_checksum = function(path_metadata, path_fastq, threads=1)
{
  # path_metadata = "/mnt/sda1/Workspace/B400_RS_002/B400_RS_002_repliseq.tsv"
  # path_fastq = "~/Workspace/B400_RS_002/fastq"
  # threads = 30

  files_df = data.frame(FASTQ_PATH=list.files(path_fastq, recursive=T, full.names=T)) %>% dplyr::mutate(FASTQ_FILE=basename(FASTQ_PATH))
  meta_df = readr::read_tsv(path_metadata) %>%
    dplyr::mutate(SAMPLE_NUMBER=match(SAMPLE_NAME, unique(SAMPLE_NAME)))

  #
  # Check for errors in arguments
  #
  error_message = ""
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

  #
  # Check MD5
  #
  if(!("MD5" %in% colnames(meta_df))) {
    error_message = paste0(error_message, "Metadata doesn't have MD5 column\n-------------------------------\n")
  }

  if(threads > 1) {
    doParallel::registerDoParallel(cores=threads)
    md5_df = foreach(f=1:nrow(files_df)) %dopar% {
      suppressMessages(library(openssl))
      data.frame(FASTQ_PATH=files_df[f,"FASTQ_PATH"], md5=as.character(as.character(openssl::md5(file(files_df[f,"FASTQ_PATH"], open = "rb")))))
    }
    md5_df = do.call(dplyr::bind_rows, md5_df)
  } else {
    md5_df = files_df %>%
      dplyr::distinct(FASTQ_PATH) %>%
      dplyr::rowwise() %>%
      dplyr::summarize(md5=as.character(openssl::md5(file(FASTQ_PATH, open = "rb")))) %>%
      dplyr::ungroup()
  }
  meta_md5_df = meta_df %>%
    inner_join(files_df, by="FASTQ_FILE") %>%
    dplyr::inner_join(md5_df, by=c("FASTQ_PATH")) %>%
    dplyr::mutate(md5_correct=MD5==md5)
  if(any(!meta_md5_df$md5_correct)) {
    error_md5 = paste0("HASH doesn't match for: ", paste(basename(meta_md5_df %>% dplyr::filter(!md5_correct) %>% .$FASTQ_FILE), collapse=", "))
    error_message = paste0(error_message, error_md5, "\n-------------------------------\n")
  }

  if(nchar(error_message)>0) {
    stop(error_message)
  } else {
    writeLines("Didn't find any problems with files...")
  }
}


pipeline_checksum_cli = function()
{
  option_list = list(
    optparse::make_option(c("-m", "--metadata"), dest="metadata", type="character", default=NULL, help="File with metadata", metavar="character"),
    optparse::make_option(c("-i", "--fastq-dir"), dest="fastq", type="character", default=".", help="Path to folder with fastq files [default= %default]", metavar="character"),
    optparse::make_option(c("-t", "--threads"), type="integer", default=1, help="Output directory [default= %default]", metavar="integer")
  )
  opt_parser = optparse::OptionParser(option_list=option_list, usage="./checksum.R [options]", description = "Use MD5 column in metadata file to validate fastq files")
  opt = optparse::parse_args(opt_parser)

  #
  # Parse CLI arguments
  #
  path_metadata = opt$metadata
  path_fastq = opt$fastq
  threads = opt$threads

  #
  # Print arguments
  #
  writeLines(paste0("metadata=", path_metadata))
  writeLines(paste0("fastq=", path_fastq))
  writeLines(paste0("threads=", threads))

  pipeline_checksum(path_metadata=path_metadata, path_fastq=path_fastq, threads=threads)
}
pipeline_checksum_cli()