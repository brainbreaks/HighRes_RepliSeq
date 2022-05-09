#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(readr))
source("utils.R")

pipeline_align = function(path_metadata, path_fastq, alignment_mapq=40, trim_quality=22, trim_minlength=40, path_database, path_output, threads=1, SAMPLE_NUMBER)
{
  path_output = R.utils::getAbsolutePath(path_output)
  if(!dir.exists(path_output)) cmd_run(paste0("mkdir ", path_output))
  if(!dir.exists(file.path(path_output, "alignments"))) cmd_run(paste0("mkdir ", file.path(path_output, "alignments")))
  if(!dir.exists(file.path(path_output, "trim"))) cmd_run(paste0("mkdir ", file.path(path_output, "trim")))

  files_df = data.frame(FASTQ_PATH=list.files(path_fastq, recursive=T, full.names=T)) %>% dplyr::mutate(FASTQ_FILE=basename(FASTQ_PATH))
  meta_df = readr::read_tsv(path_metadata) %>%
    dplyr::mutate(SAMPLE_NUMBER=match(SAMPLE_NAME, unique(SAMPLE_NAME)))
  #
  # Check external executables
  #
  cmd_is_available(c("bedtools", "samtools", "bowtie2", "picard", "trim_galore"))
  #
  # Check for errors in arguments
  #
  error_message = ""
  path_database = file.path(R.utils::getAbsolutePath(dirname(path_database)), basename(path_database))
  if (!file.exists(paste0(path_database, ".1.bt2"))){
      error_message = paste0("Bowtie index does not exist '", path_database)}

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
  arguments1_df = meta_df %>%
      dplyr::left_join(files_df, by="FASTQ_FILE") %>%
      dplyr::rowwise() %>%
      dplyr::mutate(TRIM_FASTQ_PATH=R.utils::getAbsolutePath(file.path(path_output, stringr::str_glue("trim/{file}", sample=SAMPLE_NAME, file=gsub("\\.fq.gz|\\.fastq.gz", paste0("_val_", FASTQ_PAIR, ".fq.gz"), basename(FASTQ_PATH)))))) %>%
      dplyr::ungroup()
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
  arguments_trim_df = arguments1_df %>%
    dplyr::mutate(TRIM_FASTQ_PAIR=paste0("trim_", FASTQ_PAIR)) %>% reshape2::dcast(SAMPLE_NAME~TRIM_FASTQ_PAIR, value.var="TRIM_FASTQ_PATH")
  arguments_df = arguments1_df %>%
    reshape2::dcast(SAMPLE_NAME+SAMPLE_NUMBER+SAMPLE_CONDITION+SAMPLE_FRACTION~FASTQ_PAIR, value.var="FASTQ_PATH") %>%
    dplyr::inner_join(arguments_trim_df, by="SAMPLE_NAME") %>%
    dplyr::arrange(SAMPLE_NUMBER)
  #
  # Run pipeline
  #
  writeLines("Aligning sequenced reads with reference...")


  for(i in 1:nrow(arguments_df)) {
    if(SAMPLE_NUMBER != -1 && i != SAMPLE_NUMBER){next}
    writeLines(paste0(i, "/", nrow(arguments_df), ": ", arguments_df$SAMPLE_NAME[i]))
    path_trim = R.utils::getAbsolutePath(file.path(path_output, stringr::str_glue("trim", sample=arguments_df$SAMPLE_NAME[i])))
    path_bam = R.utils::getAbsolutePath(file.path(path_output, stringr::str_glue("alignments/{sample}.bam", sample=arguments_df$SAMPLE_NAME[i])))
    path_bai = R.utils::getAbsolutePath(paste0(path_bam, ".bai"))
    path_tmpbam = R.utils::getAbsolutePath(file.path(path_output, stringr::str_glue("alignments/{sample}_tmp.bam", sample=arguments_df$SAMPLE_NAME[i])))
    path_tmpsam = R.utils::getAbsolutePath(file.path(path_output, stringr::str_glue("alignments/{sample}_tmp.sam", sample=arguments_df$SAMPLE_NAME[i])))
    path_markdup_report = R.utils::getAbsolutePath(file.path(path_output, stringr::str_glue("alignments/{sample}_markdup.txt", sample=arguments_df$SAMPLE_NAME[i])))
    if(!file.exists(path_bam)) {
      # Trim
      cmd_trim = stringr::str_glue("trim_galore --paired --length={trim_minlength} --gzip -q {trim_quality} -o '{path_trim}' -j {threads} '{input1}' '{input2}'", path_trim=path_trim, trim_minlength=trim_minlength, input1=arguments_df$`1`[i], input2=arguments_df$`2`[i], trim_quality=trim_quality, threads=threads)
      cmd_run(cmd_trim)

      # Bowtie2
      cmd_bowtie = stringr::str_glue("bowtie2 -q -p {threads} --phred33 --end-to-end --sensitive --no-mixed --no-discordant --no-unal -x '{db}' -1  '{r1}' -2  '{r2}' -S '{sam}'", bam=path_bam, sam=path_tmpsam, db=path_database, r1=arguments_df$`trim_1`[i], r2=arguments_df$`trim_2`[i], sample=arguments_df$SAMPLE_NAME[i], threads=threads)
      cmd_run(cmd_bowtie)

      # Convert to BAM
      cmd_mapq = stringr::str_glue("samtools view -bSq {alignment_mapq} -@ {threads} -o '{bam}' '{sam}'", alignment_mapq=alignment_mapq, sam=path_tmpsam, bam=path_bam, threads=threads)
      cmd_run(cmd_mapq)
      file.remove(path_tmpsam)

      # Sort
      cmd_sort = stringr::str_glue("samtools sort -@ {threads} -o '{tmp}' '{bam}'", bam=path_bam, tmp=path_tmpbam, threads=threads)
      cmd_run(cmd_sort)
      cmd_movesort = stringr::str_glue("mv '{tmpbam}' '{bam}'", tmpbam=path_tmpbam, bam=path_bam)
      cmd_run(cmd_movesort)

      # mark duplicates
      cmd_markduplicates = stringr::str_glue("picard-tools MarkDuplicates INPUT='{bam}' METRICS_FILE='{path_markdup_report}' OUTPUT='{tmpbam}' REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate", tmpbam=path_tmpbam, bam=path_bam, path_markdup_report=path_markdup_report)
      cmd_run(cmd_markduplicates)

      # Create final BAM alignment file
      cmd_movefinal = stringr::str_glue("mv '{tmpbam}' '{bam}'", tmpbam=path_tmpbam, bam=path_bam)
      cmd_run(cmd_movefinal)
    } else {
      writeLines(paste0("BAM file '", path_bam, "' exists, skipping..."))
    }
    if(!file.exists(path_bai)) {
      cmd_index = stringr::str_glue("samtools index -@ {threads} -b '{bam}'", bam=path_bam, threads=threads)
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
    optparse::make_option(c("-i", "--fastq-dir"), dest="fastq", type="character", default=".", help="Path to folder with fastq files [default= %default]", metavar="character"),
    optparse::make_option(c("-t", "--threads"), type="integer", default=1, help="Output directory [default= %default]", metavar="integer"),
    optparse::make_option('--trim-minlength', dest="trim_minlength", type="integer", default=40, help="Minimal length of the trimmed fasta sequence"),
    optparse::make_option('--alignment-mapq', dest="alignment_mapq", type="integer", default = 40, help= "MAPQ Score of Alignment"),
    optparse::make_option('--trim-quality', dest= "trim_quality", type="integer", default=22, help="Minimal quality of trimmed sequence"),
    optparse::make_option("--sample-number", dest = "SAMPLE_NUMBER", type = "integer", default=-1, help="Align a single sample from metadata table(number #SAMPLE_NUMBER)")
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
  trim_minlength = opt$trim_minlength
  alignment_mapq = opt$alignment_mapq
  trim_quality = opt$trim_quality
  SAMPLE_NUMBER = opt$SAMPLE_NUMBER

  #
  # Print arguments
  #
  writeLines(paste0("metadata=", path_metadata))
  writeLines(paste0("genome=", path_database))
  writeLines(paste0("out=", path_output))
  writeLines(paste0("fastq=", path_fastq))
  writeLines(paste0("threads=", threads))
  writeLines(paste0("trim min length=", trim_minlength))
  writeLines(paste0("MAPQ score of alignment=", alignment_mapq))
  writeLines(paste0("Trim quality=", trim_quality))
  #
  # Create missing folders
  #
  pipeline_align(path_metadata=path_metadata, SAMPLE_NUMBER = SAMPLE_NUMBER, trim_quality=trim_quality, alignment_mapq=alignment_mapq, trim_minlength=trim_minlength, path_fastq=path_fastq, path_database=path_database, path_output=path_output, threads=threads)
}
pipeline_align_cli()






