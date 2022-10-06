# helpful_commands

[![Check Links](https://github.com/paulstothard/helpful_commands/actions/workflows/links.yml/badge.svg)](https://github.com/paulstothard/helpful_commands/actions/workflows/links.yml) [![GitHub Super-Linter](https://github.com/paulstothard/helpful_commands/actions/workflows/linter.yml/badge.svg)](https://github.com/paulstothard/helpful_commands/actions/workflows/linter.yml) [![Generate TOC](https://github.com/paulstothard/helpful_commands/actions/workflows/toc.yml/badge.svg)](https://github.com/paulstothard/helpful_commands/actions/workflows/toc.yml)

Command-line tools, commands, and code snippets for performing routine data processing and bioinformatics tasks.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [helpful_commands](#helpful_commands)
  - [awk](#awk)
    - [Add a header line to a file](#add-a-header-line-to-a-file)
    - [Convert a CSV file to a FASTA file](#convert-a-csv-file-to-a-fasta-file)
    - [Print lines in file when a certain column contains a specific value](#print-lines-in-file-when-a-certain-column-contains-a-specific-value)
    - [Replace certain values in specific columns](#replace-certain-values-in-specific-columns)
    - [Add up the values in a column](#add-up-the-values-in-a-column)
    - [Create a new column from two existing columns](#create-a-new-column-from-two-existing-columns)
    - [For each category in one column, add up the values in another column](#for-each-category-in-one-column-add-up-the-values-in-another-column)
    - [Print column names and numbers](#print-column-names-and-numbers)
    - [Print the values observed in a specific column, along with the number of times each value is observed](#print-the-values-observed-in-a-specific-column-along-with-the-number-of-times-each-value-is-observed)
    - [Print the number of lines exhibiting each distinct number of fields](#print-the-number-of-lines-exhibiting-each-distinct-number-of-fields)
    - [Print lines where certain fields contain values of interest](#print-lines-where-certain-fields-contain-values-of-interest)
    - [Write each row to a separate file named after the value in a specific column](#write-each-row-to-a-separate-file-named-after-the-value-in-a-specific-column)
    - [Split a multi-FASTA file into separate files named according to the sequence title](#split-a-multi-fasta-file-into-separate-files-named-according-to-the-sequence-title)
    - [Print only specific columns, identified by name in the first row](#print-only-specific-columns-identified-by-name-in-the-first-row)
    - [Print only the lines coming after a certain starting line and before a certain ending line](#print-only-the-lines-coming-after-a-certain-starting-line-and-before-a-certain-ending-line)
    - [Print only the first non-commented line](#print-only-the-first-non-commented-line)
    - [Print the average read length for a FASTQ file](#print-the-average-read-length-for-a-fastq-file)
    - [Sort lines based on order of IDs in another file](#sort-lines-based-on-order-of-ids-in-another-file)
    - [Skip footer lines](#skip-footer-lines)
    - [Skip header lines](#skip-header-lines)
  - [Bash](#bash)
    - [View STDOUT and append it to a file](#view-stdout-and-append-it-to-a-file)
    - [Redirect STDERR to STDOUT and view both and append both to a file](#redirect-stderr-to-stdout-and-view-both-and-append-both-to-a-file)
    - [Change Bash prompt temporarily](#change-bash-prompt-temporarily)
    - [Check Bash scripts for potential issues](#check-bash-scripts-for-potential-issues)
    - [Extract a file](#extract-a-file)
    - [Save the output of a command in a variable](#save-the-output-of-a-command-in-a-variable)
    - [Set variables using values from a CSV file](#set-variables-using-values-from-a-csv-file)
    - [Perform a calculation on the command line](#perform-a-calculation-on-the-command-line)
    - [Combine FASTA files into a single file replacing each FASTA record name with the name of the input file](#combine-fasta-files-into-a-single-file-replacing-each-fasta-record-name-with-the-name-of-the-input-file)
  - [brew](#brew)
    - [List installed packages](#list-installed-packages)
    - [View available packages](#view-available-packages)
    - [Install a package](#install-a-package)
    - [Add a third-party repository](#add-a-third-party-repository)
    - [Install directly from a third-party repository](#install-directly-from-a-third-party-repository)
    - [View packages available from brewsci/bio](#view-packages-available-from-brewscibio)
    - [List installed graphical applications](#list-installed-graphical-applications)
    - [View available graphical applications](#view-available-graphical-applications)
    - [Install a graphical application](#install-a-graphical-application)
  - [Conda](#conda)
    - [Install Miniconda](#install-miniconda)
    - [Create an environment and install some packages](#create-an-environment-and-install-some-packages)
    - [Deactivate an environment](#deactivate-an-environment)
    - [Activate an environment](#activate-an-environment)
    - [List available packages](#list-available-packages)
    - [Search for a specific package](#search-for-a-specific-package)
    - [Add additional packages to an environment](#add-additional-packages-to-an-environment)
    - [List environments](#list-environments)
    - [List packages installed in the active environment](#list-packages-installed-in-the-active-environment)
    - [Remove an environment](#remove-an-environment)
  - [csvkit](#csvkit)
    - [Convert Excel to CSV](#convert-excel-to-csv)
    - [Convert JSON to CSV](#convert-json-to-csv)
    - [Print column names](#print-column-names)
    - [Select a subset of columns](#select-a-subset-of-columns)
    - [Reorder columns](#reorder-columns)
    - [Find rows with matching cells](#find-rows-with-matching-cells)
    - [Convert to JSON](#convert-to-json)
    - [Generate summary statistics](#generate-summary-statistics)
    - [Query with SQL](#query-with-sql)
    - [Merge CSV files on a specified column or columns](#merge-csv-files-on-a-specified-column-or-columns)
  - [cut](#cut)
    - [Extract columns of interest](#extract-columns-of-interest)
    - [Extract a range of columns](#extract-a-range-of-columns)
    - [Extract everything except one column](#extract-everything-except-one-column)
    - [Extract characters](#extract-characters)
    - [Change field separators](#change-field-separators)
  - [datamash](#datamash)
    - [Group records by one column and print information about each group](#group-records-by-one-column-and-print-information-about-each-group)
    - [Print statistics for a column](#print-statistics-for-a-column)
    - [Transpose](#transpose)
  - [Docker](#docker)
    - [Perform a sequence comparison using legacy BLAST](#perform-a-sequence-comparison-using-legacy-blast)
    - [Annotate sequence variants using VEP](#annotate-sequence-variants-using-vep)
    - [Annotate a bacterial genome using Prokka](#annotate-a-bacterial-genome-using-prokka)
    - [Compare sequence reads to a bacterial genome to find SNPs using Snippy](#compare-sequence-reads-to-a-bacterial-genome-to-find-snps-using-snippy)
    - [Get an interactive bash shell in a Docker container](#get-an-interactive-bash-shell-in-a-docker-container)
    - [List images](#list-images)
    - [List running containers](#list-running-containers)
    - [Stop a container](#stop-a-container)
    - [Kill all running containers](#kill-all-running-containers)
    - [Delete all containers that are not running](#delete-all-containers-that-are-not-running)
  - [File conversion](#file-conversion)
    - [Convert CSV to TSV](#convert-csv-to-tsv)
    - [Convert TSV to CSV](#convert-tsv-to-csv)
    - [Convert CSV to Markdown](#convert-csv-to-markdown)
    - [Convert CSV to HTML](#convert-csv-to-html)
    - [Convert PDF to PNG](#convert-pdf-to-png)
    - [Convert PNG to PDF](#convert-png-to-pdf)
    - [Convert DOCX to PDF](#convert-docx-to-pdf)
    - [Convert Excel to CSV with csvkit](#convert-excel-to-csv-with-csvkit)
    - [Convert CSV to Excel](#convert-csv-to-excel)
    - [Convert HTML to PDF](#convert-html-to-pdf)
    - [Convert a website to PDF](#convert-a-website-to-pdf)
    - [Convert HTML to PNG](#convert-html-to-png)
    - [Convert Markdown to PDF](#convert-markdown-to-pdf)
    - [Convert Markdown to HTML](#convert-markdown-to-html)
    - [Convert Markdown to PPTX](#convert-markdown-to-pptx)
    - [Convert PPTX to PDF](#convert-pptx-to-pdf)
    - [Convert PPTX to PNG](#convert-pptx-to-png)
  - [File downloads](#file-downloads)
    - [Download an entire website](#download-an-entire-website)
    - [Download from an FTP server](#download-from-an-ftp-server)
    - [Download from Google Drive](#download-from-google-drive)
    - [Download a GenBank file with curl](#download-a-genbank-file-with-curl)
    - [Download fastq files based on a list of SRA accessions](#download-fastq-files-based-on-a-list-of-sra-accessions)
    - [Download a reference genome FASTA file from Ensembl](#download-a-reference-genome-fasta-file-from-ensembl)
    - [Download a reference genome GTF file from Ensembl](#download-a-reference-genome-gtf-file-from-ensembl)
  - [find](#find)
    - [Perform a series of commands on files returned by find](#perform-a-series-of-commands-on-files-returned-by-find)
    - [Sort files before processing](#sort-files-before-processing)
    - [Use files as the argument list for a command](#use-files-as-the-argument-list-for-a-command)
    - [Process files in pairs](#process-files-in-pairs)
    - [Copy the files returned by find](#copy-the-files-returned-by-find)
    - [Copy the files returned by find, naming the copies after a directory in the path](#copy-the-files-returned-by-find-naming-the-copies-after-a-directory-in-the-path)
    - [Switch to the directory containing each file and execute a command](#switch-to-the-directory-containing-each-file-and-execute-a-command)
    - [Find large files](#find-large-files)
    - [Process multiple files and redirect output to one file](#process-multiple-files-and-redirect-output-to-one-file)
    - [Process multiple files and redirect output to separate files](#process-multiple-files-and-redirect-output-to-separate-files)
  - [Git](#git)
    - [Create a new Git repository](#create-a-new-git-repository)
    - [Sync a repository to your local machine](#sync-a-repository-to-your-local-machine)
    - [Clone all your repositories](#clone-all-your-repositories)
    - [Mark changed files to be included in the next commit](#mark-changed-files-to-be-included-in-the-next-commit)
    - [Undo a Git add before a commit](#undo-a-git-add-before-a-commit)
    - [Remove files from the repository](#remove-files-from-the-repository)
    - [Move or rename a file or directory](#move-or-rename-a-file-or-directory)
    - [Save the marked files to the local Git repository](#save-the-marked-files-to-the-local-git-repository)
    - [Push a commit on your local branch to a remote repository](#push-a-commit-on-your-local-branch-to-a-remote-repository)
    - [Pull a change from a remote repository to your local branch](#pull-a-change-from-a-remote-repository-to-your-local-branch)
    - [Add or edit a remote repository](#add-or-edit-a-remote-repository)
    - [Create and merge Git branches](#create-and-merge-git-branches)
    - [Specify files to ignore](#specify-files-to-ignore)
    - [Check the status of a working directory](#check-the-status-of-a-working-directory)
    - [Tag a release](#tag-a-release)
  - [grep](#grep)
    - [Count matches](#count-matches)
    - [Get the line number of a match](#get-the-line-number-of-a-match)
    - [Remove files that contain a match](#remove-files-that-contain-a-match)
    - [Remove files that do not contain a match](#remove-files-that-do-not-contain-a-match)
    - [Print non-matching lines](#print-non-matching-lines)
    - [Search using patterns from a file](#search-using-patterns-from-a-file)
  - [Image files](#image-files)
    - [Crop an image and add a white border](#crop-an-image-and-add-a-white-border)
    - [Resize an image](#resize-an-image)
    - [Take a screenshot of a window on macOS](#take-a-screenshot-of-a-window-on-macos)
    - [Take a webpage screenshot using Firefox](#take-a-webpage-screenshot-using-firefox)
    - [Record your terminal to an animated GIF](#record-your-terminal-to-an-animated-gif)
    - [Record your screen to an animated GIF](#record-your-screen-to-an-animated-gif)
    - [Create an animated GIF from a YouTube video](#create-an-animated-gif-from-a-youtube-video)
  - [join](#join)
    - [Combine rows based on shared keys with join](#combine-rows-based-on-shared-keys-with-join)
  - [Mamba](#mamba)
    - [Install Mamba](#install-mamba)
    - [Create an environment and install some packages with mamba](#create-an-environment-and-install-some-packages-with-mamba)
    - [Deactivate an environment with mamba](#deactivate-an-environment-with-mamba)
    - [Activate an environment with mamba](#activate-an-environment-with-mamba)
    - [List available packages with mamba](#list-available-packages-with-mamba)
    - [Search for a specific package with mamba](#search-for-a-specific-package-with-mamba)
    - [Add additional packages to an environment with mamba](#add-additional-packages-to-an-environment-with-mamba)
    - [List environments with mamba](#list-environments-with-mamba)
    - [List packages installed in the active environment with mamba](#list-packages-installed-in-the-active-environment-with-mamba)
    - [Remove an environment with mamba](#remove-an-environment-with-mamba)
  - [md5sum](#md5sum)
    - [Generate a file of checksums](#generate-a-file-of-checksums)
    - [Validate checksums](#validate-checksums)
  - [Miller](#miller)
    - [Extract the first 10 records of a CSV file](#extract-the-first-10-records-of-a-csv-file)
    - [Extract the last 10 records of a CSV file](#extract-the-last-10-records-of-a-csv-file)
    - [Convert formats](#convert-formats)
    - [View stats](#view-stats)
    - [Filter records](#filter-records)
    - [Sort records](#sort-records)
    - [Extract columns](#extract-columns)
    - [Edit columns](#edit-columns)
    - [Other actions](#other-actions)
    - [Combine actions](#combine-actions)
  - [Other](#other)
    - [Obtain your public IP address and network information](#obtain-your-public-ip-address-and-network-information)
    - [Copy an ssh public key to another system](#copy-an-ssh-public-key-to-another-system)
    - [Add a header to all files with a certain extension, getting the header from another file](#add-a-header-to-all-files-with-a-certain-extension-getting-the-header-from-another-file)
    - [Change the extension of multiple files](#change-the-extension-of-multiple-files)
    - [Add text to the beginning of a file](#add-text-to-the-beginning-of-a-file)
    - [Add text or a header to the beginning of all files with a particular file extension](#add-text-or-a-header-to-the-beginning-of-all-files-with-a-particular-file-extension)
    - [Find common lines between files](#find-common-lines-between-files)
    - [Format a CSV file into columns and examine its content](#format-a-csv-file-into-columns-and-examine-its-content)
    - [Format code](#format-code)
    - [Run commands at scheduled times using cron](#run-commands-at-scheduled-times-using-cron)
    - [Create a collection of MP3 files from a YouTube playlist](#create-a-collection-of-mp3-files-from-a-youtube-playlist)
    - [Perform a calculation using bc](#perform-a-calculation-using-bc)
    - [Count the bases in a fastq file](#count-the-bases-in-a-fastq-file)
    - [Count the reads in a fastq file](#count-the-reads-in-a-fastq-file)
    - [Use SQL-like queries to work with a CSV or TSV file](#use-sql-like-queries-to-work-with-a-csv-or-tsv-file)
    - [Browse, search, and edit a large CSV file](#browse-search-and-edit-a-large-csv-file)
    - [Prevent a command from stopping when you log out or exit the shell](#prevent-a-command-from-stopping-when-you-log-out-or-exit-the-shell)
    - [Edit a PDF file](#edit-a-pdf-file)
  - [parallel](#parallel)
    - [Extract files in parallel](#extract-files-in-parallel)
    - [Compress files in parallel](#compress-files-in-parallel)
    - [Process files in pairs using parallel](#process-files-in-pairs-using-parallel)
    - [Perform BLAST in parallel](#perform-blast-in-parallel)
  - [paste](#paste)
    - [Combine columns with paste](#combine-columns-with-paste)
  - [Perl](#perl)
    - [Get a random sample of lines from a text file while excluding the header line](#get-a-random-sample-of-lines-from-a-text-file-while-excluding-the-header-line)
    - [Convert a FASTA file to a CSV file with column names](#convert-a-fasta-file-to-a-csv-file-with-column-names)
    - [Count the number of lines that match a regular expression](#count-the-number-of-lines-that-match-a-regular-expression)
    - [Extract FASTA sequences from a file based on a file of sequence names of interest](#extract-fasta-sequences-from-a-file-based-on-a-file-of-sequence-names-of-interest)
    - [Add a FASTA title to the start of a sequence in RAW format](#add-a-fasta-title-to-the-start-of-a-sequence-in-raw-format)
    - [Remove commas located within quoted fields in a CSV file and create a tab-delimited file](#remove-commas-located-within-quoted-fields-in-a-csv-file-and-create-a-tab-delimited-file)
    - [Remove lines that match a pattern](#remove-lines-that-match-a-pattern)
    - [Replace commas with tabs](#replace-commas-with-tabs)
    - [Replace tabs with commas and remove quotes](#replace-tabs-with-commas-and-remove-quotes)
    - [Sort sections in a Markdown file based on headings](#sort-sections-in-a-markdown-file-based-on-headings)
    - [Search and replace text on each line](#search-and-replace-text-on-each-line)
    - [Replace ^M](#replace-m)
    - [Print matches that may span multiple lines](#print-matches-that-may-span-multiple-lines)
    - [Print matches after additional editing](#print-matches-after-additional-editing)
    - [Split a Markdown file into separate files based on a heading](#split-a-markdown-file-into-separate-files-based-on-a-heading)
    - [Format Perl code](#format-perl-code)
  - [R](#r)
    - [Compare two data sets to find differences](#compare-two-data-sets-to-find-differences)
    - [Visualize the degree of overlap among gene sets](#visualize-the-degree-of-overlap-among-gene-sets)
    - [Cluster gene lists based on overlap and identify shared genes](#cluster-gene-lists-based-on-overlap-and-identify-shared-genes)
    - [Transpose a data frame](#transpose-a-data-frame)
    - [Split two-allele genotypes into two columns](#split-two-allele-genotypes-into-two-columns)
    - [Split multi-locus genotype strings into multiple columns and decode genotypes](#split-multi-locus-genotype-strings-into-multiple-columns-and-decode-genotypes)
    - [Change values in a column based on values in another column](#change-values-in-a-column-based-on-values-in-another-column)
    - [Add comment lines to output](#add-comment-lines-to-output)
    - [Filter and sort rows](#filter-and-sort-rows)
    - [Add columns from one tibble to another](#add-columns-from-one-tibble-to-another)
    - [Combine multiple input files](#combine-multiple-input-files)
  - [rsync](#rsync)
    - [Sync a directory on local system](#sync-a-directory-on-local-system)
    - [Sync a directory to a remote system](#sync-a-directory-to-a-remote-system)
    - [Sync a directory from a remote system](#sync-a-directory-from-a-remote-system)
  - [Singularity](#singularity)
    - [Create an image from a container stored in Docker Hub](#create-an-image-from-a-container-stored-in-docker-hub)
  - [sbatch](#sbatch)
    - [Count lines in compressed fastq files](#count-lines-in-compressed-fastq-files)
  - [sed](#sed)
    - [Add a header line to a file with sed](#add-a-header-line-to-a-file-with-sed)
    - [Edit the header line with sed](#edit-the-header-line-with-sed)
    - [Print a specific line of a file](#print-a-specific-line-of-a-file)
    - [Change filenames using a regular expression](#change-filenames-using-a-regular-expression)
    - [Search and replace on lines](#search-and-replace-on-lines)
    - [Delete lines](#delete-lines)
  - [Share data with project group members](#share-data-with-project-group-members)
  - [Slurm](#slurm)
    - [View statistics related to the efficiency of resource usage of a completed job](#view-statistics-related-to-the-efficiency-of-resource-usage-of-a-completed-job)
    - [View jobs](#view-jobs)
    - [View running jobs](#view-running-jobs)
    - [View pending jobs](#view-pending-jobs)
    - [View detailed information for a specific job](#view-detailed-information-for-a-specific-job)
    - [View accounting information for completed jobs](#view-accounting-information-for-completed-jobs)
    - [Cancel a job](#cancel-a-job)
    - [Cancel all jobs](#cancel-all-jobs)
    - [Start an interactive session](#start-an-interactive-session)
  - [sort](#sort)
    - [Alphabetical sort](#alphabetical-sort)
    - [Specify the sort field](#specify-the-sort-field)
    - [Use multiple sort fields](#use-multiple-sort-fields)
    - [Sort a file with a header row](#sort-a-file-with-a-header-row)
    - [Sort by chromosome](#sort-by-chromosome)
  - [tmux](#tmux)
    - [Start a tmux session](#start-a-tmux-session)
    - [Detach a tmux session](#detach-a-tmux-session)
    - [List tmux sessions](#list-tmux-sessions)
    - [Join a tmux session](#join-a-tmux-session)
    - [Create a multi-pane tmux session](#create-a-multi-pane-tmux-session)
    - [Scroll in tmux](#scroll-in-tmux)
    - [Navigate between tmux panes](#navigate-between-tmux-panes)
    - [Kill a tmux session](#kill-a-tmux-session)
  - [tr](#tr)
    - [Translate characters](#translate-characters)
    - [Delete characters](#delete-characters)
    - [Squeeze characters](#squeeze-characters)
  - [VCF files](#vcf-files)
    - [Extract variants from a region of interest](#extract-variants-from-a-region-of-interest)
    - [Extract variants from multiple regions of interest](#extract-variants-from-multiple-regions-of-interest)
    - [Extract variants from multiple regions of interest described in a file](#extract-variants-from-multiple-regions-of-interest-described-in-a-file)
    - [Assess sex by calculating X-chromosome heterozygosity](#assess-sex-by-calculating-x-chromosome-heterozygosity)
    - [Assess sex using plink](#assess-sex-using-plink)
    - [Count genotypes](#count-genotypes)
    - [Determine the proportion of missing genotypes in each sample](#determine-the-proportion-of-missing-genotypes-in-each-sample)
    - [Determine genotype concordance](#determine-genotype-concordance)
    - [Count Mendelian errors using plink](#count-mendelian-errors-using-plink)
    - [Find runs of homozygosity using plink](#find-runs-of-homozygosity-using-plink)
    - [Add new samples](#add-new-samples)
    - [Check relatedness between samples](#check-relatedness-between-samples)
    - [Filter sites based on genotypes and other criteria](#filter-sites-based-on-genotypes-and-other-criteria)
    - [Add predicted consequences](#add-predicted-consequences)
    - [Add variant IDs](#add-variant-ids)
    - [Transfer annotations from another VCF](#transfer-annotations-from-another-vcf)
    - [Filter variants based on predicted consequences](#filter-variants-based-on-predicted-consequences)
    - [Keep only variants with a missing ID](#keep-only-variants-with-a-missing-id)
    - [Keep only variants with an assigned ID](#keep-only-variants-with-an-assigned-id)
    - [Keep only variants where FILTER is PASS](#keep-only-variants-where-filter-is-pass)
    - [Keep only biallelic SNPs](#keep-only-biallelic-snps)
    - [Keep only SNPs](#keep-only-snps)
    - [Keep only indels](#keep-only-indels)
    - [Keep only sites not found in a second VCF](#keep-only-sites-not-found-in-a-second-vcf)
    - [Count variants](#count-variants)
    - [Identify variants found in all VCF files](#identify-variants-found-in-all-vcf-files)
    - [Change sample order](#change-sample-order)
    - [Combine rows](#combine-rows)
    - [Convert a VCF file to an Excel file](#convert-a-vcf-file-to-an-excel-file)
    - [Perform case-control analysis](#perform-case-control-analysis)
    - [Remove annotations](#remove-annotations)
    - [Remove a sample](#remove-a-sample)
    - [Remove sites with any missing genotypes](#remove-sites-with-any-missing-genotypes)
    - [Remove all genotypes](#remove-all-genotypes)
    - [Print samples](#print-samples)
  - [vim](#vim)
    - [Search and replace across multiple files](#search-and-replace-across-multiple-files)
    - [Search and replace newlines](#search-and-replace-newlines)
    - [Compare two files](#compare-two-files)
    - [Copy to the clipboard](#copy-to-the-clipboard)
    - [Remove trailing whitespace](#remove-trailing-whitespace)
    - [View ^M characters](#view-m-characters)
    - [Type tab characters](#type-tab-characters)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## awk

### Add a header line to a file

```bash
awk 'BEGIN{print "my header text"}1' input.txt
```

### Convert a CSV file to a FASTA file

In this example column `1` contains the sequence title and column `3` contains the sequence:

```bash
awk -F, '{print ">"$1"\n"$3"\n"}' input.csv
```

### Print lines in file when a certain column contains a specific value

In this example lines are printed when the value in column `1` equals `9913`:

```bash
awk -F, '{if ($1 == 9913) print $0}' input.csv
```

To also print all lines starting with `#` (header lines for example), use:

```bash
awk -F, '$1 ~ /^#/ {print; next} {if ($1 == 9913) print $0}' input.csv
```

Print the first line and lines when the value in column `1` equals `9913`:

```bash
awk -F, 'NR==1; NR > 1 {if ($1 == 9913) print $0}' input.csv
```

### Replace certain values in specific columns

In this example `1` and `-1` in column `23` are replaced with `forward` and `reverse`, respectively:

```bash
awk -F\\t 'BEGIN {OFS = "\t"} {sub(/^1/, "forward", $23); sub(/^-1/, "reverse", $23); print}' input.tab
```

### Add up the values in a column

In this example values in column `5` are summed in a file with columns separated by one or more blank spaces:

```bash
awk -F' {1,}' '{sum+=$5} END {print sum}' input.txt
```

### Create a new column from two existing columns

In this example the values in columns `3` and `4` are added to create a new column:

```bash
awk -F, '{print $0,$3+$4}' input.txt
```

### For each category in one column, add up the values in another column

In this example values in column `2` are summed up for each category in column `1`:

```bash
awk -F, '{a[$1]+=$2}END{for(i in a) print i,a[i]}' input.csv
```

### Print column names and numbers

In this example the first row of the input file contains the column names:

```bash
awk -F $'\t' 'NR>1{exit};{for (i = 1; i <= NF; i++) print "column " i,"is " $i}' input.tab
```

### Print the values observed in a specific column, along with the number of times each value is observed

In this example the counts for each distinct value in column `9` are printed:

```bash
awk -F $'\t' '{count[$9]++}END{for(j in count) print j"\t"count[j]}' input.tab
```

In this example the counts for each distinct pair of values in column `1` and column `2` are printed:

```bash
awk -F $'\t' '{count[$1"\t"$2]++}END{for(j in count) print j"\t"count[j]}' input.tab
```

### Print the number of lines exhibiting each distinct number of fields

```bash
awk -F $'\t' '{count[NF]++}END{for(j in count) print "line length " j,"("count[j]" counts)"}' input.tab
```

### Print lines where certain fields contain values of interest

In this example lines where column `2` equals `7` and column `3` is between `60240145` and `60255062` are printed:

```bash
awk -F, '{ if ($2 == 7 && $3 >= 60240145 && $3 <= 60255062) print $0 }' input.csv
```

In this example the header line is printed, followed by lines where column `2` starts with `MT-` and column `3` starts with `MT-`:

```bash
awk -F $'\t' 'NR==1{print; next} $1~/^MT-/ && $2~/^MT-/' input.tab
```

### Write each row to a separate file named after the value in a specific column

In this example each file is named after the value in column `1`:

```bash
awk -F $'\t' '{ fname = $1 ".txt"; print >>fname; close(fname) }' input.tab
```

### Split a multi-FASTA file into separate files named according to the sequence title

In this example the sequences are written to a directory called `out`:

```bash
outputdir=out/
mkdir -p "$outputdir"
awk '/^>/ {OUT=substr($0,2); split(OUT, a, " "); sub(/[^A-Za-z_0-9\.\-]/, "", a[1]); OUT = "'"$outputdir"'" a[1] ".fa"}; OUT {print >>OUT; close(OUT)}' input.fasta
```

### Print only specific columns, identified by name in the first row

In this example the columns named `Affy SNP ID` and `Flank` are printed:

```bash
awk -F, 'NR==1 { for (i=1; i<=NF; i++) { ix[$i] = i } } NR>1 { print $ix["Affy SNP ID"]","$ix["Flank"] }' input.csv
```

### Print only the lines coming after a certain starting line and before a certain ending line

In this example the lines coming after a line starting with `IlmnID` and before a line starting with `[Controls]` are printed:

```bash
awk -F, '/^IlmnID/{flag=1;print;next}/^\[Controls\]/{flag=0}flag' input.csv
```

### Print only the first non-commented line

```bash
awk '/^[^#]/ { print $0;exit; }' input.txt
```

### Print the average read length for a FASTQ file

```bash
awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' input.fastq
```

### Sort lines based on order of IDs in another file

In this example the records in `file_to_sort.csv` have an identifier in column `1` that is present in `sorted_ids.txt`, which has a single column:

```bash
awk -F, 'NR==FNR {a[$1]=$0; next} ($0 in a) {print a[$0]}' file_to_sort.csv sorted_ids.txt
```

If both files have a single header line the following can be used to generate the sorted output with the restored header line:

```bash
(head -n 1 file_to_sort.csv && awk -F, 'NR==FNR {a[$1]=$0; next} ($0 in a) {print a[$0]}' <(tail -n +2 file_to_sort.csv) <(tail -n +2 sorted_ids.txt)) > sorted.csv
```

### Skip footer lines

In this example the contents of the file are printed except for the last `15` lines:

```bash
awk 'NR>c{print A[NR%c]} {A[NR%c]=$0}' c=15 input.txt
```

### Skip header lines

The following skips lines starting with `#`:

```bash
awk '/^[^#]/ { print $0 }' input.txt
```

The following skips the first `5` lines:

```bash
awk 'FNR > 5 { print $0 }' input.txt
```

## Bash

### View STDOUT and append it to a file

```bash
some_command | tee -a output.txt
```

### Redirect STDERR to STDOUT and view both and append both to a file

```bash
some_command 2>&1 | tee -a log
```

### Change Bash prompt temporarily

```bash
PS1="$ "
```

### Check Bash scripts for potential issues

Use [shellcheck](https://github.com/koalaman/shellcheck):

```bash
shellcheck some_script.sh
```

### Extract a file

The following Bash function can be used to extract a variety of file types:

```bash
extract() {
  if [ -f "$1" ]; then
    case "$1" in
    *.tar.bz2) tar xjf "$1" ;;
    *.tar.gz) tar xzf "$1" ;;
    *.bz2) bunzip2 "$1" ;;
    *.rar) unrar e "$1" ;;
    *.gz) gunzip "$1" ;;
    *.tar) tar xf "$1" ;;
    *.tbz2) tar xjf "$1" ;;
    *.tgz) tar xzf "$1" ;;
    *.zip) unzip "$1" ;;
    *.Z) uncompress "$1" ;;
    *.7z) 7z x "$1" ;;
    *) echo "'$1' cannot be extracted via extract()" ;;
    esac
  else
    echo "'$1' is not a valid file"
  fi
}
```

### Save the output of a command in a variable

Use a subshell:

```bash
result=$(echo "sqrt(16)" | bc -l)
```

### Set variables using values from a CSV file

The first `read` in the example below is used to skip a header line.

```bash
{
  read
  while IFS="," read -r column1 column2 column3
  do
    echo "$column1, $column2, $column3"
  done
} < "input.csv"
```

### Perform a calculation on the command line

Use Bash arithmetic expansion:

```bash
n=6
echo "$(( n - 1 * 2 ))"
answer="$(( n - 1 * 2 ))"
```

### Combine FASTA files into a single file replacing each FASTA record name with the name of the input file

```bash
for file in *.fasta
do
   echo ">$file" >> out.fasta
   tail -n +2 "$file" >> out.fasta
done
```

Or:

```bash
for file in *.fasta
do
   {
     echo ">$file"
     tail -n +2 "$file"
   } >> out.fasta
done
```

## brew

### List installed packages

```bash
brew list
```

### View available packages

To view packages available from the core tap via the Homebrew package manager for macOS:

- [https://formulae.brew.sh/formula/](https://formulae.brew.sh/formula/)

### Install a package

In this example `parallel`:

```bash
brew install parallel
```

### Add a third-party repository

In this example `brewsci/bio` for bioinformatics software:

```bash
brew tap brewsci/bio
```

### Install directly from a third-party repository

In this example `clustal-w` from `brewsci/bio`:

```bash
brew install brewsci/bio/clustal-w
```

### View packages available from brewsci/bio

- [https://github.com/brewsci/homebrew-bio/tree/develop/Formula](https://github.com/brewsci/homebrew-bio/tree/develop/Formula)

### List installed graphical applications

```bash
brew list --cask
```

### View available graphical applications

To view graphical applications available from the cask tap via the Homebrew package manager for macOS:

- [https://formulae.brew.sh/cask/](https://formulae.brew.sh/cask/)

### Install a graphical application

In this example the Firefox browser:

```bash
brew install firefox --cask
```

## Conda

### Install Miniconda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
source ~/miniconda3/bin/activate
conda init
source ~/.bashrc
conda update -y -n base -c defaults conda
```

### Create an environment and install some packages

In this example an environment called `ngs` is created:

```bash
conda create -y --name ngs
conda activate ngs
conda install -y -c bioconda -c conda-forge multiqc fastqc trimmomatic bowtie2 subread samtools
```

### Deactivate an environment

```bash
conda deactivate
```

### Activate an environment

```bash
conda activate ngs
```

### List available packages

```bash
conda search -c bioconda -c conda-forge
```

### Search for a specific package

```bash
conda search -c bioconda -c conda-forge blast
```

### Add additional packages to an environment

```bash
conda activate ngs
conda install -y -c bioconda -c conda-forge picard
```

### List environments

```bash
conda info --envs
```

### List packages installed in the active environment

```bash
conda list
```

### Remove an environment

In this example the environment to remove is called `my-env`:

```bash
conda deactivate
conda env remove --name my-env
```

## csvkit

The [csvkit](https://github.com/wireservice/csvkit) examples below are taken from the [csvkit documentation](https://csvkit.readthedocs.io/en/latest/).

### Convert Excel to CSV

```bash
in2csv data.xls > data.csv
```

### Convert JSON to CSV

```bash
in2csv data.json > data.csv
```

### Print column names

```bash
csvcut -n data.csv
```

### Select a subset of columns

```bash
csvcut -c column_a,column_c data.csv > new.csv
```

### Reorder columns

```bash
csvcut -c column_c,column_a data.csv > new.csv
```

### Find rows with matching cells

```bash
csvgrep -c phone_number -r "555-555-\d{4}" data.csv > new.csv
```

### Convert to JSON

```bash
csvjson data.csv > data.json
```

### Generate summary statistics

```bash
csvstat data.csv
```

### Query with SQL

```bash
csvsql --query "select name from data where age > 30" data.csv > new.csv
```

### Merge CSV files on a specified column or columns

```bash
csvjoin -c ID file1.csv file2.csv > inner_join_on_ID.csv
```

By default `csvjoin` performs an inner join. Other joins can be performed using `--outer`, `--left`, and `--right`.

## cut

### Extract columns of interest

In this example columns `1` and `2` are extracted from a from a CSV file:

```bash
cut -d, -f 1,2 sequenced_samples.csv
```

### Extract a range of columns

In this example columns `3`, `10`, `11`, and `12` are extracted from a tab-delimited file:

```bash
cut -d $'\t' -f 3,10-12 bovine_genotypes.vcf
```

### Extract everything except one column

In this example all columns except column `1` are returned (the `-f 2-` is used to mean "from field 2 to the end"):

```bash
cut -d $'\t' -f 2- bovine_genotypes.vcf
```

### Extract characters

Here `cut` is used to extract the first three characters from each line:

```bash
cut -c 1-3 sequenced_samples.csv
```

### Change field separators

`cut` can also be used to change the field separator. In this example tab-delimited values are changed to comma-separated values:

```bash
cut -d $'\t' -f 1- --output-delimiter=',' bovine_genotypes.vcf
```

To switch from comma-separated to tab-separated use:

```bash
cut -d, -f 1- --output-delimiter=$'\t' sequenced_samples.csv
```

On platforms lacking `--output-delimiter` `perl` can be used to switch tabs to commas:

```bash
perl -p -e 's/\t/,/g' sequenced_samples.csv
```

To switch commas to tabs when `--output-delimiter` isn't available the command below can be used. This  script handles cases when commas are inside of quoted fields:

```bash
perl -nle  'my @new  = (); push( @new, $+ )
while $_ =~ m{"([^\"\\]*(?:\\.[^\"\\]*)*)",?
| ([^,]+),? | ,}gx; push( @new, undef )
if substr( $text, -1, 1 ) eq '\'','\'';
for(@new){s/,/ /g} print join "\t", @new' sequenced_samples.csv
```

## datamash

[GNU datamash](https://www.gnu.org/software/datamash/) is a command-line program which performs basic numeric, textual and statistical operations on input textual data files.

### Group records by one column and print information about each group

In the following example the input CSV file has a header line. Records are grouped based on the value in column `2`, and for each group the mean value of column `5` is printed:

```bash
datamash -H -t, -g 2 mean 5 < example.csv
```

In the following example all the values in column `5` are printed for each group:

```bash
datamash -H -t, -g 2 collapse 5 < example.csv
```

### Print statistics for a column

In the following example a variety of statistics are generated for column `5`:

```bash
datamash -H -t, min 5 q1 5 median 5 q3 5 max 5 count 5 mean 5 sstdev 5 < example.csv
```

### Transpose

The following uses transpose to swap rows and columns in a CSV file with a header row:

```bash
datamash -H -t, transpose < example.csv
```

## Docker

### Perform a sequence comparison using legacy BLAST

Download a legacy BLAST Docker image:

```bash
docker pull quay.io/biocontainers/blast-legacy:2.2.26--2
```

Create a container from the image and run `formatdb` to create a formatted database. In this example the database is created from a DNA sequence file called `database.fasta`, located in the current directory:

```bash
docker run -it --rm -v "$(pwd)":/directory -u "$(id -u)":"$(id -g)" -w /directory quay.io/biocontainers/blast-legacy:2.2.26--2 formatdb -i database.fasta -p F
```

To perform a blastn search using the formatted database and a query called `query.fasta` when the file is also located in the current directory:

```bash
docker run -it --rm -v "$(pwd)":/directory -u "$(id -u)":"$(id -g)" -w /directory quay.io/biocontainers/blast-legacy:2.2.26--2 blastall -p blastn -d database.fasta -i query.fasta
```

To perform a blastn search using the formatted database and a query called `query.fasta` when the query is located in a different directory (in this example your home directory):

```bash
docker run -it --rm -v "$(pwd)":/directory/database -v "${HOME}":/directory/query -u "$(id -u)":"$(id -g)" -w /directory quay.io/biocontainers/blast-legacy:2.2.26--2 blastall -p blastn -d database/database.fasta -i query/query.fasta
```

### Annotate sequence variants using VEP

Download a VEP Docker image:

```bash
docker pull ensemblorg/ensembl-vep
```

Download cache files for the bovine genome and VEP plugins:

```bash
docker run -t -i -v "$(pwd)":/opt/vep/.vep -u "$(id -u)":"$(id -g)" ensemblorg/ensembl-vep perl INSTALL.pl -a cfp -s bos_taurus -y ARS-UCD1.2 -g all
```

Create directories for input and output files:

```bash
mkdir input
mkdir output
```

Copy the VCF files to be annotated to the newly created `input` directory and process them as follows:

```bash
find ./input -name "*.vcf" | while read f; do
    filename=$(basename -- "$f")
    filename_no_extension="${filename%.*}"
    docker run -v "$(pwd)":/opt/vep/.vep ensemblorg/ensembl-vep \
        -u "$(id -u)":"$(id -g)" \
        ./vep --cache --format vcf --tab --force_overwrite \
        --dir_cache /opt/vep/.vep/ \
        --dir_plugins /opt/vep/.vep/Plugins/ \
        --input_file /opt/vep/.vep/input/${filename} \
        --output_file /opt/vep/.vep/output/${filename_no_extension}.tab \
        --species bos_taurus --assembly ARS-UCD1.2 \
        --plugin Conservation,/opt/vep/.vep/Plugins/gerp_conservation_scores.bos_taurus.ARS-UCD1.2.bw --plugin Blosum62 --plugin Downstream --plugin Phenotypes --plugin TSSDistance --plugin miRNA \
        --variant_class --sift b --nearest gene --overlaps --gene_phenotype --regulatory --protein --symbol --ccds --uniprot --biotype --domains --check_existing --pubmed \
        --verbose
done
```

### Annotate a bacterial genome using Prokka

Download a Prokka Docker image:

```bash
docker pull staphb/prokka:latest
```

Create a container from the image and run `prokka` to annotate the sequence. In this example the genome sequence to be annotated is in a file called `sequence.fasta`, located in the current directory, and four CPUs are used:

```bash
docker run --rm -v "$(pwd)":/dir -u "$(id -u)":"$(id -g)" -w /dir staphb/prokka:latest prokka sequence.fasta --cpus 4
```

### Compare sequence reads to a bacterial genome to find SNPs using Snippy

Download a Snippy Docker image:

```bash
docker pull quay.io/biocontainers/snippy:4.6.0--hdfd78af_1
```

Create a container from the image and run `snippy` to find SNPs. In this example the reference sequence is in a file called `sequence.gbk`, and the reads are in `J10_S210_R1_001.fastq` and `J10_S210_R2_001.fastq`:

```bash
docker run -it --rm \
-u "$(id -u)":"$(id -g)" \
-v "$(pwd)":/directory \
-w /directory \
quay.io/biocontainers/snippy:4.6.0--hdfd78af_1 \
snippy --cpus 4 --outdir output --reference sequence.gbk --R1 J10_S210_R1_001.fastq --R2 J10_S210_R2_001.fastq
```

### Get an interactive bash shell in a Docker container

```bash
docker pull quay.io/biocontainers/snippy:4.6.0--hdfd78af_1
docker run -it quay.io/biocontainers/snippy:4.6.0--hdfd78af_1 /bin/bash
```

### List images

```bash
docker image ls
```

### List running containers

```bash
docker container ls
```

### Stop a container

```bash
docker container stop some_container
```

### Kill all running containers

```bash
docker container kill $(docker ps -q)
```

Or to continually kill containers:

```bash
while true; do docker container kill $(docker ps -q); sleep 2; done
```

### Delete all containers that are not running

```bash
docker container rm $(docker ps -a -q)
```

## File conversion

### Convert CSV to TSV

```bash
perl -nle  'my @new  = (); push( @new, $+ ) while $_ =~ m{"([^\"\\]*(?:\\.[^\"\\]*)*)",? | ([^,]+),? | ,}gx; push( @new, undef ) if substr( $text, -1, 1 ) eq '\'','\''; for(@new){s/,/ /g} print join "\t", @new' input.csv > output.tab
```

Or use [VisiData](https://github.com/saulpw/visidata):

```bash
vd input.csv -b -o output.tab
```

### Convert TSV to CSV

```bash
awk 'BEGIN { FS="\t"; OFS="," } {
  rebuilt=0
  for(i=1; i<=NF; ++i) {
    if ($i ~ /,/ && $i !~ /^".*"$/) {
      gsub("\"", "\"\"", $i)
      $i = "\"" $i "\""
      rebuilt=1
    }
  }
  if (!rebuilt) { $1=$1 }
  print
}' input.tsv > output.csv
```

Or use [VisiData](https://github.com/saulpw/visidata):

```bash
vd input.tsv -b -o output.tab
```

### Convert CSV to Markdown

The following uses [csv2md](https://github.com/pstaender/csv2md). The `awk` command can be used to change missing values to `.`:

```bash
awk 'BEGIN { FS = OFS = "," } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "." }; 1' input.csv > temp.csv
csv2md -p < temp.csv | sed 's/_/\\_/g' > output.md
```

Or use [VisiData](https://github.com/saulpw/visidata):

```bash
vd input.csv -b -o output.md
```

### Convert CSV to HTML

Use [VisiData](https://github.com/saulpw/visidata):

```bash
vd input.csv -b -o output.html
```

### Convert PDF to PNG

The following uses `find` and the `pdftoppm` command from the [poppler](https://poppler.freedesktop.org) package to generate a PNG image of the first page of every PDF file in the working directory:

```bash
find . -name "*.pdf" -exec pdftoppm -f 1 -l 1 -png {} {} \;
```

### Convert PNG to PDF

The following uses [ImageMagick](https://imagemagick.org):

```bash
convert *.png output.pdf
```

### Convert DOCX to PDF

The following uses [LibreOffice](https://www.libreoffice.org):

```bash
soffice --headless --convert-to pdf --outdir . word_file.docx
```

The following uses [pandoc](https://pandoc.org) and on macOS also requires [basictex](https://www.tug.org/mactex/morepackages.html):

```bash
pandoc word_file.docx --output word_file.pdf
```

### Convert Excel to CSV with csvkit

The following uses [csvkit](https://github.com/wireservice/csvkit):

```bash
in2csv data.xls > data.csv
```

Or use [VisiData](https://github.com/saulpw/visidata). Change `Sheet1` in the command below to match the name of the sheet to be converted:

```bash
vd input.xls +:Sheet1:1:1 -b -o output.csv
```

### Convert CSV to Excel

The following uses `ssconvert`, which is distributed with [Gnumeric](http://gnumeric.org):

```bash
ssconvert input.csv output.xlsx
```

Or use [VisiData](https://github.com/saulpw/visidata):

```bash
vd input.csv -b -o output.xlsx
```

### Convert HTML to PDF

The following uses [wkhtmltopdf](https://wkhtmltopdf.org):

```bash
wkhtmltopdf http://google.com google.pdf
```

### Convert a website to PDF

The following uses [wkhtmltopdf](https://wkhtmltopdf.org) and [gs](https://www.ghostscript.com/index.html):

```bash
url=https://sites.ualberta.ca/~stothard/; depth=1
wget --spider --force-html -r -l${depth} ${url} 2>&1 | grep '^--' | awk '{ print $3 }' | grep -i '\.\(html\|htm\)$' | uniq > url-list.txt
while read i; do wkhtmltopdf "$i" "$(echo "$i" | sed -e 's/https\?:\/\///' -e 's/\//-/g' ).pdf"; done < url-list.txt
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=merged-output.pdf $(ls -lrt -1 *.pdf)
```

### Convert HTML to PNG

The following uses [wkhtmltoimage](https://wkhtmltopdf.org):

```bash
wkhtmltoimage -f png http://google.com google.png
```

Another approach, which may work better for complex web sites, is to use [pageres-cli](https://github.com/sindresorhus/pageres-cli). The following creates an image that is 4485 pixels wide (5 * 897):

```bash
pageres http://google.com 897x1090 --crop --scale=5 --filename='google'
```

### Convert Markdown to PDF

The command below uses [pandoc](https://pandoc.org) and the [eisvogel.tex template](https://github.com/Wandmalfarbe/pandoc-latex-template/blob/master/eisvogel.tex).

The `head.tex` file consists of the following:

```tex
\definecolor{bgcolor}{HTML}{E0E0E0}
\let\oldtexttt\texttt

\renewcommand{\texttt}[1]{
  \colorbox{bgcolor}{\oldtexttt{#1}}
}
```

```bash
pandoc input.md -o output.pdf --pdf-engine=xelatex --from markdown --template=eisvogel.tex --highlight-style zenburn -H head.tex
```

### Convert Markdown to HTML

The command below uses [pandoc](https://pandoc.org) and the `pandoc.css` file available [here](https://gist.github.com/killercup/5917178).

```bash
pandoc -f markdown -t html -o output.html input.md --css=pandoc.css --self-contained
```

### Convert Markdown to PPTX

The command below uses [pandoc](https://pandoc.org). See the [pandoc documentation on slide shows](https://pandoc.org/MANUAL.html#slide-shows).

The `--reference-doc` parameter is optional and can be used to supply a [PowerPoint template](https://support.microsoft.com/en-us/office/create-and-save-a-powerpoint-template-ee4429ad-2a74-4100-82f7-50f8169c8aca).

```bash
pandoc input.md -o output.pptx --reference-doc template.potx
```

### Convert PPTX to PDF

The following uses [LibreOffice](https://www.libreoffice.org):

```bash
soffice --headless --convert-to pdf input.pptx
```

### Convert PPTX to PNG

The following uses [LibreOffice](https://www.libreoffice.org) and [ImageMagick](https://imagemagick.org) to create a separate file for each slide (e.g. `slide_1.png`, `slide_2.png`, etc.):

```bash
soffice --headless --convert-to pdf input.pptx
convert -density 400 -resize 3000^ -scene 1 input.pdf slide_%d.png
```

## File downloads

### Download an entire website

```bash
wget --mirror --convert-links --page-requisites --no-parent https://www.somesite.com/course/material/
```

### Download from an FTP server

Replace `host`, `account`, `password`, and `port` with their corresponding values in the following command:

```bash
wget -S -d -c -t 45 -v -r ftp://account:password@host:port/*
```

### Download from Google Drive

[Rclone](https://rclone.org) can be used to download data from many cloud storage providers.

First, follow the [configuration instructions](https://rclone.org/drive/). The commands below assume that the remote system was named `my_google_drive` during the configuration.

Note that you can use the `Add shortcut to Drive` option in Google Drive to make folders and files in `Shared with me` easier to access using rclone.

To list remote drives:

```bash
rclone listremotes
```

To list directories:

```bash
rclone lsd my_google_drve:
```

To list the contents of a directory

```bash
rclone ls my_google_drive:some_directory
```

To copy the drive to local storage:

```bash
rclone copy -P my_google_drive:some_directory ./some_directory
```

### Download a GenBank file with curl

The command below downloads a GenBank file from the Nucleotide database at NCBI:

```bash
i=NC_045512.2
curl -s  "https://eutils.ncbi.nlm.nih.gov\
/entrez/eutils/efetch.fcgi?db=nucleotide\
&id=${i}&rettype=gbwithparts&retmode=txt" \
> $i.gbk
```

### Download fastq files based on a list of SRA accessions

The following uses the [SRA Toolkit](https://github.com/ncbi/sra-tools).

Paired-end data:

```bash
cat SRR_Acc_List.txt | xargs -I{} fastq-dump -I --split-files --gzip {}
```

Single-end data:

```bash
cat SRR_Acc_List.txt | xargs -I{} fastq-dump --gzip {}
```

For better performance use `fasterq-dump` (the following works for single-end and paired-end data):

```bash
cat SRR_Acc_List.txt | xargs -I{} fasterq-dump {}
gzip *.fastq
```

Note that `fasterq-dump` will generate large cache files in `~/ncbi/public/sra`. If this directory does not have sufficient storage create a symbolic link to another directory that does, for example `/scratch`:

```bash
mv ~/ncbi/public/sra /scratch/
ln -s /scratch/sra ~/ncbi/public/sra
```

[SRA Explorer](https://sra-explorer.info/#) is an online resource that takes a list of accessions and returns a selectable list of ENA download URLs and sequencing run metadata.

### Download a reference genome FASTA file from Ensembl

To obtain a list of files available for a particular species:

```bash
rsync --list-only rsync://ftp.ensembl.org/ensembl/pub/current_fasta/oncorhynchus_mykiss/dna/
```

To download one of the files:

```bash
rsync -av --progress rsync://ftp.ensembl.org/ensembl/pub/current_fasta/oncorhynchus_mykiss/dna/Oncorhynchus_mykiss.USDA_OmykA_1.1.dna.toplevel.fa.gz .
```

### Download a reference genome GTF file from Ensembl

To obtain a list of files available for a particular species:

```bash
rsync --list-only rsync://ftp.ensembl.org/ensembl/pub/current_gtf/oncorhynchus_mykiss/
```

To download one of the files:

```bash
rsync -av --progress rsync://ftp.ensembl.org/ensembl/pub/current_gtf/oncorhynchus_mykiss/Oncorhynchus_mykiss.USDA_OmykA_1.1.106.gtf.gz .
```

## find

### Perform a series of commands on files returned by find

The command below finds `.gff` files and then each file is processed as follows:

1. `tail` is used to skip a header line.
2. `awk` is used to count the number of occurrences of each category in column 3 and print the category and counts.
3. `sort` is used to sort the categories by count from largest to smallest with ties broken by sorting on category name.
4. The results are redirected to a file named after the input file but with `.cog_counts.txt` appended.

```bash
find . -type f -name "*.gff" -print0 | xargs -0 -I{} sh -c $'tail -n +2 "$1" | awk -F $\'\t\' \'{count[$3]++}END{for(j in count) print j,count[j]}\' | sort -k 2,2nr -k 1,1> "$1.cog_counts.txt"' -- {}
```

The command below finds `.stats` files and then each file is processed as follows:

1. `count` is used to store the value in a field called `number of records`, which is obtained using `perl`.
2. `echo` is used to print the filename and `count`.
3. `perl` is used to remove `./` from the filename.
4. The filenames and counts for all the files are then sorted numerically by the count value, using `sort`.
5. `awk` is used to add a header row and `column` is used to format the output.

```bash
find . -name "*.stats" -exec sh -c $'count=$(perl -0777 -ne \'while (m/number of records:\s+?(\d+)/gm) {$out = $1; $out =~ s/\s//g; print "$out"}\' "$1"); echo "$1 $count" | perl -p -e \'s/\.\///g\'' -- {} \; | sort -k 2,2rn | awk 'BEGIN{print "file variants"}1' | column -t
```

Often it is simpler to use a loop to iterate through the files returned by `find`. The following uses `awk` to filter `.maf` files:

```bash
wd=$(pwd)
find . -name "*.maf" -type f | while IFS= read -r file; do
  dir=$(dirname -- "$file")
  fnx=$(basename -- "$file")
  fn="${fnx%.*}"

  echo "Processing file '$fnx' in directory '$dir'"
  cd "$dir"

  #always print the first two lines
  #print lines where value in column 101 is PASS
  awk -F $'\t' 'NR < 3; NR > 2 {if ($101 == "PASS") print $0}' "$fnx" > "${fn}.new"
  mv "$fnx" "${fnx}.bak"
  mv "${fn}.new" "$fnx"

  cd "$wd"
done
```

To remove the `.bak` files:

```bash
find . -name "*.bak" -type f -exec rm -rf {} \;
```

### Sort files before processing

The following used `find` to generate a list of `.vcf.gz` files, which is then sorted based on sample number using `sort`.

The files in this example are named `DNA-1.vcf.gz`, `DNA-2.vcf.gz`, etc. and the `sort` is used to sort based on the number appearing after the first `-` in the filename.

```bash
find . -name "*.vcf.gz" -type f -print0 | sort -z -k2,2n -t- | \
while IFS="" read -r -d "" file; do
  fnx=$(basename -- "$file")
  echo "'$fnx' contains these samples:"
  bcftools query -l $file
done
```

### Use files as the argument list for a command

The following finds files ending with `.vcf.gz`, sorts the files based on the number appearing after the first `-` in the filename, and prints the filenames separated by spaces by supplying all to `echo`:

```bash
find . -name "*.vcf.gz" -type f -print0 | sort -z -k2,2n -t- | xargs -r0 echo
```

To print the filenames on separate lines, use:

```bash
find . -name "*.vcf.gz" -type f -print0 | sort -z -k2,2n -t- | xargs -r0 -L 1 echo
```

### Process files in pairs

High-throughput sequencing data is often distributed as pairs of files corresponding to the two different read sets generated for each sample, e.g.:

```text
6613_S82_L001_R1_001.fastq.gz
6613_S82_L001_R2_001.fastq.gz
70532_S37_L001_R1_001.fastq.gz
70532_S37_L001_R2_001.fastq.gz
k2712_S5_L001_R1_001.fastq.gz
k2712_S5_L001_R2_001.fastq.gz
```

To analyze data from multiple samples, the following `while` loop code can be used. It iterates through the `R1` files and from each filename constructs the matching `R2` filename. Two useful variables called `fnx` and `fn` are also created for each file, storing the filename without the path to the file, and the filename without the path and without the file extension, respectively:

```bash
find . -name "*_R1_*" -type f | while IFS= read -r file; do
  fnx=$(basename -- "$file")
  fn="${fnx%.*}"

  #Construct name of other file
  file2="${file/R1_001.fastq.gz/R2_001.fastq.gz}"
  fnx2=$(basename -- "$file2")
  fn2="${fnx2%.*}"

  echo "Processing file '$fnx' and '$fnx2'"

done
```

Another option is to use [parallel](#parallel).

### Copy the files returned by find

The following finds files in the directory `output` ending in `.vcf.gz` or `.vcf.gz.tbi` and copies them to the `vcfs` directory:

```bash
mkdir vcfs
find output -type f \( -name "*.vcf.gz" -o -name "*.vcf.gz.tbi" \) -exec cp {} vcfs \;
```

### Copy the files returned by find, naming the copies after a directory in the path

The command below finds files named `star-fusion.fusion_candidates.preliminary` and parses the sample name from a directory name in the path to the file. The sample name is then used to construct a name for the copy. For example, `./231_S12_R1_001/star-fusion.fusion_candidates.preliminary` is copied to `./fusion-candidates/231_S12_R1_001.fusion_candidates.preliminary`.

```bash
find . -name star-fusion.fusion_candidates.preliminary -exec sh -c $'sample=$(perl -e \'if($ARGV[0] =~ m/^\.\/([^\/]+)/){print "$1\n"}\' $1); cp "$1" "./fusion-candidates/${sample}.fusion_candidates.preliminary"' -- {} \;
```

### Switch to the directory containing each file and execute a command

The -execdir option instructs `find` to switch to the directory containing each matching file before executing the specified command. In this example the command creates a `.zip` file for each `.vcf` file that is found:

```bash
find . -name "*.vcf" -type f -execdir zip '{}.zip' '{}' \;
```

### Find large files

The following reports the 10 largest files in the current directory or its subdirectories, sorted by size:

```bash
find . -type f -print0 | xargs -0 du -h | sort -hr | head -10
```

### Process multiple files and redirect output to one file

Print the number of lines in every `.csv` or `.tab` file in or below current directory and redirect the results to a single file:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) | while read f; do wc -l "$f" >> output.txt; done
```

Or:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -exec wc -l {} \; > output.txt
```

Or:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -print0 | xargs -0 -I{} wc -l {} > output.txt
```

### Process multiple files and redirect output to separate files

Print the number of lines in every `.csv` or `.tab` file in or below current directory and redirect the results to separate files:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) | while read f; do wc -l "$f" > "${f}.output.txt"; done
```

Or:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -exec sh -c 'wc -l "$1" > "$1.output.txt"' -- {} \;
```

Or:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -print0 | xargs -0 -I{} sh -c 'wc -l "$1" > "$1.output.txt"' -- {}
```

## Git

See [GitHub's Git documentation](https://help.github.com/en) for more information

### Create a new Git repository

```bash
git init
```

### Sync a repository to your local machine

First, copy the clone URL on the GitHub repository page by clicking `Clone or Download`. Then, enter the following command in a terminal window. The helpful\_commands repository is used as an example:

```bash
git clone https://github.com/stothard-group/helpful_commands.git
```

### Clone all your repositories

Use the official [GitHub CLI tool](https://github.com/cli/cli) `gh`.

Change `org_or_username` to your username to download all your repositories, or to your organization name to download all the repositories of your organization:

```bash
org_or_username=paulstothard
```

Authenticate with a GitHub host:

```bash
gh auth login
```

Clone up to 1000 repositories under the `$org_or_username` folder:

```bash
gh repo list "$org_or_username" --limit 1000 | while read -r repo _; do
  gh repo clone "$repo" "$repo"
done
```

Optionally, to clone new repositories and update existing ones:

```bash
gh repo list "$org_or_username" --limit 1000 | while read -r repo _; do
  gh repo clone "$repo" "$repo" -- -q 2>/dev/null || (
    cd "$repo"
    git checkout -q main 2>/dev/null || true
    git checkout -q master 2>/dev/null || true
    git pull -q
  )
done
```

### Mark changed files to be included in the next commit

To add one or more files:

```bash
git add <filename1> <filename2>
```

To add all current modifications in your project (including deletions and new files):

```bash
git add --all
```

### Undo a Git add before a commit

To undo a list of files:

```bash
git reset <filename1>
```

To undo all changes:

```bash
git reset
```

### Remove files from the repository

Note that the following instructions will remove the file/directory from both the working tree and the index.

To remove one or more files:

```bash
git rm <filename>
```

To remove a directory:

```bash
git rm -r <directory>
```

To remove a file from the index (this untracks the file, it does not delete the file itself):

```bash
git rm --cached <filename>
```

These changes must be committed with git commit.

### Move or rename a file or directory

```bash
git mv <filename-old> <filename-new>
```

This change must be committed with git commit.

### Save the marked files to the local Git repository

The commit should include a message using the -m option:

```bash
git commit -m "A concise description of the changes"
```

The following changes can be made to commits that have **not** been pushed to a remote repository. To rewrite the very last commit, with any currently staged changes:

```bash
git commit --amend -m "An updated message"
```

To commit any currently staged changes without rewriting the commit (this essentially adds the staged changes to the previous commit):

```bash
git commit --amend --no-edit
```

### Push a commit on your local branch to a remote repository

```bash
git push <remote> <branch>
```

For example, to push to the main branch:

```bash
git push -u origin main
```

### Pull a change from a remote repository to your local branch

```bash
git pull <remote> <branch>
```

For example, to pull from the main branch:

```bash
git pull origin main
```

### Add or edit a remote repository

To add a new remote:

```bash
git remote add origin <repo-url>
```

To edit an existing remote:

```bash
git remote set-url origin <new-repo-url>
```

To verify that the remote URL has changed:

```bash
git remote -v
```

### Create and merge Git branches

To view the branches in a repository:

```bash
git branch
```

To create a new branch and switch to it:

```bash
git checkout -b <new-branch>
```

To switch to a remote branch:

```bash
git fetch --all
git checkout <remote-branch>
```

After adding and committing some changes, to push this branch to remote:

```bash
git push -u origin <new-branch>
```

To merge a branch into main (local) and push the changes to remote:

```bash
git checkout main
git merge <new-branch>
git push -u origin main
```

Git merge conflicts can arise easily. For information on resolving a merge conflict, see [Resolving a merged conflict using the command line](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/resolving-a-merge-conflict-using-the-command-line).

### Specify files to ignore

Create a `.gitignore` file:

```bash
touch .gitignore
```

Add text or patterns to exclude:

```bash
echo sensitive_data.txt >> .gitignore
echo test/*.vcf >> .gitignore
git add .gitignore
git commit -m "add .gitignore file"
git push -u origin main
```

In this example, the following files will no longer be tracked: `sensitive_data.txt`, and all files with a `.vcf` extension in the directory `test`.

Note that adding a `.gitignore` file will not remove tracked files; this must be done with `git rm`. See [Remove files from the repository](#remove-files-from-the-repository).

### Check the status of a working directory

```bash
git status
```

### Tag a release

A tag allows a specific release version of code to be identified, and creates a release that can be downloaded from GitHub. A tagged version serves as a snapshot that does not change.

```bash
git tag -a v2.0.0 -m "version 2.0.0"
git push origin v2.0.0
```

Information on how to choose version numbers if available [here](https://semver.org).

## grep

### Count matches

In this example the number of lines with a match to `>` is returned:

```bash
grep -c ">" input.fasta
```

### Get the line number of a match

In this example the line numbers of lines with a match to `234829` are reported:

```bash
grep -n "234829" input.txt
```

### Remove files that contain a match

In this example `.fasta` files are removed that contain the text `complete genome` on a single line:

```bash
grep -l "complete genome" *.fasta | xargs -I{} rm -f {}
```

### Remove files that do not contain a match

In this example `.fasta` files are removed that do not contain the text `complete genome` on a single line:

```bash
grep -L "complete genome" *.fasta | xargs -I{} rm -f {}
```

### Print non-matching lines

Keep everything except lines starting with `#`:

```bash
grep -v '^#' input.txt
```

### Search using patterns from a file

This example uses patterns from the text file `ensembl_ids.txt`, which contains one gene ID per line:

```text
ENSCAFG00000018897
ENSCAFG00000006950
ENSCAFG00000013069
ENSCAFG00000013670
ENSCAFG00000003247
```

The first `grep` command is used to add the VCF header to the output file:

```bash
grep '^#' annotated_variants.vcf > annotated_variants.candidates.vcf
grep -f ensembl_ids.txt annotated_variants.vcf >> annotated_variants.candidates.vcf
```

## Image files

### Crop an image and add a white border

The following uses [ImageMagick](https://imagemagick.org) to removes any edges that are exactly the same color as the corner pixels. A 30-pixel white border is then added to the sides and the top and bottom:

```bash
convert input.png -trim -bordercolor White -border 30x30 output.png
```

### Resize an image

The following uses [ImageMagick](https://imagemagick.org) to scale the image so that its width is 4000 pixels:

```bash
convert input.png -resize 4000 output.png
```

### Take a screenshot of a window on macOS

1. Press `Command+Shift+4`.
2. Press the `space bar`.
3. Hold `Option` and click on the window.

To keep the drop shadow perform the last step without holding `Option`.

### Take a webpage screenshot using Firefox

1. Press `F12` to open the Firefox Developer Tools.
2. Enter `:screenshot` into the Web Console to download the current view as a PNG file.

To save a high-DPI webpage screenshot use `:screenshot --dpr 4`.

To save a high-DPI full-page webpage screenshot use `:screenshot --dpr 4 --fullpage`.

To delay the screenshot, to capture a menu for example, use `:screenshot --dpr 4 --fullpage --delay 5`.

### Record your terminal to an animated GIF

Use [Terminalizer](https://github.com/faressoft/terminalizer) to record a session. The following creates a file called `session.yml`:

```bash
terminalizer record session
```

To convert the session to an animated GIF:

```bash
terminalizer render session
```

For more control over rendering options create an editable configuration file at `~/.terminalizer/config.yml`:

```bash
terminalizer init
```

### Record your screen to an animated GIF

Use [LICEcap](https://www.cockos.com/licecap/).

### Create an animated GIF from a YouTube video

The following requires [youtube-dl](https://github.com/ytdl-org/youtube-dl), [mplayer](https://mplayerhq.hu/), [ImageMagick](https://imagemagick.org), and [gifsicle](https://www.lcdf.org/gifsicle/):

```bash
mkdir gif; cd gif
url=https://youtu.be/_YUAu0aP4DA
start=00:37; length=10
youtube-dl -f mp4 -o video_for_gif.mp4 $url
mplayer video_for_gif.mp4 -ao null -ss $start -endpos $length -vo png -vf scale=400:225
mogrify -format gif *.png
gifsicle --threads=2 --colors=256 --delay=4 --loopcount=0 --dither -O3 *.gif > animation.gif
```

## join

### Combine rows based on shared keys with join

`join` is used to join lines from different files when they have matching join fields. `join` expects the input files to be sorted.

Suppose `file1.txt` contains the following text:

```text
name,counts in sample 1
gene a,100
gene b,2223
gene e,575
```

Suppose `file2.txt` contains this text:

```text
counts in sample 2,name
2223,gene e
803,gene a
0,gene f
```

To join their contents based on `name`, first sort the files based on the `name` column (note that the sort field differs between the two commands):

```bash
(head -n 1 file1.txt && sort -t, -k 1 <(tail -n +2 file1.txt)) > file1.txt.sorted
(head -n 1 file2.txt && sort -t, -k 2 <(tail -n +2 file2.txt)) > file2.txt.sorted
```

Then use the `join` command to combine the rows based on field `1` in the first file and field `2` in the second file:

```bash
join -t, -1 1 -2 2 file1.txt.sorted file2.txt.sorted
```

The above command returns the following:

```text
name,counts in sample 1,counts in sample 2
gene a,100,803
gene e,575,2223
```

To format more nicely use:

```bash
join -t, -1 1 -2 2 file1.txt.sorted file2.txt.sorted | column -t -s,
```

Which gives:

```text
name    counts in sample 1  counts in sample 2
gene a  100                 803
gene e  575                 2223
```

To include genes present in just one of the files and to add a missing value character `.` for the corresponding values that could not be obtained use:

```bash
join -t, -1 1 -2 2 -a1 -a2 -e . -o auto file1.txt.sorted file2.txt.sorted | column -t -s,
```

Which gives:

```text
name    counts in sample 1  counts in sample 2
gene a  100                 803
gene b  2223                .
gene e  575                 2223
gene f  .                   0
```

Another option is to use [csvjoin](#merge-csv-files-on-a-specified-column-or-columns) from [csvkit](#csvkit).

## Mamba

[Mamba](https://github.com/mamba-org/mamba) is a reimplementation of the conda package manager in C++.

### Install Mamba

```bash
conda install mamba -n base -c conda-forge
```

### Create an environment and install some packages with mamba

In this example an environment called `ngs` is created:

```bash
mamba create -y --name ngs
conda activate ngs
mamba install -y -c bioconda -c conda-forge multiqc fastqc trimmomatic bowtie2 subread samtools
```

### Deactivate an environment with mamba

```bash
conda deactivate
```

### Activate an environment with mamba

```bash
conda activate ngs
```

### List available packages with mamba

```bash
mamba search -c bioconda -c conda-forge
```

### Search for a specific package with mamba

```bash
mamba search -c bioconda -c conda-forge blast
```

### Add additional packages to an environment with mamba

```bash
conda activate ngs
mamba install -y -c bioconda -c conda-forge picard
```

### List environments with mamba

```bash
mamba info --envs
```

### List packages installed in the active environment with mamba

```bash
mamba list
```

### Remove an environment with mamba

In this example the environment to remove is called `my-env`:

```bash
conda deactivate
mamba env remove --name my-env
```

## md5sum

### Generate a file of checksums

```bash
md5sum *.vcf.gz > md5sum.txt
```

### Validate checksums

```bash
md5sum --check md5sum.txt
```

## Miller

[Miller](https://github.com/johnkerl/miller) can be used to work with CSV, TSV, and JSON files.

### Extract the first 10 records of a CSV file

Print the column names followed by the records:

```bash
mlr --csv head -n 10 example.csv
```

Print with added formatting for readability:

```bash
mlr --icsv --opprint head -n 10 example.csv
```

Print each value with the column name in the form `column=value`:

```bash
mlr --icsv --odkvp head -n 10 example.csv
```

### Extract the last 10 records of a CSV file

Print the column names followed by the records:

```bash
mlr --csv tail -n 10 example.csv
```

Print with added formatting for readability:

```bash
mlr --icsv --opprint tail -n 10 example.csv
```

Print each field with the column name in the form `column=field`:

```bash
mlr --icsv --odkvp tail -n 10 example.csv
```

### Convert formats

From CSV to JSON:

```bash
mlr --icsv --ojson cat example.csv
```

From CSV to TSV:

```bash
mlr --icsv --otsv cat example.csv
```

### View stats

To view the `stats1` and `stats2` documentation:

```bash
mlr stats1 --help
mlr stats2 --help
```

To view several stats for the `Coverage` column:

```bash
mlr --icsv --opprint stats1 -a sum,count,min,max,mean,mode -f Coverage example.csv
```

### Filter records

To view the `filter` documentation:

```bash
mlr filter --help
```

The following filters records based on `Breed`, `Dataset`, and `Coverage`:

```bash
mlr --csv filter '$Breed == "Charolais" && $Dataset == "A" && $Coverage > 6' example.csv
```

### Sort records

To view the `sort` documentation:

```bash
mlr sort --help
```

The following first sorts alphabetically by `Breed` and then numerically by `Coverage` (from largest to smallest):

```bash
mlr --icsv --opprint sort -f Breed -nr Coverage example.csv
```

### Extract columns

To view the `cut` documentation:

```bash
mlr cut --help
```

The following extracts the `Breed` and `Coverage` columns:

```bash
mlr --csv cut -f Breed,Coverage example.csv
```

### Edit columns

To view the `put` documentation:

```bash
mlr put --help
```

The following converts `Breed` to uppercase and divides `Coverage` by 100:

```bash
mlr --csv put '$Breed = toupper($Breed); $Coverage = ($Coverage / 100)' example.csv
```

New columns can be created:

```bash
mlr --csv put '$New_coverage = ($Coverage / 100)' example.csv
```

### Other actions

To view a complete list of Miller verbs use the following:

```bash
mlr -l
```

To view documentation for a particular verb use `mlr _verb_ --help`.

### Combine actions

Perform multiple actions sequentially using `then`:

```bash
mlr --csv put '$New_coverage = ($Coverage / 100)' then sort -f Breed -nr Coverage then cut -f InterbullID,New_coverage example.csv
```

## Other

### Obtain your public IP address and network information

```bash
curl ifconfig.me/all
```

### Copy an ssh public key to another system

Generate the key pair:

```bash
ssh-keygen
```

Copy the public key to the `.ssh/authorized_keys` file on the other system using `ssh-copy-id`:

```bash
ssh-copy-id -i ~/.ssh/id_rsa.pub user@remote-host.com
```

### Add a header to all files with a certain extension, getting the header from another file

In this example the header is added to `.tab` files and comes from a file called `header.txt`. The files with the header added are saved with a `.new` extension added:

```bash
for f in *.tab; do new=`echo $f | sed 's/\(.*\)\.tab/\1.tab.new/'`; paste -sd'\n' \header.txt "$f" > "$new"; done
```

To replace the `.tab` files the `.new` files:

```bash
for f in *.new; do new=`echo $f | sed 's/\(.*\)\.new/\1/'`; mv "$f" "$new"; done
```

### Change the extension of multiple files

The following changes the `.gbff` extension to `.gbk`:

```bash
for f in *.gbff; do
    mv -- "$f" "${f%.gbff}.gbk"
done
```

If `rename` is available, this may work:

```bash
rename 's/\.gbff$/.gbk/' *.gbff
```

Or this, depending on which `rename` is installed:

```bash
rename .gbff .gbk *.gbff
```

### Add text to the beginning of a file

```bash
echo 'new header line' | cat - file.txt > temp && mv temp file.txt
```

### Add text or a header to the beginning of all files with a particular file extension

Add `my header text` to the start of all `.csv` files in the current directory (works on macOS):

```bash
find . -name "*.csv" -exec sed -i '.bak' '1s/^/my header text\'$'\n/g' {} \;
```

Add `my header text` to the start of all `.csv` files in the current directory (works on Linux):

```bash
find . -name "*.csv" -exec sed -i '1s/^/my header text\n/' {} \;
```

### Find common lines between files

Between two files, named `file1.txt` and `file2.txt`:

```bash
comm -12 <( sort file1.txt ) <( sort file2.txt )
```

Among all `.txt` files in the current directory:

```bash
number_of_files=$(find . -name "*.txt" -print | wc -l | sed 's/[^0-9]*//g')
cat *.txt | sort | uniq -c | sed -n -e "s/^ *$number_of_files \(.*\)/\1/p"
```

### Format a CSV file into columns and examine its content

The `perl` portion is used to handle empty fields:

```bash
cat data.csv | perl -pe 's/((?<=,)|(?<=^)),/ ,/g;' | column -t -s, | less -S
```

### Format code

The [Prettier](https://prettier.io) program supports many languages. The following command uses the `--write` option to reformat files in-place:

```bash
prettier --write "*html"
```

### Run commands at scheduled times using cron

The following uses `cron` to run a script to copy various files and directories to a directory backed up by Dropbox.

Create the script the `copy_to_dropbox.sh` script, editing as needed:

```bash
#!/bin/bash

PATH=/bin:/usr/bin/

home=/Users/myhome

target=${home}/Dropbox/backup

if [[ ! -e $target ]]; then
    mkdir -p $target
elif [[ ! -d $target ]]; then
    echo "$target already exists but is not a directory" 1>&2
fi

#copy files and directories of interest to $target
rsync --update -razv ${home}/.bash_profile $target/bash_profile
rsync --update -razv ${home}/lib $target
rsync --update -razv ${home}/bin $target
```

Test the script as follows:

```bash
chmod u+x copy_to_dropbox.sh
sudo env -i ./copy_to_dropbox.sh
```

Use `crontab -e` to edit the crontab (cron table) file. Add the following to the crontab to run the script everyday at noon (changing the path to `copy_to_dropbox.sh`):

```bash
0 12 * * * /path/to/copy_to_dropbox.sh >/dev/null 2>&1
```

Alternatively, use the following to run the script once per hour between 8 am and 5 pm on weekdays:

```bash
0 8-17 * * 1-5 /path/to/copy_to_dropbox.sh >/dev/null 2>&1
```

To display the crontab use `crontab -l`.

### Create a collection of MP3 files from a YouTube playlist

The following requires [youtube-dl](https://github.com/ytdl-org/youtube-dl) and [ffmpeg](https://ffmpeg.org/):

```bash
playlist_URL=https://www.youtube.com/playlist?list=PL92319EECC1754042
youtube-dl -x -i --audio-format mp3 --audio-quality 320K --embed-thumbnail --geo-bypass --rm-cache-dir --continue "$playlist_URL"
```

### Perform a calculation using bc

Use [bc](https://www.gnu.org/software/bc/):

```bash
echo "2*(42+42)" | bc
```

### Count the bases in a fastq file

```bash
zcat SRR13388732_1.fastq.gz | paste - - - - | cut -f 2 | tr -d '\n' | wc -c
```

Or:

```bash
cat SRR13388732_1.fastq | paste - - - - | cut -f 2 | tr -d '\n' | wc -c
```

### Count the reads in a fastq file

```bash
echo $(cat SRR13388732_1.fastq | wc -l) / 4 | bc
```

```bash
echo $(zcat SRR13388732_1.fastq.gz | wc -l) / 4 | bc
```

### Use SQL-like queries to work with a CSV or TSV file

The following uses [q](https://github.com/harelba/q) to count distinct values in the column named `SNP`:

```bash
q -H -d, "SELECT COUNT(DISTINCT(SNP)) FROM ./input.csv"
```

Another option is to use [csvsql](#query-with-sql) from [csvkit](#csvkit).

### Browse, search, and edit a large CSV file

Use [DB Browser for SQLite (DB4S)](https://github.com/sqlitebrowser/sqlitebrowser).

### Prevent a command from stopping when you log out or exit the shell

Use `nohup`:

```bash
nohup mycommand &
```

When using `&` the bash job ID is shown in brackets and the PID (process ID), e.g.:

```text
[1] 1963
```

The output of the command can be found in `nohup.out`.

The command can be stopped using `kill` and the PID, e.g.:

```bash
kill -9 1963
```

### Edit a PDF file

Use [LibreOffice](https://www.libreoffice.org).

## parallel

### Extract files in parallel

In the following command `{.}` is used to get the basename and remove the last extension of each input file:

```bash
parallel 'zcat {} > {.}.unpacked' ::: *.gz
```

Or:

```bash
parallel 'gunzip {}' ::: *.gz
```

### Compress files in parallel

In the following example 5 jobs are run at the same time:

```bash
parallel -j5 "gzip {}" ::: *.csv
```

### Process files in pairs using parallel

In the following example paired-end reads with names like `sampleA_1.fastq.gz` and `sampleA_2.fastq.gz` in a directory called `data` are mapped to a reference called `ref` using `bowtie2`:

```bash
parallel -j2 "bowtie2 --threads 4 -x ref -k1 -q -1 {1} -2 {2} -S {1/.}.sam >& {1/.}.log" ::: data/*_1.fastq.gz :::+ data/*_2.fastq.gz
```

The `{1/.}` removes the path and the extension from the first filename, leaving the sample name.

### Perform BLAST in parallel

Using the local system:

```bash
cat multiple_sequences.fasta | parallel --block 100k --recstart '>' --pipe blastp -evalue 0.01 -outfmt 6 -db database.fa -query - > results
```

Using the local system (denoted as `:` below) and a remote system called `server1` (connection details provided in `.ssh/config`):

```bash
cat multiple_sequences.fasta | parallel -S :,server1 --block 100k --recstart '>' --pipe blastp -evalue 0.01 -outfmt 6 -db database.fa -query - > results
```

## paste

### Combine columns with paste

In this example the columns of two files are joined. The first file is a CSV file the second is tab-delimited.

The `-d ","` specifies that the lines are to be joined with commas.

```bash
paste -d "," genotype_conversion.csv SNP_Map.tab
```

Note that the content from `SNP_Map.tab` still contains tab-delimited values.

To remove tabs from `SNP_Map.tab` you could first create a CSV version of that input file:

```bash
cut -d $'\t' -f 1- --output-delimiter=',' SNP_Map.tab > SNP_Map.csv
paste -d "," genotype_conversion.csv SNP_Map.csv
```

It is possible to do it without first creating `SNP_Map.csv` by using process substitution.

In the following the command between `<(` and `)` is first run and its output becomes input for `paste`:

```bash
paste -d "," genotype_conversion.csv \
<(cut -d $'\t' -f 1- --output-delimiter=',' SNP_Map.tab)
```

## Perl

### Get a random sample of lines from a text file while excluding the header line

In this example a random sample of 20 lines is obtained:

```bash
tail -n +2 input.txt | perl -MList::Util -e 'print List::Util::shuffle <>' | head -n 20 > output.txt
```

### Convert a FASTA file to a CSV file with column names

```bash
cat input.fasta | perl -n -0777 -e 'BEGIN{print "SNP_Name,Sequence\n"}' -e 'while ($_ =~ m/^>([^\n]+)\n([^>]+)/gm) {$name = $1; $seq = $2; $seq =~s/\s//g; print $name . "," . $seq . "\n"}' > output.csv
```

### Count the number of lines that match a regular expression

```bash
perl -lne '$a++ if /\tyes\t/; END {print $a+0}' < input.txt
```

### Extract FASTA sequences from a file based on a file of sequence names of interest

In this example the sequence names of interest are in the file `names.txt` and the FASTA sequences are in the file `input.fasta`:

```bash
cat names.txt | xargs -I{} perl -w -076 -e '$count = 0; open(SEQ, "<" . $ARGV[0]); while (<SEQ>) {if ($_ =~ m/\Q$ARGV[1]\E/) {$record = $_; $record =~ s/[\s>]+$//g; print ">$record\n"; $count = $count + 1;}} if ($count == 0) {print STDERR "No matches found for $ARGV[1]\n"} elsif ($count > 1) {print STDERR "Multiple matches found for $ARGV[1]\n"} close(SEQ);' input.fasta {} > output.fasta
```

### Add a FASTA title to the start of a sequence in RAW format

In this example the title `>KL1` is added to the beginning of the sequence in `KL1sequence.txt`:

```bash
perl -pi -e 'print ">KL1\n" if $. == 1' KL1sequence.txt
```

### Remove commas located within quoted fields in a CSV file and create a tab-delimited file

```bash
perl -nle  'my @new  = (); push( @new, $+ ) while $_ =~ m{"([^\"\\]*(?:\\.[^\"\\]*)*)",? | ([^,]+),? | ,}gx; push( @new, undef ) if substr( $text, -1, 1 ) eq '\'','\''; for(@new){s/,/ /g} print join "\t", @new' input.csv > output.tab
```

### Remove lines that match a pattern

In this example, lines from `input.txt` are written to `output.txt` unless they start with `some_text`:

```bash
perl -n -e 'print unless m/^some_text/' input.txt > output.txt
```

### Replace commas with tabs

```bash
perl -p -e 's/,/\t/g;' input.csv > output.tab
```

### Replace tabs with commas and remove quotes

```bash
perl -p -e 's/\t/,/g;' -e 's/"//g;' input.tab > output.csv
```

### Sort sections in a Markdown file based on headings

```bash
perl -0777 -ne '(undef,@paragraphs) = split /^#(?=[^#])/m; print map {"#$_"} sort { "\U$a" cmp "\U$b" } @paragraphs;' input.md
```

### Search and replace text on each line

```bash
perl -p -e 's/Red Angus/Angus/g' sequenced_samples.csv
```

Note that multiple substitutions can be performed in succession, e.g.:

```bash
echo " test pattern" | perl -pe 's/^\s+//g;' -pe 's/ /,/g;'
```

The above command produces the following:

```bash
test,pattern
```

### Replace ^M

```bash
perl -p -i -e 's/\r\n$/\n/g' file.txt
```

### Print matches that may span multiple lines

```bash
perl -0777 -ne 'while (m/^\s+\/translation="([^"]+)"/gm) {print "$1\n"}' NM_001271626.3.gbk
```

### Print matches after additional editing

```bash
perl -0777 -ne 'while (m/^\s+\/translation="([^"]+)"/gm) {$out = $1; $out =~ s/\s//g; print "$out\n"}' NM_001271626.3.gbk
```

### Split a Markdown file into separate files based on a heading

Given this input:

```text
# Answer
some answer

# Answer
some more
text

# Answer
the
final
answer
```

The command below creates three files called `output-1.md`, `output-2.md`, and `output-3.md` containing the text from the first, second, and third `# Answer` sections, respectively.

For example, `output-1.md` contains:

```text
# Answer
some answer

```

```bash
perl -ne 'BEGIN {$section_name=shift; $output_file_prefix=shift; undef $/;} $section_count = 1; $file = $ARGV; while (m/(\Q$section_name\E.*?)((?=\Q$section_name\E)|$)/sg) {open(w, ">", "$output_file_prefix-$section_count.md"); print w "$1"; $section_count++; close(w)}' "# Answer" "output" input.md
```

### Format Perl code

The following uses `perltidy` to reformat the code in `testfile.pl` and will create a file called `testfile.pl.tdy`.

```bash
perltidy testfile.pl
```

## R

### Compare two data sets to find differences

In this example SNP location information assigned by two different algorithms is compared using the `compareDF` package.

The two data sets (position information) generated by the differing algorithms are read from files. SNP name is used to match rows across the data sets, and then the 'chromosome' and 'position' are compared.

SNP records for which 'chromosome' or 'position' differ between the data sets are written to an Excel file.

```r
library(compareDF)

algorithm1_results <- read.csv2("algorithm1_results.csv", comment.char = "#", sep = ",", header = TRUE)
algorithm1_results_snp_positions <- data.frame(algorithm1_results$marker_name, algorithm1_results$chromosome, algorithm1_results$position)
colnames(algorithm1_results_snp_positions) <- c('SNP_name', 'chromosome', 'position')

algorithm2_results <- read.csv2("algorithm2_results.csv", comment.char = "#", sep = ",", header = TRUE)
algorithm2_results_snp_positions <- data.frame(algorithm2_results$SNP_name, algorithm2_results$chromosome, algorithm2_results$position)
colnames(algorithm2_results_snp_positions) <- c('SNP_name', 'chromosome', 'position')

#compare positions between data sets, matching based on SNP_name
ctable = compare_df(algorithm1_results_snp_positions, algorithm2_results_snp_positions, c('SNP_name'))
output_file <- paste("positions", "algorithm1_results", "vs", paste("algorithm2_results", ".xlsx", sep=""), sep="_")
create_output_table(ctable, output_type="xlsx", file_name=output_file, limit=1000000)
```

### Visualize the degree of overlap among gene sets

In this example, an UpSet plot is used to visualize the overlap among all combinations of gene lists in the `gene_lists` directory. In this directory each list is given as a separate `.txt` file, with a single header row and one gene name or ID per row, for example:

```text
Gene name or identifier
ENSG00000264954.2
ENSG00000224383.8
CCDS54157.1.1
```

The UpSet plot is generated using the `UpSetR` package:

```r
library(UpSetR)

setwd('/path/to/gene_lists')
filenames <- list.files(pattern = "*.txt", full.names = FALSE)

#create list of character vectors, each named after the source filename
#assumes each file has single header line (skip = 1)
gl <- sapply(filenames, scan, character(), sep="\n", skip = 1, USE.NAMES = TRUE)

#remove underscores from vector names
names(gl) <- gsub(x = names(gl), pattern = "_", replacement = " ")

#remove file extension from vector names
names(gl) <- gsub(x = names(gl), pattern = "\\..+?$", replacement = "")

upset(fromList(gl), nsets = length(gl), order.by = "freq")
```

The resulting plot displays the number of items shared among all possible combinations of overlapping sets in an easy-to-interpret and parse manner (unlike a traditional Venn diagram).

### Cluster gene lists based on overlap and identify shared genes

In this example a heatmap is used to visualize gene presence and absence for all gene lists in the `gene_lists` directory. In this directory each list is given as a separate `.txt` file, with a single header row and one gene name or ID per row, for example:

```text
Gene name or identifier
ENSG00000264954.2
ENSG00000224383.8
CCDS54157.1.1
```

The following uses the `purrr` and `RVenn` packages:

```r
library(purrr)
library(RVenn)

setwd('/path/to/gene_lists')
filenames <- list.files(pattern = "*.txt", full.names = FALSE)

#create list of character vectors, each named after the source filename
#assumes each file has single header line (skip = 1)
gl <- sapply(filenames, scan, character(), sep="\n", skip = 1, USE.NAMES = TRUE)

#remove underscores from vector names
names(gl) <- gsub(x = names(gl), pattern = "_", replacement = " ")

#remove file extension from vector names
names(gl) <- gsub(x = names(gl), pattern = "\\..+?$", replacement = "")

venn = Venn(gl)
setmap(venn, element_fontsize = 4, set_fontsize = 4)
```

The resulting heatmap displays genes and gene lists as rows and columns, respectively. The columns and rows are arranged so that genes and gene lists with similar presence / absence patterns are grouped together.

### Transpose a data frame

To convert this:

```text
snp                 sample1  sample2
ABCA12              AA       AA
APAF1               CC       CC
ARS-BFGL-BAC-10172  GG       AG
ARS-BFGL-BAC-1020   AA       GG
```

To this:

```text
snp      ABCA12  APAF1  ARS-BFGL-BAC-10172  ARS-BFGL-BAC-1020
sample1  AA      CC     GG                  AA
sample2  AA      CC     AG                  GG
```

Use this:

```r
library(dplyr)
library(tidyr)
library(janitor)

#prepare sample data frame
snp <- c('ABCA12', 'APAF1', 'ARS-BFGL-BAC-10172', 'ARS-BFGL-BAC-1020')
sample1 <- c('AA', 'CC', 'GG', 'AA')
sample2 <- c('AA', 'CC', 'AG', 'GG')
genotypes <- data.frame(snp, sample1, sample2)

genotypes %>%
  tibble::rownames_to_column() %>%
  pivot_longer(-rowname) %>%
  pivot_wider(names_from=rowname, values_from=value) ->
  genotypes_transposed

genotypes_transposed %>%
  row_to_names(row_number = 1) ->
  genotypes_transposed_with_column_names

write.table(genotypes_transposed_with_column_names, file = "genotypes_transposed_with_column_names.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " ")
```

Or use this:

```r
snp <- c('ABCA12', 'APAF1', 'ARS-BFGL-BAC-10172', 'ARS-BFGL-BAC-1020')
sample1 <- c('AA', 'CC', 'GG', 'AA')
sample2 <- c('AA', 'CC', 'AG', 'GG')
genotypes <- data.frame(snp, sample1, sample2)

genotypes_transposed <- data.frame(t(genotypes[-1]))
colnames(genotypes_transposed) <- genotypes[, 1]

write.table(genotypes_transposed, file = "genotypes_transposed_with_column_names.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " ")
```

Or use [datamash](#datamash).

### Split two-allele genotypes into two columns

To convert this:

```text
sample   ABCA12  APAF1  ARS-BFGL-BAC-10172  ARS-BFGL-BAC-1020
sample1  AA      CC     GG                  AA
sample2  AA      CC     AG                  GG
```

To this:

```text
sample   ABCA12_1  ABCA12_2  APAF1_1  APAF1_2  ARS-BFGL-BAC-10172_1  ARS-BFGL-BAC-10172_2  ARS-BFGL-BAC-1020_1  ARS-BFGL-BAC-1020_2
sample1  A         A         C        C        G                     G                     A                    A
sample2  A         A         C        C        A                     G                     G                    G
```

Use this:

```r
library(dplyr)
library(tidyr)

#prepare sample data frame
sample <- c('sample1', 'sample2')
ABCA12 <- c('AA', 'AA')
APAF1 <- c('CC', 'CC')
`ARS-BFGL-BAC-10172` <- c('GG', 'AG')
`ARS-BFGL-BAC-1020` <- c('AA', 'GG')
genotypes <- tibble(sample, ABCA12, APAF1, `ARS-BFGL-BAC-10172`, `ARS-BFGL-BAC-1020`)

#[-1] prevents first column from being split
for(column_name in names(genotypes)[-1]){
  genotypes %>%
    separate(column_name, c(paste(column_name, "1", sep = "_"), paste(column_name, "2", sep = "_")), sep = 1) ->
    genotypes
}

write.table(genotypes, file = "genotypes_split.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " ")
```

### Split multi-locus genotype strings into multiple columns and decode genotypes

To convert this:

```text
sample          alleles
HOCANF12689774  1000112112
HOCANF12689787  2011112012
HOCANF12689790  1011002122
```

To this:

```text
sample          Hapmap43437-BTA-101873  ARS-BFGL-NGS-16466  Hapmap34944-BES1_Contig627_1906  ARS-BFGL-NGS-98142  Hapmap53946-rs29015852  ARS-BFGL-NGS-114208  ARS-BFGL-NGS-66449  ARS-BFGL-BAC-32770  ARS-BFGL-NGS-65067  ARS-BFGL-BAC-32722
HOCANF12689774  AB                      BB                  BB                               BB                  AB                      AB                   AA                  AB                  AB                  AA
HOCANF12689787  AA                      BB                  AB                               AB                  AB                      AB                   AA                  BB                  AB                  AA
HOCANF12689790  AB                      BB                  AB                               AB                  BB                      BB                   AA                  AB                  AA                  AA
```

Use this:

```r
library(dplyr)
library(tidyr)

#prepare sample data frame
sample <- c('HOCANF12689774', 'HOCANF12689787', 'HOCANF12689790')
alleles <- c('1000112112', '2011112012', '1011002122')
genotypes <- data.frame(sample, alleles)

#vector of snp names
snps <- c('Hapmap43437-BTA-101873', 'ARS-BFGL-NGS-16466', 'Hapmap34944-BES1_Contig627_1906', 'ARS-BFGL-NGS-98142', 'Hapmap53946-rs29015852', 'ARS-BFGL-NGS-114208', 'ARS-BFGL-NGS-66449', 'ARS-BFGL-BAC-32770', 'ARS-BFGL-NGS-65067', 'ARS-BFGL-BAC-32722')

#create a new column for each digit in alleles
genotypes_one_column_per_snp <- separate(genotypes, col = alleles, sep = "(?<=\\d)", into = snps)

#convert each digit to textual representation
genotypes_one_column_per_snp_decoded = reshape2::dcast(
  dplyr::mutate(
    reshape2::melt(genotypes_one_column_per_snp, id.var = "sample"),
    value=plyr::mapvalues(
      value, c("0", "1", "2"), c("BB", "AB", "AA"))
  ),sample~variable)

write.table(genotypes_one_column_per_snp_decoded, file = "genotypes_one_column_per_snp_decoded.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " ")
```

### Change values in a column based on values in another column

```r
library(tidyverse)

tb <- tribble(
  ~chr, ~pos, ~sample1, ~sample2, ~A,   ~B,
  #----|-----|---------|---------|-----|-----
  "1",  2,    "AA",     "AB",     "REF","ALT",
  "1",  12,   "BB",     "AA",     "ALT","REF",
  "1",  12,   ".",      "AA",     "ALT","REF",
)

#get names of sample columns
names(tb) %>%
  str_subset(pattern = "^sample") ->
  columns_to_decode

#convert genotypes to values from A and B columns
tb %>%
  mutate_at(
    vars(one_of(columns_to_decode)),
    list(~case_when(
      . == "AA" & A == "REF" ~ "0/0",
      . == "AA" & A == "ALT" ~ "1/1",
      . == "BB" & B == "REF" ~ "0/0",
      . == "BB" & B == "ALT" ~ "1/1",
      . == "AB" ~ "0/1",
      TRUE ~ './.'))) ->
  tb_converted

print(tb_converted)

## A tibble: 3 x 6
#  chr     pos sample1 sample2 A     B
#  <chr> <dbl> <chr>   <chr>   <chr> <chr>
#1 1         2 0/0     0/1     REF   ALT
#2 1        12 0/0     1/1     ALT   REF
#3 1        12 ./.     1/1     ALT   REF
```

A different approach, using nested ifelse():

```r
library(tidyverse)

tb <- tribble(
  ~chr, ~pos, ~sample1, ~sample2, ~A,   ~B,
  #----|-----|---------|---------|-----|-----
  "1",  2,    "AA",     "AB",     "REF","ALT",
  "1",  12,   "BB",     "AA",     "ALT","REF",
  "1",  12,   ".",      "AA",     "ALT","REF",
)

#get names of sample columns
names(tb) %>%
  str_subset(pattern = "^sample") ->
  columns_to_decode

#function to convert genotypes to values from A and B columns
convert_genotypes <- function (df, col) {
  df[[col]] <-
    ifelse((df[[col]] == "AA") & (df[["A"]] == "REF"), "0/0",
           ifelse((df[[col]] == "AA") & (df[["A"]] == "ALT"), "1/1",
                  ifelse((df[[col]] == "BB") & (df[["B"]] == "REF"), "0/0",
                         ifelse((df[[col]] == "BB") & (df[["B"]] == "ALT"), "1/1",
                                ifelse(df[[col]] == "AB", "0/1", "./.")))))
  return(df)
}

for (sample in columns_to_decode) {
  tb <- convert_genotypes(tb, sample)
}

print(tb)

## A tibble: 3 x 6
#  chr     pos sample1 sample2 A     B
#  <chr> <dbl> <chr>   <chr>   <chr> <chr>
#1 1         2 0/0     0/1     REF   ALT
#2 1        12 0/0     1/1     ALT   REF
#3 1        12 ./.     1/1     ALT   REF
```

Yet another approach, by passing column names to a function that uses mutate() and case_when():

```r
library(tidyverse)

tb <- tribble(
  ~chr, ~pos, ~sample1_1, ~sample1_2, ~FORWARD_A, ~FORWARD_B,
  #----|-----|-----------|-----------|-----------|-----------
  "1",  2,    "G",        "G",        "G",        "T",
  "1",  12,   "A",        "A",        "C",        "A",
  "1",  12,   ".",        "G",        "T",        "G",
)

#get names of sample columns
names(tb) %>%
  str_subset(pattern = "^sample") ->
  columns_to_decode

#function to convert genotypes to values from A and B columns
convert <- function(df, col, format) {

  A <- paste(format, "A", sep = "_")
  B <- paste(format, "B", sep = "_")

  object <- df %>%
    mutate(
      !!col := case_when(
        eval(parse(text = col)) == eval(parse(text = A)) ~ "A",
        eval(parse(text = col)) == eval(parse(text = B)) ~ "B",
        TRUE ~ "."
      )
    )
  return(object)
}

for (sample in columns_to_decode) {
  tb <- convert(tb, sample, "FORWARD")
}

print(tb)

## A tibble: 3 x 6
# chr     pos sample1_1 sample1_2 FORWARD_A FORWARD_B
# <chr> <dbl> <chr>     <chr>     <chr>     <chr>
#1 1         2 A         A         G         T
#2 1        12 B         B         C         A
#3 1        12 .         B         T         G
```

### Add comment lines to output

```r
library(dplyr)

#vcf is a tibble
#add comment character to start of first column name
vcf <- rename(vcf, `#CHROM` = CHROM)
#write out comment line and then column names
writeLines(c("##fileformat=VCFv4.2", paste(names(vcf), collapse = "\t")), con = "genotypes.vcf")
write.table(vcf, file = "genotypes.vcf", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t", append=TRUE)
```

### Filter and sort rows

```r
library(tidyverse)

#vcf is tibble

#remove rows where POS is NA
vcf %>% drop_na(POS) ->
  vcf

#keep rows with single base in REF and ALT
vcf %>%
  filter(str_detect(REF, "^[GATCN]$")) %>%
  filter(str_detect(ALT, "^[GATCN]$")) ->
  vcf

#sort by chromosome then position
vcf %>%
  arrange(CHROM, POS) ->
  vcf
```

### Add columns from one tibble to another

```r
library(purrr)
library(tibble)

#vcf and genotypes are tibbles
#add columns from genotypes to vcf
for (column in names(genotypes)) {
  vcf %>%
    add_column(!!(column) := genotypes[[column]]) ->
    vcf
}
```

### Combine multiple input files

In the following example the rows from multiple files are combined. Additional columns are used to track the source of each row. The resulting long data is converted into two different wide formats and written to a single Excel file as separate worksheets.

The input files are structured as follows, with additional columns not shown:

```text
fusion_name, junction_read_count
Wfdc3--AABR07030443.2, 16
Wfdc3--AABR07030443.2, 2
Pdk4--Pdk2, 10
AABR07059159.1--Actr3b, 8
```

Each file is named after the source sample. For example `180_S1_R1_001.fusion_candidates.preliminary` is the data for sample `180_S1_R1_001`.

```r
library(tidyverse)
library(janitor)
library(xlsx)

#Directory containing the files
input_dir <- "data/input/"
#Pattern to use when identifying input files
file_pattern <- "*.preliminary"

#Create tibble of input files (full path and file name in separate columns)
input_files <-
  tibble(
    full_path = list.files(input_dir, pattern = file_pattern, full.names = TRUE),
    file_name = list.files(input_dir, pattern = file_pattern, full.names = FALSE)
  )

#Function to add source info to each row
read_tsv_and_add_source <- function(file_name, full_path) {
  read_tsv(full_path) %>%
    clean_names %>%
    mutate(file_name = file_name) %>%
    mutate(full_path = full_path) %>%
    #convert '180_S1_R1_001.fusion_candidates.preliminary' to '180_S1_R1_001'
    mutate(sample = str_split(file_name, ".fusion", simplify = TRUE)[[1]])
}

#Read all files into a single tibble
input_files %>%
  rowwise() %>%
  do(., read_tsv_and_add_source(.$file_name, .$full_path)) ->
  combined_data_with_source

#Group data by fusion_name and sample and for each group calculate the sum of
#junction_read_count
combined_data_with_source %>%
  group_by(fusion_name, sample) %>%
  summarise(fusion_count = sum(junction_read_count), .groups = NULL) ->
  counts_per_fusion

#Filter by fusion_name, keeping rows where fusion_name consists of two MT
#genes, for example Mt-co1--Mt-nd2
counts_per_fusion %>%
  filter(str_detect(fusion_name, "^Mt-")) %>%
  filter(str_detect(fusion_name, "--Mt-")) ->
  counts_per_MT_fusion

#Convert the data from long to wide, with values of sample becoming columns
counts_per_MT_fusion %>%
  spread(sample, fusion_count) %>%
  replace(is.na(.), 0) ->
  samples_as_columns

#Convert the data from long to wide, with values of fusion_name becoming
#columns
counts_per_MT_fusion %>%
  spread(fusion_name, fusion_count) %>%
  replace(is.na(.), 0) ->
  fusions_as_columns

#Write the data to an Excel file
write.xlsx(
  as.data.frame(samples_as_columns),
  "output.xlsx",
  sheetName = "Samples as columns",
  col.names = TRUE,
  row.names = FALSE,
  append = TRUE
)

write.xlsx(
  as.data.frame(fusions_as_columns),
  "output.xlsx",
  sheetName = "Fusions as columns",
  col.names = TRUE,
  row.names = FALSE,
  append = TRUE
)
```

## rsync

### Sync a directory on local system

```bash
rsync -avzh source_directory destination_directory
```

To sync only the contents of `source_directory` and not the directory itself use:

```bash
rsync -avzh source_directory/ destination_directory
```

### Sync a directory to a remote system

```bash
rsync -avzh source_directory user@192.168.0.101:~/destination
```

### Sync a directory from a remote system

```bash
rsync -avzh user@192.168.0.101:~/source_directory destination
```

## Singularity

### Create an image from a container stored in Docker Hub

In this example a container that can be downloaded from Docker Hub using `docker pull pstothard/cgview` is used to generate a Singularity container:

```bash
singularity build cgview.sif docker://pstothard/cgview
```

To run a command in this container, use something like the following:

```bash
singularity exec -B /scratch cgview.sif java -jar /usr/bin/cgview.jar --help
```

The `-B` is used to provide the container with access to directories.

## sbatch

### Count lines in compressed fastq files

In this example the number of lines in several `.fastq.gz` files is quickly determined by submitting jobs to Slurm using sbatch.

The naming scheme of the `.fastq.gz` files is as follows (the sample name is in the file name, for example `DG15B032198-1`):

```text
HI.5173.001.NEBNext_Index_12.DG15B032198-1_R1.fastq.gz
HI.5173.001.NEBNext_Index_12.DG15B032198-1_R2.fastq.gz
HI.5173.002.NEBNext_Index_12.DG15B032198-1_R1.fastq.gz
HI.5173.002.NEBNext_Index_12.DG15B032198-1_R2.fastq.gz
HI.5173.003.NEBNext_Index_14.DG15B032179-1_R1.fastq.gz
HI.5173.003.NEBNext_Index_14.DG15B032179-1_R2.fastq.gz
HI.5173.004.NEBNext_Index_14.DG15B032179-1_R1.fastq.gz
HI.5173.004.NEBNext_Index_14.DG15B032179-1_R2.fastq.gz
```

First create a sbatch script called `fastq.gz.lines.sbatch` to run zcat and wc (used to count lines in a compressed file):

```bash
#!/bin/bash
#SBATCH --account=def-someuser
#SBATCH --ntasks=1
#SBATCH --mem=1000M
#SBATCH --time=0-0:30

LINES="`zcat $1 | wc -l`"
echo "$1 $LINES"
```

Then use the following commands to submit a job for each `.fastq.gz` file:

```bash
files=$(find . -name "*fastq.gz" -printf '%P\n')
for f in $files
do
  command="sbatch -o ${f}.out -e ${f}.err fastq.gz.lines.sbatch $f"
  echo "Submitting a job using the following command:"
  echo "$command"
  eval "$command"
  sleep 1
done
```

Each job should create two files for each input file, for example:

```text
HI.5173.001.NEBNext_Index_12.DG15B032198-1_R1.fastq.gz.out
HI.5173.001.NEBNext_Index_12.DG15B032198-1_R1.fastq.gz.err
```

The `.out` files will contain the name of the input file and the number of lines, for example:

```text
HI.5173.001.NEBNext_Index_12.DG15B032198-1_R1.fastq.gz 229623444
```

To quickly check the `.err` files:

```bash
cat *.err | more
```

To add up the line counts for each sample using the R1 and R2 files:

```bash
cat *R1.fastq.gz.out | perl -nl -e 's/^.+?([^_\.]{10,})\S*\s*(\d+)/$1/g;' -e 'print $1 . "\t" . $2' | sort | awk '{a[$1]+=$2} END {for(i in a) print i,a[i]}' > line_counts_per_sample_R1.tab

cat *R2.fastq.gz.out | perl -nl -e 's/^.+?([^_\.]{10,})\S*\s*(\d+)/$1/g;' -e 'print $1 . "\t" . $2' | sort | awk '{a[$1]+=$2} END {for(i in a) print i,a[i]}' > line_counts_per_sample_R2.tab
```

To compare the counts obtained using the R1 and R2 files:

```bash
diff line_counts_per_sample_R1.tab line_counts_per_sample_R2.tab
```

## sed

### Add a header line to a file with sed

```bash
sed $'1s/^/my header text\\\n&/' input
```

### Edit the header line with sed

In this example `pos,reads_pp` is changed to `pos,reads` in the first line of the file. The `-i` is used to edit the file in place:

```bash
sed -i '' -e "1,1s/pos,reads_pp/pos,reads/" input.csv
```

### Print a specific line of a file

In this example line `26404` is printed:

```bash
sed -n "26404p" input.txt
```

### Change filenames using a regular expression

In this example `chr30` is replaced with `chrX`:

```bash
for f in *.fasta; do new=`echo $f | sed 's/chr30/chrX/'`; mv $f $new; done
```

### Search and replace on lines

The `sed` command can be used to find and replace text and supports its own [regular expression syntax](https://www.gnu.org/software/sed/manual/html_node/Regular-Expressions.html).

The following replaces all occurrences of `Red Angus` with `Angus`:

```bash
sed 's/Red Angus/Angus/g' sequenced_samples.csv
```

In the above the `g` indicates that all matches on a line should be replaced.

The following restricts the processing to lines 302 to 305:

```bash
sed '302,305 s/Red Angus/Angus/g' sequenced_samples.csv
```

The following replaces the first occurrence of `LIM` on each line with `---`:

```bash
sed 's/LIM/---/1' sequenced_samples.csv
```

The following adds underscores around the first three characters of each line:

```bash
sed 's/^\([a-zA-Z0-9]\{3\}\)/_\1_/g' sequenced_samples.csv
```

The above uses `\(` and `\)` to "remember" matching characters and then uses `\1` to use the stored matching text in the replacement portion. Thus the replacement text becomes `_` followed by the actual matched text, followed by `_`.

The different parts of the search part of the command have the following meanings:

* `^` match the start of a line
* `[a-zA-Z0-9]\{3\}` match three letters or numbers

### Delete lines

The following deletes the first line:

```bash
sed '1d' sequenced_samples.csv
```

The following deletes from line 500 to the end of the file (represented by `$`):

```bash
sed '500,$d' sequenced_samples.csv
```

## Share data with project group members

```bash
cd projects/some_project
chmod g+x my_dir
cd my_dir
mkdir shared_dir
chmod g+x shared_dir
chmod +t shared_dir
```

## Slurm

### View statistics related to the efficiency of resource usage of a completed job

```bash
seff <jobid>
```

### View jobs

```bash
squeue -u <username>
```

### View running jobs

```bash
squeue -u <username> -t RUNNING
```

### View pending jobs

```bash
squeue -u <username> -t PENDING
```

### View detailed information for a specific job

```bash
scontrol show job -dd <jobid>
```

### View accounting information for completed jobs

```bash
sacct -s CD --format=JobID,JobName,MaxRSS,ReqMem,Elapsed,End,State,NodeList
```

### Cancel a job

```bash
scancel <jobid>
```

### Cancel all jobs

```bash
scancel -u <username>
```

### Start an interactive session

```bash
salloc --time=2:0:0 --ntasks=2 --account=def-someuser --mem-per-cpu=8000M --mail-type=ALL --mail-user=your.email@example.com
```

## sort

### Alphabetical sort

```bash
sort -t, input.csv
```

### Specify the sort field

The `-k` option of `sort` is used to indicate that certain fields are to be used for sorting (the sorting is still applied to the entire line).

The following sorts only using column 5 (`5,5` means "starting at column 5 and ending at column 5"), numerically, from smallest to largest:

```bash
(head -n 1 sequenced_samples.csv && \
tail -n +2 sequenced_samples.csv | sort -t, -k5,5n)
```

The following sorts using column 2 to the end of line, alphabetically:

```bash
(head -n 1 sequenced_samples.csv && sort -t, -k2 \
<(tail -n +2 sequenced_samples.csv))
```

### Use multiple sort fields

The command below sorts using column 2, breaking ties using the value in column 5 (sorting numerically from largest to smallest). The input file has a single header row, which is excluded from the sort.

```bash
cat sequenced_samples.csv | awk \
'NR<2{print $0; next}{print $0| "sort -t, -k2,2 -k5,5nr"}'
```

### Sort a file with a header row

In the following examples a file with a single header line is sorted:

```bash
(head -n 1 sequenced_samples.csv && sort -t, <(tail -n +2 sequenced_samples.csv))
```

```bash
(head -n 1 sequenced_samples.csv && tail -n +2 sequenced_samples.csv | sort -t,)
```

```bash
cat sequenced_samples.csv | awk 'NR<2{print $0; next}{print $0| "sort -t','"}'
```

The above command can be modified to sort by the second column, numerically from smallest to largest:

```bash
cat sequenced_samples.csv | awk 'NR<2{print $0; next}{print $0| "sort -t',' -k2,2n"}'
```

### Sort by chromosome

This example uses the following input from the text file `input.txt`:

```text
chrX Other
chrY Information
MT !
chr1 Some
chr2 Data
chr3 Or    
chr3 Annotation
chr10 Or
chr21 Any
```

To sort by chromosome use the following:

```bash
sort -k1,1 -V -f -s input.txt
```

The above generates the following output:

```text
chr1 Some
chr2 Data
chr3 Or    
chr3 Annotation
chr10 Or
chr21 Any
chrX Other
chrY Information
MT !
```

## tmux

### Start a tmux session

```bash
tmux new -s session_name
```

### Detach a tmux session

`Ctrl`-`b` then `d`.

### List tmux sessions

```bash
tmux ls
```

### Join a tmux session

```bash
tmux attach -t session_name
```

### Create a multi-pane tmux session

The following creates three panes: one large one at the top and two smaller ones at the bottom:

```bash
tmux new-session -s multiple \; \
split-window -v -p 25 \; \
split-window -h -p 50 \; \
select-pane -t 0 \;
```

The following also creates a three-pane tmux session, and launches `vim` in the largest pane:

```bash
tmux new-session -s multiple \; \
send-keys 'vim' C-m \; \
split-window -v -p 25 \; \
split-window -h -p 50 \; \
select-pane -t 0 \;
```

The following creates six equally sized panes:

```bash
tmux new-session -s multiple \; \
split-window -h \; \
split-window -v -p 66 \; \
split-window -v \; \
select-pane -t 0 \; \
split-window -v -p 66 \; \
split-window -v \;
```

### Scroll in tmux

`Ctrl`-`b` then `PgUp`. Press `q` to return to normal mode.

Alternatively you can use `Ctrl`-`b` then `[` and then navigation keys like `Up Arrow` or `PgDn`. Press `q` to return to normal mode.

### Navigate between tmux panes

The following work with [my configuration](https://github.com/paulstothard/dotfiles/blob/master/.tmux.conf):

`Ctrl`-`j` moves up.
`Ctrl`-`k` moves down.
`Ctrl`-`h` moves left.
`Ctrl`-`l` moves right.

### Kill a tmux session

```bash
tmux kill-session -t session_name
```

## tr

### Translate characters

The following changes all `S` and `H` characters to `s` and `h`, respectively:

```bash
cat sequenced_samples.csv | tr "SH" "sh"
```

The sets of characters provided to `tr` can be the actual characters (as in the above example) or they can be one of several special characters or character groups (use `tr --help` to learn more).

For example, the following `tr` command changes uppercase text to lowercase text:

```bash
cat sequenced_samples.csv | tr "[:upper:]" "[:lower:]"
```

To change all digits to `X` use `tr` as in the following example:

```bash
cat sequenced_samples.csv | tr "[:digit:]" "X"
```

The following `tr` command changes all characters that aren't `G`, `A`, `T`, or `C` (ignoring case) or whitespace to `N`:

```bash
echo "garsxnT" | tr -c "GATCgatc[:space:]" "N"
```

The above command uses the complement option (`-c`) to convert the first set to everything that isn't in the set.

### Delete characters

Use the `-d` option with `tr` to delete characters.

The following `tr` command deletes all digits:

```bash
cat sequenced_samples.csv | tr -d "[:digit:]"
```

To delete everything except for digits using the following `tr` command:

```bash
cat sequenced_samples.csv | tr -c -d "[:digit:]"
```

### Squeeze characters

"Squeezing" is used here to mean "removing repeated instances". Typically this would be used to remove duplicated spaces or punctuation.

The following illustrates the removal extra commas by using `tr` with the `-s` option:

```bash
echo "a,b,,c,,,d" | tr -s ","
```

## VCF files

### Extract variants from a region of interest

Note that if the VCF file is gzip compressed (i.e. has a `.gz` extension), use `--gzvcf` instead of `--vcf`.

```bash
vcftools --vcf Chr5.vcf --out Chr5_filtered --chr 5 --from-bp 1 --to-bp 100000 --recode --recode-INFO-all
```

### Extract variants from multiple regions of interest

```bash
bgzip Chr5.vcf
tabix -p vcf Chr5.vcf.gz
bcftools view -r 5:1-10000,5:200000-210000 -o output.vcf Chr5.vcf.gz
```

### Extract variants from multiple regions of interest described in a file

In this example the regions of interest are stored in a text file called `regions.txt`. Each line describes the chromosome, start, and end of a region:

```text
3 62148416 62200719
4 54643953 54720351
4 63732381 63795159
5 10163746 10218801
5 10784272 10841310
```

```bash
bgzip input.vcf
tabix -p vcf input.vcf.gz
tabix --print-header -R regions.txt input.vcf.gz > regions_of_interest.vcf
```

### Assess sex by calculating X-chromosome heterozygosity

```bash
#vcftools adds .het extension to output automatically, so output becomes output.vcf.het
vcftools --vcf input.vcf --chr X --het --out output.vcf
#add column to vcftools output reporting percent heterozygous genotypes
awk -F$'\t' 'BEGIN{OFS="\t"}; {if(NR==1){print $0,"Percent HET"} else {print $0, ($4 - $2) / $4 * 100}}' output.vcf.het > output.vcf.percent_het
```

### Assess sex using plink

Use [--check-sex](https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex) in `plink`.

First convert the VCF file to a plink binary fileset:

```bash
species=dog
pseudoX_start=340475
pseudoX_end=6642728
mkdir plink_bed
plink --vcf input.vcf --make-bed --allow-extra-chr \
--$species --split-x $pseudoX_start $pseudoX_end \
--out plink_bed/input
```

Next edit the `.fam` file in `plink_bed` to describe family relationships, sex, and case / control status.

The following is from the `plink` documentation:

> A text file with no header line, and one line per sample with the following six fields:
>
> Family ID ('FID')
> Within-family ID ('IID'; cannot be '0')
> Within-family ID of father ('0' if father isn't in dataset)
> Within-family ID of mother ('0' if mother isn't in dataset)
> Sex code ('1' = male, '2' = female, '0' = unknown)
> Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

Use `--check-sex`:

```bash
mkdir plink_check_sex
plink --bfile plink_bed/input --check-sex --allow-extra-chr \
--$species --out plink_check_sex/input
cat plink_check_sex/input.sexcheck
```

In the output `1` = male, `2` = female, `other` = unknown.

### Count genotypes

Use [--geno-counts](https://www.cog-genomics.org/plink/2.0/basic_stats#geno_counts) in `plink2`:

```bash
plink2 --vcf snps.vcf --geno-counts 'cols=chrom,pos,ref,alt,homref,refalt,altxy,hapref,hapalt,missing' --allow-extra-chr --chr-set 95
```

Or use `bcftools`. First split multiallelic sites into biallelic records using `bcftools norm`:

```bash
bgzip snps.vcf
tabix -fp vcf snps.vcf.gz
bcftools norm -m-any snps.vcf.gz -Oz > snps.norm.vcf.gz
```

Then generate a CSV report providing genotype counts for each variant:

```bash
INPUT=snps.norm.vcf.gz
OUTPUT=snps.norm.vcf.genotype.counts.csv

paste \
<(bcftools view "$INPUT" | awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
<(bcftools query -f '[\t%SAMPLE=%GT]\n' "$INPUT" | awk 'BEGIN {print "nHomRef"} {print gsub(/0\|0|0\/0/, "")}') \
<(bcftools query -f '[\t%SAMPLE=%GT]\n' "$INPUT" | awk 'BEGIN {print "nHet"} {print gsub(/0\|1|1\|0|0\/1|1\/0/, "")}') \
<(bcftools query -f '[\t%SAMPLE=%GT]\n' "$INPUT" | awk 'BEGIN {print "nHomAlt"} {print gsub(/1\|1|1\/1/, "")}') \
| sed 's/,\t/\t/g' | sed 's/,$//g' > "$OUTPUT"
```

In some cases a good alternative is to use `bcftools` to calculate additional `INFO` tags from the genotypes. The following adds `AC_Hom` and `AC_Het` tags to the VCF file.

```bash
bgzip source.vcf
tabix -fp vcf source.vcf.gz
bcftools +fill-tags source.vcf.gz -Oz -o source.additional-fields.vcf.gz -- -t AC_Hom,AC_Het
```

### Determine the proportion of missing genotypes in each sample

```bash
INPUT=SNPs.vcf
paste \
<(bcftools query -f '[%SAMPLE\t]\n' "$INPUT" | head -1 | tr '\t' '\n') \
<(bcftools query -f '[%GT\t]\n' "$INPUT" | awk -v OFS="\t" '{for (i=1;i<=NF;i++) if ($i == "./.") sum[i]+=1 } END {for (i in sum) print i, sum[i] / NR }' | sort -k1,1n | cut -f 2)
```

The above produces output like:

```text
sample1 0.956705
sample2 0.124076
sample3 0.0281456
```

### Determine genotype concordance

Use `bcftools`:

```bash
bgzip A.vcf
tabix -fp vcf A.vcf.gz

bgzip B.vcf
tabix -fp vcf B.vcf.gz

bcftools gtcheck -g A.vcf.gz B.vcf.gz
```

Or use SnpSift:

```bash
java -jar SnpSift.jar concordance A.vcf B.vcf
```

Or use Picard. In this example the comparison is limited to the sample named `HOCANF12689774`:

```bash
java -jar ~/bin/picard.jar GenotypeConcordance \
CALL_VCF=A.vcf \
CALL_SAMPLE=HOCANF12689774 \
O=gc_concordance.vcf \
TRUTH_VCF=B.vcf \
TRUTH_SAMPLE=HOCANF12689774
```

### Count Mendelian errors using plink

Use [--mendel](https://www.cog-genomics.org/plink/1.9/basic_stats#mendel) in `plink`.

First convert the VCF file to a plink binary fileset:

```bash
species=dog
pseudoX_start=340475
pseudoX_end=6642728
mkdir plink_bed
plink --vcf input.vcf --make-bed --allow-extra-chr \
--$species --split-x $pseudoX_start $pseudoX_end \
--out plink_bed/input
```

Next edit the `.fam` file in `plink_bed` to describe family relationships, sex, and case / control status.

The following is from the `plink` documentation:

> A text file with no header line, and one line per sample with the following six fields:
>
> Family ID ('FID')
> Within-family ID ('IID'; cannot be '0')
> Within-family ID of father ('0' if father isn't in dataset)
> Within-family ID of mother ('0' if mother isn't in dataset)
> Sex code ('1' = male, '2' = female, '0' = unknown)
> Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

Use `--mendel`:

```bash
mkdir plink_mendelian_errors
plink --bfile plink_bed/input --mendel --geno 0 --allow-extra-chr \
--$species --out plink_mendelian_errors/input
cat plink_mendelian_errors/input.imendel
```

### Find runs of homozygosity using plink

Use [--homozyg](https://www.cog-genomics.org/plink/1.9/ibd#homozyg) in `plink`.

First convert the VCF file to a plink binary fileset:

```bash
species=dog
pseudoX_start=340475
pseudoX_end=6642728
mkdir plink_bed
plink --vcf input.vcf --make-bed --allow-extra-chr \
--$species --split-x $pseudoX_start $pseudoX_end \
--out plink_bed/input
```

Next edit the `.fam` file in `plink_bed` to describe family relationships, sex, and case / control status.

The following is from the `plink` documentation:

> A text file with no header line, and one line per sample with the following six fields:
>
> Family ID ('FID')
> Within-family ID ('IID'; cannot be '0')
> Within-family ID of father ('0' if father isn't in dataset)
> Within-family ID of mother ('0' if mother isn't in dataset)
> Sex code ('1' = male, '2' = female, '0' = unknown)
> Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

Next, create a text file containing the Family ID and Within-family ID of each sample to be included in the analysis.

In this example the file `roh_samples.txt` consists of:

```text
family1 M-2
family4 M-9
family3 M-10
family2 M-15
family2 M-16
family1 M-3
family1 M-4
family1 M-5
```

Use `--homozyg`:

```bash
mkdir plink_roh
cd plink_roh
plink --bfile ../plink_bed/input \
--homozyg group extend \
--$species \
--keep roh_samples.txt \
--allow-extra-chr \
--homozyg-snp 50 \
--homozyg-kb 100 \
--homozyg-density 2 \
--homozyg-gap 100 \
--homozyg-window-snp 50 \
--homozyg-window-het 2 \
--homozyg-window-missing 10 \
--homozyg-window-threshold 0.05
```

The `plink.hom.overlap` file in `plink_roh` can be filtered to, for example, obtain regions of homozygosity that are found in `5` cases and `0` controls:

```bash
awk '{ if ($2 == "CON" && $4 == "5:0") print $5" "$8" "$9 }' plink.hom.overlap > plink.hom.overlap.filtered
```

Variants in the filtered regions of homozygosity can be extracted from the VCF file:

```bash
bgzip input.vcf
tabix -p vcf input.vcf.gz
tabix --print-header -R plink_roh/plink.hom.overlap.filtered \
input.vcf.gz > input.hom.vcf
```

### Add new samples

From the `bcftools` `merge` documentation:

> Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file. For example, when merging file A.vcf.gz containing samples S1, S2 and S3 and file B.vcf.gz containing samples S3 and S4, the output file will contain four samples named S1, S2, S3, 2:S3 and S4.

```bash
bcftools merge *.vcf.gz -Oz -o merged.vcf.gz
```

### Check relatedness between samples

Use the [KING algorithm](https://academic.oup.com/bioinformatics/article/26/22/2867/228512), which is implemented in `vcftools` and is accessed using the `--relatedness2`:

```bash
mkdir relatedness
vcftools --vcf input.vcf --not-chr X --max-missing-count 0 --relatedness2 \
--out relatedness/input.vcf
cat relatedness/input.vcf.relatedness2 | column -t
```

### Filter sites based on genotypes and other criteria

Use [bcftools view](https://samtools.github.io/bcftools/bcftools.html#view) to filter sites based on an [expression](https://samtools.github.io/bcftools/bcftools.html#expressions). Sites for which the expression is true can be kept using the `-i` option or excluded using the `-e` option.

In the following example sites are excluded if any sample is both homozygous for an alternate allele and has a genotype quality greater than 30:

```bash
bcftools view -e 'GT[*]="AA" & GQ[*]>30' SNPs.vcf
```

[SnpSift](http://pcingola.github.io/SnpEff/ss_introduction/) can also be used to filter VCF files based on sample genotypes.

The following approach can be used to exclude sites where any sample meets the following criteria: is homozygous and the genotype quality is greater than `30` and the genotype is not `0/0`. Worded another way, a site is kept if no sample exhibits a good-quality homozygous alternate genotype.

In this example there are 6 samples in the VCF file.

First generate the filter string:

```bash
FILTER=$(perl -e '@array = (); foreach(0..5) {push @array, "( isHom( GEN[$_] ) & GEN[$_].GQ > 30 & isVariant( GEN[$_] ) )";} print " !( " . join(" | ", @array) . " )\n";')
```

Examine the filter string:

```bash
echo $FILTER
```

This produces:

```text
!( ( isHom( GEN[0] ) & GEN[0].GQ > 30 & isVariant( GEN[0] ) ) | ( isHom( GEN[1] ) & GEN[1].GQ > 30 & isVariant( GEN[1] ) ) | ( isHom( GEN[2] ) & GEN[2].GQ > 30 & isVariant( GEN[2] ) ) | ( isHom( GEN[3] ) & GEN[3].GQ > 30 & isVariant( GEN[3] ) ) | ( isHom( GEN[4] ) & GEN[4].GQ > 30 & isVariant( GEN[4] ) ) | ( isHom( GEN[5] ) & GEN[5].GQ > 30 & isVariant( GEN[5] ) ) )
```

Use the filter string and `SnpSift.jar` to complete the filtering step (the `set +H` is used so that the `!` in the filter string doesn't activate Bash history expansion):

```bash
set +H
cat snps.vcf | java -jar SnpSift.jar filter "$FILTER" > snps_with_no_homozygous_alt_genotypes.vcf
set -H
```

### Add predicted consequences

Use [SnpEff](http://pcingola.github.io/SnpEff/se_introduction/) to predict variant effects.

List available pre-built databases for annotation:

```bash
java -jar snpEff.jar databases
```

Download a database:

```bash
java -jar snpEff.jar snpEff download -v CanFam3.1.99
```

Annotate a VCF file:

```bash
java -jar snpEff.jar -Xmx8g CanFam3.1.99 input.vcf > input.ann.vcf
```

Alternatively, use Ensembl VEP to predict variant effects. VEP can be used to annotate structural variants.

```bash
SPECIES=canis_lupus
ASSEMBLY=CanFam3.1
INPUT=input.vcf
VEP_CACHE=vep
OUTDIR=vep_annotated
WD="$(pwd)"
export PERL5LIB="$PERL5LIB:$WD/$VEP_CACHE"
export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$WD/$VEP_CACHE/htslib"

#create cache directory and output directory
mkdir -p $OUTDIR
mkdir -p $VEP_CACHE

#build cache
vep_install -a cfp -s $SPECIES -y $ASSEMBLY \
-c $VEP_CACHE \
-d $VEP_CACHE \
--PLUGINS all --CONVERT

#run vep
SPECIES=canis_lupus_familiaris
vep --cache --format vcf --vcf \
--dir_cache $VEP_CACHE \
--dir_plugins $VEP_CACHE/Plugins \
--input_file $INPUT \
--output_file $OUTDIR/$INPUT \
--species $SPECIES --assembly $ASSEMBLY \
--max_sv_size 1000000000 \
--force_overwrite \
--plugin Blosum62 --plugin Downstream --plugin Phenotypes --plugin TSSDistance --plugin miRNA \
--variant_class --sift b --nearest gene --overlaps --gene_phenotype --regulatory --protein \
--symbol --ccds --uniprot --biotype --domains --check_existing --no_check_alleles --pubmed \
--verbose
```

### Add variant IDs

Use [SnpSift](http://pcingola.github.io/SnpEff/ss_introduction/) `annotate` to add variant IDs.

In this example the file `canis_lupus_familiaris.sorted.vcf` has variant IDs to be transferred to `input.vcf`.

```bash
bgzip canis_lupus_familiaris.sorted.vcf
tabix -p vcf canis_lupus_familiaris.sorted.vcf.gz
java -jar SnpSift.jar annotate \
canis_lupus_familiaris.sorted.vcf.gz \
input.vcf > input.rsID.vcf
```

### Transfer annotations from another VCF

In this example some annotations are added to the file `source.vcf` (calculated from genotypes) and then annotations are transferred to the file `input.vcf`. The results are written to `output.vcf`.

The transfer of annotations is done using [vcfanno](https://github.com/brentp/vcfanno).

This procedure is useful for annotating variants in one file with information obtained, for example, from a much larger population of samples. The resulting tags can be used in downstream filtering operations.

First fill some additional annotation fields to `source.vcf`:

```bash
bgzip source.vcf
tabix -fp vcf source.vcf.gz
bcftools +fill-tags source.vcf.gz -Oz -o source.additional-fields.vcf.gz -- -t AC_Hom,AC_Het
```

Create a `config.toml` file for vcfanno:

```text
[[annotation]]
file="source.additional-fields.vcf.gz"
fields=["AF", "AC_Hom", "AC_Het"]
names=["source.AF", "source.AC_Hom", "source.AC_Het"]
ops=["self", "self", "self"]
```

Prepare the `input.vcf` file:

```bash
bgzip input.vcf
tabix -fp vcf input.vcf.gz
```

Create an index for `source.additional-fields.vcf.gz`:

```bash
tabix -fp vcf source.additional-fields.vcf.gz
```

Perform the annotation transfer:

```bash
vcfanno config.toml input.vcf.gz > out.vcf
```

The file `out.vcf` should now include `INFO` tags `source.AF`, `source.AC_Hom`, and `source.AC_Het` with values calculated from the samples in `source.vcf`.

### Filter variants based on predicted consequences

Use [SnpSift](http://pcingola.github.io/SnpEff/ss_introduction/) to filter VCF files that have been annotated using [SnpEff](http://pcingola.github.io/SnpEff/se_introduction/).

The following keeps variants that are predicted to have `HIGH` or `MODERATE` impacts:

```bash
cat input.ann.vcf | SnpSift filter "((ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE'))" > input.ann.high_or_moderate.vcf
```

### Keep only variants with a missing ID

```bash
awk -F $'\t' '$1 ~ /^#/ {print; next} $3~/^\./' input.vcf > input.noID.vcf
```

### Keep only variants with an assigned ID

```bash
awk -F $'\t' '$1 ~ /^#/ {print; next} $3~/^\./ {next} {print}' input.vcf > input.ID.vcf
```

### Keep only variants where FILTER is PASS

Use [SnpSift](http://pcingola.github.io/SnpEff/ss_introduction/) to filter VCF files.

The following keeps variants that have a `FILTER` value of `PASS`:

```bash
cat input.vcf | SnpSift filter "( FILTER = 'PASS' )" > input.PASS.vcf
```

Or, use `bcftools`:

```bash
bcftools view -f PASS input.vcf > input.PASS.vcf
```

### Keep only biallelic SNPs

```bash
bgzip input.vcf
tabix -fp vcf input.vcf.gz
bcftools view --min-alleles 2 --max-alleles 2 --types snps --output-file biallelic-snp.vcf.gz input.vcf.gz
```

### Keep only SNPs

```bash
bgzip input.vcf
tabix -fp vcf input.vcf.gz
bcftools view --types snps --output-file snp.vcf.gz input.vcf.gz
```

### Keep only indels

```bash
bgzip input.vcf
tabix -fp vcf input.vcf.gz
bcftools view --types indels --output-file indels.vcf.gz input.vcf.gz
```

### Keep only sites not found in a second VCF

Use the [bcftools isec command](https://samtools.github.io/bcftools/bcftools.html#isec).

```bash
bcftools isec --complement -c some \
file1.vcf.gz \
file2.vcf.gz \
-p sites-in-file-1-not-in-file-2
```

The sites in `file1.vcf.gz` that are not in `file2.vcf.gz` can be found in `sites-in-file-1-not-in-file-2/0000.vcf`.

The `-c some` causes sites to be considered equivalent when some subset of ALT alleles match. To require identical REF and ALT alleles use `-c none`. To match based only on position, use `-c all`.

### Count variants

```bash
grep -c -v '^#' input.ann.vcf
```

Or:

```bash
zgrep -c -v '^#' input.ann.vcf.gz
```

### Identify variants found in all VCF files

In this example there are `5` VCF files in the directory `vcfs`.

Prepare index files:

```bash
cd vcfs
find . -name "*.vcf" -exec bgzip {} \;
find . -name "*.vcf.gz" -exec tabix -p vcf {} \;
```

Determine the intersection (change `5` to match the number of input files):

```bash
mkdir overlaps
bcftools isec -p overlaps \
-n=5 -c all -w1 \
*.vcf.gz
mv overlaps/0000.vcf overlaps/intersection.vcf
```

Count the number of variants in `overlaps/intersection.vcf`:

```bash
grep -c -v '^#' overlaps/intersection.vcf
```

### Change sample order

First determine current order of samples:

```bash
bcftools query -l input.vcf > sample_order.txt
```

Next edit `sample_order.txt` to reflect the desired sample order.

For example, change the contents from this:

```text
M-9
M-10
M-15
M-16
M-2
```

To this:

```text
M-2
M-9
M-10
M-15
M-16
```

Generate a VCF file with the new sample order:

```bash
bgzip input.vcf
tabix -p vcf input.vcf.gz
bcftools view -S sample_order.txt \
input.vcf.gz > input.revised_sample_order.vcf
```

### Combine rows

Note that when using this approach the source files must have the same sample columns appearing in the same order.

In this example the contents of `snps.vcf` and `indels.vcf` are combined.

```bash
bgzip snps.vcf
bgzip indels.vcf
tabix -p vcf snps.vcf.gz
tabix -p vcf indels.vcf.gz
bcftools concat --allow-overlaps \
snps.vcf.gz \
indels.vcf.gz \
-Oz -o snps_and_indels.vcf.gz
```

### Convert a VCF file to an Excel file

First remove header content:

```bash
grep -v '^##' input.vcf > input.tsv
```

Then use [VisiData](https://github.com/saulpw/visidata) to convert the TSV file to an Excel file:

```bash
vd input.csv -b -o input.xlsx
```

### Perform case-control analysis

[SnpSift](http://pcingola.github.io/SnpEff/ss_introduction/) can generate p-values for different models.

In this example the `+++++` specifies that the first five samples are cases. The `-------` specifies that the next seven samples are controls. The `000` specifies that the last three samples should be ignored.

```bash
cat input.vcf | java -jar SnpSift.jar caseControl "+++++-------000" > input.case-control.vcf
```

The results can be filtered based on p-value. The following keeps variants where the p-value under the recessive model is less than `0.001`:

```bash
cat input.case-control.vcf | java -jar SnpSift.jar filter "CC_REC[0] < 0.001" > input.case-control.sig.vcf
```

### Remove annotations

Use `bcftools` `annotate`.

The following removes all `INFO` fields and all `FORMAT` fields except for `GT` and `GQ`:

```bash
bcftools annotate -x INFO,^FORMAT/GT,FORMAT/GQ input.vcf > input_GT_GQ.vcf
```

### Remove a sample

In this example sample `GM-2` is removed:

```bash
vcftools --remove-indv GM-2 --vcf input.vcf --recode --out output.vcf
```

### Remove sites with any missing genotypes

```bash
bcftools view -e 'GT[*]="mis"' input.vcf > output.vcf
```

### Remove all genotypes

```bash
bgzip input.vcf
tabix -p vcf input.vcf.gz
bcftools view -G input.vcf.gz -Oz -o output.vcf.gz
```

Or:

```bash
bcftools view -G input.vcf > output.vcf
```

### Print samples

```bash
bcftools query -l input.vcf.gz
```

Or:

```bash
bcftools query -l input.vcf
```

## vim

### Search and replace across multiple files

In this example search and replace operations are performed on all the `.html` files in a directory. First, open the files in multiple buffers in vim:

```bash
vim *.html
```

Then use `argdo` to perform a search and replace across all the files. In this example blank lines are removed:

```text
:argdo %s/^$//ge
```

In this example the text between `<p class="last-updated">` and `</p>` are replaced with the current date. Note the use of `\zs` and `\ze` so that the text between those tags is replaced and not the tags themselves:

```text
:argdo %s/<p class="last-updated">\zs[^<]*\ze<\/p>/\=strftime("%c")/ge
```

### Search and replace newlines

In replacement syntax use `\r` instead of `\n` to represent newlines. For example, to replace commas with newlines:

```text
:%s/,/\r/g
```

### Compare two files

```bash
vimdiff file1 file2
```

### Copy to the clipboard

```text
"+y
```

### Remove trailing whitespace

```text
:%s/\s\+$//e
```

### View ^M characters

```text
:e ++ff=unix
```

### Type tab characters

In insert mode type `Ctrl`-`v` then `tab`.
