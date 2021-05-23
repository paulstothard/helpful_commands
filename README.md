# Table of Contents

<!-- toc -->

- [awk](#awk)
  * [Add a header line to a file](#add-a-header-line-to-a-file)
  * [Convert a CSV file to a FASTA file](#convert-a-csv-file-to-a-fasta-file)
  * [Print lines in file when a certain column contains a specific value](#print-lines-in-file-when-a-certain-column-contains-a-specific-value)
  * [Replace certain values in specific columns](#replace-certain-values-in-specific-columns)
  * [Add up the values in a column](#add-up-the-values-in-a-column)
  * [Create a new column from two existing columns](#create-a-new-column-from-two-existing-columns)
  * [For each category in one column, add up the values in another column](#for-each-category-in-one-column-add-up-the-values-in-another-column)
  * [Print column names and numbers](#print-column-names-and-numbers)
  * [Print the values observed in a specific column, along with the number of times each value is observed](#print-the-values-observed-in-a-specific-column-along-with-the-number-of-times-each-value-is-observed)
  * [Print the number of lines exhibiting each distinct number of fields](#print-the-number-of-lines-exhibiting-each-distinct-number-of-fields)
  * [Print lines where certain fields contain values of interest](#print-lines-where-certain-fields-contain-values-of-interest)
  * [Write each row to a separate file named after the value in a specific column](#write-each-row-to-a-separate-file-named-after-the-value-in-a-specific-column)
  * [Split a multi-FASTA file into separate files named according to the sequence title](#split-a-multi-fasta-file-into-separate-files-named-according-to-the-sequence-title)
  * [Print only specific columns, identified by name in the first row](#print-only-specific-columns-identified-by-name-in-the-first-row)
  * [Print only the lines coming after a certain starting line and before a certain ending line](#print-only-the-lines-coming-after-a-certain-starting-line-and-before-a-certain-ending-line)
  * [Print only the first non-commented line](#print-only-the-first-non-commented-line)
  * [Print the average read length for a FASTQ file](#print-the-average-read-length-for-a-fastq-file)
  * [Sort lines based on order of IDs in another file](#sort-lines-based-on-order-of-ids-in-another-file)
- [brew](#brew)
  * [List installed packages](#list-installed-packages)
  * [View available packages](#view-available-packages)
  * [Install a package](#install-a-package)
  * [Add a third-party repository](#add-a-third-party-repository)
  * [Install directly from a third-party repository](#install-directly-from-a-third-party-repository)
  * [View packages available from brewsci/bio](#view-packages-available-from-brewscibio)
  * [List installed graphical applications](#list-installed-graphical-applications)
  * [View available graphical applications](#view-available-graphical-applications)
  * [Install a graphical application](#install-a-graphical-application)
- [Conda](#conda)
  * [Install Miniconda](#install-miniconda)
  * [Create an environment and install some packages](#create-an-environment-and-install-some-packages)
  * [Deactivate an environment](#deactivate-an-environment)
  * [Activate an environment](#activate-an-environment)
  * [List available packages](#list-available-packages)
  * [Search for a specific package](#search-for-a-specific-package)
  * [Add additional packages to an environment](#add-additional-packages-to-an-environment)
  * [List environments](#list-environments)
  * [List packages installed in the active environment](#list-packages-installed-in-the-active-environment)
- [csvkit](#csvkit)
  * [Convert Excel to CSV](#convert-excel-to-csv)
  * [Convert JSON to CSV](#convert-json-to-csv)
  * [Print column names](#print-column-names)
  * [Select a subset of columns](#select-a-subset-of-columns)
  * [Reorder columns](#reorder-columns)
  * [Find rows with matching cells](#find-rows-with-matching-cells)
  * [Convert to JSON](#convert-to-json)
  * [Generate summary statistics](#generate-summary-statistics)
  * [Query with SQL](#query-with-sql)
- [cut](#cut)
  * [Extract columns of interest](#extract-columns-of-interest)
  * [Extract a range of columns](#extract-a-range-of-columns)
  * [Extract everything except one column](#extract-everything-except-one-column)
  * [Extract characters](#extract-characters)
  * [Change field separators](#change-field-separators)
- [datamash](#datamash)
  * [Group records by one column and print information about each group](#group-records-by-one-column-and-print-information-about-each-group)
  * [Print statistics for a column](#print-statistics-for-a-column)
  * [Transpose](#transpose)
- [Docker](#docker)
  * [Perform a sequence comparison using legacy BLAST](#perform-a-sequence-comparison-using-legacy-blast)
  * [Annotate sequence variants using VEP](#annotate-sequence-variants-using-vep)
  * [Annotate a bacterial genome using Prokka](#annotate-a-bacterial-genome-using-prokka)
  * [List running containers](#list-running-containers)
  * [Stop a container](#stop-a-container)
  * [Kill all running containers](#kill-all-running-containers)
  * [Delete all containers that are not running](#delete-all-containers-that-are-not-running)
- [find](#find)
  * [Perform a series of commands on files returned by find](#perform-a-series-of-commands-on-files-returned-by-find)
  * [Copy the files returned by find, naming the copies after a directory in the path](#copy-the-files-returned-by-find-naming-the-copies-after-a-directory-in-the-path)
  * [Switch to the directory containing each file and execute a command](#switch-to-the-directory-containing-each-file-and-execute-a-command)
  * [Find large files](#find-large-files)
- [Git](#git)
  * [Create a new Git repository](#create-a-new-git-repository)
  * [Sync a repository to your local machine](#sync-a-repository-to-your-local-machine)
  * [Mark changed files to be included in the next commit](#mark-changed-files-to-be-included-in-the-next-commit)
  * [Undo a Git add before a commit](#undo-a-git-add-before-a-commit)
  * [Remove files from the repository](#remove-files-from-the-repository)
  * [Move or rename a file or directory](#move-or-rename-a-file-or-directory)
  * [Save the marked files to the local Git repository](#save-the-marked-files-to-the-local-git-repository)
  * [Push a commit on your local branch to a remote repository](#push-a-commit-on-your-local-branch-to-a-remote-repository)
  * [Pull a change from a remote repository to your local branch](#pull-a-change-from-a-remote-repository-to-your-local-branch)
  * [Add or edit a remote repository](#add-or-edit-a-remote-repository)
  * [Create and merge Git branches](#create-and-merge-git-branches)
  * [Specify files to ignore](#specify-files-to-ignore)
  * [Check the status of a working directory](#check-the-status-of-a-working-directory)
  * [Tag a release](#tag-a-release)
- [grep](#grep)
  * [Count matches](#count-matches)
  * [Get the line number of a match](#get-the-line-number-of-a-match)
  * [Remove files that contain a match](#remove-files-that-contain-a-match)
  * [Remove files that do not contain a match](#remove-files-that-do-not-contain-a-match)
  * [Remove lines that match](#remove-lines-that-match)
- [Miller](#miller)
  * [Extract the first 10 records of a CSV file](#extract-the-first-10-records-of-a-csv-file)
  * [Extract the last 10 records of a CSV file](#extract-the-last-10-records-of-a-csv-file)
  * [Convert formats](#convert-formats)
  * [View stats](#view-stats)
  * [Filter records](#filter-records)
  * [Sort records](#sort-records)
  * [Extract columns](#extract-columns)
  * [Edit columns](#edit-columns)
  * [Other actions](#other-actions)
  * [Combine actions](#combine-actions)
- [Other](#other)
  * [Obtain your public IP address and network information](#obtain-your-public-ip-address-and-network-information)
  * [Copy an ssh public key to another system](#copy-an-ssh-public-key-to-another-system)
  * [Download files from an FTP server](#download-files-from-an-ftp-server)
  * [Download files from Google Drive](#download-files-from-google-drive)
  * [Extract a file](#extract-a-file)
  * [Add a header to all files with a certain extension, getting the header from another file](#add-a-header-to-all-files-with-a-certain-extension-getting-the-header-from-another-file)
  * [View STDOUT and append it to a file](#view-stdout-and-append-it-to-a-file)
  * [Redirect STDERR to STDOUT and view both and append both to a file](#redirect-stderr-to-stdout-and-view-both-and-append-both-to-a-file)
  * [Change the extension of multiple files](#change-the-extension-of-multiple-files)
  * [Add text to the beginning of a file](#add-text-to-the-beginning-of-a-file)
  * [Add text or a header to the beginning of all files with a particular file extension](#add-text-or-a-header-to-the-beginning-of-all-files-with-a-particular-file-extension)
  * [Find common lines between files](#find-common-lines-between-files)
  * [Convert a CSV file to a Markdown table](#convert-a-csv-file-to-a-markdown-table)
  * [Convert PDF files to PNG files](#convert-pdf-files-to-png-files)
  * [Convert PNG files to a single PDF file](#convert-png-files-to-a-single-pdf-file)
  * [Convert a DOCX file to a PDF file](#convert-a-docx-file-to-a-pdf-file)
  * [Convert an Excel file to a CSV file](#convert-an-excel-file-to-a-csv-file)
  * [Convert a CSV file to an Excel file](#convert-a-csv-file-to-an-excel-file)
  * [Convert a TSV file to an Excel file](#convert-a-tsv-file-to-an-excel-file)
  * [Convert an HTML file to a PDF file](#convert-an-html-file-to-a-pdf-file)
  * [Convert a website to a PDF file](#convert-a-website-to-a-pdf-file)
  * [Convert an HTML file to a PNG file](#convert-an-html-file-to-a-png-file)
  * [Convert a Markdown file to a PDF file](#convert-a-markdown-file-to-a-pdf-file)
  * [Convert a Markdown file to an HTML file](#convert-a-markdown-file-to-an-html-file)
  * [Crop an image and add a white border](#crop-an-image-and-add-a-white-border)
  * [Resize an image](#resize-an-image)
  * [Format a CSV file into columns and examine its content](#format-a-csv-file-into-columns-and-examine-its-content)
  * [Format code](#format-code)
  * [Change Bash prompt temporarily](#change-bash-prompt-temporarily)
  * [Check bash/sh shell scripts for potential issues](#check-bashsh-shell-scripts-for-potential-issues)
  * [Take a screenshot of a window on macOS](#take-a-screenshot-of-a-window-on-macos)
  * [Take a webpage screenshot using Firefox](#take-a-webpage-screenshot-using-firefox)
  * [Create PowerPoint slides from a Markdown file](#create-powerpoint-slides-from-a-markdown-file)
  * [Run commands at scheduled times using cron](#run-commands-at-scheduled-times-using-cron)
  * [Record your terminal to an animated GIF](#record-your-terminal-to-an-animated-gif)
  * [Create an animated GIF from a YouTube video](#create-an-animated-gif-from-a-youtube-video)
  * [Create a collection of MP3 files from a YouTube playlist](#create-a-collection-of-mp3-files-from-a-youtube-playlist)
  * [Download a GenBank file with curl](#download-a-genbank-file-with-curl)
  * [Perform a calculation on the command line](#perform-a-calculation-on-the-command-line)
  * [Save the output of a command in a variable](#save-the-output-of-a-command-in-a-variable)
- [parallel](#parallel)
  * [Extract files in parallel](#extract-files-in-parallel)
  * [Compress files in parallel](#compress-files-in-parallel)
  * [Process files in pairs](#process-files-in-pairs)
  * [Perform BLAST in parallel](#perform-blast-in-parallel)
- [paste](#paste)
  * [Combine columns with paste](#combine-columns-with-paste)
- [Perl](#perl)
  * [Get a random sample of lines from a text file while excluding the header line](#get-a-random-sample-of-lines-from-a-text-file-while-excluding-the-header-line)
  * [Convert a FASTA file to a CSV file with column names](#convert-a-fasta-file-to-a-csv-file-with-column-names)
  * [Count the number of lines that match a regular expression](#count-the-number-of-lines-that-match-a-regular-expression)
  * [Extract FASTA sequences from a file based on a file of sequence names of interest](#extract-fasta-sequences-from-a-file-based-on-a-file-of-sequence-names-of-interest)
  * [Add a FASTA title to the start of a sequence in RAW format](#add-a-fasta-title-to-the-start-of-a-sequence-in-raw-format)
  * [Remove commas located within quoted fields in a CSV file and create a tab-delimited file](#remove-commas-located-within-quoted-fields-in-a-csv-file-and-create-a-tab-delimited-file)
  * [Replace tabs with commas and remove quotes](#replace-tabs-with-commas-and-remove-quotes)
  * [Sort sections in a Markdown file based on headings](#sort-sections-in-a-markdown-file-based-on-headings)
  * [Search and replace text on each line](#search-and-replace-text-on-each-line)
  * [Print matches that may span multiple lines](#print-matches-that-may-span-multiple-lines)
  * [Print matches after additional editing](#print-matches-after-additional-editing)
  * [Format Perl code](#format-perl-code)
- [Process multiple files](#process-multiple-files)
  * [for loop](#for-loop)
  * [while loop](#while-loop)
  * [find with -exec](#find-with--exec)
  * [find with xargs](#find-with-xargs)
  * [parallel](#parallel-1)
- [Process multiple files in pairs](#process-multiple-files-in-pairs)
- [R](#r)
  * [Compare two data sets to find differences](#compare-two-data-sets-to-find-differences)
  * [Visualize the degree of overlap among gene sets](#visualize-the-degree-of-overlap-among-gene-sets)
  * [Cluster gene lists based on overlap and identify shared genes](#cluster-gene-lists-based-on-overlap-and-identify-shared-genes)
  * [Transpose a data frame](#transpose-a-data-frame)
  * [Split two-allele genotypes into two columns](#split-two-allele-genotypes-into-two-columns)
  * [Split multi-locus genotype strings into multiple columns and decode genotypes](#split-multi-locus-genotype-strings-into-multiple-columns-and-decode-genotypes)
  * [Change values in a column based on values in another column](#change-values-in-a-column-based-on-values-in-another-column)
  * [Add comment lines to output](#add-comment-lines-to-output)
  * [Filter and sort rows](#filter-and-sort-rows)
  * [Add columns from one tibble to another](#add-columns-from-one-tibble-to-another)
  * [Combine multiple input files](#combine-multiple-input-files)
- [rsync](#rsync)
  * [Sync a directory on local system](#sync-a-directory-on-local-system)
  * [Sync a directory to a remote system](#sync-a-directory-to-a-remote-system)
  * [Sync a directory from a remote system](#sync-a-directory-from-a-remote-system)
- [Singularity](#singularity)
  * [Creating an image from a container stored in Docker Hub](#creating-an-image-from-a-container-stored-in-docker-hub)
- [sbatch](#sbatch)
  * [Count lines in compressed fastq files](#count-lines-in-compressed-fastq-files)
- [sed](#sed)
  * [Add a header line to a file](#add-a-header-line-to-a-file-1)
  * [Print a specific line of a file](#print-a-specific-line-of-a-file)
  * [Change filenames using a regular expression](#change-filenames-using-a-regular-expression)
  * [Search and replace on lines](#search-and-replace-on-lines)
  * [Delete lines](#delete-lines)
- [Share data with project group members](#share-data-with-project-group-members)
- [Slurm](#slurm)
  * [View statistics related to the efficiency of resource usage of a completed job](#view-statistics-related-to-the-efficiency-of-resource-usage-of-a-completed-job)
  * [View jobs](#view-jobs)
  * [View running jobs](#view-running-jobs)
  * [View pending jobs](#view-pending-jobs)
  * [View detailed information for a specific job](#view-detailed-information-for-a-specific-job)
  * [View accounting information for completed jobs](#view-accounting-information-for-completed-jobs)
  * [Cancel a job](#cancel-a-job)
  * [Cancel all jobs](#cancel-all-jobs)
  * [Start an interactive session](#start-an-interactive-session)
- [sort](#sort)
  * [Alphabetical sort](#alphabetical-sort)
  * [Specify the sort field](#specify-the-sort-field)
  * [Use multiple sort fields](#use-multiple-sort-fields)
  * [Sort a file with a header row](#sort-a-file-with-a-header-row)
- [tmux](#tmux)
  * [Start a tmux session](#start-a-tmux-session)
  * [Detach a tmux session](#detach-a-tmux-session)
  * [List tmux sessions](#list-tmux-sessions)
  * [Join a tmux session](#join-a-tmux-session)
  * [Create a multi-pane tmux session](#create-a-multi-pane-tmux-session)
  * [Navigate between tmux panes](#navigate-between-tmux-panes)
  * [Kill a tmux session](#kill-a-tmux-session)
- [tr](#tr)
  * [Translate characters](#translate-characters)
  * [Delete characters](#delete-characters)
  * [Squeeze characters](#squeeze-characters)
- [vcftools and bcftools](#vcftools-and-bcftools)
  * [Extract variants from a region of interest and write to a new vcf file](#extract-variants-from-a-region-of-interest-and-write-to-a-new-vcf-file)
  * [Extract variants from multiple regions of interest and write to a new vcf file](#extract-variants-from-multiple-regions-of-interest-and-write-to-a-new-vcf-file)
- [vim](#vim)
  * [Search and replace across multiple files](#search-and-replace-across-multiple-files)
  * [Search and replace newlines](#search-and-replace-newlines)
  * [Compare two files](#compare-two-files)
  * [Copy to the clipboard](#copy-to-the-clipboard)

<!-- tocstop -->

## awk

### Add a header line to a file

```bash
awk 'BEGIN{print "my header text"}1' input
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
awk -F $'\t' '{count[$9]++}END{for(j in count) print j,"("count[j]" counts)"}' input.tab
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

### Write each row to a separate file named after the value in a specific column

In this example each file is named after the value in column `1`:

```bash
awk -F '\t' '{ fname = $1 ".txt"; print >>fname; close(fname) }' input.tab
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

In this example, the records in `file_to_sort.csv` have an identifier in column `1` that is present in `sorted_ids.txt`, which has a single column:

```bash
awk -F, 'NR==FNR {a[$1]=$0; next} ($0 in a) {print a[$0]}' file_to_sort.csv sorted_ids.txt
```

If both files have a single header line the following can be used to generate the sorted output with the restored header line:

```bash
(head -n 1 file_to_sort.csv && awk -F, 'NR==FNR {a[$1]=$0; next} ($0 in a) {print a[$0]}' <(tail -n +2 file_to_sort.csv) <(tail -n +2 sorted_ids.txt)) > sorted.csv
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

Download the legacy BLAST Docker image:

```bash
docker pull quay.io/biocontainers/blast-legacy:2.2.26--2
```

Create a container from the image and run `formatdb` to create a formatted database. In this example the database is created from a DNA sequence file called `database.fasta`, located in the current directory:

```bash
docker run -it --rm -v $(pwd):/directory -w /directory quay.io/biocontainers/blast-legacy:2.2.26--2 formatdb -i database.fasta -p F
```

To perform a blastn search using the formatted database and a query called `query.fasta` when the file is also located in the current directory:

```bash
docker run -it --rm -v $(pwd):/directory -w /directory quay.io/biocontainers/blast-legacy:2.2.26--2 blastall -p blastn -d database.fasta -i query.fasta
```

To perform a blastn search using the formatted database and a query called `query.fasta` when the query is located in a different directory (in this example your home directory):

```bash
docker run -it --rm -v $(pwd):/directory/database -v ${HOME}:/directory/query -w /directory quay.io/biocontainers/blast-legacy:2.2.26--2 blastall -p blastn -d database/database.fasta -i query/query.fasta
```

### Annotate sequence variants using VEP

Download the VEP Docker image:

```bash
docker pull ensemblorg/ensembl-vep
```

Download cache files for the bovine genome and VEP plugins:

```bash
docker run -t -i -v $(pwd):/opt/vep/.vep ensemblorg/ensembl-vep perl INSTALL.pl -a cfp -s bos_taurus -y ARS-UCD1.2 -g all
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
    docker run -v $(pwd):/opt/vep/.vep ensemblorg/ensembl-vep \
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

Download the Prokka Docker image:

```bash
docker pull staphb/prokka:latest
```

Create a container from the image and run `prokka` to annotate the sequence. In this example the genome sequence to be annotated is in a file called `sequence.fasta`, located in the current directory, and four CPUs are used:

```bash
docker run --rm -v $(pwd):/dir -w /dir staphb/prokka:latest prokka sequence.fasta --cpus 4
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

### Delete all containers that are not running

```bash
docker container rm $(docker ps -a -q)
```

## find

### Perform a series of commands on files returned by find

In this example `$'...'` is used for quoting, as it can contain escaped single quotes, and `tail` is used to skip a header line, `awk` is used to count the number of occurrences of each category in column 3 and print the category and counts, and `sort` is used to sort the categories by count from largest to smallest with ties broken by sorting on category name:

```bash
find . -type f -name "*.gff" -print0 | xargs -0 -I{} sh -c $'tail -n +2 "$1" | awk -F $\'\t\' \'{count[$3]++}END{for(j in count) print j,count[j]}\' | sort -k 2,2nr -k 1,1> "$1.cog_counts.txt"' -- {}
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

```
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

Note that adding a `.gitignore` file will not remove tracked files; this must be done with `git rm`. See [Removing files from the repository](#removing-files-from-the-repository).

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

### Remove lines that match

Keep everything except lines starting with `#`:

```bash
grep -v '^#' input.txt
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

### Download files from an FTP server

Replace `host`, `account`, `password`, and `port` with their corresponding values in the following command:

```bash
wget -S -d -c -t 45 -v -r ftp://account:password@host:port/*
```

### Download files from Google Drive

[Rclone](https://rclone.org) can be used to download data from many cloud storage providers.

First, follow the [configuration instructions](https://rclone.org/drivei). The commands below assume that the remote system was named `my_google_drive` during the configuration.

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

### Extract a file

The following bash function can be used to extract a variety of file types.

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

### Add a header to all files with a certain extension, getting the header from another file

In this example the header is added to `.tab` files and comes from a file called `header.txt`. The files with the header added are saved with a `.new` extension added:

```bash
for f in *.tab; do new=`echo $f | sed 's/\(.*\)\.tab/\1.tab.new/'`; paste -sd'\n' \header.txt "$f" > "$new"; done
```

To replace the `.tab` files the `.new` files:

```bash
for f in *.new; do new=`echo $f | sed 's/\(.*\)\.new/\1/'`; mv "$f" "$new"; done
```

### View STDOUT and append it to a file

```bash
some_command | tee -a output.txt
```

### Redirect STDERR to STDOUT and view both and append both to a file

```bash
some_command 2>&1 | tee -a log
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

### Convert a CSV file to a Markdown table

The following uses [csv2md](https://github.com/pstaender/csv2md). The `awk` command can be used if some rows of the input have missing fields on the end:

```bash
awk -F, -v OFS="," 'NR==1 {cols=NF} {$1=$1; for (i=NF+1; i <= cols; i++) $i = "."} 1' input.csv > temp.csv
csv2md -p < temp.csv > output.md
```

### Convert PDF files to PNG files

The following uses `find` and the `pdftoppm` command from the [poppler](https://poppler.freedesktop.org) package to generate a PNG image of the first page of every PDF file in the working directory:

```bash
find . -name "*.pdf" -exec pdftoppm -f 1 -l 1 -png {} {} \;
```

### Convert PNG files to a single PDF file

The following uses [ImageMagick](https://imagemagick.org):

```bash
convert *.png output.pdf
```

### Convert a DOCX file to a PDF file

The following uses [LibreOffice](https://www.libreoffice.org):

```bash
soffice --headless --convert-to pdf --outdir . word_file.docx
```

The following uses [pandoc](https://pandoc.org) and on macOS also requires [basictex](https://www.tug.org/mactex/morepackages.html):

```bash
pandoc word_file.docx --output word_file.pdf
```

### Convert an Excel file to a CSV file

The following uses [csvkit](https://github.com/wireservice/csvkit):

```bash
in2csv data.xls > data.csv
```

### Convert a CSV file to an Excel file

The following uses `ssconvert`, which is distributed with Gnumeric:

```bash
ssconvert input.csv output.xlsx
```

### Convert a TSV file to an Excel file

The following uses `ssconvert`, which is distributed with Gnumeric:

```bash
ssconvert input.tsv output.xls
```

### Convert an HTML file to a PDF file

The following uses [wkhtmltopdf](https://wkhtmltopdf.org):

```bash
wkhtmltopdf http://google.com google.pdf
```

### Convert a website to a PDF file

The following uses [wkhtmltopdf](https://wkhtmltopdf.org) and [gs](https://www.ghostscript.com/index.html):

```bash
url=http://www.3rs-reduction.co.uk/html/main_menu.html; depth=1
wget --spider --force-html -r -l${depth} ${url} 2>&1 | grep '^--' | awk '{ print $3 }' | grep -i '\.\(html\|htm\)$' | uniq > url-list.txt
while read i; do wkhtmltopdf "$i" "$(echo "$i" | sed -e 's/https\?:\/\///' -e 's/\//-/g' ).pdf"; done < url-list.txt
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=merged-output.pdf $(ls -lrt -1 *.pdf)
```

### Convert an HTML file to a PNG file

The following uses [wkhtmltoimage](https://wkhtmltopdf.org):

```bash
wkhtmltoimage -f png http://google.com google.png
```

Another approach, which may work better for complex web sites, is to use [pageres-cli](https://github.com/sindresorhus/pageres-cli). The following creates an image that is 4485 pixels wide (5 * 897):

```bash
pageres http://google.com 897x1090 --crop --scale=5 --filename='google'
```

### Convert a Markdown file to a PDF file

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

### Convert a Markdown file to an HTML file

The commmand below uses [pandoc](https://pandoc.org) and the `pandoc.css` file available [here](https://gist.github.com/killercup/5917178).

```bash
pandoc -f markdown -t html -o output.html input.md --css=pandoc.css --self-contained
```

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

### Change Bash prompt temporarily

```bash
PS1="$ "
```

### Check bash/sh shell scripts for potential issues

Use [shellcheck](https://github.com/koalaman/shellcheck):

```bash
shellcheck some_script.sh
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

### Create PowerPoint slides from a Markdown file

[Pandoc](https://pandoc.org) can be used to generate PowerPoint slides. The following Markdown text describes several slides with notes:

````
% Presentation title
% Name
% Date

# Section title

## Slide title

Single bulleted list:

- list item
- list item
- list item

::: notes

Speaker notes

:::

## Slide title

Single bulleted list:

- list item
  - list item
  - list item
    - list item
- list item

::: notes

Speaker notes

:::

## Slide title

Single ordered list:

1. list item
   1. list item
   1. list item
      1. list item
   1. list item
1. list item

::: notes

Speaker notes

:::

## Slide title

:::::::::::::: {.columns}

::: {.column width="50%"}

Left column:

- list item
  - list item
  - list item
    - list item
- list item

:::

::: {.column width="50%"}

Right column:

1. list item
   1. list item
   1. list item
      1. list item
   1. list item
1. list item

:::

::::::::::::::

::: notes

Speaker notes

:::

## Slide title

:::::::::::::: {.columns}

::: {.column width="50%"}

Left column:

- Bullet
- Bullet
- Bullet

:::

::: {.column width="50%"}

![image caption](image.png)

:::

::::::::::::::

::: notes

Speaker notes

:::

## Slide title

| Program | Version | Purpose                                   |
|---------|---------|-------------------------------------------|
| Prokka  | 1.14.5  | Genome annotation                         |
| FastQC  | 0.11.9  | Illumina read quality assessment          |
| Snippy  | 4.6.0   | Genome comparisons                        |
| CGView  | 2.0.2   | Genome visualization                      |
| SPAdes  | 3.12.0  | Genome assembly                           |
| Quast   | 5.0.2   | Genome assembly quality assessment        |
| NUCmer  | 3.1     | Genome alignment                          |
| Circos  | 0.69-8  | Genome visualization                      |
| Infoseq | 6.6.0.0 | Sequence summary statistics               |
| R       | 3.6.2   | Visualizing overlaps between variant sets |

::: notes

Speaker notes

:::

## Slide title

```r
snp <- c('ABCA12', 'APAF1', 'ARS-BFGL-BAC-10172', 'ARS-BFGL-BAC-1020')
sample1 <- c('AA', 'CC', 'GG', 'AA')
sample2 <- c('AA', 'CC', 'AG', 'GG')
genotypes <- data.frame(snp, sample1, sample2)
```

::: notes

Speaker notes

:::
````

To convert the Markdown file to PowerPoint slides:

```bash
pandoc input.md -o slides.pptx
```

Alternatively, to create the slides using an existing PowerPoint template or presentation for formatting:

```bash
pandoc input.md -o slides.pptx --reference-doc some_template.potx
```

The formatting of individual slides can then be adjusted in PowerPoint, using the `Design` tab and the `Design Ideas` button. Slide numbers and headers and footers can be added using `View->Slide Master` followed by `Insert`, and then `Header & Footer`.

To reduce the size of the file, use `File->Compress Pictures...`.

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

### Create a collection of MP3 files from a YouTube playlist

The following requires [youtube-dl](https://github.com/ytdl-org/youtube-dl) and [ffmpeg](https://ffmpeg.org/):

```bash
youtube-dl -x -i --audio-format mp3 --audio-quality 320K --embed-thumbnail --geo-bypass https://www.youtube.com/playlist?list=PL92319EECC1754042
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

### Perform a calculation on the command line

Use [bc](https://www.gnu.org/software/bc/):

```bash
echo "2*(42+42)" | bc
```

Or use Bash arithmetic expansion:

```bash
n=6
echo "$(( n - 1 * 2 ))"
answer="$(( n - 1 * 2 ))"
```

### Save the output of a command in a variable

Use a subshell:

```bash
result=$(echo "sqrt(16)" | bc -l)
```

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

### Process files in pairs

In the following example paired-end reads with names like `sampleA_1.fastq.gz` and `sampleA_2.fastq.gz` in a directory called `data` are mapped to a reference called `ref` using `bowtie2`:

```bash
parallel -j2 "bowtie2 --threads 4 -x ref -k1 -q -1 {1} -2 {2} -S {1/.}.sam >& {1/.}.log" ::: data/*_1.fastq.gz :::+ data/*_2.fastq.gz
```

The `{1/.}` and `{2/.}` remove the path and the extension from the two files in each pair.

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
$ paste -d "," genotype_conversion.csv SNP_Map.tab
```

Note that the content from `SNP_Map.tab` still contains tab-delimited values.

To remove tabs from `SNP_Map.tab` you could first create a CSV version of that input file:

```bash
$ cut -d $'\t' -f 1- --output-delimiter=',' SNP_Map.tab > SNP_Map.csv
$ paste -d "," genotype_conversion.csv SNP_Map.csv
```

It is possible to do it without first creating `SNP_Map.csv` by using process substitution.

In the following the command between `<(` and `)` is first run and its output becomes input for `paste`:

```bash
$ paste -d "," genotype_conversion.csv \
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

### Replace tabs with commas and remove quotes

```bash
perl -p -e 's/\t/,/g;' -e 's/"//g' input.tab > output.csv
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

### Print matches that may span multiple lines

```bash
perl -0777 -ne 'while (m/^\s+\/translation="([^"]+)"/gm) {print "$1\n"}' NM_001271626.3.gbk
```

### Print matches after additional editing

```bash
perl -0777 -ne 'while (m/^\s+\/translation="([^"]+)"/gm) {$out = $1; $out =~ s/\s//g; print "$out\n"}' NM_001271626.3.gbk
```

### Format Perl code

The following uses `perltidy` to reformat the code in `testfile.pl` and will create a file called `testfile.pl.tdy`.

```bash
perltidy testfile.pl
```

## Process multiple files

### for loop

Change all `.fasta` files in the current directory to `.fna` files:

```bash
for f in *.fasta; do new=`echo $f | sed 's/\(.*\)\.fasta/\1.fna/'`; mv "$f" "$new"; done
```

### while loop

Print the number of lines in every `.csv` or `.tab` file in or below current directory:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) | while read f; do wc -l "$f"; done
```

Print the number of lines in every `.csv` or `.tab` file in or below current directory and redirect the results to a single file:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) | while read f; do wc -l "$f" >> output.txt; done
```

Print the number of lines in every `.csv` or `.tab` file in or below current directory and redirect the results to separate files:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) | while read f; do wc -l "$f" > "${f}.output.txt"; done
```

### find with -exec

Change all `.fasta` files in current directory to `.fna` files by appending a `.fna` extension:

```bash
find . -type f -name "*.fasta" -exec mv {} {}.fna \;
```

Print the number of lines in every `.csv` or `.tab` file in or below current directory:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -exec wc -l {} \;
```

Print the number of lines in every `.csv` or `.tab` file in or below current directory and redirect the results to a single file:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -exec wc -l {} \; > output.txt
```

Print the number of lines in every `.csv` or `.tab` file in or below current directory and redirect the results to separate files:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -exec sh -c 'wc -l "$1" > "$1.output.txt"' -- {} \;
```

### find with xargs

Change all `.fasta` files in current directory to `.fna` files by appending a `.fna` extension:

```bash
find . -type f -name "*.fasta" -print0 | xargs -0 -I{} mv {} {}.fna
```

Print the number of lines in every `.csv` or `.tab` file in or below current directory:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -print0 | xargs -0 -I{} wc -l {}
```

Print the number of lines in every `.csv` or `.tab` file in or below current directory and redirect the results to a single file:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -print0 | xargs -0 -I{} wc -l {} > output.txt
```

Print the number of lines in every `.csv` or `.tab` file in or below current directory and redirect the results to separate files:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -print0 | xargs -0 -I{} sh -c 'wc -l "$1" > "$1.output.txt"' -- {}
```

Print the number of lines in every `.csv` or `.tab` file in or below current directory and redirect the results to separate files. Process up to `4` files in parallel:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -print0 | xargs -n1 -P4 -0 -I{} sh -c 'wc -l "$1" > "$1.output.txt"' -- {}
```

### parallel

See the [parallel examples](#parallel).

## Process multiple files in pairs

High-throughput sequencing data is often distributed as pairs of files corresponding to the two different read sets generated for each sample, e.g.:

```
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

## R

### Compare two data sets to find differences

In this example SNP location information assigned by two different algorithms is compared using the `compareDF` package. The two data sets (position information) generated by the differing algorithms are read from files. SNP name is used to match rows across the data sets, and then the 'chromosome' and 'position' are compared. SNP records for which 'chromosome' or 'position' differ between the data sets are written to an Excel file:

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

```
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

```
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

```
snp                 sample1  sample2
ABCA12              AA       AA
APAF1               CC       CC
ARS-BFGL-BAC-10172  GG       AG
ARS-BFGL-BAC-1020   AA       GG
```

To this:

```
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

```
sample   ABCA12  APAF1  ARS-BFGL-BAC-10172  ARS-BFGL-BAC-1020
sample1  AA      CC     GG                  AA
sample2  AA      CC     AG                  GG
```

To this:

```
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

```
sample          alleles
HOCANF12689774  1000112112
HOCANF12689787  2011112012
HOCANF12689790  1011002122
```

To this:

```
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

```
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

### Creating an image from a container stored in Docker Hub

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

```
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

```
HI.5173.001.NEBNext_Index_12.DG15B032198-1_R1.fastq.gz.out
HI.5173.001.NEBNext_Index_12.DG15B032198-1_R1.fastq.gz.err
```

The `.out` files will contain the name of the input file and the number of lines, for example:

```
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

### Add a header line to a file

```bash
sed $'1s/^/my header text\\\n&/' input
```

### Print a specific line of a file

In this example line `26404`:

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
salloc --time=3:00:00 --nodes=1 --ntasks-per-node=1 --cpus-per-task=8 --mem=64000M --account=def-someuser
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
(head -n 1 sequenced_samples.csv && tail -n +2 sequenced_samples.csv | sort)
```

```bash
cat sequenced_samples.csv | awk 'NR<2{print $0; next}{print $0| "sort"}'
```

## tmux

### Start a tmux session

```bash
tmux new -s session_name
```

### Detach a tmux session

`C-b` refers to Control+B. The following means press Control+B, release, and then press "d":

```bash
C-b d
```

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

### Navigate between tmux panes

The commands that work will depend on how tmux is configured:

`C-j` moves up.
`C-k` moves down.
`C-h` moves left.
`C-l` moves right.

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

## vcftools and bcftools

### Extract variants from a region of interest and write to a new vcf file

Note that if the vcf file is gzip compressed (i.e. has a `.gz` extension), use `--gzvcf` instead of `--vcf`.

```bash
vcftools --vcf Chr5.vcf --out Chr5_filtered --chr 5 --from-bp 1 --to-bp 100000 --recode --recode-INFO-all
```

### Extract variants from multiple regions of interest and write to a new vcf file

```bash
bgzip Chr5.vcf
tabix -fp vcf Chr5.vcf.gz 
bcftools view -r 5:1-10000,5:200000-210000 -o output.vcf Chr5.vcf.gz
``` 

## vim

### Search and replace across multiple files

In this example search and replace operations are performed on all the `.html` files in a directory. First, open the files in multiple buffers in vim:

```bash
vim *.html
```

Then use `argdo` to perform a search and replace across all the files. In this example blank lines are removed:

```
:argdo %s/^$//ge
```

In this example the text between `<p class="lastupdated">` and `</p>` are replaced with the current date. Note the use of `\zs` and `\ze` so that the text between those tags is replaced and not the tags themselves:

```
:argdo %s/<p class="lastupdated">\zs[^<]*\ze<\/p>/\=strftime("%c")/ge
```

### Search and replace newlines

In replacement syntax use `\r` instead of `\n` to represent newlines. For example, to replace commas with newlines:

```
:%s/,/\r/g 
```

### Compare two files

```bash
vimdiff file1 file2 
```

### Copy to the clipboard

```
"+y
```
