# Table of Contents

<!-- toc -->

- [awk](#awk)
  * [Convert a CSV file to a FASTA file](#convert-a-csv-file-to-a-fasta-file)
  * [Print lines in file when a certain column contains a specific value](#print-lines-in-file-when-a-certain-column-contains-a-specific-value)
  * [Replace certain values in specific columns](#replace-certain-values-in-specific-columns)
  * [Sum values in one column based on categories given in another column](#sum-values-in-one-column-based-on-categories-given-in-another-column)
  * [Print column names and numbers](#print-column-names-and-numbers)
  * [Print the values observed in a specific column, along with the number of times each value is observed](#print-the-values-observed-in-a-specific-column-along-with-the-number-of-times-each-value-is-observed)
  * [Print the number of lines exhibiting each distinct number of fields](#print-the-number-of-lines-exhibiting-each-distinct-number-of-fields)
  * [Print lines where certain fields contain values of interest](#print-lines-where-certain-fields-contain-values-of-interest)
  * [Write each row to a separate file named after the value in a specific column](#write-each-row-to-a-separate-file-named-after-the-value-in-a-specific-column)
  * [Split a multi-FASTA file into separate files named according to the sequence title](#split-a-multi-fasta-file-into-separate-files-named-according-to-the-sequence-title)
  * [Print only specific columns, identified by name in the first row](#print-only-specific-columns-identified-by-name-in-the-first-row)
  * [Print only the lines coming after a certain starting line and before a certain ending line](#print-only-the-lines-coming-after-a-certain-starting-line-and-before-a-certain-ending-line)
  * [Print only the first non-commented line](#print-only-the-first-non-commented-line)
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
  * [Create an environment and install some programs](#create-an-environment-and-install-some-programs)
  * [Deactivate an environment](#deactivate-an-environment)
  * [Activate an environment](#activate-an-environment)
  * [Add additional programs to an environment](#add-additional-programs-to-an-environment)
  * [List environments](#list-environments)
- [Docker](#docker)
  * [Perform a sequence comparison using legacy BLAST](#perform-a-sequence-comparison-using-legacy-blast)
  * [Annotate sequence variants using VEP](#annotate-sequence-variants-using-vep)
  * [Annotate a bacterial genome using Prokka](#annotate-a-bacterial-genome-using-prokka)
- [find](#find)
  * [Perform a series of commands on files returned by find](#perform-a-series-of-commands-on-files-returned-by-find)
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
- [Other](#other)
  * [Obtain your public IP address and network information](#obtain-your-public-ip-address-and-network-information)
  * [Copy an ssh public key to another system](#copy-an-ssh-public-key-to-another-system)
  * [Download files from an FTP server](#download-files-from-an-ftp-server)
  * [Download files from Google Drive](#download-files-from-google-drive)
  * [Combine the columns in two tab-delimited files](#combine-the-columns-in-two-tab-delimited-files)
  * [Add a header to all files with a certain extension, getting the header from another file](#add-a-header-to-all-files-with-a-certain-extension-getting-the-header-from-another-file)
  * [View STDOUT and append it to a file](#view-stdout-and-append-it-to-a-file)
  * [Redirect STDERR to STDOUT and view both and append both to a file](#redirect-stderr-to-stdout-and-view-both-and-append-both-to-a-file)
  * [Change the extension of multiple files](#change-the-extension-of-multiple-files)
  * [Add text or a header to the beginning of all files with a particular file extension](#add-text-or-a-header-to-the-beginning-of-all-files-with-a-particular-file-extension)
  * [Find common lines between files](#find-common-lines-between-files)
  * [Convert a CSV file to a Markdown table](#convert-a-csv-file-to-a-markdown-table)
  * [Convert PDF files to PNG files](#convert-pdf-files-to-png-files)
  * [Convert PNG files to a single PDF file](#convert-png-files-to-a-single-pdf-file)
  * [Convert a DOCX file to a PDF file](#convert-a-docx-file-to-a-pdf-file)
  * [Convert an HTML file to a PDF file](#convert-an-html-file-to-a-pdf-file)
  * [Convert a website to a PDF file](#convert-a-website-to-a-pdf-file)
  * [Convert an HTML file to a PNG file](#convert-an-html-file-to-a-png-file)
  * [Crop an image and add a white border](#crop-an-image-and-add-a-white-border)
  * [Resize an image](#resize-an-image)
  * [Format a CSV file into columns and examine its content](#format-a-csv-file-into-columns-and-examine-its-content)
  * [Format code](#format-code)
  * [Change Bash prompt temporarily](#change-bash-prompt-temporarily)
  * [Check bash/sh shell scripts for potential issues](#check-bashsh-shell-scripts-for-potential-issues)
  * [Take a screenshot of a window on macOS](#take-a-screenshot-of-a-window-on-macos)
  * [Create PowerPoint slides from a Markdown file](#create-powerpoint-slides-from-a-markdown-file)
  * [Run commands at scheduled times using cron](#run-commands-at-scheduled-times-using-cron)
  * [Record your terminal to an animated GIF](#record-your-terminal-to-an-animated-gif)
  * [Create an animated GIF from a YouTube video](#create-an-animated-gif-from-a-youtube-video)
  * [Create a collection of MP3 files from a YouTube playlist](#create-a-collection-of-mp3-files-from-a-youtube-playlist)
- [Perl](#perl)
  * [Get a random sample of lines from a text file while excluding the header line](#get-a-random-sample-of-lines-from-a-text-file-while-excluding-the-header-line)
  * [Convert a FASTA file to a CSV file with column names](#convert-a-fasta-file-to-a-csv-file-with-column-names)
  * [Count the number of lines that match a regular expression](#count-the-number-of-lines-that-match-a-regular-expression)
  * [Extract FASTA sequences from a file based on a file of sequence names of interest](#extract-fasta-sequences-from-a-file-based-on-a-file-of-sequence-names-of-interest)
  * [Add a FASTA title to the start of a sequence in RAW format](#add-a-fasta-title-to-the-start-of-a-sequence-in-raw-format)
  * [Remove commas located within quoted fields in a CSV file and create a tab-delimited file](#remove-commas-located-within-quoted-fields-in-a-csv-file-and-create-a-tab-delimited-file)
  * [Replace tabs with commas and remove quotes in a CSV file](#replace-tabs-with-commas-and-remove-quotes-in-a-csv-file)
  * [Sort sections in a Markdown file based on headings](#sort-sections-in-a-markdown-file-based-on-headings)
  * [Format Perl code](#format-perl-code)
- [Process multiple files](#process-multiple-files)
  * [for loop](#for-loop)
  * [while loop](#while-loop)
  * [find with -exec](#find-with--exec)
  * [find with xargs](#find-with-xargs)
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
- [sbatch](#sbatch)
  * [Count lines in compressed fastq files](#count-lines-in-compressed-fastq-files)
- [sed](#sed)
  * [Print a specific line of a file](#print-a-specific-line-of-a-file)
  * [Change filenames using a regular expression](#change-filenames-using-a-regular-expression)
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
- [tmux](#tmux)
  * [Start a tmux session](#start-a-tmux-session)
  * [Detach a tmux session](#detach-a-tmux-session)
  * [List tmux sessions](#list-tmux-sessions)
  * [Join a tmux session](#join-a-tmux-session)
  * [Create a multi-pane tmux session](#create-a-multi-pane-tmux-session)
  * [Navigate between tmux panes](#navigate-between-tmux-panes)
  * [Kill a tmux session](#kill-a-tmux-session)
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

### Convert a CSV file to a FASTA file

In this example column **1** contains the sequence title and column **3** contains the sequence:

```bash
awk -F, '{print ">"$1"\n"$3"\n"}' input.csv > output.fasta
```

### Print lines in file when a certain column contains a specific value

In this example lines are printed when the value in column **1** equals **9913**:

```bash
awk -F, '{if ($1 == 9913) print $0}' input.csv > output.csv
```

### Replace certain values in specific columns

In this example **1** and **-1** in column **23** are replaced with **forward** and **reverse**, respectively:

```bash
awk -F\\t 'BEGIN {OFS = "\t"} {sub(/^1/, "forward", $23); sub(/^-1/, "reverse", $23); print}' input.tab > output.tab
```

### Sum values in one column based on categories given in another column

In this example values in column **2** are added up for each category in column **1**:

```bash
awk -F, '{a[$1]+=$2}END{for(i in a) print i,a[i]}' input.csv
```

### Print column names and numbers

In this example the first row of the input file contains the column names:

```bash
awk -F $'\t' 'NR>1{exit};{for (i = 1; i <= NF; i++) print "column " i,"is " $i}' input.tab
```

### Print the values observed in a specific column, along with the number of times each value is observed 

In this example the counts for each distinct value in column **9** are printed:

```bash
awk -F $'\t' '{count[$9]++}END{for(j in count) print j,"("count[j]" counts)"}' input.tab
```

### Print the number of lines exhibiting each distinct number of fields

```bash
awk -F $'\t' '{count[NF]++}END{for(j in count) print "line length " j,"("count[j]" counts)"}' input.tab
```

### Print lines where certain fields contain values of interest

In this example lines where column **2** equals **7** and column **3** is between **60240145** and **60255062** are printed:

```bash
awk -F, '{ if ($2 == 7 && $3 >= 60240145 && $3 <= 60255062) print $0 }' input.csv
```

### Write each row to a separate file named after the value in a specific column

In this example each file is named after the value in column **1**:

```bash
awk -F '\t' '{ fname = $1 ".txt"; print >>fname; close(fname) }' input.tab
```

### Split a multi-FASTA file into separate files named according to the sequence title

In this example the sequences are written to a directory called **out**:

```bash
outputdir=out/
mkdir -p "$outputdir"
awk '/^>/ {OUT=substr($0,2); split(OUT, a, " "); sub(/[^A-Za-z_0-9\.\-]/, "", a[1]); OUT = "'"$outputdir"'" a[1] ".fa"}; OUT {print >>OUT; close(OUT)}' input.fasta
```

### Print only specific columns, identified by name in the first row

In this example the columns named **Affy SNP ID** and **Flank** are printed:

```bash
awk -F, 'NR==1 { for (i=1; i<=NF; i++) { ix[$i] = i } } NR>1 { print $ix["Affy SNP ID"]","$ix["Flank"] }' input.csv > output.csv
```

### Print only the lines coming after a certain starting line and before a certain ending line

In this example the lines coming after a line starting with **IlmnID** and before a line starting with **[Controls]** are printed:

```bash
awk -F, '/^IlmnID/{flag=1;print;next}/^\[Controls\]/{flag=0}flag' input.csv > output.csv
```

### Print only the first non-commented line

```bash
awk '/^[^#]/ { print $0;exit; }' input
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

In this example **parallel**:

```bash
brew install parallel
```

### Add a third-party repository

In this example **brewsci/bio** for bioinformatics software:

```bash
brew tap brewsci/bio
```

### Install directly from a third-party repository

In this example **clustal-w** from **brewsci/bio**:

```bash
brew install brewsci/bio/clustal-w
```

### View packages available from brewsci/bio

- [https://github.com/brewsci/homebrew-bio/tree/develop/Formula](https://github.com/brewsci/homebrew-bio/tree/develop/Formula)

### List installed graphical applications

```bash
brew cask list
```

### View available graphical applications

To view graphical applications available from the cask tap via the Homebrew package manager for macOS:

- [https://formulae.brew.sh/cask/](https://formulae.brew.sh/cask/)

### Install a graphical application 

In this example the Firefox browser:

```bash
brew cask install firefox
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

### Create an environment and install some programs

In this example an environment called **ngs** is created:

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

### Add additional programs to an environment

```bash
conda activate ngs
conda install -y -c bioconda -c conda-forge picard
```

### List environments

```bash
conda info --envs
```

## Docker

### Perform a sequence comparison using legacy BLAST

Download the legacy BLAST Docker image:

```bash
docker pull quay.io/biocontainers/blast-legacy:2.2.26--2
```

Create a container from the image and run **formatdb** to create a formatted database. In this example the database is created from a DNA sequence file called `sequence.fasta`, located in the current directory:

```bash
docker run -it --rm -v $(pwd):/directory -w /directory quay.io/biocontainers/blast-legacy:2.2.26--2 formatdb -i sequence.fasta -p F
```

To perform a blastn search using the formatted database and a query called `query.fasta` when the file is also located in the current directory:

```bash
docker run -it --rm -v $(pwd):/directory -w /directory quay.io/biocontainers/blast-legacy:2.2.26--2 blastall -p blastn -d sequence.fasta -i query.fasta
```

To perform a blastn search using the formatted database and a query called `query.fasta` when the query is located in a different directory (in this example your home directory):

```bash
docker run -it --rm -v $(pwd):/directory/database -v ${HOME}:/directory/query -w /directory quay.io/biocontainers/blast-legacy:2.2.26--2 blastall -p blastn -d database/sequence.fasta -i query/query.fasta
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

Create a container from the image and run **prokka** to annotate the sequence. In this example the genome sequence to be annotated is in a file called `sequence.fasta`, located in the current directory, and four CPUs are used:

```bash
docker run --rm -v $(pwd):/dir -w /dir staphb/prokka:latest prokka sequence.fasta --genus --cpus 4
```

## find

### Perform a series of commands on files returned by find

In this example `$'...'` is used for quoting, as it can contain escaped single quotes, and **tail** is used to skip a header line, **awk** is used to count the number of occurrences of each category in column 3 and print the category and counts, and **sort** is used to sort the categories by count from largest to smallest with ties broken by sorting on category name:

```bash
find . -type f -name "*.gff" -print0 | xargs -0 -I{} sh -c $'tail -n +2 "$1" | awk -F $\'\t\' \'{count[$3]++}END{for(j in count) print j,count[j]}\' | sort -k 2,2nr -k 1,1> "$1.cog_counts.txt"' -- {}
```

### Switch to the directory containing each file and execute a command

The -execdir option instructs **find** to switch to the directory containing each matching file before executing the specified command. In this example the command creates a **.zip** file for each **.vcf** file that is found:

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

First, copy the clone URL on the GitHub repository page by clicking **Clone or Download**. Then, enter the following command in a terminal window. The helpful\_commands repository is used as an example:

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

The following changes can be made to commits that have **not** been pushed to a remote repository:
To rewrite the very last commit, with any currently staged changes:

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

For example, to push to the master branch:

```bash
git push -u origin master
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

To merge a branch into Master (local) and push the changes to remote:

```bash
git checkout master
git merge <new-branch>
git push -u origin master
```

Git merge conflicts can arise easily. For information on resolving a merge conflict, see [Resolving a merged conflict using the command line](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/resolving-a-merge-conflict-using-the-command-line)

### Specify files to ignore

Create a .gitignore file:

```bash
touch .gitignore
```

Add text or patterns to exclude:

```bash
echo sensitive_data.txt >> .gitignore
echo test/*.vcf >> .gitignore
git add .gitignore
git commit -m "add .gitignore file"
git push -u origin master
```

In this example, the following files will no longer be tracked: `sensitive_data.txt`, and all files with a **.vcf** extension in the directory `test`.

Note that adding a .gitignore file will not remove tracked files; this must be done with `git rm`. See [Removing files from the repository](#removing-files-from-the-repository)

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

In this example the number of lines with a match to **>** is returned:

```bash
grep -c ">" input.fasta
```

### Get the line number of a match

In this example the line numbers of lines with a match to **234829** are reported:

```bash
grep -n "234829" input.txt
```

### Remove files that contain a match

In this example **.fasta** files are removed that contain the text **complete genome** on a single line:

```bash
grep -l "complete genome" *.fasta | xargs -I{} rm -f {}
```

### Remove files that do not contain a match

In this example **.fasta** files are removed that do not contain the text **complete genome** on a single line:

```bash
grep -L "complete genome" *.fasta | xargs -I{} rm -f {}
```

### Remove lines that match

Keep everything except lines starting with **#**:

```bash
grep -v '^#' input.txt > output.txt
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

Copy the public key to the `.ssh/authorized_keys` file on the other system using **ssh-copy-id**:

```bash
ssh-copy-id -i ~/.ssh/id_rsa.pub user@remote-host.com
```

### Download files from an FTP server

Replace `host`, `account`, `password`, and `port` with their corresponding values in the following command:

```bash
wget -S -d -c -t 45 -v -r ftp://account:password@host:port/*
```

### Download files from Google Drive

Rclone can be used to download data from many cloud storage providers.

First, follow the [configuration instructions](https://rclone.org/drivei). The commands below assume that the remote system was named **my_google_drive** during the configuration.

Note that you can use the **Add shortcut to Drive** option in Google Drive to make folders and files in **Shared with me** easier to access using rclone.

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

### Combine the columns in two tab-delimited files

```bash
paste -d"\t" input1.tab input2.tab > output.tab
```

### Add a header to all files with a certain extension, getting the header from another file

In this example the header is added to **.tab** files and comes from a file called `header.txt`. The files with the header added are saved with a **.new** extension added:

```bash
for f in *.tab; do new=`echo $f | sed 's/\(.*\)\.tab/\1.tab.new/'`; paste -sd'\n' \header.txt "$f" > "$new"; done
```

To replace the **.tab** files the **.new** files:

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

The following changes the **.gbff** extension to **.gbk**:

```bash
for f in *.gbff; do 
    mv -- "$f" "${f%.gbff}.gbk"
done
```

If **rename** is available, this may work:

```bash
rename 's/\.gbff$/.gbk/' *.gbff 
```

Or this, depending on which **rename** is installed:

```bash
rename .gbff .gbk *.gbff 
```

### Add text or a header to the beginning of all files with a particular file extension

Add **my header text** to the start of all **.csv** files in the current directory (works on macOS):

```bash
find . -name "*.csv" -exec sed -i '.bak' '1s/^/my header text\'$'\n/g' {} \;
```

Add **my header text** to the start of all **.csv** files in the current directory (works on Linux):

```bash
find . -name "*.csv" -exec sed -i '1s/^/my header text\n/' {} \;
```

### Find common lines between files

Between two files, named **file1.txt** and **file2.txt**:

```bash
comm -12 <( sort file1.txt ) <( sort file2.txt )
```

Among all **.txt** files in the current directory:

```bash
number_of_files=$(find . -name "*.txt" -print | wc -l | sed 's/[^0-9]*//g')
cat *.txt | sort | uniq -c | sed -n -e "s/^ *$number_of_files \(.*\)/\1/p"
```

### Convert a CSV file to a Markdown table

The following uses **csv2md**. The **awk** command can be used if some rows of the input have missing fields on the end:

```bash
awk -F, -v OFS="," 'NR==1 {cols=NF} {$1=$1; for (i=NF+1; i <= cols; i++) $i = "."} 1' input.csv > temp.csv
csv2md -p < temp.csv > output.md
```

### Convert PDF files to PNG files

The following uses **find** and the **pdftoppm** command from the **poppler** package to generate a PNG image of the first page of every PDF file in the working directory:

```bash
find . -name "*.pdf" -exec pdftoppm -f 1 -l 1 -png {} {} \;
```

### Convert PNG files to a single PDF file

The following uses **imagemagick**:

```bash
convert *.png output.pdf
```

### Convert a DOCX file to a PDF file

The following uses **LibreOffice**:

```bash
soffice --headless --convert-to pdf --outdir . word_file.docx
```

The following uses **pandoc** and on macOS also requires **basictex**:

```bash
pandoc word_file.docx --output word_file.pdf
```

### Convert an HTML file to a PDF file

The following uses **wkhtmltopdf**:

```bash
wkhtmltopdf http://google.com google.pdf
```

### Convert a website to a PDF file

The following uses **wkhtmltopdf** and **gs**:

```bash
url=http://www.3rs-reduction.co.uk/html/main_menu.html; depth=1
wget --spider --force-html -r -l${depth} ${url} 2>&1 | grep '^--' | awk '{ print $3 }' | grep -i '\.\(html\|htm\)$' | uniq > url-list.txt
while read i; do wkhtmltopdf "$i" "$(echo "$i" | sed -e 's/https\?:\/\///' -e 's/\//-/g' ).pdf"; done < url-list.txt
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=merged-output.pdf $(ls -lrt -1 *.pdf)
```

### Convert an HTML file to a PNG file

The following uses **wkhtmltoimage**:

```bash
wkhtmltoimage -f png http://google.com google.png
```

Another approach, which may work better for complex web sites, is to use **pageres**. The following creates an image that is 4485 pixels wide (5 * 897):

```bash
pageres http://google.com 897x1090 --crop --scale=5 --filename='google'
```

### Crop an image and add a white border

The following uses **imagemagick** to removes any edges that are exactly the same color as the corner pixels. A 30-pixel white border is then added to the sides and the top and bottom:

```bash
convert input.png -trim -bordercolor White -border 30x30 output.png
```

### Resize an image 

The following uses **imagemagick** to scale the image so that its width is 4000 pixels:

```bash
convert input.png -resize 4000 output.png
```

### Format a CSV file into columns and examine its content

```bash
cat data.csv | perl -pe 's/((?<=,)|(?<=^)),/ ,/g;' | column -t -s, | less -S
```

### Format code

The Prettier program supports many languages. The following command uses the `--write` option to reformat files in-place:

```bash
prettier --write "*html"
```

### Change Bash prompt temporarily

```bash
PS1="$ "
```

### Check bash/sh shell scripts for potential issues

```bash
shellcheck some_script.sh
```

### Take a screenshot of a window on macOS

1. Press **Command+Shift+4**.
2. Press the **space bar**.
3. Hold **Option** and click on the window.

To keep the drop shadow perform the last step without holding **Option**. 

### Create PowerPoint slides from a Markdown file

Pandoc can be used to generate PowerPoint slides. The following Markdown text describes several slides with notes:

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

The formatting of individual slides can then be adjusted in PowerPoint, using the **Design** tab and the **Design Ideas** button. Slide numbers and headers and footers can be added using **View->Slide Master** followed by **Insert**, and then **Header & Footer**.

To reduce the size of the file, use **File->Compress Pictures...**.

### Run commands at scheduled times using cron

The following uses **cron** to run a script to copy various files and directories to a directory backed up by Dropbox.

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

Use Terminalizer to record a session. The following creates a file called `session.yml`:

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

The following requires **youtube-dl**, **mplayer**, **imagemagick**, and **gifsicle**:

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

The following requires **youtube-dl** and **ffmpeg**:

```bash
youtube-dl -x -i --audio-format mp3 --audio-quality 320K --embed-thumbnail --geo-bypass https://www.youtube.com/playlist?list=PL92319EECC1754042
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

In this example the title **>KL1** is added to the beginning of the sequence in `KL1sequence.txt`:

```bash
perl -pi -e 'print ">KL1\n" if $. == 1' KL1sequence.txt
```

### Remove commas located within quoted fields in a CSV file and create a tab-delimited file

```bash
perl -nle  'my @new  = (); push( @new, $+ ) while $_ =~ m{"([^\"\\]*(?:\\.[^\"\\]*)*)",? | ([^,]+),? | ,}gx; push( @new, undef ) if substr( $text, -1, 1 ) eq '\'','\''; for(@new){s/,/ /g} print join "\t", @new' input.csv > output.tab
```

### Replace tabs with commas and remove quotes in a CSV file

```bash
perl -p -e 's/\t/,/g;' -e 's/"//g' input.csv > output.csv
```

### Sort sections in a Markdown file based on headings

```bash
perl -0777 -ne '(undef,@paragraphs) = split /^#(?=[^#])/m; print map {"#$_"} sort { "\U$a" cmp "\U$b" } @paragraphs;' input.md
```

### Format Perl code

The following uses **perltidy** to reformat the code in `testfile.pl` and will create a file called `testfile.pl.tdy`.

```bash
perltidy testfile.pl
```

## Process multiple files

### for loop

Change all **.fasta** files in the current directory to **.fna** files:

```bash
for f in *.fasta; do new=`echo $f | sed 's/\(.*\)\.fasta/\1.fna/'`; mv "$f" "$new"; done
```

### while loop

Print the number of lines in every **.csv** or **.tab** file in or below current directory:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) | while read f; do wc -l "$f"; done
```

Print the number of lines in every **.csv** or **.tab** file in or below current directory and redirect the results to a single file:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) | while read f; do wc -l "$f" >> output.txt; done
```

Print the number of lines in every **.csv** or **.tab** file in or below current directory and redirect the results to separate files:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) | while read f; do wc -l "$f" > "${f}.output.txt"; done
```

### find with -exec

Change all **.fasta** files in current directory to **.fna** files by appending a **.fna** extension:

```bash
find . -type f -name "*.fasta" -exec mv {} {}.fna \;
```

Print the number of lines in every **.csv** or **.tab** file in or below current directory:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -exec wc -l {} \;
```

Print the number of lines in every **.csv** or **.tab** file in or below current directory and redirect the results to a single file:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -exec wc -l {} \; > output.txt
```

Print the number of lines in every **.csv** or **.tab** file in or below current directory and redirect the results to separate files:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -exec sh -c 'wc -l "$1" > "$1.output.txt"' -- {} \;
```

### find with xargs

Change all **.fasta** files in current directory to **.fna** files by appending a **.fna** extension:

```bash
find . -type f -name "*.fasta" -print0 | xargs -0 -I{} mv {} {}.fna
```

Print the number of lines in every **.csv** or **.tab** file in or below current directory:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -print0 | xargs -0 -I{} wc -l {}
```

Print the number of lines in every **.csv** or **.tab** file in or below current directory and redirect the results to a single file:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -print0 | xargs -0 -I{} wc -l {} > output.txt
```

Print the number of lines in every **.csv** or **.tab** file in or below current directory and redirect the results to separate files:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -print0 | xargs -0 -I{} sh -c 'wc -l "$1" > "$1.output.txt"' -- {}
```

Print the number of lines in every **.csv** or **.tab** file in or below current directory and redirect the results to separate files. Process up to **4** files in parallel:

```bash
find . -type f \( -name "*.csv" -o -name "*.tab" \) -print0 | xargs -n1 -P4 -0 -I{} sh -c 'wc -l "$1" > "$1.output.txt"' -- {}
```

## R

### Compare two data sets to find differences

In this example SNP location information assigned by two different algorithms is compared using the **compareDF** package. The two data sets (position information) generated by the differing algorithms are read from files. SNP name is used to match rows across the data sets, and then the 'chromosome' and 'position' are compared. SNP records for which 'chromosome' or 'position' differ between the data sets are written to an Excel file:

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

In this example, an UpSet plot is used to visualize the overlap among all combinations of gene lists in the **gene_lists** directory. In this directory each list is given as a separate **.txt** file, with a single header row and one gene name or ID per row, for example:

```
Gene name or identifier
ENSG00000264954.2
ENSG00000224383.8
CCDS54157.1.1
```

The UpSet plot is generated using the **UpSetR** package:

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

In this example a heatmap is used to visualize gene presence and absence for all gene lists in the **gene_lists** directory. In this directory each list is given as a separate **.txt** file, with a single header row and one gene name or ID per row, for example:

```
Gene name or identifier
ENSG00000264954.2
ENSG00000224383.8
CCDS54157.1.1
```

The following uses the **purrr** and **RVenn** packages:

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

genotypes_one_column_per_snp <- separate(genotypes, col = alleles, sep = "(?<=\\d)", into = snps)

#change both instances of 'sample' to match name of first column
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

## sbatch

### Count lines in compressed fastq files

In this example the number of lines in several **.fastq.gz** files is quickly determined by submitting jobs to Slurm using sbatch.

The naming scheme of the **.fastq.gz** files is as follows (the sample name is in the file name, for example **DG15B032198-1**):

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

Then use the following commands to submit a job for each **.fastq.gz** file:

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

The **.out** files will contain the name of the input file and the number of lines, for example:

```
HI.5173.001.NEBNext_Index_12.DG15B032198-1_R1.fastq.gz 229623444
```

To quickly check the **.err** files:

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

### Print a specific line of a file

In this example line **26404**:

```bash
sed -n "26404p" input.txt
```

### Change filenames using a regular expression

In this example **chr30** is replaced with **chrX**:

```bash
for f in *.fasta; do new=`echo $f | sed 's/chr30/chrX/'`; mv $f $new; done
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
salloc --time=2:0:0 --ntasks=1 --mem-per-cpu=2000M --account=def-someuser
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

The following also creates a three-pane tmux session, and launches **vim** in the largest pane:

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

## vcftools and bcftools

### Extract variants from a region of interest and write to a new vcf file

Note that if the vcf file is gzip compressed (i.e. has a **.gz** extension), use `--gzvcf` instead of `--vcf`.

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

In this example search and replace operations are performed on all the **.html** files in a directory. First, open the files in multiple buffers in vim:

```bash
vim *.html
```

Then use **argdo** to perform a search and replace across all the files. In this example blank lines are removed:

```
:argdo %s/^$//ge
```

In this example the text between `<p class="lastupdated">` and `</p>` are replaced with the current date. Note the use of **\zs** and **\ze** so that the text between those tags is replaced and not the tags themselves:

```
:argdo %s/<p class="lastupdated">\zs[^<]*\ze<\/p>/\=strftime("%c")/ge
```

### Search and replace newlines

In replacement syntax use **\r** instead of **\n** to represent newlines. For example, to replace commas with newlines:

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
