# Table of Contents

<!-- toc -->

- [Process multiple files](#process-multiple-files)
  * [for loop](#for-loop)
  * [while loop](#while-loop)
  * [find with -exec](#find-with--exec)
  * [find with xargs](#find-with-xargs)
- [grep](#grep)
  * [Count matches](#count-matches)
  * [Get the line number of a match](#get-the-line-number-of-a-match)
  * [Remove files that contain a match](#remove-files-that-contain-a-match)
  * [Remove files that do not contain a match](#remove-files-that-do-not-contain-a-match)
  * [Remove lines that match](#remove-lines-that-match)
- [awk](#awk)
  * [Convert a CSV file to a FASTA file](#convert-a-csv-file-to-a-fasta-file)
  * [Print lines in file when a certain column contains a specific value](#print-lines-in-file-when-a-certain-column-contains-a-specific-value)
  * [Replace certain values in specific columns](#replace-certain-values-in-specific-columns)
  * [Sum values in one column based on categories given in another column](#sum-values-in-one-column-based-on-categories-given-in-another-column)
  * [Print column names and numbers](#print-column-names-and-numbers)
  * [Print the types of values observed in a specific column, along with the number of times each type is observed](#print-the-types-of-values-observed-in-a-specific-column-along-with-the-number-of-times-each-type-is-observed)
  * [Print the number of lines exhibiting each distinct number of fields](#print-the-number-of-lines-exhibiting-each-distinct-number-of-fields)
  * [Print lines where certain fields contain values of interest](#print-lines-where-certain-fields-contain-values-of-interest)
  * [Write each row to a separate file named after the value in a specific column](#write-each-row-to-a-separate-file-named-after-the-value-in-a-specific-column)
  * [Split a multi-FASTA file into separate files named according to the sequence title](#split-a-multi-fasta-file-into-separate-files-named-according-to-the-sequence-title)
  * [Print only specific columns, identified by name in the first row](#print-only-specific-columns-identified-by-name-in-the-first-row)
  * [Print only the lines coming after a certain starting line and before a certain ending line](#print-only-the-lines-coming-after-a-certain-starting-line-and-before-a-certain-ending-line)
- [sed](#sed)
  * [Print a specific line of a file](#print-a-specific-line-of-a-file)
  * [Change filenames using a regular expression](#change-filenames-using-a-regular-expression)
- [Perl](#perl)
  * [Get a random sample of lines from a text file while excluding the header line](#get-a-random-sample-of-lines-from-a-text-file-while-excluding-the-header-line)
  * [Convert a FASTA file to a CSV file with column names](#convert-a-fasta-file-to-a-csv-file-with-column-names)
  * [Count the number of lines that match a regular expression](#count-the-number-of-lines-that-match-a-regular-expression)
  * [Extract FASTA sequences from a file based on a file of sequence names of interest](#extract-fasta-sequences-from-a-file-based-on-a-file-of-sequence-names-of-interest)
  * [Add a FASTA title to the start of a sequence in RAW format](#add-a-fasta-title-to-the-start-of-a-sequence-in-raw-format)
  * [Remove commas located within quoted fields in a CSV file and create a tab-delimited file](#remove-commas-located-within-quoted-fields-in-a-csv-file-and-create-a-tab-delimited-file)
  * [Replace tabs with commas and remove quotes in a CSV file](#replace-tabs-with-commas-and-remove-quotes-in-a-csv-file)
- [find](#find)
  * [Perform a series of commands on files returned by find](#perform-a-series-of-commands-on-files-returned-by-find)
- [Other](#other)
  * [Combine the columns in two tab-delimited files](#combine-the-columns-in-two-tab-delimited-files)
  * [Add a header to all files with a certain extension, getting the header from another file](#add-a-header-to-all-files-with-a-certain-extension-getting-the-header-from-another-file)
  * [View STDOUT and append it to a file](#view-stdout-and-append-it-to-a-file)
  * [Redirect STDERR to STDOUT and view both and append both to a file](#redirect-stderr-to-stdout-and-view-both-and-append-both-to-a-file)
  * [Convert pdf files to png](#convert-pdf-files-to-png)
- [sbatch](#sbatch)
  * [Count lines in compressed fastq files](#count-lines-in-compressed-fastq-files)
- [General Slurm commands](#general-slurm-commands)
  * [View statistics related to the efficiency of resource usage of a completed job](#view-statistics-related-to-the-efficiency-of-resource-usage-of-a-completed-job)
  * [View jobs](#view-jobs)
  * [View running jobs](#view-running-jobs)
  * [View pending jobs](#view-pending-jobs)
  * [View detailed information for a specific job](#view-detailed-information-for-a-specific-job)
  * [View accounting information for completed jobs](#view-accounting-information-for-completed-jobs)
  * [Cancel a job](#cancel-a-job)
  * [Cancel all jobs](#cancel-all-jobs)
  * [Start an interactive session](#start-an-interactive-session)
- [Share data with project group members](#share-data-with-project-group-members)
- [Use Conda to install NGS tools](#use-conda-to-install-ngs-tools)
  * [Install Miniconda](#install-miniconda)
  * [Create an environment and install some programs](#create-an-environment-and-install-some-programs)
  * [Deactivate an environment](#deactivate-an-environment)
  * [Activate an environment](#activate-an-environment)
- [Add additional programs to an environment](#add-additional-programs-to-an-environment)
- [Run a program using Docker](#run-a-program-using-docker)
- [Use brew to install software](#use-brew-to-install-software)
  * [List installed packages](#list-installed-packages)
  * [View available packages](#view-available-packages)
  * [Install a package](#install-a-package)
  * [Add a third-party repository](#add-a-third-party-repository)
  * [Install directly from a third-party repository](#install-directly-from-a-third-party-repository)
  * [View packages available from brewsci/bio](#view-packages-available-from-brewscibio)
  * [List installed graphical applications](#list-installed-graphical-applications)
  * [View available graphical applications](#view-available-graphical-applications)
  * [Install a graphical application](#install-a-graphical-application)
- [Version Control with Git](#version-control-with-git)
  * [Create a new Git repository](#create-a-new-git-repository)
  * [Sync a repository to your local machine](#sync-a-repository-to-your-local-machine)
  * [Mark changed files to be included in the next commit](#mark-changed-files-to-be-included-in-the-next-commit)
  * [Undo a Git add before a commit](#undo-a-git-add-before-a-commit)
  * [Remove files from the repository](#remove-files-from-the-repository)
  * [Move or rename a file or directory](#move-or-rename-a-file-or-directory)
  * [Save the marked files to the local Git repository](#save-the-marked-files-to-the-local-git-repository)
  * [Push a commit on your local branch to a remote repository](#push-a-commit-on-your-local-branch-to-a-remote-repository)
  * [Add or edit a remote repository](#add-or-edit-a-remote-repository)
  * [Creat and merge Git branches](#creat-and-merge-git-branches)
  * [Specify files to ignore](#specify-files-to-ignore)
  * [Check the status of a working directory](#check-the-status-of-a-working-directory)
  * [Tag a release](#tag-a-release)
- [vim](#vim)
  * [Search and replace across multiple files](#search-and-replace-across-multiple-files)
  * [Search and replace newlines](#search-and-replace-newlines)

<!-- tocstop -->

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

### Print the types of values observed in a specific column, along with the number of times each type is observed 

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

## find

### Perform a series of commands on files returned by find

In this example `$'...'` is used for quoting, as it can contain escaped single quotes, and **tail** is used to skip a header line, **awk** is used to count the number of occurrences of each category in column 3 and print the category and counts, and **sort** is used to sort the categories by count from largest to smallest with ties broken by sorting on category name:

```bash
find . -type f -name "*.gff" -print0 | xargs -0 -I{} sh -c $'tail -n +2 "$1" | awk -F $\'\t\' \'{count[$3]++}END{for(j in count) print j,count[j]}\' | sort -k 2,2nr -k 1,1> "$1.cog_counts.txt"' -- {}
```

## Other

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

### Convert pdf files to png

The following uses **find** and the **pdftoppm** command from the **poppler** package to generate a png image of the first page of every pdf file in the working directory:

```bash
find . -name "*.pdf" -exec pdftoppm -f 1 -l 1 -png {} {} \;
```

## sbatch

### Count lines in compressed fastq files

In this example, the number of lines in several **.fastq.gz** files is quickly determined by submitting jobs to Slurm using sbatch.

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

## General Slurm commands

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

## Share data with project group members

```bash
cd projects/some_project
chmod g+x my_dir
cd my_dir
mkdir shared_dir
chmod g+x shared_dir
chmod +t shared_dir
```

## Use Conda to install NGS tools

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

## Add additional programs to an environment

```bash
conda activate ngs
conda install -y -c bioconda -c conda-forge picard
```

## Run a program using Docker

In this example a Docker container is used to run legacy BLAST.

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

## Use brew to install software

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

## Version Control with Git

See [Github's Git documentation](https://help.github.com/en) for more information

### Create a new Git repository

```bash
git init
```

### Sync a repository to your local machine

First, copy the clone URL on the Github repository page by clicking **Clone or Download**. Then, enter the following command in a terminal window. The helpful_commands repository is used as an example:

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

### Creat and merge Git branches

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

Information on how to choose version numbers if availble [here](https://semver.org).

## vim

### Search and replace across multiple files

In this example search and replace operations are perfomed on all the **.html** files in a directory. First, open the files in multiple buffers in vim:

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
