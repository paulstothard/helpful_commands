# helpful_commands

Command-line tools and commands for performing a variety of tasks, several related to bioinformatics.

## Processing multiple files

### Using loops

Change all **.fasta** files in the current directory to **.fna** files:

```bash
for f in *.fasta; do new=`echo $f | sed 's/\(.*\)\.fasta/\1.fna/'`; mv "$f" "$new"; done
```

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

### Using find with -exec

Change all **.fasta** files in current directory to **.fna** files:

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

### Using find with xargs

Change all **.fasta** files in current directory to **.fna** files:

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

Count matches. In this example the number of lines with a match to **>** is returned:

```bash
grep -c ">" input.fasta
```

Get the line number of a match. In this example the line numbers of lines with a match to **234829** are reported:

```bash
grep -n "234829" input.txt
```

Remove files that contain a match. In this example **.fasta** files are removed that contain the text **complete genome** on a single line:

```bash
grep -l "complete genome" *.fasta | xargs -I{} rm -f {}
```

Remove files that do not contain a match. In this example **.fasta** files are removed that do not contain the text **complete genome** on a single line:

```bash
grep -L "complete genome" *.fasta | xargs -I{} rm -f {}
```

Keep everything except lines starting with **#**:

```bash
grep -v '^#' input.txt > output.txt
```

## awk

Convert a CSV file to a FASTA file. In this example column **1** contains the sequence title and column **3** contains the sequence:

```bash
awk -F, '{print ">"$1"\n"$3"\n"}' input.csv > output.fasta
```

Print lines in file when a certain column contains a specific value. In this example lines are printed when the value in column **1** equals **9913**:

```bash
awk -F, '{if ($1 == 9913) print $0}' input.csv > output.csv
```

Replace certain values in specific columns. In this example **1** and **-1** in column **23** are replaced with **forward** and **reverse**, respectively:

```bash
awk -F\\t 'BEGIN {OFS = "\t"} {sub(/^1/, "forward", $23); sub(/^-1/, "reverse", $23); print}' input.tab > output.tab
```

Sum values in one column based on categories given in another column. In this example values in column **2** are added up for each category in column **1**:

```bash
awk -F, '{a[$1]+=$2}END{for(i in a) print i,a[i]}' input.csv
```

Print column names and numbers. In this example the first row of the input file contains the column names:

```bash
awk -F $'\t' 'NR>1{exit};{for (i = 1; i <= NF; i++) print "column " i,"is " $i}' input.tab
```

Print the types of values observed in a specific column, along with the number of times each type is observed. In this example the counts for each distinct value in column **9** are printed:

```bash
awk -F $'\t' '{count[$9]++}END{for(j in count) print j,"("count[j]" counts)"}' input.tab
```

Print the number of lines exhibiting each distinct number of fields:

```bash
awk -F $'\t' '{count[NF]++}END{for(j in count) print "line length " j,"("count[j]" counts)"}' input.tab
```

Print lines where certain fields contain values of interest. In this example lines where column **2** equals **7** and column **3** is between **60240145** and **60255062** are printed:

```bash
awk -F, '{ if ($2 == 7 && $3 >= 60240145 && $3 <= 60255062) print $0 }' input.csv
```

Write each row to a separate file named after the value in a specific column. In this example each file is named after the value in column **1**:

```bash
awk -F '\t' '{ fname = $1 ".txt"; print >>fname; close(fname) }' input.tab
```

Split a multi-FASTA file into separate files, one per sequence named according to the sequence title. In this example the sequences are written to a directory called **out**:

```bash
outputdir=out/
mkdir -p "$outputdir"
awk '/^>/ {OUT=substr($0,2); split(OUT, a, " "); sub(/[^A-Za-z_0-9\.\-]/, "", a[1]); OUT = "'"$outputdir"'" a[1] ".fa"}; OUT {print >>OUT; close(OUT)}' input.fasta
```

Print only specific columns, identified by name in the first row. In this example the columns named **Affy SNP ID** and **Flank** are printed:

```bash
awk -F, 'NR==1 { for (i=1; i<=NF; i++) { ix[$i] = i } } NR>1 { print $ix["Affy SNP ID"]","$ix["Flank"] }' input.csv > output.csv
```

Print only the lines coming after a certain starting line and before a certain ending line. In this example the lines coming after a line starting with **IlmnID** and before a line starting with **[Controls]** are printed:

```bash
awk -F, '/^IlmnID/{flag=1;print;next}/^\[Controls\]/{flag=0}flag' input.csv > output.csv
```

## sed

Print a specific line of a file. In this example line **26404**:

```bash
sed -n "26404p" input.txt
```

Change filenames using a regular expression. In this example **chr30** is replaced with **chrX**:

```bash
for f in *.fasta; do new=`echo $f | sed 's/chr30/chrX/'`; mv $f $new; done
```

## Perl

Get a random sample of lines from a text file when you don't want to include a header line. In this example a random sample of 20 lines is obtained:

```bash
tail -n +2 input.txt | perl -MList::Util -e 'print List::Util::shuffle <>' | head -n 20 > output.txt
```

Convert a FASTA file to a CSV file with column names:

```bash
cat input.fasta | perl -n -0777 -e 'BEGIN{print "SNP_Name,Sequence\n"}' -e 'while ($_ =~ m/^>([^\n]+)\n([^>]+)/gm) {$name = $1; $seq = $2; $seq =~s/\s//g; print $name . "," . $seq . "\n"}' > output.csv
```

Count the number of lines that match a regular expression:

```bash
perl -lne '$a++ if /\tyes\t/; END {print $a+0}' < input.txt
```

Extract FASTA sequences from a file based on a file of sequence names of interest. In this example the sequence names of interest are in the file **names.txt** and the FASTA sequences are in the file **input.fasta**:

```bash
cat names.txt | xargs -I{} perl -w -076 -e '$count = 0; open(SEQ, "<" . $ARGV[0]); while (<SEQ>) {if ($_ =~ m/\Q$ARGV[1]\E/) {$record = $_; $record =~ s/[\s>]+$//g; print ">$record\n"; $count = $count + 1;}} if ($count == 0) {print STDERR "No matches found for $ARGV[1]\n"} elsif ($count > 1) {print STDERR "Multiple matches found for $ARGV[1]\n"} close(SEQ);' input.fasta {} > output.fasta
```

Add a FASTA title to the start of a sequence in RAW format. In this example the title **>KL1** is added to the beginning of the sequence in **KL1sequence.txt**:

```bash
perl -pi -e 'print ">KL1\n" if $. == 1' KL1sequence.txt
```

Remove commas located within quoted fields in a CSV file and create a tab-delimited file:

```bash
perl -nle  'my @new  = (); push( @new, $+ ) while $_ =~ m{"([^\"\\]*(?:\\.[^\"\\]*)*)",? | ([^,]+),? | ,}gx; push( @new, undef ) if substr( $text, -1, 1 ) eq '\'','\''; for(@new){s/,/ /g} print join "\t", @new' input.csv > output.tab
```

Replace tabs with commas and remove quotes in a CSV file:

```bash
perl -p -e 's/\t/,/g;' -e 's/"//g' input.csv > output.csv
```

## Other

Combine the columns in two tab-delimited files:

```bash
paste -d"\t" input1.tab input2.tab > output.tab
```

Add a header to all files with a certain extension, getting the header from another file. In this example the header is added to **.tab** files and comes from a file called **header.txt**. The files with the header added are saved with a **.new** extension added:

```bash
for f in *.tab; do new=`echo $f | sed 's/\(.*\)\.tab/\1.tab.new/'`; paste -sd'\n' \header.txt "$f" > "$new"; done
```

View STDOUT and append it to a file:

```bash
some_command | tee -a output.txt
```

Redirect STDERR to STDOUT and view both and append both to a file:

```bash
some_command 2>&1 | tee -a log
```

## sbatch

### Counting lines in compressed fastq files

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

First create a sbatch script called **fastq.gz.lines.sbatch** to run zcat and wc (used to count lines in a compressed file):

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

The **out** file will contain the name of the input file and the number of lines, for example:

```
HI.5173.001.NEBNext_Index_12.DG15B032198-1_R1.fastq.gz 229623444
```

To quickly check the **err** files:

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

### General Slurm commands

View statistics related to the efficiency of resource usage of a completed job:

```bash
seff <jobid>
```

View jobs:

```bash
squeue -u <username>
```

View running jobs:

```bash
squeue -u <username> -t RUNNING
```

View pending jobs:

```bash
squeue -u <username> -t PENDING
```

View detailed information for a specific job:

```bash
scontrol show job -dd <jobid>
```

View accounting information for completed jobs:

```bash
sacct -s CD --format=JobID,JobName,MaxRSS,ReqMem,Elapsed,End,State,NodeList
```

Cancel a job:

```bash
scancel <jobid>
```

Cancel all jobs:

```bash
scancel -u <username>
```

Start an interactive session:

```bash
salloc --time=2:0:0 --ntasks=1 --mem-per-cpu=2000M --account=def-someuser
```

## Sharing data with project group members

```bash
cd projects/some_project
chmod g+x my_dir
cd my_dir
mkdir shared_dir
chmod g+x shared_dir
chmod +t shared_dir
```

## Using Conda to install NGS tools

Install Miniconda and create an environment called **ngs** with several NGS-related tools:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
source ~/miniconda3/bin/activate
conda init
source ~/.bashrc
conda update -y -n base -c defaults conda
conda create -y --name ngs
conda activate ngs
conda install -y -c bioconda -c conda-forge multiqc fastqc trimmomatic bowtie2 subread samtools
```

To deactivate the environment:

```bash
conda deactivate
```

To activate the environment:

```bash
conda activate ngs
```

To add additional programs to the ngs environment (e.g. picard):

```bash
conda activate ngs
conda install -y -c bioconda -c conda-forge picard
```

## Running a program using Docker

In this example a Docker container is used to run legacy BLAST.

Download the legacy BLAST Docker image:

```bash
docker pull quay.io/biocontainers/blast-legacy:2.2.26--2
```

Create a container from the image and run formatdb to create a formatted database. In this example the database is created from a DNA sequence file called **sequence.fasta**, located in the current directory:

```bash
docker run -it --rm -v $(pwd):/directory -w /directory quay.io/biocontainers/blast-legacy:2.2.26--2 formatdb -i sequence.fasta -p F
```

To perform a blastn search using the formatted database and a query called **query.fasta** when the file is also located in the current directory:

```bash
docker run -it --rm -v $(pwd):/directory -w /directory quay.io/biocontainers/blast-legacy:2.2.26--2 blastall -p blastn -d sequence.fasta -i query.fasta
```

To perform a blastn search using the formatted database and a query called **query.fasta** when the query is located in a different directory (in this example your home directory):

```bash
docker run -it --rm -v $(pwd):/directory/database -v ${HOME}:/directory/query -w /directory quay.io/biocontainers/blast-legacy:2.2.26--2 blastall -p blastn -d database/sequence.fasta -i query/query.fasta
```

## Using brew to install software

For a listing of installed packages:

```bash
brew list
```

For a listing of all packages available from the core tap via the Homebrew package manager for macOS:

- [https://formulae.brew.sh/formula/](https://formulae.brew.sh/formula/)

To install a package, in this example **parallel**:

```bash
brew install parallel
```

To add a third-party repository, in this example **brewsci/bio** for bioinformatics software:

```bash
brew tap brewsci/bio
```

To install directly from a third-party repository, in this example **clustal-w** from **brewsci/bio**:

```bash
brew install brewsci/bio/clustal-w
```

For a listing of packages available from **brewsci/bio**:

- [https://github.com/brewsci/homebrew-bio/tree/develop/Formula](https://github.com/brewsci/homebrew-bio/tree/develop/Formula)

For a listing of installed graphical applications:

```bash
brew cask list
```

For a listing of all graphical applications available from the cask tap via the Homebrew package manager for macOS:

- [https://formulae.brew.sh/cask/](https://formulae.brew.sh/cask/)

To install a graphical application, in this example the Firefox browser:

```bash
brew cask install firefox
```

## Version Control with Git

See [Github's Git documentation](https://help.github.com/en) for more information

### Creating a new Git repository

```
git init
```

### Syncing a repository to your local machine

First, copy the clone URL on the Github repository page by clicking **Clone or Download**. Then, enter the following command in a terminal window. The helpful_commands repository is used as an example:

```
git clone https://github.com/stothard-group/helpful_commands.git
```

### Marking changed files to be included in the next commit

To add one or more files:

```
git add <filename1> <filename2>
```

To add all current modifications in your project (including deletions and new files):

```
git add --all
```

### Undoing a Git add **before** a commit

To undo a list of files:

```
git reset <filename1>
```

To undo all changes:

```
git reset
```

### Removing files from the repository

Note that the following instructions will remove the file/directory from both the working tree and the index.

To remove one or more files:

```
git rm <filename>
```

To remove a directory:

```
git rm -r <directory>
```

To remove a file from the index (this untracks the file, it does not delete the file itself):

```
git rm --cached <filename>
```

These changes must be committed with git commit.

### Moving or renaming a file or directory

```
git mv <filename-old> <filename-new>
```

This change must be committed with git commit.

### Saving the marked files to the local Git repository

The commit should include a message using the -m option:

```
git commit -m "A concise description of the changes"
```

The following changes can be made to commits that have **not** been pushed to a remote repository:
To rewrite the very last commit, with any currently staged changes:

```
git commit --amend -m "An updated message"
```

To commit any currently staged changes without rewriting the commit (this essentially adds the staged changes to the previous commit):

```
git commit --amend --no-edit
```

### To push a commit on your local branch to a remote repository

```
git push <remote> <branch>
```

For example, to push to the master branch:

```
git push -u origin master
```

### Adding or editing a remote repository

To add a new remote:

```
git remote add origin <repo-url>
```

To edit an existing remote:

```
git remote set-url origin <new-repo-url>
```

To verify that the remote URL has changed:

```
git remote -v
```

### Creating and merging Git branches

To view the branches in a repository:

```
git branch
```

To create a new branch and switch to it:

```
git checkout -b <new-branch>
```

To switch to a remote branch:

```
git fetch --all
git checkout <remote-branch>
```

After adding and committing some changes, to push this branch to remote:

```
git push -u origin <new-branch>
```

To merge a branch into Master (local) and push the changes to remote:

```
git checkout master
git merge <new-branch>
git push -u origin master
```

Git merge conflicts can arise easily. For information on resolving a merge conflict, see [Resolving a merged conflict using the command line](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/resolving-a-merge-conflict-using-the-command-line)

### Specifying intentionally untracked files to ignore

Create a .gitignore file:

```
touch .gitignore
```

Add text or patterns to exclude:

```
echo sensitive_data.txt >> .gitignore
echo test/*.vcf >> .gitignore
git add .gitignore
git commit -m "add .gitignore file"
git push -u origin master
```

In this example, the following files will no longer be tracked: `sensitive_data.txt`, and all files with a .vcf extension in the directory `test`.

Note that adding a .gitignore file will not remove tracked files; this must be done with `git rm`. See [Removing files from the repository](#removing-files-from-the-repository)

### Checking the status of a working directory

```
git status
```
