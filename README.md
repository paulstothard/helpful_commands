# helpful_commands
Useful Bash commands for performing a variety of tasks, mostly related to bioinformatics

## Processing multiple files

### Using loops

Change all **.fasta** files in the current directory to **.fna** files:

```bash
for f in *.fasta; do new=`echo $f | sed 's/\(.*\)\.fasta/\1.fna/'`; mv "$f" "$new"; done
```

Print the number of lines in every **.html** or **.htm** file in or below current directory:

```bash
find . -type f \( -name "*.html" -o -name "*.htm" \) | while read f; do wc -l "$f"; done
```

Print the number of lines in every **.html** or **.htm** file in or below current directory and redirect the results to a single file:

```bash
find . -type f \( -name "*.html" -o -name "*.htm" \) | while read f; do wc -l "$f" >> output.txt; done
```

Print the number of lines in every **.html** or **.htm** file in or below current directory and redirect the results to separate files:

```bash
find . -type f \( -name "*.html" -o -name "*.htm" \) | while read f; do wc -l "$f" > "${f}.output.txt"; done
```

### Using find with -exec

Change all **.fasta** files in current directory to **.fna** files:

```bash
find . -type f -name "*.fasta" -exec mv {} {}.fna \;
```

Print the number of lines in every **.html** or **.htm** file in or below current directory:

```bash
find . -type f \( -name "*.html" -o -name "*.htm" \) -exec wc -l {} \;
```

Print the number of lines in every **.html** or **.htm** file in or below current directory and redirect the results to a single file:

```bash
find . -type f \( -name "*.html" -o -name "*.htm" \) -exec wc -l {} \; > output.txt
```

Print the number of lines in every **.html** or **.htm** file in or below current directory and redirect the results to separate files:

```bash
find . -type f \( -name "*.html" -o -name "*.htm" \) -exec sh -c 'wc -l "$1" > "$1.output.txt"' -- {} \;
```

### Using find with xargs

Change all **.fasta** files in current directory to **.fna** files:

```bash
find . -type f -name "*.fasta" -print0 | xargs -0 -I{} mv {} {}.fna
```

Print the number of lines in every **.html** or **.htm** file in or below current directory:

```bash
find . -type f \( -name "*.html" -o -name "*.htm" \) -print0 | xargs -0 -I{} wc -l {}
```

Print the number of lines in every **.html** or **.htm** file in or below current directory and redirect the results to a single file:

```bash
find . -type f \( -name "*.html" -o -name "*.htm" \) -print0 | xargs -0 -I{} wc -l {} > output.txt
```

Print the number of lines in every **.html** or **.htm** file in or below current directory and redirect the results to separate files:

```bash
find . -type f \( -name "*.html" -o -name "*.htm" \) -print0 | xargs -0 -I{} sh -c 'wc -l "$1" > "$1.output.txt"' -- {}
```

Print the number of lines in every **.html** or **.htm** file in or below current directory and redirect the results to separate files; process up to 4 files in parallel:

```bash
find . -type f \( -name "*.html" -o -name "*.htm" \) -print0 | xargs -n1 -P4 -0 -I{} sh -c 'wc -l "$1" > "$1.output.txt"' -- {}
```

## Using grep

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

## Using awk

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

Split a mulit-FASTA file into separate files, one per sequence named according to the sequence title. In this example the sequences are written to a directory called **output\_directory**:

```bash
outputdir=output_directory/
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

## Using sed

Print a specific line of a file. In this example line **26404**:

```bash
sed -n "26404p" input.txt
```

Change filenames using a regular expression. In this example **chr30** is replaced with **chrX**:

```bash
for f in *.fasta; do new=`echo $f | sed 's/chr30/chrX/'`; mv $f $new; done
```

## Using perl

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

Remove commas located within quoted fields in a CSV file:

```bash
perl -nle  'my @new  = (); push( @new, $+ ) while $_ =~ m{"([^\"\\]*(?:\\.[^\"\\]*)*)",? | ([^,]+),? | ,}gx; push( @new, undef ) if substr( $text, -1, 1 ) eq '\'','\'';
 for(@new){s/,/ /g} print join "\t", @new' input.csv > output.csv
```

Replace tabs with commas and remove quotes in a CSV file:

```bash
perl -p -e 's/\t/,/g;' -e 's/"//g' input.csv > output.csv
```

## Other

Combine the columns in two tab-delimited files.

```bash
paste -d"\t" input1.tab input2.tab > output.tab
```

Add a header to all files with a certain extension, getting the header from another file. In this example the header is added to **.tab** files and comes from a file called **header.txt**. The files with the header added are saved with a **.new** extension added:

```bash
for f in *.tab; do new=`echo $f | sed 's/\(.*\)\.tab/\1.tab.new/'`; paste -sd'\n' \header.txt "$f" > "$new"; done
```