# geneturtle
A python script to draw gene structure map using turtle

# Introduction
This is a tool for draw gene structures (for prokaryote, no module for draw intron). Every gene will be marked by an arrow, genes belongs to same group (assign by users) will have same color. 

# Required Packages
turtle, math, biopython

# Input file
Now geneturtle could process two type of files, GenBank format, and the five columns table could both be recognized.
A five columns table should be prepared as input before start, consist by 'gene name', 'start', 'end', 'direction', 'group', for example:

```
dnaA 1 1375 + DNA
dnaN 1529 2665 + DNA
RS00015 2757 3638 + NA
RS00020 3648 3851 - NA
RS00025 3955 4314 + NA
ychF 4400 5515 + enzyme
cps2A 5600 6700 - CPS
cps2B 7000 9000 - CPS
```

Notice: Long sequence could be processed, but I strongly recommend the sequence not longer than 40, 000bp, if the sequence is too long, there may have some unpredictable error.

# Usage
Put the python script and database folder into the folder contains your sequence file

```
geneturtle [-i] [-t] [-o] [--color_mode] [--force ] [-v]
Input and Output:
  -i, --input             The file include the detail to draw the figure
  -t, --type              The file type of your input file, "gbk" and "list" are available
  -o, --output            Output file
Parameters:
  -c, --color_mode        We have three mode for color assignment for every gene group, you can choose "random", "preset", or upload the color file set by yourself. [default: preset]
  -f, --force             In case the sequence is too long and you still want to run it.
  -v, --version           Show version number and exit
```


# Quick Start   
``` Python
python geneturtle.py -i inputfile.txt -t list
```


# Colormode
Now we provide three colormode, "preset", "random", and you can also upload a file with group and it corresponding color. "preset" mode will provide 10 prepared color, so if you have more than 10 groups, then "preset" mode cannot be used; "random" mode means the color will be generated randomly, so there are no limited for group number. However, I recommend that do not make too many groups, that will make the figure very hard to see.
The third mode allows the users set the color by themselves, provide an extra table in this format:

```
Transcriptional regulator 221,160,221
Transposase 79,79,79
```

It should be noticed that if you are using GenBank format data and you want to set the color yourself, you should edit the gene product in GenBank file, or the concise_gbk_feature funtion in geneturtle.


# Example output
10,000bp for one line, if the sequence length is more than 10,000bp, geneturtle will auto change to a new line. (if you added command "-f" or "-force", one line will be more than 10,000bp)

One line:
![gene_map_2](https://user-images.githubusercontent.com/108860907/228293743-ad35e723-41d7-4671-80f3-9c7fb43f06c2.jpg)

Multiple lines:
![gene_map_1](https://user-images.githubusercontent.com/108860907/228293664-1a3f8951-15fb-45e0-8c5d-f093a7c24c6a.jpg)
