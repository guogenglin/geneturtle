# geneturtle
A python script to draw gene structure map using turtle

# Introduction
This is a tool for draw gene structures (for prokaryote, no module for draw intron). Every gene will be marked by an arrow, genes belongs to same group (assign by users) will have same color. 

# Required Packages
turtle

# Input file
A five columns table should be prepared as input before start, consist by 'gene name', 'start', 'end', 'direction', 'group', for example:

'''
dnaA 1 1375 + DNA
dnaN 1529 2665 + DNA
RS00015 2757 3638 + NA
RS00020 3648 3851 - NA
RS00025 3955 4314 + NA
ychF 4400 5515 + enzyme
cps2A 5600 6700 - CPS
cps2B 7000 9000 - CPS
'''

Notice: Please ensure that your sequence not longer than 40000bp, and the group of genes no more than 10.

# Quick Start   
``` Python
python geneturtle.py inputfile.txt
```

# Example output
10000bp for one line, if the sequence length is more than 10000bp, geneturtle will auto change to a new line.

One line:
![gene_map_1](https://user-images.githubusercontent.com/108860907/228293664-1a3f8951-15fb-45e0-8c5d-f093a7c24c6a.jpg)

Multiple lines:
![gene_map_2](https://user-images.githubusercontent.com/108860907/228293743-ad35e723-41d7-4671-80f3-9c7fb43f06c2.jpg)
