This code enables the user to choose his preferred linker type- 2A/fusion and it presents provideprovides himthemhim  with the best linkers in of this type.

**Function inputs:**

•	Coding sequence (target + conjugate), amino acids alphabet.

•	Type of linker – 2A/fusion linker.

•	Ouput path – the path that the user would like to save the disorder profile plots in.

•	Linkers.csv (downloaded from https://www.ibi.vu.nl/programs/linkerdbwww/)

**Output: 10 options for best fusion linkers / 4 options for best 2A linkers (with descending order) based on the following principles:
2A linkers**

2A linkers:

•	There are four common 2A sequences.

•	The conventional way to rank the 2A peptides is from the most efficient P2A, followed by T2A, E2A and F2A.


Fusion linkers:

•	Essential gene's score- a Euclidian distance between the essential gene disorder profile before and after fusion.

•	Target gene's score- a Euclidian distance between the target gene disorder profile before and after fusion.

•	Final score = 1/2 * Essential gene's score + 1/2 * Target gene's score

•	Chosen linker= linkers [argmin(final score)]


**Required packages:**

•	IUPPED2A (https://iupred2a.elte.hu/download_new)
