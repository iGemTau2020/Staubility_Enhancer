The code optimzes a given sequence, composed of target gene, linker and conjugates (essential) gene. 

**optimization inputs:** 
1. coding sequence (target + linker + essential), nucleotides alphabet. make sure that the target and essential are divisible by three. 
2. optimization parameter: tAI, NTE, TDR or codon usage fraction. 
3. optimization method: `use_best_codon`, `match_codon_usage`, or `harmonize_rca`
4. GC content - give a range (minimum and maximum). window size is default.  
5. restriction enzyme pattern to avoid (from the list of `rest_dict.keys()` from rest_dict in bioPython)
6. linker indices - for the extraction of target gene, essential gene and linker.  
7. output path - to save the report. 
8. number of optimized EFM sites

**NOTE:** another input for the funtion is the linker indices as a tuple. It in not provided by the user but by us. 

**output: optimized sequence based on the following steps:** 
1. optimization of initiation of translation - reduce the mRNA folding energy in the start of the sequence (about 15 codons). 
2. optimization of the translation efficiency - combine low GC content with fast codons - report included (**genbank and pdf formats)**.  
3. Unsupervised optimization of stability - EFM optimizer 

**Required packages:**  
* Vienna-RNA (https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.14.tar.gz) 
* forgi (pip install)
* dnachisel (pip install)
* biopython (pip install)
* sequenticon (pip install)
