This code uses disorder prediction results for all yeast's genes and creates relevant features from this data. For further reading about the data visit our wiki – model page.

**Input:**

•	Q0045_COX1_SGDID_S000007260.merged. This is one example file out of 5500 files used. We uploaded one file for example due to copyrights.
Each file is prediction result of one gene. In each file there is per-residue prediction results from disorder prediction algorithms and ANCHOR disordered binding sites

**Code:** The data has been proceeded with MATLAB 'disorder.m' script.

**Output: 'disorder_output.xlsx' (columns=features, rows=genes).**

The file contains the following features:

1.	Gene_name- the gene's name taken from the file's name.

2.	Mean_score- mean on all the numerical columns in the original file on all the positions.

3.	 Disorder_percentage1- for every gene this score is calculated using all the gene's positions by consensus approach #1 (Explained on our wiki). We counted the number of disorder sites in the gene (3/4 in column 7 in the original file) and divided by the number of positions.

4.	Disorder_percentage2- for every gene this score is calculated using all the gene's positions by consensus approach #2 (Explained on our wiki). We counted the number of disorder sites in the gene ("DISORDERED" in column 8 in the original file) and divided by the number of positions.

5.	Disorder_percentage_both_algorithms- for every gene this score is calculated using all the gene's positions. We counted the number of disorder sites in the gene when both 
consensus approaches determined the site as disordered ("DISORDERED" in column 8 and 3/4 value in column 7 in the original file) and divided by the number of positions.

6.	Anchor_disorder_percentage- for every gene this score is calculated using all the gene's positions. We counted the number of disorder sites in the gene (1 in column 13 in the original file) and divided by the number of positions.

7.	Anchor_percentage-. We used mean function on all the gene's positions in column 12 in the original file (Anchor score, explained on our wiki).

8.	IUPRED_percentage- We used mean function on all the gene's positions in column 3 in the original file (IUPRED score, explained on our wiki).

9.	VSL2B_percentage- We used mean function on all the gene's positions in column 4 in the original file (VSL2B score, explained on our wiki).

10.	DisEMBL_percentage- for every gene this score is calculated using all the gene's positions. We counted the number of disorder sites in the gene ("DISORDERED" in column 5 in the original file,DisEMBL, , explained on our wiki) and divided by the number of positions.

11.	MoreRONN_percentage- We used mean function on all the gene's positions in column 6 in the original file (MoreRONN score, explained on our wiki).

12.	Disorder30_percentage-- for every gene this score is calculated using all the gene's positions. We counted the number of disorder30 sites in the gene ("DIS_30" in column 9 in the original file, explained on our wiki) and divided by the number of positions.

13.	Disorder_percentage1_window1-Disorder_percentage1_window30- each window contains 50 gene's positions and the starting point for each window is the number of the window. For example, window 17 is from position 17 until position 66. In each window we counted the number for disorders according to consensus #1 and divided by 50. If the gene's length is smaller than 79 than some windows aren't full and therefore these window's value will be NaN.

14.	Disorder_std1_window1-Disorder_std1_window30-in each window detailed in 13, we calculated the standard deviation (If not possible due to incomplete window the value assigned is NaN).

15.	Disorder_percentage2_window1-Disorder_percentage2_window30- same as described in 13 using consensus #2 (If not possible due to incomplete window the value assigned is NaN).

16.	Disorder_std2_window1-Disorder_std2_window30- same as described in 14 using consensus #2 (If not possible due to incomplete window the value assigned is NaN).

17.	Disorder30_percentage_window1- Disorder30_percentage_window30- the windows are same as described before. In each window we counted the number of "DIS_30" values and divided by 50 (If not possible due to incomplete window the value assigned is NaN).

18.	Disorder30_std_window1- Disorder30_std_window30- we created a vector of 0 and 1 when 1="DIS_30" and 0="ORD_30" and for each window and calculated the vector's std (If not possible due to incomplete window the value assigned is NaN).

19.	Disorder_anchor_percentage_window1- Disorder_anchor_percentage_window30- the windows are same as describe before. The calculation is same as described 6 but for every window (If not possible due to incomplete window the value assigned is NaN).

20.	Disorder_anchor_std_window1- Disorder_anchor_std_window30- the anchor's values are 0/1, therefore std was easily computed for every window (if not possible due to incomplete window the value assigned is NaN).

21.	Max_disorder_percentage_window1- Max_disorder_percentage_window30- for every window the value is the maximum score between all the positions and all the score columns in the original file (3,4,6,12) (if not possible due to incomplete window the value assigned is NaN).

22.	Disorder_percentage1_lastwindow, Disorder_std1_lastwindow, Disorder_percentage2_lastwindow, Disorder_std2_lastwindow,  Disorder30_percentage_lastwindow,Disorder30_std_lastwindow, Disorder_anchor_percentage_lastwindow, Disorder_anchor_std_lastwindow, Max_disorder_percentage_lastwindow-  the following were calculate the same way described in 13-21 but the window was different. Here the features were calculated only to the last window- 50 positions from the end and backward and therefore there aren't NaN values.

**NOTE:If the sequence is smaller than 129 the last window share's some of its positions with the first 30 windows.**

