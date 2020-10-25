filename = "BIOGRID-ALL-3.5.184.tab3 - Copy.txt"
import csv
f = open('firstrun.csv', 'w')
writer = csv.writer(f)

a = 0;
with open(filename, newline = '') as files:                                                                                          
	file_reader = csv.reader(files, delimiter='\t')
	for file in file_reader:	
		if file[15] == "559292" or file[16] == "559292":
			print('writing line', a)
			a+=1
			writer.writerows([game])
