
import csv
a=0;
with open('firstrun.csv', 'r', newline ='') as f,open('finaldata.csv', 'w',newline = '') as f_out:
	reader = csv.reader(f)
	writer = csv.writer(f_out)
	for row in reader:
		a+=1
		if row != []:
			print('im in line', a)
			writer.writerow(row)