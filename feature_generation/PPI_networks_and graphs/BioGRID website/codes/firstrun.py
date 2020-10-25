filename = "BIOGRID-ALL-3.5.184.tab3 - Copy.txt"
import csv
f = open('firstrun.csv', 'w')
writer = csv.writer(f)

a = 0;
with open(filename, newline = '') as games:                                                                                          
	game_reader = csv.reader(games, delimiter='\t')
	for game in game_reader:	
		if game[15] == "559292" or game[16] == "559292":
			print('writing line', a)
			a+=1
			writer.writerows([game])