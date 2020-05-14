import csv

in_counts_file = "infected_number.txt"
out_counts_file = "case_counts.csv"

in_contacts_file = "contact_data_compact.txt"

data = []
t = t0 = 0
c = c0 = 0
with open(in_file, 'r') as f:
    reader = csv.reader(f, dialect="excel")
    next(reader)
    for row in reader:
    	t = int(row[0])
    	c = int(row[4])
    	if t0 < t:
    		data.append([t0, c0])
    	t0 = t
    	c0 = c
    	if 120 < t:
    		break

data.append([t0, c0])
	

with open(out_counts_file, "w") as f:
    writer = csv.writer(f, dialect="excel")    
    writer.writerow(["Time", "Count"])
    for row in data:
        writer.writerow(row)
