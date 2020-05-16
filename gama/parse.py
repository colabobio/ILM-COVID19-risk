import csv

in_counts_file = "infected_number.txt"
out_counts_file = "case_counts.csv"

debug_contacts = True
in_contacts_file = "contact_data_compact.txt"

data = []
t = t0 = 0
r = r0 = lastr = 0
i = i0 = lasti = 0
with open(in_counts_file, 'r') as f:
    reader = csv.reader(f, dialect="excel")
    next(reader)
    for row in reader:
        t = int(row[0])
        i = int(row[4])
        r = int(row[5]) + int(row[6])
        if t0 < t:
            dr = max(0, r0 - lastr)
            di = max(0, i0 - lasti + dr)
            data.append([t0, di])
            lastr = r0
            lasti = i0
        t0 = t
        i0 = i
        r0 = r
        if 120 < t:
            break

dr = max(0, r0 - lastr)
di = max(0, i0 - lasti + dr)
data.append([t0, di])

with open(out_counts_file, "w") as f:
    writer = csv.writer(f, dialect="excel")    
    writer.writerow(["Time", "Count"])
    for row in data:
        writer.writerow(row)

# ------------------------------------------

history = []
infected = {}
t0 = 0
with open(in_contacts_file, 'r') as f:
    reader = csv.reader(f, dialect="excel")
    next(reader)
    for row in reader:
        t = int(row[0])
        sind = row[2]
        scov = float(row[3])
        iind = row[4]
        icov = float(row[5])

        if iind in infected:
            cdat = infected[iind]
            clst = cdat["contacts"]
        else:
            clst = []
            cdat = {"infectivity":icov, "contacts":clst}
            infected[iind] = cdat

        clst += [scov]
         
        if t0 != t:
            history += [{"time":t, "infected":infected}]
            infected = {}
        t0 = t

values = []
indices =[-1] * (t + 1)
nval = 0
for cdat in history:
   t = cdat["time"]
   if debug_contacts: print("=======> time")
   infected = cdat["infected"]
   indices[t] = nval
   for inf in infected:
       infectivity = infected[inf]["infectivity"]
       contacts = infected[inf]["contacts"]
       ilen = len(contacts)
       nval += 2 + ilen
       values += [ilen, infectivity]
       for v in contacts:
           values += [v]
       if debug_contacts: print(inf, infectivity, contacts)
   values += [-1]
   nval += 1
   # print(len(values), nval)

if debug_contacts:
    for t in range(0, len(indices)):
       idx = indices[t]
       if idx < 0: continue
       print("=======> time", t)
       while -1 < values[idx]:
           ncont = int(values[idx])
           idx += 1
           inf = values[idx]
           idx += 1       
           print(inf, end =" ")
           susc = values[idx:idx+ncont]
           idx += ncont
           print(susc)

# Adding missing indices
for i in range(len(indices), 0, -1):
    if indices[i-1] < 0: indices[i-1] = indices[i]
print(indices)

with open('indices', 'w') as f:
    for idx in indices:
        f.write("%i\n" % idx)

with open('contacts', 'w') as f:
    for v in values:
        f.write("%.3f\n" % v)

