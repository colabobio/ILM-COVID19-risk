import csv, os, argparse

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--main_dir', nargs=1, default=['./gama'],
                    help='Directory containing GAMA data')
parser.add_argument('-o', '--output_dir', nargs=1, default=['./output'],
                    help='Sub-directory containing output')
args = parser.parse_args()

main_dir = args.main_dir[0]
output_dir = args.output_dir[0]

param_file = os.path.join(main_dir, "parameter_values.txt")

point_mle_file = os.path.join(main_dir, output_dir, "param_point_estimates.csv")
mean_mle_file = os.path.join(main_dir, output_dir, "param_mean_sdev.csv")

true_a0 = 0
true_a1 = 0
true_b0 = 0
true_b1 = 0
with open(param_file, 'r') as f:
    reader = csv.reader(f, dialect="excel")
    next(reader)
    for row in reader:
        if row[0] == 'A0': true_a0 = float(row[1])
        if row[0] == 'A1': true_a1 = float(row[1])
        if row[0] == 'B0': true_b0 = float(row[1])
        if row[0] == 'B1': true_b1 = float(row[1])

#print(true_a0, true_a1, true_b0, true_b1)

point_a00 = 0
point_a01 = 0
point_a10 = 0
point_a11 = 0
with open(point_mle_file, 'r') as f:
    reader = csv.reader(f, dialect="excel")
    next(reader)
    for row in reader:
        if row[0] == 'a00': point_a00 = float(row[1])
        if row[0] == 'a01': point_a01 = float(row[1])
        if row[0] == 'a10': point_a10 = float(row[1])
        if row[0] == 'a11': point_a11 = float(row[1])

#print(point_a00, point_a01, point_a10, point_a11)

mean_a00 = 0
mean_a01 = 0
mean_a10 = 0
mean_a11 = 0
with open(mean_mle_file, 'r') as f:
    reader = csv.reader(f, dialect="excel")
    next(reader)
    for row in reader:
        if row[0] == 'a00': mean_a00 = float(row[1])
        if row[0] == 'a01': mean_a01 = float(row[1])
        if row[0] == 'a10': mean_a10 = float(row[1])
        if row[0] == 'a11': mean_a11 = float(row[1])

#print(mean_a00, mean_a01, mean_a10, mean_a11)

print("T=%.3f P=%.3f M=%.3f" % (true_b0/true_b1, point_a00/point_a01, mean_a00/mean_a01))
print("T=%.3f P=%.3f M=%.3f" % (true_a0/true_a1, point_a01/point_a11, mean_a01/mean_a11))