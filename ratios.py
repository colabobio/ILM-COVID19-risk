import csv, os, argparse

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--main_dir', nargs=1, default=['./gama'],
                    help='Directory containing GAMA data')
parser.add_argument('-o', '--output_dir', nargs=1, default=['./output'],
                    help='Sub-directory containing output')
parser.add_argument("-n", "--num_values", type=int, nargs=1, default=[10],
                    help="Numbr of likelihood values to use in mean calculation")
args = parser.parse_args()

main_dir = args.main_dir[0]
output_dir = args.output_dir[0]
num_values = args.num_values[0]

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

tr_a0a1 = true_a0/true_a1
tr_b0b1 = true_b0/true_b1

p1_a0a1 = point_a01/point_a11
p2_a0a1 = point_a00/point_a10

m1_a0a1 = mean_a01/mean_a11
m2_a0a1 = mean_a00/mean_a10

p1_b0b1 = point_a00/point_a01
p2_b0b1 = point_a10/point_a11

m1_b0b1 = mean_a00/mean_a01
m2_b0b1 = mean_a10/mean_a11

print("R1) a0/a1=%.3f point=%.3f mean=%.3f" % (tr_a0a1, p1_a0a1, m1_a0a1))
print("R2) a0/a1=%.3f point=%.3f mean=%.3f" % (tr_a0a1, p2_a0a1, m2_a0a1))
print("R1) b0/b1=%.3f point=%.3f mean=%.3f" % (tr_b0b1, p1_b0b1, m1_b0b1))
print("R2) b0/b1=%.3f point=%.3f mean=%.3f" % (tr_b0b1, p2_b0b1, m2_b0b1))

pm_a0a1 = 0.5 * (p1_a0a1 + p2_a0a1)
mm_a0a1 = 0.5 * (m1_a0a1 + m2_a0a1)

pm_b0b1 = 0.5 * (p1_b0b1 + p2_b0b1)
mm_b0b1 = 0.5 * (m1_b0b1 + m2_b0b1)

print("R1+R2) a0/a1=%.3f point=%.3f mean=%.3f" % (tr_a0a1, pm_a0a1, mm_a0a1))
print("R1+R2) b0/b1=%.3f point=%.3f mean=%.3f" % (tr_b0b1, pm_b0b1, mm_b0b1))

print("============================")

import pandas as pd
import numpy as np

search_file = os.path.join(main_dir, output_dir, "param_mle_global_search.csv")
df = pd.read_csv(search_file)
df = df.sort_values(['loglik'], ascending=[False]).head(num_values)
df['a0a1_1'] = df['a01'] / df['a11']
df['a0a1_2'] = df['a00'] / df['a10']
df['b0b1_1'] = df['a00'] / df['a01']
df['b0b1_2'] = df['a10'] / df['a11']

raw_a0a1 = np.concatenate((np.array(df['a0a1_1']), np.array(df['a0a1_2'])))
raw_b0b1 = np.concatenate((np.array(df['b0b1_1']), np.array(df['b0b1_2'])))

print("a0/a1=%.3f - MLE mean=%.3f std=%.3f" % (tr_a0a1, np.mean(raw_a0a1), np.std(raw_a0a1)))
print("b0/b1=%.3f - MLE mean=%.3f std=%.3f" % (tr_b0b1, np.mean(raw_b0b1), np.std(raw_b0b1)))

