#1/usr/bin/python3

import subprocess
import csv
from collections import defaultdict

# Extracts columns of data from an excel worksheet.

# Jon Parsons
# 8-5-19

################################################################################

# Get workbook directory

print("Input workbook name with relative directory.\n")
book = input("Book name: ")

# Decide which data you want
row_start = int(input("Input the row the data begins on:\n"))
col1 = int(input("Input first Column number\n"))
col2 = int(input("Input second Column number\n"))

# Match indices
row_start -= 1
col1 -= 1
col2 -= 1

# Read in data
reader = csv.reader(open(book), delimiter=',')
data = [var for var in [row for row in reader]]

# output
file = open("data_out.dat", 'w')

counter = 0

for row in data:
    if counter >= row_start:
        file.write("{} {}\n".format(data[counter][col1], data[counter][col2]))

    counter += 1


file.close()
