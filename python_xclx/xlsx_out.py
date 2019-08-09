#1/usr/bin/python3

import subprocess
import csv
from collections import defaultdict

# Extracts columns of data from an excel worksheet.

# Jon Parsons
# 8-5-19

################################################################################

# Convert xlsx to readable without external library
def xlsx(fname,sheet):
    import zipfile
    from xml.etree.ElementTree import iterparse
    z = zipfile.ZipFile(fname)
    sheet_name = 'xl/worksheets/sheet'+sheet+'.xml'
    strings = [el.text for e, el in iterparse(z.open('xl/sharedStrings.xml')) if el.tag.endswith('}t')]
    rows = []
    row = {}
    value = ''
    for e, el in iterparse(z.open(sheet_name)):
        if el.tag.endswith('}v'):
            value = el.text
        if el.tag.endswith('}c'):
            if el.attrib.get('t') == 's':
                value = strings[int(value)]
            letter = el.attrib['r']
            while letter[-1].isdigit():
                letter = letter[:-1]
            row[letter] = value
            value = ''
        if el.tag.endswith('}row'):
            rows.append(row)
            row = {}
    return rows


# Get workbook directory

print("Input workbook name with relative directory.\n")
book = input("Book name: ")
page = input("\nName of sheet:")

# Decide which data you want
row_start = int(input("Input the row the data begins on:\n"))
col1 = int(input("Input first Column number\n"))
col2 = int(input("Input second Column number\n"))

# Match indices
row_start -= 1
col1 -= 1
col2 -= 1

data = xlsx(book,page)

# output
file = open("data_out_xls.dat", 'w')

counter = 0

for row in data:
    if counter >= row_start:
        file.write("{} {}\n".format(data[counter][col1], data[counter][col2]))

    counter += 1


file.close()
