"""
This script adds the openmc XML input files into one file for version control. 

Type this into command line: 
python xml_to_input.py --o (name of output file)

"""

import argparse

parser = argparse.ArgumentParser(description='xml parser')
parser.add_argument("--o")

args = parser.parse_args()
out = args.o

with open("settings.xml") as f1:
    lines1 = f1.readlines()
with open("materials.xml") as f2: 
    lines2 = f2.readlines()
with open("geometry.xml") as f3: 
    lines3 = f3.readlines()
    with open(out, "w") as fout:
        fout.writelines(lines1)
        fout.write('\n')
        fout.writelines(lines2)
        fout.write('\n')
        fout.writelines(lines3)