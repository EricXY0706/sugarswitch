protein = "NB1"
pos = "60"

lines = open(f"/sdata2/WORK/Xuyi/{protein}/evc/_ECs.txt").readlines()
for line in lines:
    line = line.strip("\n").split(" ")
    if line[0] == pos or line[2] == pos:
        if abs(float(line[-1])) > 0.5:
            print(line)