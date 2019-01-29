with open("nbt3536-S3.tsv", "r") as inp, open("nbt3536.fa", "w") as oup:
    for line in inp:
        line = line.strip().split()
        oup.write(">{}\n{}\n".format(line[0], line[1]))
