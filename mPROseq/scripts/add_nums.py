import sys

with open(sys.argv[1]) as f:
    c = 1
    for inline in f:
        ln = inline.strip()
        print(str(c)+'\t'+ln)
        c += 1
