import glob
import gzip
import pickle

gmts = glob.glob("gmts/*.gmt")

gmt_lists: list[list[set[str]]] = []
for gmt in gmts:
    gmt_list = []
    with open(gmt) as r:
        data = r.read().split("\n")
    for line in data:
        if "\t" in line:
            vals = line.split("\t")
            gmt_list.append(set(vals[2:]))
    gmt_lists.append(gmt_list)

with gzip.open("gmt.pkl.gz", "wb") as w:
    pickle.dump(gmt_lists, w)
