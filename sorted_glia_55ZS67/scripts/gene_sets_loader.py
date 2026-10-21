import os
def load_sets():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out = {}
    for line in open(os.path.join(root, "data", "modules.gmt")):
        p = line.rstrip("\n").split("\t"); out[p[0]] = p[2:]
    return out
