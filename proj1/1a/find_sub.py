from collections import Counter, defaultdict
from statistics import mode
import pickle

chars = "ATCG"
char_to_ind = {c: i for i, c in enumerate(chars)}

with open("project1a_reference_genome.fasta") as gf:
    gf.readline()  # first line is irrelevant
    genome = "".join(l.strip() for l in gf.readlines())

queries = []
with open("project1a_with_error_paired_reads.fasta") as rf:
    for v, i in enumerate(rf):
        if v % 2 == 1:
            queries.append(i.strip())

candidates = defaultdict(list)
for v, q in enumerate(queries):
    for i in range(len(genome) - len(q) + 1):
        mismatches = [j for j, c in enumerate(q) if c != genome[i + j]]

        if len(mismatches) <= 3:
            for j in mismatches:
                candidates[i + j].append(q[j])

for pos, goto in candidates.items():
    if len(goto) >= 10:
        print(f">S{pos} {genome[pos]} {mode(goto)}")
