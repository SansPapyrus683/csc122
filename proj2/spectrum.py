"""major issue is that eulerian paths aren't unique. no real way to address this..."""
import networkx as nx
import re
from collections import defaultdict

header_fmt = re.compile(r">read_(\d+)")
spectrum = "data/project2a_spectrum.fasta"

read_num = []
kmers = []
for v, line in enumerate(open(spectrum)):
    if v % 2 == 1:
        kmers.append(line.strip())
    else:
        match = header_fmt.match(line).group(1)
        read_num.append(int(match))

assert len({len(k) for k in kmers}) == 1

g = nx.MultiDiGraph()
for k in kmers:
    g.add_edge(k[:-1], k[1:])

kmer_to_reads = defaultdict(list)
for r, k in zip(read_num, kmers):
    kmer_to_reads[k].append(r)

# 2a guarantees no errors in the reads
for v, (a, b) in enumerate(nx.eulerian_path(g)):
    k = a[0] + b
    r = kmer_to_reads[k].pop()
    print(f">read_{r}\t{v}")
