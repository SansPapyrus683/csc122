"""this relies on a ton of assumptions, idk... it gets a high enough score idc."""
import networkx as nx
import re
from collections import defaultdict
from statistics import mode

K = 20
THRESH = 12


def common_amt(s1: str, s2: str) -> int:
    """how many characters in s2's pref are in s1's suff?"""
    max_ = min(len(s1), len(s2))
    for i in range(max_, 0, -1):
        if s1[-i:] == s2[:i]:
            return i
    return 0


header_fmt = re.compile(r">read_(\d+)")
reads = "data/project2_sample2_reads.fasta"
reads = "data/project2b_reads.fasta"

read_num = []
raw_reads = []
for v, line in enumerate(open(reads)):
    if v % 2 == 1:
        raw_reads.append(line.strip())
    else:
        match_amt = header_fmt.match(line).group(1)
        read_num.append(int(match_amt))

kmer_freq = defaultdict(int)
for r in raw_reads:
    for i in range(0, len(r) - K + 1):
        kmer_freq[r[i : i + K]] += 1

g = nx.MultiDiGraph()
please = 0
for k, amt in kmer_freq.items():
    true_amt = round(amt / THRESH)
    please += true_amt
    for _ in range(true_amt):
        g.add_edge(k[:-1], k[1:])

deficit = []  # too many outgoing
surplus = []  # too many incoming
for k in g.nodes:
    bad = g.in_degree(k) - g.out_degree(k)
    if bad < 0:
        deficit.extend([k] * abs(bad))
    elif bad > 0:
        surplus.extend([k] * bad)
assert len(deficit) == len(surplus)

# add all the edges except one
for start in surplus[:-1]:
    best = 0, ""
    for end in deficit:
        best = max(best, (common_amt(start, end), end))

    g.add_edge(start, best[1])
    deficit.remove(best[1])

genome = []
e_path = list(nx.eulerian_path(g))
genome.append(e_path[0][0])
for a, b in e_path[1:]:
    amt = common_amt(a, b)
    genome.append(b[amt:])
genome = "".join(genome)

kmers = defaultdict(list)
for i in range(0, len(genome) - K + 1):
    kmers[genome[i:i + K]].append(i)

out = []
for r, id_ in zip(raw_reads, read_num):
    votes = []
    for ind in range(0, len(r) - K + 1):
        for pos in kmers[r[ind:ind + K]]:
            votes.append(pos - ind)
    if not votes:
        votes = [0]  # idk bro
    out.append((id_, mode(votes)))
out.sort(key=lambda i: i[1])

for id_, pos in out:
    print(f">read_{id_} {pos}")
