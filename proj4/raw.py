"""
this was code i initially wrote for the textbook leetcode question
it's horribly unoptimized, but at least you can tell to a degree what it's doing
the numpy version in the notebook is genuinely awful
"""

import numpy as np


class LP:  # short for log prob
    def __init__(self, val: float):
        self.val = np.log(max(val, 0))

    @classmethod
    def from_raw(cls, val: float):
        ret = cls(1)
        ret.val = val  # can't think of a better way HAHA
        return ret

    # these operations are the only ones i should have to support :clueless:
    def __add__(self, other: LP):
        return self.from_raw(np.logaddexp(self.val, other.val))

    def __mul__(self, other: LP):
        return self.from_raw(self.val + other.val)

    def __truediv__(self, other: LP):
        return self.from_raw(self.val - other.val)

    def __repr__(self):
        return str(self.val)


STATES = "ABCD"
EMIT_CHARS = "xyz"

marks = "zxxzyyyzxxxxxxzyzxzxzyzzzxzxxxyxyyyzxyzyxzxxyyxxxzyxxyzxzzzxzyxyzzyxzxyyxzxyyzzyyxzyxyyxxzyzzxzzzzyx"

TRANSITION = {
    "A": {"A": 0.077, "B": 0.351, "C": 0.286, "D": 0.287},
    "B": {"A": 0.337, "B": 0.297, "C": 0.237, "D": 0.129},
    "C": {"A": 0.344, "B": 0.409, "C": 0.068, "D": 0.179},
    "D": {"A": 0.214, "B": 0.44, "C": 0.203, "D": 0.143},
}
E_PROB = {
    "A": {"x": 0.601, "y": 0.264, "z": 0.134},
    "B": {"x": 0.209, "y": 0.665, "z": 0.126},
    "C": {"x": 0.361, "y": 0.547, "z": 0.092},
    "D": {"x": 0.343, "y": 0.159, "z": 0.498},
}

ep = E_PROB
prob = TRANSITION
for s in STATES:
    for ns in STATES:
        prob[s][ns] = LP(prob[s][ns])
    for e in EMIT_CHARS:
        ep[s][e] = LP(ep[s][e])

for _ in range(100):
    fwd = [{c: LP(1 / len(STATES)) * ep[c][marks[0]] for c in STATES}]
    for i in marks[1:]:
        new_c = {c: LP(0) for c in STATES}
        for c in STATES:
            for prev_c in STATES:
                new_c[c] += fwd[-1][prev_c] * prob[prev_c][c]
            new_c[c] *= ep[c][i]
        fwd.append(new_c)

    sink = sum(fwd[-1].values(), LP(0))

    bkwd = [{c: LP(1) for c in STATES}]
    for i in reversed(marks[1:]):
        new_c = {c: LP(0) for c in STATES}
        for c in STATES:
            for next_c in STATES:
                new_c[c] += bkwd[-1][next_c] * prob[c][next_c] * ep[next_c][i]
        bkwd.append(new_c)
    bkwd.reverse()

    star = []
    for i in range(len(fwd)):
        ans = {}
        for c in STATES:
            ans[c] = fwd[i][c] * bkwd[i][c] / sink
        star.append(ans)

    sstar = []
    for v, i in enumerate(star[:-1]):
        ans = {}
        for c1 in STATES:
            for c2 in STATES:
                weight = prob[c1][c2] * ep[c2][marks[v + 1]]
                ans[c1 + c2] = fwd[v][c1] * weight * bkwd[v + 1][c2] / sink
        sstar.append(ans)

    new_prob = {c: {} for c in STATES}
    new_ep = {c: {} for c in STATES}
    for c1 in STATES:
        for c2 in STATES:
            new_prob[c1][c2] = sum((s[c1 + c2] for s in sstar), LP(0))

        for e in EMIT_CHARS:
            tot = LP(0)
            for i in range(len(marks)):
                if marks[i] == e:
                    tot += star[i][c1]
            new_ep[c1][e] = tot

    for c in STATES:
        tot = sum(new_prob[c].values(), LP(0))
        new_prob[c] = {k: v / tot for k, v in new_prob[c].items()}
        tot = sum(new_ep[c].values(), LP(0))
        new_ep[c] = {k: v / tot for k, v in new_ep[c].items()}

    prob = new_prob
    ep = new_ep

print("\t" + "\t".join(STATES))
for a, b in prob.items():
    print(a + "\t" + "\t".join(str(round(np.exp(b[s].val), 3)) for s in STATES))

print("\t" + "\t".join(EMIT_CHARS))
for a, b in ep.items():
    print(a + "\t" + "\t".join(str(round(np.exp(b[s].val), 3)) for s in EMIT_CHARS))
