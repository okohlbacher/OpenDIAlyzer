#!/usr/bin/env python3
"""S1 gate: does the fragment-sharing graph over non-specific peptides percolate?

Family-level FDR only works if peptides cluster into many small families. If the
fragment-sharing graph forms one giant component, "family-level FDR" is
proteome-level FDR and the whole idea collapses.

Link criterion: two peptides share >= k b-ions iff they share their first k
residues; they share >= k y-ions iff they share their last k residues. So
"shares >= k fragments" reduces exactly to a common prefix or suffix of length
k -- no mass arithmetic needed.

That reduction also makes the graph cheap to build. Every peptide joins exactly
one prefix-group and one suffix-group, so the peptide graph's components are the
components of the bipartite prefix/suffix group graph. Union-find over groups is
O(n) instead of O(n^2) over peptide pairs.

Sequence-based, so it ignores mass degeneracy (I/L, and different residue
combinations of equal mass). That makes it a LOWER bound on connectivity: the
real graph is at least this connected, never less.
"""
import sys
from collections import defaultdict


def read_fasta(path):
    seqs, cur = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


class UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        p = self.p
        if x not in p:
            p[x] = x
            return x
        root = x
        while p[root] != root:
            root = p[root]
        while p[x] != root:            # path compression
            p[x], x = root, p[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def percolate(seqs, k, lo=8, hi=12):
    """Return (n_peptides, n_components, largest_component_fraction)."""
    uf = UF()
    weight = defaultdict(int)          # key-pair -> peptide count
    n_pep = 0
    for s in seqs:
        n = len(s)
        for i in range(n - lo + 1):
            for L in range(lo, hi + 1):
                if i + L > n:
                    break
                pkey = ("p", s[i:i + k])
                skey = ("s", s[i + L - k:i + L])
                uf.union(pkey, skey)
                weight[pkey] += 1
                n_pep += 1
    if not n_pep:
        return 0, 0, 0.0
    comp = defaultdict(int)
    for key, w in weight.items():
        comp[uf.find(key)] += w
    largest = max(comp.values())
    return n_pep, len(comp), largest / n_pep


def main():
    fasta = sys.argv[1]
    seqs = read_fasta(fasta)
    print(f"{len(seqs)} proteins, {sum(map(len, seqs)):,} residues\n")

    subsets = [int(x) for x in (sys.argv[2].split(",") if len(sys.argv) > 2
                                else ["200", "1000", "5000", "20594"])]
    print(f"{'k':>3} {'proteins':>9} {'peptides':>12} {'components':>12} "
          f"{'largest comp':>13}  verdict")
    for k in (3, 4, 5, 6, 7, 8):
        for n in subsets:
            npep, ncomp, frac = percolate(seqs[:n], k)
            verdict = ("PERCOLATED" if frac > 0.5 else
                       "large" if frac > 0.1 else "fragmented")
            print(f"{k:>3} {n:>9} {npep:>12,} {ncomp:>12,} "
                  f"{frac:>12.1%}  {verdict}")
        print()


if __name__ == "__main__":
    main()
