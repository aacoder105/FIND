"""Revision analyses for the FIND manuscript (PONE-D-26-34008).

Input: gnomAD v4.1 browser CSV exports (one per gene, "Export variants to CSV",
all consequence boxes and "Filtered variants" ticked), placed in a folder as
FLNC.csv, TMEM127.csv, MYH7.csv, BRCA1.csv, BRCA2.csv.

Usage:  python revision_analyses.py <data_folder>

Outputs to the screen:
 1. FIND reproduction check (should list 26 variants)
 2. Fisher's exact tests for the 20 Table 1 variants (S1 Table)
 3. Threshold sweep and extra variants at relaxed settings (S2 Table)
 4. Permutation test, 1,000 permutations, seed 20261015 (S3 Table)
Requires: Python 3, pandas, numpy, scipy.
"""
import sys, warnings, numpy as np, pandas as pd
from scipy.stats import fisher_exact
import find_rules as f
warnings.filterwarnings('ignore')
data = sys.argv[1] if len(sys.argv) > 1 else '.'
G = ['FLNC','TMEM127','MYH7','BRCA1','BRCA2']
T1 = {'FLNC':['c.3791-1G>C'],'TMEM127':['c.245-1G>C'],
 'MYH7':['p.Arg1053Gln','p.Arg869His','p.Glu497Asp'],
 'BRCA1':['p.Gln1756ProfsTer74','p.Glu23ValfsTer17','p.Leu63Ter','p.Glu797ThrfsTer3'],
 'BRCA2':['p.Asn257LysfsTer17','p.Gln1089SerfsTer10','p.Val3079PhefsTer4','p.Ser1982ArgfsTer22','p.Glu1308Ter','p.Ile1859LysfsTer3','p.Leu2092ProfsTer7','p.Arg2318Ter','p.Arg2494Ter','c.9118-2A>G','p.Cys3233TrpfsTer15']}
A = pd.concat([f.load(f'{data}/{g}.csv').assign(gene=g) for g in G], ignore_index=True)
A['t1'] = [h in T1[g] for g, h in zip(A.gene, A['HGVS Consequence'])]
ac = A[[f'Allele Count {p}' for p in f.POPS]].to_numpy(float)
an = A[[f'Allele Number {p}' for p in f.POPS]].to_numpy(float)
print(f'{len(A)} P/LP or canonical-transcript pLoF variants screened; Bonferroni cutoff {0.05/len(A):.2e}')

def fisher(i, j):
    a, n = ac[i, j], an[i, j]; ao, no = ac[i].sum() - a, an[i].sum() - n
    return fisher_exact([[a, n - a], [ao, no - ao]], alternative='greater')[1]

ok, pop, taf, oth = f.apply(A)
print('\n1. FIND reproduction at defaults (5 alleles, tenfold, AF > 0.00008):', ok.sum(), 'variants')
for i in np.where(ok)[0]:
    print(f'   {A.gene[i]:8s} {A["HGVS Consequence"][i]:22s} {pop[i]:26s} in Table 1: {A.t1[i]}')

print('\n2. Fisher tests, Table 1 variants')
for i in np.where(A.t1)[0]:
    j = f.POPS.index(pop[i]); fold = 'inf' if oth[i] == 0 else f'{taf[i]/oth[i]:.1f}'
    print(f'   {A.gene[i]:8s} {A["HGVS Consequence"][i]:22s} fold(FIND) {fold:>6s}  p = {fisher(i, j):.2e}')

print('\n3. Threshold sweep (Table 1 variants retained of 20; total returned)')
for mc, fo in [(3,10),(4,10),(5,10),(6,10),(8,10),(10,10),(5,5),(5,15),(5,20),(5,50)]:
    o = f.apply(A, mc, fo)[0]
    print(f'   min count {mc:2d}, fold {fo:2d}: total {o.sum():2d}, Table 1 retained {(o & A.t1).sum():2d}')
for lab, (mc, fo) in {'4 alleles, tenfold': (4,10), '5 alleles, fivefold': (5,5)}.items():
    o, pp, _, _ = f.apply(A, mc, fo)
    print(f'   Extra variants at {lab}:')
    for i in np.where(o & ~ok)[0]:
        print(f'     {A.gene[i]:6s} {A["HGVS Consequence"][i]:22s} {pp[i]:26s} ClinVar: {A["ClinVar Germline Classification"][i]}  p = {fisher(i, f.POPS.index(pp[i])):.1e}')

print('\n4. Permutation test (1,000 permutations, seed 20261015)')
rng = np.random.default_rng(20261015)
tot = ac.sum(1).astype(int); prob = an / an.sum(1, keepdims=True)
S = {(5,10): [], (4,10): [], (5,5): []}
for it in range(1000):
    perm = np.array([rng.multinomial(t, p) for t, p in zip(tot, prob)], float)
    P = A.copy()
    for k, pp in enumerate(f.POPS): P[f'Allele Count {pp}'] = perm[:, k]
    for key in S: S[key].append(f.apply(P, *key)[0].sum())
for (mc, fo), v in S.items():
    v = np.array(v)
    print(f'   {mc} alleles, fold {fo}: mean {v.mean():.3f} hits; {(v>0).sum()} of 1000 permutations with >=1 hit; distribution {np.bincount(v).tolist()}')
