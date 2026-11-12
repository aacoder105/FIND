FIND revision analyses (PLOS ONE PONE-D-26-34008)

What this does
  Reproduces the numbers in S1-S3 Tables of the revised manuscript:
  Fisher's exact tests for the 20 Table 1 variants, the threshold sweep,
  and the permutation test.

Files
  find_rules.py          FIND's filtering rules, re-implemented from the Methods
  revision_analyses.py   runs all analyses and prints the results
  gnomad_data/           gnomAD browser CSV exports used (downloaded Oct 15, 2026;
                         browser v4.1.2, variant frequencies unchanged from v4.1.1)

How to run
  Requires Python 3 with pandas, numpy and scipy.
  From this folder:  python revision_analyses.py gnomad_data
  Takes about one minute (the permutations are the slow part).

pLoF rule
  As in FIND's app.py, pLoF variants (stop gained, frameshift, splice acceptor, splice
  donor) are screened only when the pLoF consequence is on the gene's canonical
  transcript, whatever their ClinVar classification. 2,057 variants are screened.

Check built in
  Step 1 of the output must list the same 26 variants that FIND itself returns
  for FLNC, TMEM127, MYH7, BRCA1 and BRCA2. The remaining analyses rely on that match.

Written with AI assistance (Claude, Anthropic); see the AI-use statement in Methods.
