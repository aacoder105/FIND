# FIND (Founder candidates hidden IN Data)

FIND is a lightweight web-based tool for identifying population-enriched pathogenic and loss-of-function variants in gnomAD that may represent candidate founder mutations.

The tool applies simple frequency-based filtering criteria to highlight variants enriched in a single population relative to others.

## Features

- Browser-based interface (no installation required)
- Query up to 10 genes simultaneously
- Screens pathogenic, likely pathogenic, and loss-of-function variants
- Identifies variants with strong population enrichment
- Highlights variants with a unique predominant population
- Export results to Excel format

## Filtering Criteria

FIND identifies variants that:

- are classified as pathogenic, likely pathogenic, or loss-of-function
- have allele frequency ≥ 0.00008 in one population
- are at least 10-fold more common in one population than any other
- treat populations with ≤4 alleles as frequency zero

## Web Interface

Users enter gene symbols and receive a table displaying:

- variant identifiers
- allele frequencies across populations
- predominant enriched population
- downloadable results

## Example Use Cases

FIND recovered known founder mutations in genes including:

- BRCA2
- MYH7
- FLNC
- TMEM127

and identified additional candidate founder variants for further investigation.

## Access

Web application:
https://ethnic-variant-mutation-finder.onrender.com/

When searching multiple genes, use the format "gene1, gene2," exactly.

Preprint:
https://www.biorxiv.org/content/10.64898/2026.06.05.730273v1

## Citation

Horowitz, Aaron; Liebman, Adam; Liebman, Susan. (2026). FIND: a software tool for identifying population-enriched pathogenic variants in gnomAD. bioRxiv, https://doi.org/10.64898/2026.06.05.730273. 

## Contact

Aaron L. Horowitz
aaron.l.horowitz@gmail.com<mailto:aaron.l.horowitz@gmail.com>


