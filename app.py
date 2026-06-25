#!/usr/bin/env python3
"""
Flask web server for GnomAD variant analysis
Provides API endpoint for fetching pathogenic/pLoF variants
ONLY returns variants with a unique population
"""

from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
import requests
from collections import defaultdict
from typing import List, Dict, Any

app = Flask(__name__)
CORS(app)

# ---------------- CONFIG ----------------
API_URL = "https://gnomad.broadinstitute.org/api"
REFERENCE_GENOME = "GRCh38"
DATASET = "gnomad_r4"

HEADERS = {"Content-Type": "application/json", "User-Agent": "GnomAD-Variant-Analyzer/1.0"}

# pLoF terms
PLOF_TERMS = {
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
}

# Put in here the populations you dont want
EXCLUDED_POPULATIONS = {"remaining"}

# population id -> human-readable name
POP_NAME_OVERRIDES = {
    "afr": "African/African American",
    "amr": "Admixed American (Latino)",
    "eas": "East Asian",
    "nfe": "Non-Finnish European",
    "fin": "Finnish",
    "sas": "South Asian",
    "oth": "Other",
    "asj": "Ashkenazi Jewish",
    "mid": "Middle Eastern",
    "ami": "Amish",
}

# GraphQL query
GRAPHQL_QUERY = """
query($gene_symbol: String!, $reference_genome: ReferenceGenomeId!, $dataset: DatasetId!) {
  gene(gene_symbol: $gene_symbol, reference_genome: $reference_genome) {
    gene_id
    canonical_transcript_id
    variants(dataset: $dataset) {
      variant_id
      rsid
      chrom
      pos
      transcript_consequence {
        consequence_terms
        gene_id
        transcript_id
      }
      joint { populations { id ac an } }
      genome { populations { id ac an } }
      exome { populations { id ac an } }
    }
    clinvar_variants {
      variant_id
      clinical_significance
    }
  }
}
"""


# ---------------- HELPERS ----------------
def run_graphql(query: str, variables: Dict[str, Any]) -> Dict[str, Any]:
    payload = {"query": query, "variables": variables}
    resp = requests.post(API_URL, json=payload, headers=HEADERS, timeout=60)
    try:
        return resp.json()
    except Exception:
        return {"errors": [{"message": "invalid json response", "text": resp.text}]}


def safe_int(x):
    try:
        return int(x)
    except Exception:
        return 0


def compute_af(ac: int, an: int) -> float:
    return (ac / an) if (an and an > 0) else 0.0


def is_excluded_population(pid: str) -> bool:
    """Return True for populations that should be dropped (e.g. 'remaining').

    Also covers the sex-stratified splits such as 'remaining_XX' / 'remaining_XY'
    so the bucket is removed in every form it can appear in.
    """
    if not pid:
        return True
    base = pid.split("_")[0]
    return base in EXCLUDED_POPULATIONS


def is_plof(node: Dict[str, Any], canonical_transcript_id: str = None) -> bool:
    """Check if variant is LoF based on consequence terms."""
    tc = node.get("transcript_consequence") or node.get("transcript_consequences")
    if not tc:
        return False
    entries = tc if isinstance(tc, list) else [tc]

    for entry in entries:
        if not isinstance(entry, dict):
            continue

        if canonical_transcript_id:
            transcript_id = entry.get("transcript_id")
            if transcript_id != canonical_transcript_id:
                continue

        terms = entry.get("consequence_terms")
        if not terms:
            continue
        if isinstance(terms, list):
            for t in terms:
                if t in PLOF_TERMS:
                    return True
        else:
            if terms in PLOF_TERMS:
                return True

    return False


def is_pathogenic(sig):
    """Check if ClinVar says pathogenic/likely pathogenic"""
    if not sig:
        return False
    s = str(sig).lower()
    if "benign" in s or "conflicting" in s or "uncertain" in s:
        return False
    return "pathogenic" in s


def gather_population_ac_an(node: Dict[str, Any]) -> Dict[str, Dict[str, int]]:
    pop_ac = defaultdict(int)
    pop_an = defaultdict(int)
    for block in ("joint",):
        blk = node.get(block)
        if not isinstance(blk, dict):
            continue
        pops = blk.get("populations") or []
        for p in pops:
            pid = p.get("id")
            if not pid:
                continue
            # Skip the "remaining" bucket so it cannot interfere.
            if is_excluded_population(pid):
                continue
            ac = safe_int(p.get("ac"))
            an = safe_int(p.get("an"))
            pop_ac[pid] += ac
            pop_an[pid] += an
    out = {}
    for pid in pop_ac:
        out[pid] = {"ac": pop_ac[pid], "an": pop_an.get(pid, 0)}
    return out


def process_gene(gene: str) -> List[Dict[str, Any]]:
    """Process a single gene and return variant results - ONLY variants with unique population."""
    rows: List[Dict[str, Any]] = []

    print(f"Fetching gene {gene} ...")
    variables = {"gene_symbol": gene, "reference_genome": REFERENCE_GENOME, "dataset": DATASET}
    resp = run_graphql(GRAPHQL_QUERY, variables)

    if "errors" in resp:
        print(f"GraphQL error for {gene}: {resp.get('errors')}")
        return rows

    data = resp.get("data", {}) or {}
    gene_obj = data.get("gene") or {}
    canonical_transcript_id = gene_obj.get("canonical_transcript_id")

    clinvar_list = gene_obj.get("clinvar_variants") or []
    clinvar_map = {c.get("variant_id"): c.get("clinical_significance") for c in clinvar_list}

    for var in gene_obj.get("variants") or []:
        vid = var.get("variant_id")
        if not vid:
            continue

        plof_flag = is_plof(var, canonical_transcript_id)
        path_flag = is_pathogenic(clinvar_map.get(vid))
        if not plof_flag and not path_flag:
            continue

        pop_map = gather_population_ac_an(var)
        pop_af_named: Dict[str, float] = {}
        pop_ac_named: Dict[str, int] = {}
        all_pop_names = set()

        for pid, vals in pop_map.items():
            ac = vals.get("ac", 0)
            an = vals.get("an", 0)
            pname = POP_NAME_OVERRIDES.get(pid, pid)
            all_pop_names.add(pname)

            if an and an > 0 and ac >= 5:
                af = compute_af(ac, an)
                pop_af_named[pname] = af
                pop_ac_named[pname] = ac

        # Determine unique_population
        unique_pop = ""
        if len(pop_af_named) > 0:
            sorted_pops = sorted(pop_af_named.items(), key=lambda x: x[1], reverse=True)
            max_pop, max_af = sorted_pops[0]
            if max_af > 0.00008:
                if len(sorted_pops) == 1:
                    unique_pop = max_pop
                elif len(sorted_pops) > 1:
                    second_af = sorted_pops[1][1]
                    if max_af >= 10 * second_af:
                        unique_pop = max_pop

        # FILTER: Only include variants with a unique population
        if not unique_pop:
            continue

        # Build row
        row: Dict[str, Any] = {"gene": gene, "variant_id": vid, "unique_population": unique_pop}
        for pname in sorted(all_pop_names):
            if pname in pop_af_named:
                row[pname] = f"{pop_af_named[pname]:.6f}"
            else:
                row[pname] = "0.000000"
        rows.append(row)

    return rows


# ---------------- API ENDPOINTS ----------------
@app.route('/')
def home():
    return send_file('index.html')


@app.route('/api/analyze', methods=['POST'])
def analyze():
    """
    API endpoint to analyze genes
    Expects JSON: {"genes": ["BRCA1", "BRCA2"]}
    Returns JSON: {"results": [...], "total": 10}
    """
    try:
        data = request.get_json()
        genes = data.get('genes', [])

        if not genes:
            return jsonify({"error": "No genes provided"}), 400

        if len(genes) > 10:
            return jsonify({"error": "Maximum 10 genes allowed"}), 400

        all_results = []
        for gene in genes:
            gene = gene.strip().upper()
            if gene:
                gene_results = process_gene(gene)
                all_results.extend(gene_results)

        return jsonify({
            "results": all_results,
            "total": len(all_results)
        })

    except Exception as e:
        print(f"Error in analyze endpoint: {str(e)}")
        return jsonify({"error": str(e)}), 500


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
