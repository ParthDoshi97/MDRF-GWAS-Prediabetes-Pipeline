#!/usr/bin/env python3
"""
Step 5: Parse and synthesize multi-database biological & clinical annotations
from Ensembl VEP, GTEx v8/v10, Open Targets Platform (L2G & Druggability),
Reactome Pathway Analysis, and ClinVar.
"""

import os
import json
import glob

def summarize_vep():
    vep_files = sorted(glob.glob("results/vep_*.json"))
    vep_summary = {}
    for vf in vep_files:
        rsid = os.path.basename(vf).replace("vep_", "").replace(".json", "")
        with open(vf) as f:
            data = json.load(f)
        if not data:
            continue
        rec = data[0]
        most_severe = rec.get("most_severe_consequence")
        tcs = rec.get("transcript_consequences", [])
        genes = sorted(list({tc.get("gene_symbol") for tc in tcs if tc.get("gene_symbol")}))
        # Find MANE Select or canonical transcript
        primary_tc = []
        for tc in tcs:
            if tc.get("mane_select") or tc.get("canonical") == 1:
                primary_tc.append({
                    "transcript_id": tc.get("transcript_id"),
                    "gene_symbol": tc.get("gene_symbol"),
                    "gene_id": tc.get("gene_id"),
                    "biotype": tc.get("biotype"),
                    "consequence_terms": tc.get("consequence_terms"),
                    "impact": tc.get("impact"),
                    "conservation": tc.get("conservation"),
                    "loeuf": tc.get("loeuf"),
                    "dosage_haplo": tc.get("dosage_sensitivity_haplo"),
                    "dosage_triplo": tc.get("dosage_sensitivity_triplo"),
                    "mane_select": tc.get("mane_select")
                })
        if not primary_tc and tcs:
            tc = tcs[0]
            primary_tc.append({
                "transcript_id": tc.get("transcript_id"),
                "gene_symbol": tc.get("gene_symbol"),
                "gene_id": tc.get("gene_id"),
                "biotype": tc.get("biotype"),
                "consequence_terms": tc.get("consequence_terms"),
                "impact": tc.get("impact"),
                "conservation": tc.get("conservation"),
                "loeuf": tc.get("loeuf"),
                "dosage_haplo": tc.get("dosage_sensitivity_haplo"),
                "dosage_triplo": tc.get("dosage_sensitivity_triplo")
            })
        vep_summary[rsid] = {
            "most_severe_consequence": most_severe,
            "overlapping_genes": genes,
            "primary_transcripts": primary_tc
        }
    return vep_summary

def summarize_reactome():
    if not os.path.exists("results/reactome_enrichment.json"):
        return []
    with open("results/reactome_enrichment.json") as f:
        data = json.load(f)
    pathways = []
    for p in data.get("pathways", [])[:15]:
        pathways.append({
            "stId": p.get("stId"),
            "name": p.get("name"),
            "pValue": p.get("entities", {}).get("pValue"),
            "fdr": p.get("entities", {}).get("fdr"),
            "found": p.get("entities", {}).get("found"),
            "total": p.get("entities", {}).get("total")
        })
    return pathways

def summarize_gtex():
    gtex_top = {}
    for gf in sorted(glob.glob("results/gtex_top_*.json")):
        gene = os.path.basename(gf).replace("gtex_top_", "").replace(".json", "")
        with open(gf) as f:
            gtex_top[gene] = json.load(f)
    return gtex_top

def summarize_opentargets():
    ot_summary = {"druggability": {}, "l2g": {}}
    for df in sorted(glob.glob("results/ot_*_druggability.json")):
        gene = os.path.basename(df).replace("ot_", "").replace("_druggability.json", "")
        with open(df) as f:
            ot_summary["druggability"][gene] = json.load(f)
    for lf in sorted(glob.glob("results/ot_l2g_*.json")):
        rsid = os.path.basename(lf).replace("ot_l2g_", "").replace(".json", "")
        with open(lf) as f:
            ot_summary["l2g"][rsid] = json.load(f)
    return ot_summary

if __name__ == "__main__":
    out = {
        "vep": summarize_vep(),
        "reactome_top_pathways": summarize_reactome(),
        "gtex_top_tissues": summarize_gtex(),
        "opentargets": summarize_opentargets()
    }
    with open("results/task5_biological_interpretation.json", "w") as f:
        json.dump(out, f, indent=2)
    print("Saved results/task5_biological_interpretation.json")
