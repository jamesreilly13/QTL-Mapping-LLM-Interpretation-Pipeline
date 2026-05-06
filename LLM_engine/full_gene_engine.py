# full_gene_engine.py

import pickle
import pandas as pd
from sentence_transformers import SentenceTransformer
import faiss
import numpy as np
import pickle
import faiss
from sentence_transformers import SentenceTransformer
from llama_cpp import Llama
import requests
import time
from collections import defaultdict
import re


MODEL_NAME = "pritamdeka/BioBERT-mnli-snli-scinli-scitail-mednli-stsb"

PHENOTYPE_CONTEXT = {
    "NPP_Exposure_ActivityCounts_GreyZone_Total":
        "locomotor activity in a neutral zone during novelty exposure",

    "NPP_Exposure_ExplorationCounts_GreyZone_Total":
        "exploratory investigation of a neutral zone during novelty exposure",

    "NPP_Exposure_EntranceCounts_BlackZone_Total":
        "entries into the dark zone during novelty exposure",

    "NPP_Exposure_EntranceCounts_WhiteZone_Total":
        "entries into the illuminated zone during novelty exposure",

    "NPP_Exposure_EntranceCounts_GreyZone_Total":
        "entries into the neutral zone during novelty exposure",

    "NPP_Exposure_ZoneTime_BlackZone_Total":
        "time spent in the dark zone during novelty exposure",

    "NPP_Exposure_ZoneTime_GreyZone_Total":
        "time spent in a neutral zone during novelty exposure",

    "NPP_NoveltyPreference_EntranceCounts_GreyWhiteBlack_Total":
        "zone entry counts during novelty preference testing",

    "NPP_NoveltyPreference_EntranceCounts_WhiteVsBlack_Total":
        "relative preference for illuminated versus dark zones",

    "NPP_NoveltyPreference_ExplorationCounts_GreyWhiteBlack_Total":
        "exploratory investigation across zones during novelty preference testing",

    "NPP_Test_ActivityCounts_GreyZone_Total":
        "locomotor activity in a neutral zone during novelty preference testing",

    "NPP_Test_ActivityCounts_WhiteZone_Total":
        "locomotor activity in an illuminated zone during novelty preference testing",

    "NPP_Test_ExplorationCounts_BlackZone_Total":
        "exploration of the dark zone during novelty preference testing",

    "NPP_Test_ExplorationCounts_WhiteZone_Total":
        "exploration of the illuminated zone during novelty preference testing",

    "LD_pct_time_in_light":
        "percentage of time spent in the illuminated compartment of the light-dark box",

    "LD_pct_ambulatory_Counts_in_light":
        "proportion of locomotor activity occurring in the illuminated compartment",

    "LD_pct_distance_traveled_in_light":
        "percentage of locomotor distance traveled in the illuminated compartment",

    "LD_pct_resting_time_in_light":
        "proportion of resting behavior in the illuminated compartment",

    "LD_total_transitions":
        "number of transitions between light and dark compartments",
    
    "OFA_total_distance_traveled_in_center":
        "total distance traveled in the center of the open field",
    
    "OFA_total_distance_traveled_in_corner":
        "total distance traveled in the corners of the open field",
    
    "OFA_total_distance_traveled_in_perimeter":
        "total distance traveled along arena walls during open field testing",
    
    "OFA_total_resting_time_in_corner":
        "resting behavior in the corners of the open field",
    
    "OFA_total_resting_time_in_perimeter":
        "resting behavior along arena walls during open field testing",

    "OFA_pct_time_center":
        "percentage of time spent in the center of the open field",

    "OFA_total_time_in_center":
        "time spent in the center region of the open field",

    "OFA_total_resting_time_in_center":
        "resting behavior in the center of the open field",

    "OFA_pct_time_corner":
        "percentage of time spent in the corners of the open field",

    "OFA_total_time_in_corner":
        "time spent in the corners of the open field",

    "OFA_pct_resting_time_corner":
        "resting behavior in the corners of the open field",

    "OFA_pct_time_perimeter":
        "percentage of time spent along arena walls",

    "OFA_total_time_in_perimeter":
        "time spent along arena walls",

    "OFA_total_distance_traveled":
        "total locomotor distance traveled during open field testing",

    "OFA_total_ambulatory_time":
        "time spent in active locomotion",

    "OFA_total_resting_time":
        "total time spent immobile",

    "OFA_distance_traveled_first_five":
        "locomotor activity during initial exposure to the open field",

    "OFA_distance_traveled_last_five":
        "locomotor activity during the later phase of the session",

    "OFA_distance_traveled_slope":
        "rate of change in locomotor activity across the open field session",

    "OFA_pct_dist_center":
        "percentage of locomotor distance traveled in the center",

    "OFA_pct_dist_corner":
        "percentage of locomotor distance traveled in corner regions",

    "OFA_pct_dist_perimeter":
        "percentage of locomotor distance traveled along arena walls",

    "HB_Total_Entries":
        "total number of exploratory head-dip events in the hole board test",

    "HB_Novel_Entries":
        "number of exploratory head dips into previously unvisited holes",

    "HB_Repeat_Entries":
        "number of repeated head dips into previously visited holes"
}

class FullGeneEngine:

    def __init__(self):

        MODEL_NAME = "pritamdeka/BioBERT-mnli-snli-scinli-scitail-mednli-stsb"

        with open("data/qtl_data.pkl", "rb") as f:
            self.qtl_data = pickle.load(f)

        self.model = SentenceTransformer(MODEL_NAME)

        self.pubmed_index = faiss.read_index("data/pubmed_index.faiss")

        with open("data/pubmed_metadata.pkl", "rb") as f:
            self.pubmed_records = pickle.load(f)

        self.llm = Llama(
            model_path="models/Meta-Llama-3.1-8B-Instruct-Q4_K_M.gguf",
            n_ctx=8192,
            n_threads=8,
            n_gpu_layers=50,
            verbose=False
        )

    # phenotype helpers
    def describe_phenotype(self, code):
        return PHENOTYPE_CONTEXT.get(code, code)
    
    def phenotype_reverse_map(self):
        return {v: k for k, v in PHENOTYPE_CONTEXT.items()}

    # --------------------------------------------------
    # GET FULL GENE TABLE
    # --------------------------------------------------

    def get_gene_table(self, group=None, phenotype=None):

        dfs = []

        for g, tables in self.qtl_data.items():

            if group and g != group:
                continue

            genes = tables.get("genes")

            if genes is None:
                continue

            df = pd.DataFrame(genes)

            if phenotype and "phenotype" in df.columns:
                df = df[df["phenotype"] == phenotype]

            df["group"] = g

            dfs.append(df)

        if not dfs:
            return pd.DataFrame()

        return pd.concat(dfs, ignore_index=True)
    

    def get_go_table(self, group=None, phenotype=None):

        dfs = []

        for g, tables in self.qtl_data.items():

            if group and g != group:
                continue

            go = tables.get("go")

            if go is None:
                continue

            df = pd.DataFrame(go)

            if phenotype and "phenotype" in df.columns:
                df = df[df["phenotype"] == phenotype]

            dfs.append(df)

        return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
    

    def filter_go_by_genes(self, go_df, gene_list, min_overlap=2):

        if go_df.empty:
            return pd.DataFrame()

        gene_set = set(gene_list)

        rows = []

        for _, r in go_df.iterrows():

            raw = r.get("geneID")

            if not isinstance(raw, str):
                continue

            go_genes = set(raw.split("/"))

            overlap = gene_set & go_genes

            if len(overlap) >= min_overlap:
                rows.append({
                    "Description": r.get("Description"),
                    "gene_overlap": list(overlap),
                    "overlap_count": len(overlap),
                    "FoldEnrichment": r.get("FoldEnrichment"),
                    "qvalue": r.get("qvalue")
                })

        df = pd.DataFrame(rows)

        if df.empty:
            return df

        return df.sort_values(
            ["overlap_count", "qvalue"],
            ascending=[False, True]
        )
    

    def build_go_context(self, go_df, max_terms=10):

        if go_df.empty:
            return "No enriched biological processes identified."

        lines = []

        for _, r in go_df.head(max_terms).iterrows():

            genes = ", ".join(r["gene_overlap"])

            line = (
                f"{r['Description']} "
                f"(genes: {genes}; "
                f"enrichment: {r['FoldEnrichment']}, "
                f"qvalue={r['qvalue']})"
            )

            lines.append(line)

        return "\n".join(lines)


    def classify_variant_impact(self, consequence):

        if pd.isna(consequence):
            return "Unknown"

        terms = {t.strip() for t in str(consequence).split(",")}

        HIGH = {
            "stop_gained","stop_lost","frameshift_variant",
            "splice_donor_variant","splice_acceptor_variant",
            "initiator_codon_variant"
        }

        MODERATE = {"missense_variant","inframe_insertion"}

        LOW = {"synonymous_variant","stop_retained_variant"}

        REG = {
            "upstream_gene_variant","downstream_gene_variant",
            "5_prime_UTR_variant","3_prime_UTR_variant",
            "intron_variant","splice_region_variant"
        }

        if terms & HIGH:
            return "High impact"

        if terms & MODERATE:
            return "Moderate impact"

        if terms & LOW:
            return "Low impact"

        if terms & REG:
            return "Regulatory"

        if "intergenic_variant" in terms:
            return "Intergenic"

        return "Other"


    def get_variant_table(self, group=None, phenotype=None):

        dfs = []

        for g, tables in self.qtl_data.items():

            if group and g != group:
                continue

            snps = tables.get("snps")

            if snps is None:
                continue

            df = pd.DataFrame(snps)

            if phenotype and "phenotype" in df.columns:
                df = df[df["phenotype"] == phenotype]

            dfs.append(df)

        return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


    def build_variant_summary_table(self, df, gene_list):

        if df.empty:
            return pd.DataFrame()

        df = df.copy()

        # classify impact
        df["impact_category"] = df["consequence"].apply(self.classify_variant_impact)

        df = df[df["gene_name"].isin(gene_list)]

        if df.empty:
            return pd.DataFrame()

        summaries = []

        for gene, g in df.groupby("gene_name"):

            g = g.dropna(subset=["lod"])

            if g.empty:
                continue

            g_sorted = g.sort_values("lod", ascending=False)

            summaries.append({
                "gene": gene,
                "variant_count": len(g),
                "mean_lod": g["lod"].mean(),
                "max_lod": g["lod"].max(),

                "high_impact": (g["impact_category"] == "High impact").sum(),
                "moderate_impact": (g["impact_category"] == "Moderate impact").sum(),
                "low_impact": (g["impact_category"] == "Low impact").sum(),

                "top_snps": ", ".join(
                    g_sorted["snp"].dropna().astype(str).head(5)
                )
            })

        summary_df = pd.DataFrame(summaries)

        if summary_df.empty:
            return summary_df

        return summary_df.sort_values("max_lod", ascending=False)


    def get_coef_table(self, group=None, phenotype=None):

        dfs = []

        for g, tables in self.qtl_data.items():

            if group and g != group:
                continue

            coef = tables.get("coef")

            if coef is None:
                continue

            df = pd.DataFrame(coef)

            if phenotype and "phenotype" in df.columns:
                df = df[df["phenotype"] == phenotype]

            dfs.append(df)

        return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

    # --------------------------------------------------
    # COMPRESSED GENE CONTEXT (NAME + DESCRIPTION)
    # --------------------------------------------------

    def build_context(self, df, max_genes=300):

        # ---- SAFETY CHECKS ----
        if "Name" not in df.columns:
            raise ValueError("Gene table missing 'Name' column")

        # Drop bad rows
        df = df.dropna(subset=["Name"])

        if df.empty:
            return [], "No valid gene names found."

        # ---- GENES ----
        genes = df["Name"].astype(str).unique().tolist()
        gene_subset = genes[:max_genes]

        gene_text = ", ".join(gene_subset)

        # ---- DESCRIPTIONS (SAFE) ----
        if "description" in df.columns:
            desc_df = (
                df[["Name", "description"]]
                .dropna()
                .drop_duplicates()
                .head(50)
            )

            desc_text = "\n".join(
                f"{str(r['Name'])}: {str(r['description'])}"
                for _, r in desc_df.iterrows()
            )
        else:
            desc_text = "No gene descriptions available."

        context = f"""
    Total genes: {len(genes)}

    Gene list (subset):
    {gene_text}

    Gene descriptions:
    {desc_text}
    """

        return genes, context


    # get available phenotypes for a group or overall
    def get_available_phenotypes(self, group=None):

        phenotypes = set()

        for g, tables in self.qtl_data.items():

            if group and g != group:
                continue

            genes = tables.get("genes")
            if genes is None:
                continue

            df = pd.DataFrame(genes)

            if "phenotype" in df.columns:
                phenotypes.update(df["phenotype"].dropna().unique())

        # Convert to human-readable
        return sorted([
            self.describe_phenotype(p) for p in phenotypes
        ])

    # --------------------------------------------------
    # PUBMED RETRIEVAL
    # --------------------------------------------------

    def aggregate_gene_literature(self, pubmed_results):
        """
        Converts flat PubMed results into gene-centered evidence summaries.
        """

        ADDICTION_TERMS = {
            "addiction", "reward", "dopamine", "cocaine", "opioid",
            "alcohol", "reinforcement", "dependence", "self-administration",
            "behavior", "locomotor"
        }

        gene_map = defaultdict(lambda: {
            "papers": [],
            "pmids": set(),
            "addiction_hits": 0,
            "total_hits": 0,
            "themes": defaultdict(int)
        })

        for r in pubmed_results:

            gene = r.get("gene")
            if not gene:
                continue

            title = (r.get("title") or "").lower()
            abstract = (r.get("abstract") or "").lower()
            text = title + " " + abstract

            gene_entry = gene_map[gene]

            gene_entry["papers"].append(r)
            gene_entry["pmids"].add(r["pmid"])
            gene_entry["total_hits"] += 1

            # addiction relevance
            if any(term in text for term in ADDICTION_TERMS):
                gene_entry["addiction_hits"] += 1

            # simple theme extraction (cheap but effective)
            for term in ADDICTION_TERMS:
                if term in text:
                    gene_entry["themes"][term] += 1

        # ---- finalize structure ----
        output = {}

        for gene, data in gene_map.items():

            paper_count = len(data["pmids"])
            addiction_hits = data["addiction_hits"]

            # simple evidence score (you can refine later)
            score = (
                0.5 * paper_count +
                2.0 * addiction_hits
            )

            # normalize (soft scaling)
            score = min(score / 20.0, 1.0)

            top_pmids = list(data["pmids"])[:5]

            top_themes = sorted(
                data["themes"].items(),
                key=lambda x: x[1],
                reverse=True
            )[:5]

            output[gene] = {
                "paper_count": paper_count,
                "addiction_hits": addiction_hits,
                "top_pmids": top_pmids,
                "top_themes": [t[0] for t in top_themes],
                "evidence_score": round(score, 3)
            }

        return output


    def build_gene_query(self, gene, phenotype_text):
        return f"""
        ({gene}[Title/Abstract])
        AND (mouse OR murine)
        AND (brain OR behavior OR neuroscience OR {phenotype_text})
        """


    def pubmed_esearch(self, query, retmax=10):
        base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

        params = {
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retmax": retmax
        }

        r = requests.get(base, params=params)
        r.raise_for_status()

        return r.json()["esearchresult"]["idlist"]


    def pubmed_efetch(self, pmids):
        if not pmids:
            return []

        base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

        params = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "xml"
        }

        r = requests.get(base, params=params)
        r.raise_for_status()

        # lightweight XML parsing
        import xml.etree.ElementTree as ET

        root = ET.fromstring(r.text)

        results = []

        for article in root.findall(".//PubmedArticle"):
            pmid = article.findtext(".//PMID")

            title = article.findtext(".//ArticleTitle")
            abstract = article.findtext(".//AbstractText")

            results.append({
                "pmid": pmid,
                "title": title,
                "abstract": abstract
            })

        return results


    def retrieve_pubmed(self, genes, phenotype):
        
        ADDIC_KEYWORDS = [
            "addiction", "reward", "dopamine", "substance",
            "cocaine", "alcohol", "opioid", "behavior",
            "reinforcement", "dependence", "morphine",
            "metamphetamine"
        ]

        def score_paper(title, abstract, gene):
            text = (title or "") + " " + (abstract or "")
            text = text.lower()

            score = 0

            if gene.lower() in text:
                score += 2

            score += sum(k in text for k in ADDIC_KEYWORDS)

            return score
        
        results = {}
        phenotype_text = self.describe_phenotype(phenotype)

        for gene in genes[:10]:

            query = self.build_gene_query(gene, phenotype_text)

            try:
                pmids = self.pubmed_esearch(query, retmax=5)
                papers = self.pubmed_efetch(pmids)

                for p in papers:
                    pmid = p.get("pmid")

                    if not pmid:
                        continue

                    # avoid duplicates
                    if pmid not in results:
                        results[pmid] = {
                            "pmid": pmid,
                            "title": p.get("title"),
                            "abstract": p.get("abstract"),
                            "gene": gene,
                            "score": score_paper(p.get("title"), p.get("abstract"), gene)
                        }

                time.sleep(0.34)  # be polite to NCBI servers

            except Exception as e:
                print(f"PubMed error for {gene}: {e}")
                continue

        return sorted(results.values(), key=lambda x: x["score"], reverse=True)[:10]

    # --------------------------------------------------
    # ANSWER
    # --------------------------------------------------

    def answer(self, query, group=None, phenotype=None):

        # -------- GENES --------
        df = self.get_gene_table(group, phenotype)

        if df.empty:
            return {"answer": "No gene data found."}

        genes, context = self.build_context(df)

        # -------- GO --------
        go_df = self.get_go_table(group, phenotype)

        filtered_go = self.filter_go_by_genes(go_df, genes)

        go_context = self.build_go_context(filtered_go)

        # -------- VARIANTS --------
        raw_variants = self.get_variant_table(group, phenotype)
        variant_df = self.build_variant_summary_table(raw_variants, genes)

        # -------- COEFFICIENTS --------
        coef_df = self.get_coef_table(group, phenotype)

        # -------- PUBMED --------
        pubmed = self.retrieve_pubmed(genes, phenotype)
        gene_lit_summary = self.aggregate_gene_literature(pubmed)

        pubmed_text = "\n\n".join(
            f"""GENE: {gene}
        Paper count: {v['paper_count']}
        Addiction-related hits: {v['addiction_hits']}
        Evidence score: {v['evidence_score']}
        Top themes: {", ".join(v['top_themes'])}
        Top PMIDs: {", ".join(v['top_pmids'])}
        """
            for gene, v in gene_lit_summary.items()
        )

        # -------- LLM PROMPT --------
        prompt = f"""
You are an expert in addiction genetics.

Your task is to produce a biologically rigorous interpretation of genetic and pathway data related to a behavioral phenotype.

This is NOT a question-answering task.
This is NOT a numerical task.
Do NOT produce a "final answer" or any numeric output.

--------------------------------------------------
QUESTION
{query}

PHENOTYPE
{phenotype}

--------------------------------------------------
GENE DATA SUMMARY
{context}

ENRICHED BIOLOGICAL PROCESSES
{go_context}

LITERATURE EVIDENCE (GENE-LEVEL AGGREGATED)
{pubmed_text}

--------------------------------------------------
HOW TO USE THE LITERATURE DATA

Each gene includes:

- Paper count → overall study depth
- Addiction hits → direct relevance to addiction biology
- Evidence score (0–1) → strength of support
- Top themes → dominant biological functions

MANDATORY RULES:

1. You MUST use the literature data when interpreting genes
2. If evidence_score > 0.3 → treat as WELL-SUPPORTED and explicitly state this
3. If addiction_hits > 0 → explicitly describe addiction-related relevance
4. You MUST NOT claim "no evidence" if addiction_hits > 0
5. You MUST distinguish:
   - well-supported genes
   - moderately supported genes
   - poorly studied genes
6. Base conclusions ONLY on provided data (do NOT rely on generic assumptions)

--------------------------------------------------
OUTPUT FORMAT (STRICT)

- Do NOT include titles like "Answer" or "Final Answer"
- Do NOT include any numeric-only lines (e.g., "0")
- Do NOT include placeholders or meta commentary
- Do NOT repeat the question
- Do NOT list all genes

Your response MUST begin EXACTLY with:

1. Key biological themes

Then continue with EXACTLY these sections:

2. Mechanisms
3. Literature integration
4. Interpretation

--------------------------------------------------
WRITING STYLE

- Be concise but biologically precise
- Prioritize mechanism over description
- Synthesize pathways, not gene lists
- Only mention specific genes when they illustrate a mechanism or are supported by literature

Begin now.
"""

        response = self.llm(
            prompt,
            max_tokens=800,
            temperature=0.2
        )

        answer = response["choices"][0]["text"]

        return {
            "answer": answer,
            "genes": df,
            "variants": variant_df,
            "coef": coef_df,
            "go": go_df,
            "gene_lit_summary": gene_lit_summary,
            "pubmed_text": pubmed_text
        }
