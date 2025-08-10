#!/usr/bin/env python3
import csv
import re
from rapidfuzz import fuzz

# ------------------------------------------------------------------------------
# 1. FILE PATHS
# ------------------------------------------------------------------------------
FDA_FILE = "ANALYSIS/00-05-PyPGx/fda_pharmacogenomics_associations.csv"
ALL_GENES_FILE = "ANALYSIS/00-05-PyPGx/all_genes_results.tsv"

# ------------------------------------------------------------------------------
# 2. CATEGORY PATTERNS
#    We handle recognized synonyms for poor, intermediate, normal (extensive),
#    rapid, and ultrarapid metabolizers.
# ------------------------------------------------------------------------------
FDA_CATEGORY_PATTERNS = [
    (re.compile(r"\bpoor metabolizer(s)?\b", re.IGNORECASE), "PM"),
    (re.compile(r"\bintermediate metabolizer(s)?\b", re.IGNORECASE), "IM"),
    (re.compile(r"\bnormal metabolizer(s)?\b", re.IGNORECASE), "NM"),
    (re.compile(r"\bextensive metabolizer(s)?\b", re.IGNORECASE), "NM"),
    (re.compile(r"\brapid metabolizer(s)?\b", re.IGNORECASE), "RM"),
    (re.compile(r"\bultra[-\s]?rapid metabolizer(s)?\b", re.IGNORECASE), "UM"),
]

PGX_CATEGORY_PATTERNS = [
    (re.compile(r"\bpoor metabolizer(s)?\b", re.IGNORECASE), "PM"),
    (re.compile(r"\bintermediate metabolizer(s)?\b", re.IGNORECASE), "IM"),
    (re.compile(r"\bnormal metabolizer(s)?\b", re.IGNORECASE), "NM"),
    (re.compile(r"\bextensive metabolizer(s)?\b", re.IGNORECASE), "NM"),
    (re.compile(r"\brapid metabolizer(s)?\b", re.IGNORECASE), "RM"),
    (re.compile(r"\bultra[-\s]?rapid metabolizer(s)?\b", re.IGNORECASE), "UM"),
]

# ------------------------------------------------------------------------------
# 3. FUZZY MATCH CONFIG
#    We do partial-ratio on Phenotype vs. the FDA's Affected Subgroups text,
#    if category-based logic doesn't yield a match.
# ------------------------------------------------------------------------------
FUZZY_THRESHOLD = 85

# ------------------------------------------------------------------------------
# 4. CATEGORY PARSING FUNCTIONS
# ------------------------------------------------------------------------------
def extract_fda_categories(affected_subgroups_text: str):
    """
    Identify recognized categories (PM,IM,NM,RM,UM) from FDA 'Affected Subgroups+'.
    Returns a set of zero or more categories.
    """
    cats = set()
    text_lc = affected_subgroups_text.lower()
    for (pattern, label) in FDA_CATEGORY_PATTERNS:
        if pattern.search(text_lc):
            cats.add(label)
    return cats

def extract_pgx_category(phenotype_text: str):
    """
    Identify recognized categories (PM,IM,NM,RM,UM) from pypgx 'Phenotype'.
    Returns a set of zero or more categories (though typically one).
    """
    cats = set()
    text_lc = phenotype_text.lower()
    for (pattern, label) in PGX_CATEGORY_PATTERNS:
        if pattern.search(text_lc):
            cats.add(label)
    return cats

# ------------------------------------------------------------------------------
# 5. LOAD FDA CSV
# ------------------------------------------------------------------------------
def load_fda_csv(csv_path):
    """
    The CSV is expected to have columns:
      [Drug, Gene, Affected Subgroups+, Description of Gene-Drug Interaction].
    """
    rows = []
    with open(csv_path, mode="r", encoding="utf-8") as f:
        reader = csv.reader(f)
        # skip header
        next(reader, None)

        for line in reader:
            if len(line) < 4:
                continue
            drug = line[0].strip()
            gene = line[1].strip().upper()
            affected = line[2].strip()
            recommendation = line[3].strip()

            fda_cats = extract_fda_categories(affected)

            rows.append({
                "Drug": drug,
                "Gene": gene,
                "AffectedSubgroupsText": affected,
                "Categories": fda_cats,   # e.g. {'PM','IM'} etc.
                "Recommendation": recommendation
            })
    return rows

# ------------------------------------------------------------------------------
# 6. LOAD PYPgx TSV
# ------------------------------------------------------------------------------
def parse_population_and_number(sample_id: str):
    """
    e.g. "CHOPCCAS-250" => ("CHOPCCAS","250").
    If no dash, fallback to (sample_id, "").
    """
    if "-" in sample_id:
        pop, num = sample_id.split("-", 1)
        return pop, num
    else:
        return sample_id, ""

def load_pgx_tsv(tsv_path):
    """
    Returns list of dicts with:
     FullSampleID, PopulationName, IndividualNumber, Gene, Phenotype, PGxCategories (set)
    """
    data = []
    with open(tsv_path, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample_id = row["Unnamed: 0"]
            pop_name, indiv_num = parse_population_and_number(sample_id)

            gene = (row.get("Gene","") or "").upper()
            phenotype = row.get("Phenotype","") or ""
            cat = extract_pgx_category(phenotype)

            data.append({
                "FullSampleID": sample_id,
                "PopulationName": pop_name,
                "IndividualNumber": indiv_num,
                "Gene": gene,
                "Phenotype": phenotype,
                "Category": cat
            })
    return data

# ------------------------------------------------------------------------------
# 7. "No Contradiction" Check
#    e.g., if the user is "PM" but the FDA text is explicitly "normal metabolizers only"
#    with no mention of "poor," that's contradictory => skip.
# ------------------------------------------------------------------------------
def is_contradictory(user_cats: set, fda_cats: set):
    """
    If the FDA categories are a strict subset that does NOT include the user's category,
    it might be contradictory.
    - Example: user_cats = {'PM'} but fda_cats = {'NM'} => contradiction
    - Example: user_cats = {'PM'} but fda_cats = {'PM','IM','NM'} => not contradictory
      because PM is included in the FDA set.
    - If either side is empty, we don't treat it as contradictory (since we can't be sure).
    """
    if user_cats and fda_cats:
        # If there's no intersection => contradiction
        # (only if fda_cats is not "broad" enough to include user_cats)
        if user_cats.isdisjoint(fda_cats):
            return True
    return False

# ------------------------------------------------------------------------------
# 8. MATCHING LOGIC
#    We attempt:
#      (1) If user & FDA categories overlap => match
#      (2) Else fallback to fuzzy match if not contradictory
# ------------------------------------------------------------------------------
from collections import defaultdict

def match_fda_recommendations(fda_list, pgx_list):
    fda_by_gene = defaultdict(list)
    for f in fda_list:
        fda_by_gene[f["Gene"]].append(f)

    matches = []

    for p in pgx_list:
        gene = p["Gene"]
        if not gene:
            continue

        possible_rows = fda_by_gene.get(gene, [])
        if not possible_rows:
            continue

        user_cats = p["Category"]   # set of e.g. {"PM"}
        user_pheno = p["Phenotype"]

        for frow in possible_rows:
            fda_cats = frow["Categories"]  # e.g. {"PM","IM"}
            fda_text = frow["AffectedSubgroupsText"]

            # 1) Check if they share categories (non-empty intersection)
            intersection = user_cats.intersection(fda_cats)
            if intersection:
                # no contradiction => direct match
                # (If you'd like to also do fuzzy check for a tie-break, you can)
                matches.append(build_match_row(p, frow, method="Category-Overlap"))
                continue  # done with this FDA row

            # 2) If categories don't overlap, check if that's a contradiction
            #    e.g. user = PM, FDA = NM only => skip
            if is_contradictory(user_cats, fda_cats):
                # skip
                continue

            # 3) Fallback fuzzy match of phenotype vs. FDA text
            #    If high enough => consider it a match
            score = fuzz.partial_ratio(user_pheno.lower(), fda_text.lower())
            if score >= FUZZY_THRESHOLD:
                # We consider it a fallback match
                matches.append(build_match_row(p, frow, method=f"Fuzzy={score:.1f}"))

    return matches

def build_match_row(pgx_row, fda_row, method=""):
    return {
        "PopulationName": pgx_row["PopulationName"],
        "IndividualNumber": pgx_row["IndividualNumber"],
        "FullSampleID": pgx_row["FullSampleID"],
        "Gene": pgx_row["Gene"],
        "Phenotype": pgx_row["Phenotype"],
        "PGxCategory": ", ".join(sorted(pgx_row["Category"])),
        "FDA_Drug": fda_row["Drug"],
        "FDA_AffectedSubgroups": fda_row["AffectedSubgroupsText"],
        "FDA_Categories": ", ".join(sorted(fda_row["Categories"])),
        "FDA_Recommendation": fda_row["Recommendation"],
        "MatchMethod": method
    }

# ------------------------------------------------------------------------------
# 9. MAIN
# ------------------------------------------------------------------------------
def main():
    fda_data = load_fda_csv(FDA_FILE)
    pgx_data = load_pgx_tsv(ALL_GENES_FILE)

    matched = match_fda_recommendations(fda_data, pgx_data)

    # Output tab-delimited
    print("Population\tIndividual\tGene\tPhenotype\tPGxCategory\tFDA_Drug\tAffectedSubgroups\tFDA_Categories\tFDA_Recommendation\tMatchMethod")
    for m in matched:
        print(
            f"{m['PopulationName']}\t"
            f"{m['IndividualNumber']}\t"
            f"{m['Gene']}\t"
            f"{m['Phenotype']}\t"
            f"{m['PGxCategory']}\t"
            f"{m['FDA_Drug']}\t"
            f"{m['FDA_AffectedSubgroups']}\t"
            f"{m['FDA_Categories']}\t"
            f"{m['FDA_Recommendation']}\t"
            f"{m['MatchMethod']}"
        )

if __name__ == "__main__":
    main()

