# fetch_pubmed.py (저널명 추가)
from Bio import Entrez, Medline
import pandas as pd

Entrez.email = "hanq719@gmail.com"
Entrez.api_key = "213cb12dcb378a3f33f98e78a4a59d3a1c09"

journals = [
    "British Journal of Anaesthesia",
    "Anesthesiology",
    "Anaesthesia",
    "European Journal of Anaesthesiology",
    "Pain",
    "Journal of Clinical Anesthesia",
    "Regional Anesthesia and Pain Medicine",
    "Anesthesia and Analgesia",
    "European Journal of Pain",
    "Pain Medicine"
]
query = " OR ".join(f'"{j}"[Journal]' for j in journals) + " AND 2025[DP]"

# 1) UID 검색
handle = Entrez.esearch(db="pubmed", term=query, retmax=500, sort="pub date")
uids = Entrez.read(handle)["IdList"]

# 2) 상세 정보
handle = Entrez.efetch(db="pubmed", id=uids, rettype="medline", retmode="text")
records = Medline.parse(handle)

rows = []
for r in records:
    pmid = r.get("PMID", "")
    journal = r.get("JT", "")   # ✅ 저널명 추가
    rows.append({
        "pmid": pmid,
        "journal": journal,     # ✅ CSV에 journal 필드 포함
        "title": r.get("TI", ""),
        "authors": "; ".join(r.get("AU", [])),
        "abstract": r.get("AB", ""),
        "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
    })

df = pd.DataFrame(rows)
df.to_csv("pubmed_basic.csv", index=False)
print("Saved pubmed_basic.csv with journal column")
