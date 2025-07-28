# fetch_pubmed.py (보안 개선 버전)
from Bio import Entrez, Medline
import pandas as pd
import os
from dotenv import load_dotenv

# .env 파일에서 환경변수 로드
load_dotenv()

# 환경변수에서 설정 가져오기
ENTREZ_EMAIL = os.getenv("ENTREZ_EMAIL")
ENTREZ_API_KEY = os.getenv("ENTREZ_API_KEY")

if not ENTREZ_EMAIL:
    raise ValueError("""
    ❌ ENTREZ_EMAIL이 설정되지 않았습니다.
    
    .env 파일에 다음과 같이 추가하세요:
    ENTREZ_EMAIL=your_email@example.com
    """)

if not ENTREZ_API_KEY:
    print("⚠️ ENTREZ_API_KEY가 설정되지 않았습니다.")
    print("💡 API 키가 없어도 동작하지만, 요청 제한이 있을 수 있습니다.")
    print("   NCBI API 키 발급: https://www.ncbi.nlm.nih.gov/account/settings/")

# Entrez 설정
Entrez.email = ENTREZ_EMAIL
if ENTREZ_API_KEY:
    Entrez.api_key = ENTREZ_API_KEY

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

# 검색 쿼리 생성
query = " OR ".join(f'"{j}"[Journal]' for j in journals) + " AND 2025[DP]"
print(f"🔍 검색 쿼리: {query}")

try:
    # 1) UID 검색
    print("📡 PubMed에서 논문 ID 검색 중...")
    handle = Entrez.esearch(db="pubmed", term=query, retmax=500, sort="pub date")
    search_results = Entrez.read(handle)
    uids = search_results["IdList"]
    handle.close()
    
    print(f"✅ {len(uids)}개의 논문 ID 발견")
    
    if not uids:
        print("❌ 검색 결과가 없습니다. 검색 조건을 확인하세요.")
        exit(1)

    # 2) 상세 정보 가져오기
    print("📄 논문 상세 정보 가져오는 중...")
    handle = Entrez.efetch(db="pubmed", id=uids, rettype="medline", retmode="text")
    records = list(Medline.parse(handle))
    handle.close()
    
    print(f"✅ {len(records)}개의 논문 정보 수집 완료")

except Exception as e:
    print(f"❌ PubMed API 호출 실패: {e}")
    print("💡 네트워크 연결이나 API 키를 확인하세요.")
    exit(1)

# 3) 데이터 처리
rows = []
for r in records:
    pmid = r.get("PMID", "")
    journal = r.get("JT", "")   # 저널명
    title = r.get("TI", "")
    authors = r.get("AU", [])
    abstract = r.get("AB", "")
    
    rows.append({
        "pmid": pmid,
        "journal": journal,
        "title": title,
        "authors": "; ".join(authors) if authors else "",
        "abstract": abstract,
        "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""
    })

# 4) DataFrame 생성 및 저장
df = pd.DataFrame(rows)

# 기본 통계 출력
print(f"\n📊 수집된 데이터 통계:")
print(f"   - 총 논문 수: {len(df)}")
print(f"   - Abstract가 있는 논문: {len(df[df['abstract'].notna() & (df['abstract'] != '')])}")
print(f"   - 저널별 분포:")
for journal, count in df['journal'].value_counts().head(10).items():
    print(f"     • {journal}: {count}개")

# CSV 저장
output_file = "pubmed_basic.csv"
df.to_csv(output_file, index=False, encoding='utf-8')
print(f"\n✅ 데이터 저장 완료 → {output_file}")
print(f"💡 다음 단계: python prepare_for_gemini.py")