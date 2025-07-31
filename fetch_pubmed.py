# fetch_pubmed.py (보안 개선 버전 + 날짜 추출)
from Bio import Entrez, Medline
import pandas as pd
import os
from dotenv import load_dotenv
from datetime import datetime
import re

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

def parse_publication_date(dp_field):
    """
    DP 필드에서 발행 날짜를 파싱하여 datetime 객체로 변환
    DP 필드 형식: "2025 Jan", "2025 Jan 15", "2025" 등
    """
    if not dp_field:
        return None
    
    try:
        # 다양한 날짜 형식 처리
        dp_field = dp_field.strip()
        
        # "2025 Jan 15" 형식
        match = re.match(r'(\d{4})\s+([A-Za-z]{3,})\s+(\d{1,2})', dp_field)
        if match:
            year, month, day = match.groups()
            return datetime.strptime(f"{year} {month} {day}", "%Y %b %d")
        
        # "2025 Jan" 형식
        match = re.match(r'(\d{4})\s+([A-Za-z]{3,})', dp_field)
        if match:
            year, month = match.groups()
            return datetime.strptime(f"{year} {month}", "%Y %b")
        
        # "2025" 형식
        match = re.match(r'(\d{4})', dp_field)
        if match:
            year = match.group(1)
            return datetime.strptime(year, "%Y")
            
        return None
    except Exception as e:
        print(f"⚠️ 날짜 파싱 오류 ({dp_field}): {e}")
        return None

try:
    # 1) UID 검색
    print("📡 PubMed에서 논문 ID 검색 중...")
    handle = Entrez.esearch(db="pubmed", term=query, retmax=300, sort="pub date")
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

# 3) 데이터 처리 (날짜 정보 포함)
rows = []
publication_dates = []

for r in records:
    pmid = r.get("PMID", "")
    journal = r.get("JT", "")   # 저널명
    title = r.get("TI", "")
    authors = r.get("AU", [])
    abstract = r.get("AB", "")
    
    # 발행 날짜 추출 (DP 필드)
    dp_field = r.get("DP", "")
    pub_date = parse_publication_date(dp_field)
    
    # 날짜 통계를 위해 파싱된 날짜 저장
    if pub_date:
        publication_dates.append(pub_date)
    
    rows.append({
        "pmid": pmid,
        "journal": journal,
        "title": title,
        "authors": "; ".join(authors) if authors else "",
        "abstract": abstract,
        "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
        "publication_date_raw": dp_field,
        "publication_date": pub_date.strftime("%Y-%m-%d") if pub_date else "",
        "publication_year": pub_date.year if pub_date else "",
        "publication_month": pub_date.strftime("%Y-%m") if pub_date else ""
    })

# 4) DataFrame 생성 및 저장
df = pd.DataFrame(rows)

# 날짜 통계 계산
if publication_dates:
    newest_date = max(publication_dates)
    oldest_date = min(publication_dates)
    
    print(f"\n📅 발행 날짜 통계:")
    print(f"   - 최신 논문: {newest_date.strftime('%Y년 %m월 %d일')}")
    print(f"   - 가장 오래된 논문: {oldest_date.strftime('%Y년 %m월 %d일')}")
    print(f"   - 날짜 정보가 있는 논문: {len(publication_dates)}개")
    
    # 날짜 통계를 별도 파일로 저장
    date_stats = {
        "newest_date": newest_date.strftime("%Y-%m-%d"),
        "oldest_date": oldest_date.strftime("%Y-%m-%d"),
        "newest_date_formatted": newest_date.strftime("%Y년 %m월 %d일"),
        "oldest_date_formatted": oldest_date.strftime("%Y년 %m월 %d일"),
        "total_papers_with_dates": len(publication_dates),
        "date_range_days": (newest_date - oldest_date).days
    }
    
    import json
    with open("publication_date_stats.json", "w", encoding="utf-8") as f:
        json.dump(date_stats, f, ensure_ascii=False, indent=2)
    
    print(f"✅ 날짜 통계 저장 완료 → publication_date_stats.json")
else:
    print("⚠️ 발행 날짜 정보를 찾을 수 없습니다.")

# 기본 통계 출력
print(f"\n📊 수집된 데이터 통계:")
print(f"   - 총 논문 수: {len(df)}")
print(f"   - Abstract가 있는 논문: {len(df[df['abstract'].notna() & (df['abstract'] != '')])}")
print(f"   - 저널별 분포:")
for journal, count in df['journal'].value_counts().head(10).items():
    print(f"     • {journal}: {count}개")

# 월별 논문 분포 (상위 5개월)
if len(df[df['publication_month'] != '']) > 0:
    print(f"   - 월별 분포 (상위 5개월):")
    for month, count in df[df['publication_month'] != '']['publication_month'].value_counts().head(5).items():
        print(f"     • {month}: {count}개")

# CSV 저장
output_file = "pubmed_basic.csv"
df.to_csv(output_file, index=False, encoding='utf-8')
print(f"\n✅ 데이터 저장 완료 → {output_file}")
print(f"💡 다음 단계: python prepare_for_gemini.py")