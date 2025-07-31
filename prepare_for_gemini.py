# prepare_for_gemini.py (날짜 정보 포함)
import pandas as pd
import json
from datetime import datetime

# CSV 파일 로드
df = pd.read_csv("pubmed_basic.csv").dropna(subset=["abstract"])

# abstract 길이 필터 (예: 100자 미만 제외)
df = df[df["abstract"].str.len() > 100].reset_index(drop=True)

print(f"📊 데이터 전처리 완료:")
print(f"   - 총 논문 수: {len(df)}")
print(f"   - Abstract 평균 길이: {df['abstract'].str.len().mean():.0f}자")

# 날짜 정보가 있는 논문 확인
date_info_count = len(df[df['publication_date'] != ''])
print(f"   - 날짜 정보가 있는 논문: {date_info_count}개 ({date_info_count/len(df)*100:.1f}%)")

if date_info_count > 0:
    # 최신/최오래된 논문 정보 출력
    df_with_dates = df[df['publication_date'] != ''].copy()
    df_with_dates['pub_date_parsed'] = pd.to_datetime(df_with_dates['publication_date'])
    
    newest_paper = df_with_dates.loc[df_with_dates['pub_date_parsed'].idxmax()]
    oldest_paper = df_with_dates.loc[df_with_dates['pub_date_parsed'].idxmin()]
    
    print(f"   - 최신 논문: {newest_paper['publication_date']} ({newest_paper['journal']})")
    print(f"   - 가장 오래된 논문: {oldest_paper['publication_date']} ({oldest_paper['journal']})")

# JSON으로 저장 (날짜 정보 포함)
output_data = []
for _, row in df.iterrows():
    paper_data = {
        "pmid": row["pmid"],
        "journal": row["journal"],
        "title": row["title"],
        "authors": row.get("authors", ""),   # ✅ 저자 필드 추가
        "abstract": row["abstract"],
        "link": row["link"],
        "publication_date": row.get("publication_date", ""),
        "publication_date_raw": row.get("publication_date_raw", ""),
        "publication_year": row.get("publication_year", ""),
        "publication_month": row.get("publication_month", "")
    }
    output_data.append(paper_data)

# JSON 파일로 저장
with open("abstracts_for_gemini.json", "w", encoding="utf-8") as f:
    json.dump(output_data, f, ensure_ascii=False, indent=2)

print(f"✅ {len(output_data)}개의 abstract 준비 완료 (저자 포함) → abstracts_for_gemini.json")
print("💡 다음 단계: python analyze_with_gemini.py")