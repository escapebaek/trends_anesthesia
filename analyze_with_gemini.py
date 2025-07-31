import json
import os
import re
import google.generativeai as genai
from dotenv import load_dotenv
from datetime import datetime

# .env 파일에서 환경변수 로드
load_dotenv()

# 환경변수에서 API 키 가져오기
API_KEY = os.getenv("GEMINI_API_KEY")

if not API_KEY:
    raise ValueError("""
    ❌ GEMINI_API_KEY가 설정되지 않았습니다.
    
    다음 중 하나의 방법으로 설정하세요:
    
    1. .env 파일 생성 (권장):
       GEMINI_API_KEY=your_api_key_here
    
    2. 시스템 환경변수 설정:
       export GEMINI_API_KEY=your_api_key_here  # Linux/Mac
       set GEMINI_API_KEY=your_api_key_here     # Windows
    
    3. Python 실행 시 환경변수 설정:
       GEMINI_API_KEY=your_api_key python analyze_with_gemini.py
    """)

genai.configure(api_key=API_KEY)

# 1. 데이터 로드
try:
    with open("abstracts_for_gemini.json", "r", encoding="utf-8") as f:
        abstracts = json.load(f)
except FileNotFoundError:
    print("❌ abstracts_for_gemini.json 파일을 찾을 수 없습니다.")
    print("💡 먼저 fetch_pubmed.py와 prepare_for_gemini.py를 실행하세요.")
    exit(1)

# 날짜 정보 분석
papers_with_dates = [paper for paper in abstracts if paper.get("publication_date")]
if papers_with_dates:
    dates = [datetime.strptime(paper["publication_date"], "%Y-%m-%d") for paper in papers_with_dates]
    newest_date = max(dates)
    oldest_date = min(dates)
    
    print(f"📅 논문 발행 날짜 범위:")
    print(f"   - 최신: {newest_date.strftime('%Y년 %m월 %d일')}")
    print(f"   - 최오래된: {oldest_date.strftime('%Y년 %m월 %d일')}")
    print(f"   - 날짜 정보가 있는 논문: {len(papers_with_dates)}개")

# 2. 마취학 분류 카테고리 정의
anesthesia_categories = [
    "마취전 관리 (Pre-op Evaluation)",
    "마취 약리(Pharmacology of Anesthetics)",
    "법의학 및 윤리(Forensic and Ethical Considerations in Anesthesia)",
    "마취장비 및 감시(Anesthesia Equipment and Monitoring)",
    "기도관리(Airway Management)",
    "흡입마취(Inhalation Anesthesia)",
    "정맥마취(Intravenous Anesthesia)",
    "신경근차단(Neuromuscular Blockade)",
    "부위마취(Regional Anesthesia)",
    "수액 및 수혈(Fluid Management and Transfusion)",
    "산과마취(Obstetric Anesthesia)",
    "소아마취(Pediatric Anesthesia)",
    "심장마취(Cardiac Anesthesia)",
    "폐마취(Thoracic Anesthesia)",
    "뇌신경마취(Neuroanesthesia)",
    "수술장 밖 진정 및 마취(Sedation and Anesthesia Outside the Operating Room)",
    "수술 후 통증관리(Postoperative Pain Management)",
    "통증관리(Pain Management)",
    "노인마취(Geriatric Anesthesia)",
    "외래마취(Outpatient Anesthesia)",
    "심폐소생술(CPR)",
    "중환자관리(Critical Care Management)",
    "장기이식(Transplantation Anesthesia)"
]

# 3. 프롬프트 생성 (새로운 분류 방식)
date_context = ""
if papers_with_dates:
    date_context = f"""
    
    IMPORTANT DATE CONTEXT:
    - This analysis covers papers published from {oldest_date.strftime('%Y-%m-%d')} to {newest_date.strftime('%Y-%m-%d')}
    - Total papers with date information: {len(papers_with_dates)}
    - This represents the most current research trends in anesthesia for 2025
    """

categories_text = "\n".join(f"- {cat}" for cat in anesthesia_categories)

prompt = f"""You are an expert anesthesiologist and research analyst. Analyze the following anesthesia-related paper abstracts and classify them into specific categories.

{date_context}

CLASSIFICATION CATEGORIES:
{categories_text}

TASK:
1. For each abstract, determine the most appropriate category from the list above
2. Within each category, identify specific subtopics (e.g., "Kidney transplantation", "Liver transplantation" under "장기이식")
3. Provide a concise summary of each abstract (2-3 sentences)
4. Return the results in the following JSON structure:

{{
  "마취전 관리 (Pre-op Evaluation)": {{
    "Preoperative Risk Assessment": [
      {{
        "pmid": "12345678",
        "title": "Risk factors for postoperative complications",
        "author": "Kim HS, Lee JW",
        "journal": "Anesthesiology",
        "link": "https://pubmed.ncbi.nlm.nih.gov/12345678/",
        "issue_date": "2025-07",
        "abstract_summary": "This study investigated preoperative risk factors..."
      }}
    ]
  }},
  "마취 약리(Pharmacology of Anesthetics)": {{
    "Propofol Pharmacokinetics": [
      {{
        "pmid": "...",
        "title": "...",
        "author": "...",
        "journal": "...",
        "link": "...",
        "issue_date": "...",
        "abstract_summary": "..."
      }}
    ]
  }}
}}

INSTRUCTIONS:
- Create specific subtopic names based on the content (avoid generic terms)
- Each subtopic should contain an array of papers
- Abstract summaries should be concise but informative (2-3 sentences)
- Use the exact PMID, title, author, journal, and link from the provided data
- Format issue_date as "YYYY-MM" if available
- If a paper doesn't clearly fit any category, classify it as the closest match
- Do not include markdown, code fences, or explanations - return only valid JSON

ABSTRACTS TO ANALYZE:
"""

# 논문 데이터를 프롬프트에 추가
for i, item in enumerate(abstracts):
    prompt += f"\n{i+1}. PMID: {item.get('pmid', 'N/A')} | Journal: {item['journal']} | Date: {item.get('publication_date', 'Unknown')} | Link: {item['link']}\n"
    prompt += f"Title: {item['title']}\n"
    prompt += f"Authors: {item.get('authors', 'N/A')}\n"
    prompt += f"Abstract: {item['abstract']}\n"

# 4. 모델 호출
try:
    print("🤖 Gemini API 호출 중...")
    print(f"📊 분석할 논문 수: {len(abstracts)}개")
    model = genai.GenerativeModel("gemini-2.5-pro")
    response = model.generate_content(prompt)
    print("✅ API 호출 성공")
except Exception as e:
    print(f"❌ Gemini API 호출 실패: {e}")
    print("💡 API 키가 올바른지, 할당량이 남아있는지 확인하세요.")
    exit(1)

# 5. JSON 추출 및 파싱
raw_text = response.text.strip()
print("🔍 JSON 추출 중...")

# JSON 블록 찾기
json_match = re.search(r'\{.*\}', raw_text, re.DOTALL)
if json_match:
    json_str = json_match.group(0)
else:
    # 백틱이나 다른 마크다운 요소 제거
    json_str = re.sub(r'```json\s*', '', raw_text)
    json_str = re.sub(r'```\s*$', '', json_str)
    json_str = json_str.strip()

# 6. JSON 파싱
try:
    classified_data = json.loads(json_str)
    print(f"✅ JSON 파싱 성공")
    
    # 분류 결과 통계
    total_papers = 0
    category_counts = {}
    
    for category, subtopics in classified_data.items():
        category_count = 0
        for subtopic, papers in subtopics.items():
            category_count += len(papers)
        category_counts[category] = category_count
        total_papers += category_count
    
    print(f"📊 분류 결과:")
    print(f"   - 총 분류된 논문: {total_papers}개")
    print(f"   - 활성 카테고리: {len([c for c in category_counts.values() if c > 0])}개")
    
    # 상위 5개 카테고리 출력
    sorted_categories = sorted(category_counts.items(), key=lambda x: x[1], reverse=True)
    print(f"   - 상위 카테고리:")
    for cat, count in sorted_categories[:5]:
        if count > 0:
            print(f"     • {cat}: {count}개")

except json.JSONDecodeError as e:
    print("❌ JSON 파싱 실패:", e)
    print("Raw 출력 (처음 1000자):\n", raw_text[:1000])
    print("\n마지막 1000자:\n", raw_text[-1000:])
    
    # 빈 구조로 초기화
    classified_data = {}
    for category in anesthesia_categories:
        classified_data[category] = {}

# 7. 메타데이터와 함께 저장
output_data = {
    "metadata": {
        "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "total_papers_analyzed": len(abstracts),
        "total_papers_classified": sum(
            len(papers) for subtopics in classified_data.values() 
            for papers in subtopics.values()
        ) if classified_data else 0,
        "papers_with_dates": len(papers_with_dates),
        "date_range": {
            "oldest": oldest_date.strftime("%Y-%m-%d") if papers_with_dates else None,
            "newest": newest_date.strftime("%Y-%m-%d") if papers_with_dates else None,
            "oldest_formatted": oldest_date.strftime("%Y년 %m월 %d일") if papers_with_dates else None,
            "newest_formatted": newest_date.strftime("%Y년 %m월 %d일") if papers_with_dates else None
        },
        "categories_used": len([cat for cat, subtopics in classified_data.items() 
                              if any(len(papers) > 0 for papers in subtopics.values())]) if classified_data else 0,
        "category_distribution": {
            category: sum(len(papers) for papers in subtopics.values())
            for category, subtopics in classified_data.items()
        } if classified_data else {}
    },
    "classified_abstracts": classified_data
}

# 8. 파일 저장
# 기본 분류 결과
output_file = "anesthesia_classified_abstracts.json"
with open(output_file, "w", encoding="utf-8") as f:
    json.dump(classified_data, f, ensure_ascii=False, indent=2)

# 메타데이터 포함 버전
output_file_with_meta = "anesthesia_classified_with_metadata.json"
with open(output_file_with_meta, "w", encoding="utf-8") as f:
    json.dump(output_data, f, ensure_ascii=False, indent=2)

print(f"\n✅ 분류 결과 저장 완료:")
print(f"   → {output_file} (분류 결과만)")
print(f"   → {output_file_with_meta} (메타데이터 포함)")

if classified_data:
    print(f"📈 분류 통계:")
    print(f"   - 총 카테고리: {len(anesthesia_categories)}개")
    print(f"   - 사용된 카테고리: {len([cat for cat, subtopics in classified_data.items() if any(len(papers) > 0 for papers in subtopics.values())])}개")
    print(f"   - 총 세부주제: {sum(len(subtopics) for subtopics in classified_data.values())}개")

if papers_with_dates:
    print(f"📅 분석 기간: {oldest_date.strftime('%Y년 %m월 %d일')} ~ {newest_date.strftime('%Y년 %m월 %d일')}")

print(f"💡 다음 단계: python visualize_classified_trends.py")