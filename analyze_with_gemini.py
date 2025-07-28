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

# 2. 프롬프트 (날짜 정보 포함, 세부 논문 링크 포함, 세분화된 키워드 그룹화)
date_context = ""
if papers_with_dates:
    date_context = f"""
    
    IMPORTANT DATE CONTEXT:
    - This analysis covers papers published from {oldest_date.strftime('%Y-%m-%d')} to {newest_date.strftime('%Y-%m-%d')}
    - Total papers with date information: {len(papers_with_dates)}
    - This represents the most current research trends in anesthesia for 2025
    """

prompt = (
    "You are an expert research assistant. Analyze the following anesthesia-related paper abstracts grouped by journal.\n\n"
    f"{date_context}\n\n"
    "Task:\n"
    "1. For each journal, identify 12-15 research topic clusters with specific focus areas.\n"
    "2. Do NOT create overly broad clusters (e.g., avoid 'Regional Anesthesia'); instead, split them into detailed subtopics (e.g., 'Adductor Canal Block', 'Erector Spinae Plane Block').\n"
    "3. For each cluster, provide:\n"
    "   - topic (string): a specific representative name (avoid generic terms)\n"
    "   - count (integer): number of abstracts mentioning any keyword in this cluster\n"
    "   - related_keywords (array): 3-5 specific keywords contained in the cluster\n"
    "   - article_links (array): PubMed links of 3-5 representative articles related to this cluster\n"
    "   - description (string): concise summary (5-10 words)\n"
    "4. Use only the provided article links when filling the 'article_links' field.\n"
    "5. Avoid generic words like 'patients', 'surgery', 'pain', 'human'.\n"
    "6. Include drug names, specific techniques, biomarkers, and precise clinical outcomes where applicable.\n"
    "7. Consider the temporal context - these are recent 2025 publications representing current research trends.\n"
    "8. Return strictly a JSON object where each key is a journal name, and its value is an array of topic clusters.\n"
    "9. Do not include markdown, code fences, or explanations.\n\n"
    "Example cluster entry: \n"
    "{ \"topic\": \"Adductor Canal Block\", \"count\": 7, \"related_keywords\": [\"ACB\", \"nerve block\", \"postoperative analgesia\"], \"article_links\": [\"https://pubmed.ncbi.nlm.nih.gov/12345678/\", \"https://pubmed.ncbi.nlm.nih.gov/23456789/\"], \"description\": \"Specific nerve block for knee surgery pain\" }\n\n"
    "Abstracts (with journal, publication date, and article link):\n" +
    "\n\n".join(f"{i+1}. Journal: {item['journal']} | Date: {item.get('publication_date', 'Unknown')} | Link: {item['link']} | Abstract: {item['abstract']}" for i, item in enumerate(abstracts))
)

# 3. 모델 호출
try:
    print("🤖 Gemini API 호출 중...")
    model = genai.GenerativeModel("gemini-2.5-pro")
    response = model.generate_content(prompt)
    print("✅ API 호출 성공")
except Exception as e:
    print(f"❌ Gemini API 호출 실패: {e}")
    print("💡 API 키가 올바른지, 할당량이 남아있는지 확인하세요.")
    exit(1)

# 4. JSON 추출
raw_text = response.text.strip()
match = re.search(r"\{.*\}", raw_text, re.S)
json_str = match.group(0) if match else raw_text

# 5. JSON 파싱
try:
    trends_by_journal = json.loads(json_str)
    print(f"✅ JSON 파싱 성공 - {len(trends_by_journal)}개 저널 데이터")
except json.JSONDecodeError as e:
    print("❌ JSON 파싱 실패:", e)
    print("Raw 출력 (처음 500자):\n", raw_text[:500])
    trends_by_journal = {}

# 6. 저장 (날짜 메타데이터 포함)
output_data = {
    "metadata": {
        "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "total_papers": len(abstracts),
        "papers_with_dates": len(papers_with_dates),
        "date_range": {
            "oldest": oldest_date.strftime("%Y-%m-%d") if papers_with_dates else None,
            "newest": newest_date.strftime("%Y-%m-%d") if papers_with_dates else None,
            "oldest_formatted": oldest_date.strftime("%Y년 %m월 %d일") if papers_with_dates else None,
            "newest_formatted": newest_date.strftime("%Y년 %m월 %d일") if papers_with_dates else None
        }
    },
    "trends_by_journal": trends_by_journal
}

# 기존 형식으로도 저장 (하위 호환성)
output_file = "anesthesia_trends_by_journal_with_article_links.json"
with open(output_file, "w", encoding="utf-8") as f:
    json.dump(trends_by_journal, f, ensure_ascii=False, indent=2)

# 메타데이터 포함 버전 저장
output_file_with_meta = "anesthesia_trends_with_metadata.json"
with open(output_file_with_meta, "w", encoding="utf-8") as f:
    json.dump(output_data, f, ensure_ascii=False, indent=2)

print(f"✅ 저장 완료:")
print(f"   → {output_file} (기존 형식)")
print(f"   → {output_file_with_meta} (메타데이터 포함)")
print(f"📊 분석된 총 토픽 수: {sum(len(topics) for topics in trends_by_journal.values())}")

if papers_with_dates:
    print(f"📅 분석 기간: {oldest_date.strftime('%Y년 %m월 %d일')} ~ {newest_date.strftime('%Y년 %m월 %d일')}")
print(f"💡 다음 단계: python visualize_trends.py")