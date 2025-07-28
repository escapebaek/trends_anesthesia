import json
import os
import re
import google.generativeai as genai

API_KEY = "AIzaSyB-TP4cwbW3a1qPezmDEal7z2Q-9NAQvMM"
genai.configure(api_key=API_KEY)

# 1. 데이터 로드
with open("abstracts_for_gemini.json", "r", encoding="utf-8") as f:
    abstracts = json.load(f)

# 2. 프롬프트 (세부 논문 링크 포함, 세분화된 키워드 그룹화)
prompt = (
    "You are an expert research assistant. Analyze the following anesthesia-related paper abstracts grouped by journal.\n\n"
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
    "7. Return strictly a JSON object where each key is a journal name, and its value is an array of topic clusters.\n"
    "8. Do not include markdown, code fences, or explanations.\n\n"
    "Example cluster entry: \n"
    "{ \"topic\": \"Adductor Canal Block\", \"count\": 7, \"related_keywords\": [\"ACB\", \"nerve block\", \"postoperative analgesia\"], \"article_links\": [\"https://pubmed.ncbi.nlm.nih.gov/12345678/\", \"https://pubmed.ncbi.nlm.nih.gov/23456789/\"], \"description\": \"Specific nerve block for knee surgery pain\" }\n\n"
    "Abstracts (with journal and article link):\n" +
    "\n\n".join(f"{i+1}. Journal: {item['journal']} | Link: {item['link']} | Abstract: {item['abstract']}" for i, item in enumerate(abstracts))
)

# 3. 모델 호출
model = genai.GenerativeModel("gemini-2.5-pro")
response = model.generate_content(prompt)

# 4. JSON 추출
raw_text = response.text.strip()
match = re.search(r"\{.*\}", raw_text, re.S)
json_str = match.group(0) if match else raw_text

# 5. JSON 파싱
try:
    trends_by_journal = json.loads(json_str)
except json.JSONDecodeError as e:
    print("❌ JSON 파싱 실패:", e)
    print("raw 출력:\n", raw_text[:500])
    trends_by_journal = {}

# 6. 저장
with open("anesthesia_trends_by_journal_with_article_links.json", "w", encoding="utf-8") as f:
    json.dump(trends_by_journal, f, ensure_ascii=False, indent=2)

print(f"✅ 저장 완료 → anesthesia_trends_by_journal_with_article_links.json (세부 키워드별 논문 링크 포함)")
