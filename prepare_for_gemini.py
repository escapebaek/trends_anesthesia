# prepare_for_gemini.py
import pandas as pd

df = pd.read_csv("pubmed_basic.csv").dropna(subset=["abstract"])
# abstract 길이 필터 (예: 100자 미만 제외)
df = df[df["abstract"].str.len() > 100].reset_index(drop=True)
df.to_json("abstracts_for_gemini.json", orient="records", force_ascii=False)
print(f"{len(df)}개의 abstract 준비 완료 → abstracts_for_gemini.json")
