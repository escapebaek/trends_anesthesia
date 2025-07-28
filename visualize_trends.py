import json
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.offline as pyo
from yattag import Doc
import webbrowser
import os
import subprocess
import sys
from collections import Counter
import re
from datetime import datetime

# GitHub 설정 (사용자가 수정해야 할 부분)
GITHUB_REPO_PATH = "."  # 현재 디렉토리가 git 레포지토리라고 가정
GITHUB_REPO_URL = "https://github.com/escapebaek/trends_anesthesia.git"
# 예: "https://github.com/username/repo-name"
AUTO_DEPLOY = True      # 자동 배포 여부

def setup_git_repo():
    """Git 레포지토리 초기 설정"""
    if not os.path.exists(".git"):
        print("📁 Git 레포지토리를 초기화합니다...")
        subprocess.run(["git", "init"], check=True)
        
        # .gitignore 생성
        with open(".gitignore", "w") as f:
            f.write("""
# Python
__pycache__/
*.pyc
*.pyo
*.pyd
.Python
*.so
.coverage
.pytest_cache/

# Data files (optional - 보안상 민감한 데이터는 제외)
# *.json

# OS
.DS_Store
Thumbs.db
""")
        
        print("✅ Git 레포지토리가 초기화되었습니다.")
        print("🔗 GitHub에서 레포지토리를 생성하고 다음 명령어를 실행하세요:")
        print("   git remote add origin https://github.com/escapebaek/trends_anesthesia.git")
        return False
    return True

def deploy_to_github():
    """GitHub Pages로 자동 배포"""
    try:
        # Git 상태 확인
        result = subprocess.run(["git", "status", "--porcelain"], 
                              capture_output=True, text=True, check=True)
        
        if result.stdout.strip():
            print("📤 변경사항을 GitHub에 업로드합니다...")
            
            # 현재 시간으로 커밋 메시지 생성
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            commit_message = f"Update dashboard - {timestamp}"
            
            # Git 명령어 실행
            subprocess.run(["git", "add", "."], check=True)
            subprocess.run(["git", "commit", "-m", commit_message], check=True)
            subprocess.run(["git", "push"], check=True)
            
            print("✅ GitHub에 업로드 완료!")
            
            # GitHub Pages URL 추정
            try:
                # 원격 URL 가져오기
                result = subprocess.run(["git", "remote", "get-url", "origin"], 
                                      capture_output=True, text=True, check=True)
                remote_url = result.stdout.strip()
                
                # GitHub Pages URL 생성
                if "github.com" in remote_url:
                    # https://github.com/user/repo.git -> user/repo
                    repo_path = remote_url.split("github.com/")[1].replace(".git", "")
                    username, repo_name = repo_path.split("/")
                    pages_url = f"https://{username}.github.io/{repo_name}/"
                    
                    print(f"🌐 GitHub Pages URL: {pages_url}")
                    print("⏳ 배포까지 5-10분 정도 소요될 수 있습니다.")
                    
                    # 브라우저에서 GitHub Pages 열기 (선택사항)
                    open_browser = input("GitHub Pages를 브라우저에서 열까요? (y/n): ")
                    if open_browser.lower() == 'y':
                        webbrowser.open(pages_url)
                        
            except Exception as e:
                print(f"⚠️ GitHub Pages URL을 자동으로 확인할 수 없습니다: {e}")
                
        else:
            print("ℹ️ 변경사항이 없습니다.")
            
    except subprocess.CalledProcessError as e:
        print(f"❌ Git 명령어 실행 실패: {e}")
        print("🔧 해결방법:")
        print("   1. Git이 설치되어 있는지 확인")
        print("   2. GitHub 레포지토리가 연결되어 있는지 확인")
        print("   3. 인증 정보가 올바른지 확인")
        return False
    except Exception as e:
        print(f"❌ 배포 중 오류 발생: {e}")
        return False
    
    return True

# 기존 시각화 코드 (이전과 동일)
# 1. JSON 로드
json_path = "anesthesia_trends_by_journal_with_article_links.json"
if not os.path.exists(json_path):
    print(f"❌ {json_path} 파일을 찾을 수 없습니다.")
    sys.exit(1)

with open(json_path, "r", encoding="utf-8") as f:
    data = json.load(f)

# 2. DataFrame 변환 및 데이터 전처리
rows = []
all_keywords = []
for journal, clusters in data.items():
    for cluster in clusters:
        keywords = cluster.get("related_keywords", [])
        all_keywords.extend(keywords)
        rows.append({
            "journal": journal,
            "topic": cluster["topic"],
            "count": cluster["count"],
            "description": cluster.get("description", ""),
            "keywords": ", ".join(keywords),
            "num_keywords": len(keywords),
            "links": cluster.get("article_links", []),
            "num_links": len(cluster.get("article_links", []))
        })

df = pd.DataFrame(rows)

# 키워드 빈도 분석
keyword_counts = Counter(all_keywords)
top_keywords = dict(keyword_counts.most_common(20))

# 저널별 통계
journal_stats = df.groupby('journal').agg({
    'count': ['sum', 'mean', 'max'],
    'topic': 'count'
}).round(2)
journal_stats.columns = ['total_articles', 'avg_per_topic', 'max_topic', 'num_topics']
journal_stats = journal_stats.reset_index()

# 3. 색상 팔레트 정의 (현대적인 색상)
colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8', '#F7DC6F', '#BB8FCE', '#85C1E9']
journal_colors = {journal: colors[i % len(colors)] for i, journal in enumerate(df['journal'].unique())}

# 4. HTML 문서 생성
doc, tag, text = Doc().tagtext()

def create_modern_css():
    return """
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            color: #333;
        }
        
        .container {
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }
        
        .header {
            text-align: center;
            color: white;
            margin-bottom: 40px;
            padding: 40px 0;
        }
        
        .header h1 {
            font-size: 3em;
            font-weight: 700;
            margin-bottom: 10px;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
        }
        
        .header p {
            font-size: 1.2em;
            opacity: 0.9;
        }
        
        .dashboard-grid {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 30px;
            margin-bottom: 40px;
        }
        
        .chart-container {
            background: white;
            border-radius: 20px;
            padding: 30px;
            box-shadow: 0 20px 40px rgba(0,0,0,0.1);
            backdrop-filter: blur(10px);
            border: 1px solid rgba(255,255,255,0.2);
        }
        
        .full-width {
            grid-column: 1 / -1;
        }
        
        .chart-title {
            font-size: 1.5em;
            font-weight: 600;
            margin-bottom: 20px;
            color: #2c3e50;
            text-align: center;
        }
        
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }
        
        .stat-card {
            background: rgba(255,255,255,0.95);
            border-radius: 15px;
            padding: 25px;
            text-align: center;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
            border: 1px solid rgba(255,255,255,0.3);
        }
        
        .stat-number {
            font-size: 2.5em;
            font-weight: 700;
            color: #667eea;
            margin-bottom: 10px;
        }
        
        .stat-label {
            color: #666;
            font-size: 1.1em;
            font-weight: 500;
        }
        
        .journal-section {
            background: white;
            border-radius: 20px;
            margin: 30px 0;
            padding: 30px;
            box-shadow: 0 20px 40px rgba(0,0,0,0.1);
        }
        
        .journal-header {
            background: linear-gradient(135deg, #667eea, #764ba2);
            color: white;
            padding: 20px 30px;
            border-radius: 15px;
            margin-bottom: 25px;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        
        .journal-title {
            font-size: 1.8em;
            font-weight: 600;
        }
        
        .journal-count {
            background: rgba(255,255,255,0.2);
            padding: 8px 16px;
            border-radius: 20px;
            font-size: 0.9em;
        }
        
        .topics-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(350px, 1fr));
            gap: 20px;
        }
        
        .topic-card {
            background: #f8f9fa;
            border-radius: 15px;
            padding: 20px;
            border-left: 5px solid #667eea;
            transition: all 0.3s ease;
            cursor: pointer;
        }
        
        .topic-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 15px 35px rgba(0,0,0,0.1);
            background: #fff;
        }
        
        .topic-title {
            font-size: 1.2em;
            font-weight: 600;
            color: #2c3e50;
            margin-bottom: 10px;
        }
        
        .topic-count {
            display: inline-block;
            background: #667eea;
            color: white;
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 0.8em;
            font-weight: 600;
            margin-bottom: 10px;
        }
        
        .topic-description {
            color: #666;
            font-size: 0.95em;
            margin-bottom: 15px;
            line-height: 1.4;
        }
        
        .topic-keywords {
            margin-bottom: 15px;
        }
        
        .keyword-tag {
            display: inline-block;
            background: #e3f2fd;
            color: #1976d2;
            padding: 3px 8px;
            border-radius: 12px;
            font-size: 0.8em;
            margin: 2px;
        }
        
        .topic-links a {
            display: inline-block;
            background: #28a745;
            color: white;
            text-decoration: none;
            padding: 6px 12px;
            border-radius: 15px;
            font-size: 0.8em;
            margin: 2px;
            transition: all 0.3s ease;
        }
        
        .topic-links a:hover {
            background: #218838;
            transform: scale(1.05);
        }
        
        .footer {
            text-align: center;
            color: white;
            margin-top: 50px;
            padding: 20px;
            opacity: 0.8;
        }
        
        @media (max-width: 768px) {
            .dashboard-grid {
                grid-template-columns: 1fr;
            }
            .header h1 {
                font-size: 2em;
            }
            .topics-grid {
                grid-template-columns: 1fr;
            }
        }
    </style>
    """

# HTML 구조 생성 (기존과 동일하지만 푸터 추가)
doc.asis("<!DOCTYPE html>")
with tag("html", lang="en"):
    with tag("head"):
        doc.asis('<meta charset="UTF-8">')
        doc.asis('<meta name="viewport" content="width=device-width, initial-scale=1.0">')
        with tag("title"):
            text("Anesthesia Research Trends - Interactive Dashboard")
        doc.asis(create_modern_css())
        doc.asis('<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>')
    
    with tag("body"):
        with tag("div", klass="container"):
            # 헤더
            with tag("div", klass="header"):
                with tag("h1"):
                    text("🏥 Anesthesia Research Trends")
                with tag("p"):
                    text("Interactive Analysis of Current Research Topics Across Major Journals")
                with tag("p", style="font-size: 0.9em; margin-top: 10px; opacity: 0.7;"):
                    text(f"Last updated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
            
            # 통계 카드
            with tag("div", klass="stats-grid"):
                with tag("div", klass="stat-card"):
                    with tag("div", klass="stat-number"):
                        text(str(len(df)))
                    with tag("div", klass="stat-label"):
                        text("Research Topics")
                
                with tag("div", klass="stat-card"):
                    with tag("div", klass="stat-number"):
                        text(str(df['count'].sum()))
                    with tag("div", klass="stat-label"):
                        text("Total Articles")
                
                with tag("div", klass="stat-card"):
                    with tag("div", klass="stat-number"):
                        text(str(len(df['journal'].unique())))
                    with tag("div", klass="stat-label"):
                        text("Journals")
                
                with tag("div", klass="stat-card"):
                    with tag("div", klass="stat-number"):
                        text(str(len(all_keywords)))
                    with tag("div", klass="stat-label"):
                        text("Keywords")

# 차트 생성 (기존과 동일)
# 메인 바 차트
fig1 = px.bar(
    df.sort_values('count', ascending=True).tail(20),
    x="count",
    y="topic",
    color="journal",
    color_discrete_map=journal_colors,
    orientation="h",
    hover_data={"description": True, "keywords": True},
    title="🔝 Top 20 Research Topics by Article Count",
    labels={"count": "Number of Articles", "topic": "Research Topic"}
)
fig1.update_layout(
    height=900,
    font=dict(family="Arial, sans-serif", size=11),
    plot_bgcolor='rgba(0,0,0,0)',
    paper_bgcolor='rgba(0,0,0,0)',
    title_font_size=18,
    title_x=0.5,
    title_y=0.98,
    showlegend=True,
    legend=dict(
        orientation="h", 
        yanchor="top", 
        y=-0.05, 
        xanchor="center", 
        x=0.5,
        bgcolor="rgba(255,255,255,0.8)",
        bordercolor="rgba(0,0,0,0.1)",
        borderwidth=1
    ),
    margin=dict(l=200, r=50, t=80, b=120),
    yaxis=dict(
        tickfont=dict(size=10),
        automargin=True
    ),
    xaxis=dict(
        tickfont=dict(size=11),
        title_font=dict(size=12)
    )
)

# 도넛 차트
journal_totals = df.groupby('journal')['count'].sum().reset_index()
fig2 = px.pie(
    journal_totals,
    values='count',
    names='journal',
    color='journal',
    color_discrete_map=journal_colors,
    title="📊 Article Distribution by Journal",
    hole=0.4
)
fig2.update_traces(
    textposition='inside',
    textinfo='percent+label',
    hovertemplate='<b>%{label}</b><br>Articles: %{value}<br>Percentage: %{percent}<extra></extra>'
)
fig2.update_layout(
    font=dict(family="Arial, sans-serif", size=12),
    plot_bgcolor='rgba(0,0,0,0)',
    paper_bgcolor='rgba(0,0,0,0)',
    title_font_size=20,
    title_x=0.5
)

# 키워드 차트
keyword_df = pd.DataFrame(list(top_keywords.items()), columns=['keyword', 'frequency'])
fig3 = px.bar(
    keyword_df.sort_values('frequency', ascending=True),
    x='frequency',
    y='keyword',
    orientation='h',
    title="🏷️ Most Frequent Keywords",
    color='frequency',
    color_continuous_scale='Viridis'
)
fig3.update_layout(
    height=600,
    font=dict(family="Arial, sans-serif", size=12),
    plot_bgcolor='rgba(0,0,0,0)',
    paper_bgcolor='rgba(0,0,0,0)',
    title_font_size=20,
    title_x=0.5,
    showlegend=False
)

# HTML에 차트 추가
with tag("div", klass="dashboard-grid"):
    with tag("div", klass="chart-container full-width"):
        doc.asis(fig1.to_html(full_html=False, include_plotlyjs=False, div_id="main-chart"))
    
    with tag("div", klass="chart-container"):
        doc.asis(fig2.to_html(full_html=False, include_plotlyjs=False, div_id="journal-pie"))
    
    with tag("div", klass="chart-container"):
        doc.asis(fig3.to_html(full_html=False, include_plotlyjs=False, div_id="keywords-chart"))

# 저널별 상세 섹션 (기존과 동일)
for journal in sorted(df["journal"].unique()):
    journal_data = df[df["journal"] == journal].sort_values('count', ascending=False)
    
    with tag("div", klass="journal-section"):
        with tag("div", klass="journal-header"):
            with tag("div", klass="journal-title"):
                text(f"📖 {journal}")
            with tag("div", klass="journal-count"):
                text(f"{len(journal_data)} topics • {journal_data['count'].sum()} articles")
        
        with tag("div", klass="topics-grid"):
            for _, row in journal_data.iterrows():
                with tag("div", klass="topic-card"):
                    with tag("div", klass="topic-title"):
                        text(row["topic"])
                    
                    with tag("div", klass="topic-count"):
                        text(f"{row['count']} articles")
                    
                    with tag("div", klass="topic-description"):
                        text(row["description"])
                    
                    if row["keywords"]:
                        with tag("div", klass="topic-keywords"):
                            for keyword in row["keywords"].split(", "):
                                with tag("span", klass="keyword-tag"):
                                    text(keyword)
                    
                    if row["links"]:
                        with tag("div", klass="topic-links"):
                            for i, url in enumerate(row["links"]):
                                with tag("a", href=url, target="_blank"):
                                    text(f"📄 Article {i+1}")

# 푸터 추가
with tag("div", klass="footer"):
    with tag("p"):
        text("Generated with Python, Plotly & GitHub Pages")
    with tag("p", style="font-size: 0.9em; margin-top: 5px;"):
        text("Auto-deployed research trends dashboard")

# JavaScript (기존과 동일)
with tag("script"):
    doc.asis("""
    document.querySelectorAll('.topic-card').forEach(card => {
        card.addEventListener('mouseenter', function() {
            this.style.borderLeftWidth = '8px';
        });
        card.addEventListener('mouseleave', function() {
            this.style.borderLeftWidth = '5px';
        });
    });
    
    window.addEventListener('load', function() {
        document.querySelectorAll('.chart-container, .topic-card, .stat-card').forEach((el, index) => {
            el.style.opacity = '0';
            el.style.transform = 'translateY(20px)';
            el.style.transition = 'all 0.6s ease';
            
            setTimeout(() => {
                el.style.opacity = '1';
                el.style.transform = 'translateY(0)';
            }, index * 100);
        });
    });
    """)

# HTML 저장
output_html = "index.html"  # GitHub Pages를 위해 index.html로 저장
with open(output_html, "w", encoding="utf-8") as f:
    f.write(doc.getvalue())

print(f"✅ 대시보드 생성 완료 → {output_html}")

# 자동 배포 실행
if AUTO_DEPLOY:
    print("\n🚀 GitHub Pages 자동 배포를 시작합니다...")
    
    # Git 레포지토리 확인/설정
    if setup_git_repo():
        # 배포 실행
        if deploy_to_github():
            print("🎉 배포가 완료되었습니다!")
        else:
            print("⚠️ 배포 중 문제가 발생했습니다.")
            # 로컬에서라도 열기
            webbrowser.open("file://" + os.path.abspath(output_html))
    else:
        print("📝 Git 설정을 완료한 후 다시 실행해주세요.")
        webbrowser.open("file://" + os.path.abspath(output_html))
else:
    # 로컬에서만 열기
    webbrowser.open("file://" + os.path.abspath(output_html))
    print("💡 자동 배포를 원하시면 스크립트 상단의 AUTO_DEPLOY = True로 설정하세요.")