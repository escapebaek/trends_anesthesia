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
import signal
import threading

# GitHub 설정 (사용자가 수정해야 할 부분)
GITHUB_REPO_PATH = "."  # 현재 디렉토리가 git 레포지토리라고 가정
GITHUB_REPO_URL = "https://github.com/escapebaek/trends_anesthesia.git"
AUTO_DEPLOY = True      # 자동 배포 여부
AUTO_OPEN_BROWSER = True  # 자동으로 브라우저 열기 여부

def safe_input(prompt, timeout=10, default='n'):
    """타임아웃이 있는 안전한 입력 함수"""
    def timeout_handler():
        print(f"\n⏰ {timeout}초 타임아웃 - 기본값 '{default}' 사용")
        return default
    
    try:
        print(prompt, end='', flush=True)
        timer = threading.Timer(timeout, timeout_handler)
        timer.start()
        
        try:
            result = input().strip().lower()
            timer.cancel()
            return result if result else default
        except (EOFError, KeyboardInterrupt):
            timer.cancel()
            print(f"\n⚠️ 입력 취소됨 - 기본값 '{default}' 사용")
            return default
        except Exception:
            timer.cancel()
            print(f"\n❌ 입력 오류 - 기본값 '{default}' 사용")
            return default
            
    except Exception:
        print(f"\n🔧 안전한 입력 모드 - 기본값 '{default}' 사용")
        return default

def setup_git_repo():
    """Git 레포지토리 초기 설정"""
    if not os.path.exists(".git"):
        print("📁 Git 레포지토리를 초기화합니다...")
        try:
            subprocess.run(["git", "init"], check=True)
            
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

# 환경변수 파일 (중요: API 키 보호)
.env
.env.local
.env.production
.env.staging

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
        except subprocess.CalledProcessError as e:
            print(f"❌ Git 초기화 실패: {e}")
            return False
    return True

def deploy_to_github():
    """GitHub Pages로 자동 배포"""
    try:
        result = subprocess.run(["git", "status", "--porcelain"], 
                              capture_output=True, text=True, check=True)
        
        if result.stdout.strip():
            print("📤 변경사항을 GitHub에 업로드합니다...")
            
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            commit_message = f"Update anesthesia classification dashboard - {timestamp}"
            
            subprocess.run(["git", "add", "."], check=True, timeout=30)
            subprocess.run(["git", "commit", "-m", commit_message], check=True, timeout=30)
            subprocess.run(["git", "push"], check=True, timeout=60)
            
            print("✅ GitHub에 업로드 완료!")
            
            try:
                result = subprocess.run(["git", "remote", "get-url", "origin"], 
                                      capture_output=True, text=True, check=True, timeout=10)
                remote_url = result.stdout.strip()
                
                if "github.com" in remote_url:
                    repo_path = remote_url.split("github.com/")[1].replace(".git", "")
                    username, repo_name = repo_path.split("/")
                    pages_url = f"https://{username}.github.io/{repo_name}/"
                    
                    print(f"🌐 GitHub Pages URL: {pages_url}")
                    print("⏳ 배포까지 5-10분 정도 소요될 수 있습니다.")
                    
                    if AUTO_OPEN_BROWSER:
                        print("🚀 자동으로 GitHub Pages를 브라우저에서 엽니다...")
                        try:
                            webbrowser.open(pages_url)
                            print("✅ 브라우저에서 열었습니다!")
                        except Exception as e:
                            print(f"⚠️ 브라우저 열기 실패: {e}")
                    else:
                        open_browser = safe_input(
                            "GitHub Pages를 브라우저에서 열까요? (y/n, 10초 후 자동으로 'n'): ", 
                            timeout=10, 
                            default='n'
                        )
                        
                        if open_browser == 'y':
                            try:
                                webbrowser.open(pages_url)
                                print("✅ 브라우저에서 열었습니다!")
                            except Exception as e:
                                print(f"⚠️ 브라우저 열기 실패: {e}")
                        else:
                            print("📝 수동으로 URL을 복사해서 브라우저에서 확인하세요.")
                    
                    return pages_url
                        
            except Exception as e:
                print(f"⚠️ GitHub Pages URL을 자동으로 확인할 수 없습니다: {e}")
                
        else:
            print("ℹ️ 변경사항이 없습니다.")
            
    except subprocess.TimeoutExpired:
        print("❌ Git 명령어 실행 시간 초과")
        return False
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

# 1. JSON 로드
json_path = "anesthesia_classified_abstracts.json"
if not os.path.exists(json_path):
    print(f"❌ {json_path} 파일을 찾을 수 없습니다.")
    print("💡 먼저 analyze_with_gemini.py를 실행하세요.")
    sys.exit(1)

print("📊 분류된 데이터 로드 중...")
with open(json_path, "r", encoding="utf-8") as f:
    classified_data = json.load(f)

# 메타데이터 로드 (있다면)
metadata = {}
meta_path = "anesthesia_classified_with_metadata.json"
if os.path.exists(meta_path):
    with open(meta_path, "r", encoding="utf-8") as f:
        full_data = json.load(f)
        metadata = full_data.get("metadata", {})

# 2. 데이터 전처리
category_stats = []
subtopic_stats = []
all_papers = []

for category, subtopics in classified_data.items():
    category_count = 0
    category_subtopics = 0
    
    for subtopic, papers in subtopics.items():
        if papers:  # 빈 리스트가 아닌 경우만
            category_count += len(papers)
            category_subtopics += 1
            
            # 세부주제 통계
            subtopic_stats.append({
                "category": category,
                "subtopic": subtopic,
                "count": len(papers),
                "category_short": category.split("(")[0].strip()
            })
            
            # 개별 논문 데이터
            for paper in papers:
                paper_data = paper.copy()
                paper_data["category"] = category
                paper_data["subtopic"] = subtopic
                paper_data["category_short"] = category.split("(")[0].strip()
                all_papers.append(paper_data)
    
    if category_count > 0:  # 논문이 있는 카테고리만
        category_stats.append({
            "category": category,
            "category_short": category.split("(")[0].strip(),
            "total_papers": category_count,
            "subtopics": category_subtopics
        })

# DataFrame 생성
df_categories = pd.DataFrame(category_stats)
df_subtopics = pd.DataFrame(subtopic_stats)
df_papers = pd.DataFrame(all_papers)

print(f"✅ 데이터 처리 완료:")
print(f"   - 활성 카테고리: {len(df_categories)}개")
print(f"   - 총 세부주제: {len(df_subtopics)}개")
print(f"   - 총 논문: {len(df_papers)}개")

# 3. 색상 팔레트 정의
colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8', '#F7DC6F', 
          '#BB8FCE', '#85C1E9', '#F8C471', '#82E0AA', '#F1948A', '#AED6F1']

category_colors = {}
if len(df_categories) > 0:
    for i, category in enumerate(df_categories['category'].unique()):
        category_colors[category] = colors[i % len(colors)]

print("📈 차트 생성 중...")

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
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
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
        
        .category-section {
            background: white;
            border-radius: 20px;
            margin: 30px 0;
            padding: 30px;
            box-shadow: 0 20px 40px rgba(0,0,0,0.1);
        }
        
        .category-header {
            background: linear-gradient(135deg, #667eea, #764ba2);
            color: white;
            padding: 20px 30px;
            border-radius: 15px;
            margin-bottom: 25px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            flex-wrap: wrap;
        }
        
        .category-title {
            font-size: 1.5em;
            font-weight: 600;
        }
        
        .category-stats {
            display: flex;
            gap: 15px;
            flex-wrap: wrap;
        }
        
        .category-stat {
            background: rgba(255,255,255,0.2);
            padding: 8px 16px;
            border-radius: 20px;
            font-size: 0.9em;
        }
        
        .subtopics-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(400px, 1fr));
            gap: 25px;
        }
        
        .subtopic-card {
            background: #f8f9fa;
            border-radius: 15px;
            padding: 25px;
            border-left: 5px solid #667eea;
            transition: all 0.3s ease;
            cursor: pointer;
        }
        
        .subtopic-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 15px 35px rgba(0,0,0,0.1);
            background: #fff;
        }
        
        .subtopic-title {
            font-size: 1.3em;
            font-weight: 600;
            color: #2c3e50;
            margin-bottom: 15px;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        
        .paper-count {
            background: #667eea;
            color: white;
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 0.8em;
            font-weight: 600;
        }
        
        .papers-list {
            max-height: 300px;
            overflow-y: auto;
            padding-right: 10px;
        }
        
        .paper-item {
            background: white;
            border-radius: 10px;
            padding: 15px;
            margin-bottom: 10px;
            border: 1px solid #e9ecef;
            transition: all 0.3s ease;
        }
        
        .paper-item:hover {
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
            border-color: #667eea;
        }
        
        .paper-title {
            font-weight: 600;
            color: #2c3e50;
            margin-bottom: 8px;
            font-size: 1.05em;
            line-height: 1.3;
        }
        
        .paper-details {
            display: grid;
            grid-template-columns: auto 1fr auto;
            gap: 10px;
            align-items: center;
            margin-bottom: 10px;
            font-size: 0.9em;
            color: #666;
        }
        
        .paper-journal {
            font-weight: 500;
            color: #495057;
        }
        
        .paper-date {
            background: #e3f2fd;
            color: #1976d2;
            padding: 2px 8px;
            border-radius: 10px;
            font-size: 0.8em;
        }
        
        .paper-summary {
            color: #666;
            font-size: 0.95em;
            line-height: 1.4;
            margin-bottom: 10px;
        }
        
        .paper-link {
            display: inline-block;
            background: #28a745;
            color: white;
            text-decoration: none;
            padding: 6px 12px;
            border-radius: 15px;
            font-size: 0.85em;
            transition: all 0.3s ease;
        }
        
        .paper-link:hover {
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
            .subtopics-grid {
                grid-template-columns: 1fr;
            }
            .category-header {
                flex-direction: column;
                align-items: flex-start;
                gap: 15px;
            }
        }
        
        /* 스크롤바 스타일링 */
        .papers-list::-webkit-scrollbar {
            width: 6px;
        }
        
        .papers-list::-webkit-scrollbar-track {
            background: #f1f1f1;
            border-radius: 3px;
        }
        
        .papers-list::-webkit-scrollbar-thumb {
            background: #c1c1c1;
            border-radius: 3px;
        }
        
        .papers-list::-webkit-scrollbar-thumb:hover {
            background: #a8a8a8;
        }
    </style>
    """

# HTML 구조 생성
doc.asis("<!DOCTYPE html>")
with tag("html", lang="en"):
    with tag("head"):
        doc.asis('<meta charset="UTF-8">')
        doc.asis('<meta name="viewport" content="width=device-width, initial-scale=1.0">')
        with tag("title"):
            text("Anesthesia Research Classification - Interactive Dashboard")
        doc.asis(create_modern_css())
        doc.asis('<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>')
    
    with tag("body"):
        with tag("div", klass="container"):
            # 헤더
            with tag("div", klass="header"):
                with tag("h1"):
                    text("🏥 마취학 연구 분류 대시보드")
                with tag("p"):
                    text("Anesthesia Research Classification Dashboard")
                with tag("p", style="font-size: 0.9em; margin-top: 10px; opacity: 0.7;"):
                    if metadata.get("analysis_date"):
                        text(f"Last updated: {metadata['analysis_date']}")
                    else:
                        text(f"Last updated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
            
            # 통계 카드
            with tag("div", klass="stats-grid"):
                with tag("div", klass="stat-card"):
                    with tag("div", klass="stat-number"):
                        text(str(len(df_categories)))
                    with tag("div", klass="stat-label"):
                        text("활성 카테고리")
                
                with tag("div", klass="stat-card"):
                    with tag("div", klass="stat-number"):
                        text(str(len(df_subtopics)))
                    with tag("div", klass="stat-label"):
                        text("세부 주제")
                
                with tag("div", klass="stat-card"):
                    with tag("div", klass="stat-number"):
                        text(str(len(df_papers)))
                    with tag("div", klass="stat-label"):
                        text("분류된 논문")
                
                with tag("div", klass="stat-card"):
                    with tag("div", klass="stat-number"):
                        text(str(metadata.get("total_papers_analyzed", len(df_papers))))
                    with tag("div", klass="stat-label"):
                        text("분석된 총 논문")

# 차트 생성 (데이터가 있는 경우에만)
if len(df_categories) > 0:
    # 카테고리별 논문 수 바 차트
    fig1 = px.bar(
        df_categories.sort_values('total_papers', ascending=True),
        x="total_papers",
        y="category_short",
        color="category_short",
        orientation="h",
        title="📊 카테고리별 논문 분포",
        labels={"total_papers": "논문 수", "category_short": "카테고리"},
        color_discrete_sequence=colors
    )
    fig1.update_layout(
        height=max(400, len(df_categories) * 40),
        font=dict(family="Arial, sans-serif", size=11),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        title_font_size=18,
        title_x=0.5,
        showlegend=False,
        margin=dict(l=200, r=50, t=80, b=60),
        yaxis=dict(tickfont=dict(size=10)),
        xaxis=dict(tickfont=dict(size=11))
    )

    # 도넛 차트
    fig2 = px.pie(
        df_categories,
        values='total_papers',
        names='category_short',
        title="🥧 카테고리별 비율",
        hole=0.4,
        color_discrete_sequence=colors
    )
    fig2.update_traces(
        textposition='inside',
        textinfo='percent+label',
        hovertemplate='<b>%{label}</b><br>논문 수: %{value}<br>비율: %{percent}<extra></extra>'
    )
    fig2.update_layout(
        font=dict(family="Arial, sans-serif", size=12),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        title_font_size=18,
        title_x=0.5
    )

    # 세부주제 상위 20개 차트
    top_subtopics = df_subtopics.sort_values('count', ascending=True).tail(20)
    fig3 = px.bar(
        top_subtopics,
        x='count',
        y='subtopic',
        color='category_short',
        orientation='h',
        title="🔍 상위 세부주제 (Top 20)",
        labels={"count": "논문 수", "subtopic": "세부주제"},
        color_discrete_sequence=colors
    )
    fig3.update_layout(
        height=800,
        font=dict(family="Arial, sans-serif", size=10),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        title_font_size=18,
        title_x=0.5,
        margin=dict(l=250, r=50, t=80, b=60),
        yaxis=dict(tickfont=dict(size=9)),
        legend=dict(title="카테고리")
    )

    # HTML에 차트 추가
    with tag("div", klass="dashboard-grid"):
        with tag("div", klass="chart-container"):
            doc.asis(fig1.to_html(full_html=False, include_plotlyjs=False, div_id="category-chart"))
        
        with tag("div", klass="chart-container"):
            doc.asis(fig2.to_html(full_html=False, include_plotlyjs=False, div_id="category-pie"))
        
        with tag("div", klass="chart-container full-width"):
            doc.asis(fig3.to_html(full_html=False, include_plotlyjs=False, div_id="subtopic-chart"))

# 카테고리별 상세 섹션
for _, cat_row in df_categories.sort_values('total_papers', ascending=False).iterrows():
    category = cat_row['category']
    category_subtopics = df_subtopics[df_subtopics['category'] == category].sort_values('count', ascending=False)
    
    with tag("div", klass="category-section"):
        with tag("div", klass="category-header"):
            with tag("div", klass="category-title"):
                text(f"📚 {category}")
            with tag("div", klass="category-stats"):
                with tag("div", klass="category-stat"):
                    text(f"{cat_row['total_papers']}편")
                with tag("div", klass="category-stat"):
                    text(f"{cat_row['subtopics']}개 주제")
        
        with tag("div", klass="subtopics-grid"):
            for _, subtopic_row in category_subtopics.iterrows():
                subtopic = subtopic_row['subtopic']
                papers = [p for p in all_papers if p['category'] == category and p['subtopic'] == subtopic]
                
                with tag("div", klass="subtopic-card"):
                    with tag("div", klass="subtopic-title"):
                        with tag("span"):
                            text(subtopic)
                        with tag("span", klass="paper-count"):
                            text(f"{len(papers)}편")
                    
                    with tag("div", klass="papers-list"):
                        for paper in papers:
                            with tag("div", klass="paper-item"):
                                with tag("div", klass="paper-title"):
                                    text(paper.get('title', 'No title'))
                                
                                with tag("div", klass="paper-details"):
                                    with tag("span", klass="paper-journal"):
                                        text(paper.get('journal', 'Unknown journal'))
                                    with tag("span"):
                                        text(f"by {paper.get('author', 'Unknown author')}")
                                    if paper.get('issue_date'):
                                        with tag("span", klass="paper-date"):
                                            text(paper['issue_date'])
                                
                                if paper.get('abstract_summary'):
                                    with tag("div", klass="paper-summary"):
                                        text(paper['abstract_summary'])
                                
                                if paper.get('link'):
                                    with tag("a", href=paper['link'], target="_blank", klass="paper-link"):
                                        text("📄 PubMed에서 보기")

# 푸터 추가
with tag("div", klass="footer"):
    with tag("p"):
        text("Generated with Python, Gemini AI & GitHub Pages")
    if metadata.get("date_range", {}).get("oldest_formatted"):
        with tag("p", style="font-size: 0.9em; margin-top: 5px;"):
            date_range = metadata["date_range"]
            text(f"논문 발행 기간: {date_range['oldest_formatted']} ~ {date_range['newest_formatted']}")

# JavaScript 추가
with tag("script"):
    doc.asis("""
    document.querySelectorAll('.subtopic-card').forEach(card => {
        card.addEventListener('mouseenter', function() {
            this.style.borderLeftWidth = '8px';
        });
        card.addEventListener('mouseleave', function() {
            this.style.borderLeftWidth = '5px';
        });
    });
    
    // 부드러운 등장 애니메이션
    window.addEventListener('load', function() {
        document.querySelectorAll('.chart-container, .subtopic-card, .stat-card, .category-section').forEach((el, index) => {
            el.style.opacity = '0';
            el.style.transform = 'translateY(20px)';
            el.style.transition = 'all 0.6s ease';
            
            setTimeout(() => {
                el.style.opacity = '1';
                el.style.transform = 'translateY(0)';
            }, index * 50);
        });
    });
    
    // 논문 아이템 호버 효과
    document.querySelectorAll('.paper-item').forEach(item => {
        item.addEventListener('mouseenter', function() {
            this.style.transform = 'translateX(5px)';
        });
        item.addEventListener('mouseleave', function() {
            this.style.transform = 'translateX(0)';
        });
    });
    """)

# HTML 저장
output_html = "index.html"  # GitHub Pages를 위해 index.html로 저장
print("💾 HTML 파일 생성 중...")
with open(output_html, "w", encoding="utf-8") as f:
    f.write(doc.getvalue())

print(f"✅ 마취학 분류 대시보드 생성 완료 → {output_html}")

# 자동 배포 실행
if AUTO_DEPLOY:
    print("\n🚀 GitHub Pages 자동 배포를 시작합니다...")
    
    # Git 레포지토리 확인/설정
    if setup_git_repo():
        # 배포 실행
        pages_url = deploy_to_github()
        if pages_url:
            print("🎉 배포가 완료되었습니다!")
        else:
            print("⚠️ 배포 중 문제가 발생했습니다. 로컬에서 확인합니다.")
            try:
                webbrowser.open("file://" + os.path.abspath(output_html))
            except Exception:
                print(f"📁 수동으로 파일을 열어주세요: {os.path.abspath(output_html)}")
    else:
        print("📝 Git 설정을 완료한 후 다시 실행해주세요.")
        try:
            webbrowser.open("file://" + os.path.abspath(output_html))
        except Exception:
            print(f"📁 수동으로 파일을 열어주세요: {os.path.abspath(output_html)}")
else:
    # 로컬에서만 열기
    try:
        webbrowser.open("file://" + os.path.abspath(output_html))
        print("🌐 로컬 브라우저에서 대시보드를 열었습니다.")
    except Exception:
        print(f"📁 수동으로 파일을 열어주세요: {os.path.abspath(output_html)}")
    print("💡 자동 배포를 원하시면 스크립트 상단의 AUTO_DEPLOY = True로 설정하세요.")

print("\n🏁 마취학 연구 분류 대시보드 생성이 완료되었습니다!")
print("📊 대시보드 주요 기능:")
print("   - 카테고리별 논문 분포 차트")
print("   - 세부주제별 상세 정보")
print("   - 각 논문의 요약 및 PubMed 링크")
print("   - 반응형 디자인으로 모바일 지원")