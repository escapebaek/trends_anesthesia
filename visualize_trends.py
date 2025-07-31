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

# 전체 마취학 대주제 구조 정의
ANESTHESIA_CATEGORIES = [
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

# 2. 전체 마취학 구조를 기반으로 데이터 전처리
def create_complete_category_structure(data):
    """전체 마취학 카테고리 구조 생성"""
    complete_structure = []
    
    for category in ANESTHESIA_CATEGORIES:
        category_short = category.split("(")[0].strip()
        
        if category in data and data[category]:
            # 데이터가 있는 카테고리
            subtopic_count = 0
            paper_count = 0
            subtopics = []
            
            for subtopic, papers in data[category].items():
                if papers:  # 빈 리스트가 아닌 경우만
                    subtopic_count += 1
                    paper_count += len(papers)
                    subtopics.append({
                        'name': subtopic,
                        'count': len(papers)
                    })
            
            complete_structure.append({
                'category': category,
                'category_short': category_short,
                'status': 'active',
                'paper_count': paper_count,
                'subtopic_count': subtopic_count,
                'subtopics': subtopics
            })
        else:
            # 데이터가 없는 카테고리
            complete_structure.append({
                'category': category,
                'category_short': category_short,
                'status': 'empty',
                'paper_count': 0,
                'subtopic_count': 0,
                'subtopics': []
            })
    
    return complete_structure

category_structure = create_complete_category_structure(classified_data)

# 활성 카테고리와 전체 카테고리 분리
active_categories = [cat for cat in category_structure if cat['status'] == 'active']
empty_categories = [cat for cat in category_structure if cat['status'] == 'empty']

# DataFrame 생성
df_all_categories = pd.DataFrame(category_structure)
df_active = pd.DataFrame(active_categories)

# 개별 논문 데이터 처리
all_papers = []
subtopic_stats = []

for category, subtopics in classified_data.items():
    for subtopic, papers in subtopics.items():
        if papers:  # 빈 리스트가 아닌 경우만
            subtopic_stats.append({
                "category": category,
                "subtopic": subtopic,
                "count": len(papers),
                "category_short": category.split("(")[0].strip()
            })
            
            for paper in papers:
                paper_data = paper.copy()
                paper_data["category"] = category
                paper_data["subtopic"] = subtopic
                paper_data["category_short"] = category.split("(")[0].strip()
                all_papers.append(paper_data)

df_subtopics = pd.DataFrame(subtopic_stats)
df_papers = pd.DataFrame(all_papers)

print(f"✅ 데이터 처리 완료:")
print(f"   - 전체 마취학 카테고리: {len(ANESTHESIA_CATEGORIES)}개")
print(f"   - 활성 카테고리: {len(active_categories)}개")
print(f"   - 빈 카테고리: {len(empty_categories)}개")
print(f"   - 총 세부주제: {len(df_subtopics)}개")
print(f"   - 총 논문: {len(df_papers)}개")

# 3. 향상된 색상 팔레트 정의
def create_color_palette():
    """카테고리별 색상 팔레트 생성"""
    colors = [
        '#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8', '#F7DC6F', 
        '#BB8FCE', '#85C1E9', '#F8C471', '#82E0AA', '#F1948A', '#AED6F1',
        '#D7BDE2', '#A9DFBF', '#F9E79F', '#F8D7DA', '#D4EDDA', '#D1ECF1',
        '#E2E3E5', '#F8F9FA', '#FFF3CD', '#FCF8E3', '#E7F3FF'
    ]
    
    category_colors = {}
    for i, category in enumerate(ANESTHESIA_CATEGORIES):
        if i < len(colors):
            category_colors[category] = colors[i]
        else:
            # 추가 색상이 필요한 경우 색상 재사용
            category_colors[category] = colors[i % len(colors)]
            
    return category_colors

category_colors = create_color_palette()

print("📈 향상된 차트 생성 중...")

# 4. HTML 문서 생성
doc, tag, text = Doc().tagtext()

def create_enhanced_css():
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
            max-width: 1600px;
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
            font-size: 3.5em;
            font-weight: 700;
            margin-bottom: 15px;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
        }
        
        .header p {
            font-size: 1.3em;
            opacity: 0.9;
            margin-bottom: 10px;
        }
        
        .overview-stats {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 25px;
            margin: 40px 0;
        }
        
        .overview-card {
            background: rgba(255,255,255,0.95);
            border-radius: 20px;
            padding: 30px;
            text-align: center;
            box-shadow: 0 15px 35px rgba(0,0,0,0.15);
            border: 1px solid rgba(255,255,255,0.3);
            transition: transform 0.3s ease;
        }
        
        .overview-card:hover {
            transform: translateY(-5px);
        }
        
        .overview-number {
            font-size: 3em;
            font-weight: 700;
            margin-bottom: 10px;
        }
        
        .overview-number.active { color: #28a745; }
        .overview-number.empty { color: #6c757d; }
        .overview-number.subtopics { color: #17a2b8; }
        .overview-number.papers { color: #dc3545; }
        
        .overview-label {
            color: #666;
            font-size: 1.2em;
            font-weight: 500;
        }
        
        .charts-section {
            background: white;
            border-radius: 25px;
            margin: 40px 0;
            padding: 40px;
            box-shadow: 0 20px 40px rgba(0,0,0,0.1);
        }
        
        .charts-title {
            text-align: center;
            font-size: 2.2em;
            font-weight: 600;
            color: #2c3e50;
            margin-bottom: 40px;
        }
        
        .chart-grid {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 40px;
            margin-bottom: 40px;
        }
        
        .chart-container {
            background: #f8f9fa;
            border-radius: 20px;
            padding: 30px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.05);
        }
        
        .full-width-chart {
            grid-column: 1 / -1;
            background: #f8f9fa;
            border-radius: 20px;
            padding: 30px;
            margin-top: 30px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.05);
        }
        
        .chart-title {
            font-size: 1.6em;
            font-weight: 600;
            margin-bottom: 25px;
            color: #2c3e50;
            text-align: center;
        }
        
        .category-overview {
            background: white;
            border-radius: 25px;
            margin: 40px 0;
            padding: 40px;
            box-shadow: 0 20px 40px rgba(0,0,0,0.1);
        }
        
        .category-overview-title {
            text-align: center;
            font-size: 2.2em;
            font-weight: 600;
            color: #2c3e50;
            margin-bottom: 40px;
        }
        
        .category-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(320px, 1fr));
            gap: 25px;
        }
        
        .category-item {
            border-radius: 15px;
            padding: 25px;
            transition: all 0.3s ease;
            cursor: pointer;
            border: 2px solid transparent;
        }
        
        .category-item.active {
            background: linear-gradient(135deg, #e8f5e8, #f0f8f0);
            border-color: #28a745;
        }
        
        .category-item.empty {
            background: linear-gradient(135deg, #f8f9fa, #e9ecef);
            border-color: #dee2e6;
            opacity: 0.7;
        }
        
        .category-item:hover {
            transform: translateY(-3px);
            box-shadow: 0 10px 25px rgba(0,0,0,0.1);
        }
        
        .category-name {
            font-size: 1.3em;
            font-weight: 600;
            margin-bottom: 15px;
            color: #2c3e50;
        }
        
        .category-stats-inline {
            display: flex;
            gap: 15px;
            flex-wrap: wrap;
        }
        
        .category-stat-badge {
            padding: 6px 12px;
            border-radius: 20px;
            font-size: 0.9em;
            font-weight: 500;
        }
        
        .active .category-stat-badge {
            background: #28a745;
            color: white;
        }
        
        .empty .category-stat-badge {
            background: #6c757d;
            color: white;
        }
        
        .detailed-categories {
            margin-top: 50px;
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
            padding: 25px 35px;
            border-radius: 15px;
            margin-bottom: 30px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            flex-wrap: wrap;
        }
        
        .category-title-main {
            font-size: 1.6em;
            font-weight: 600;
        }
        
        .category-stats-header {
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
        }
        
        .category-stat-header {
            background: rgba(255,255,255,0.2);
            padding: 10px 18px;
            border-radius: 25px;
            font-size: 1em;
            font-weight: 500;
        }
        
        .subtopics-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(450px, 1fr));
            gap: 30px;
        }
        
        .subtopic-card {
            background: #f8f9fa;
            border-radius: 15px;
            padding: 25px;
            border-left: 6px solid #667eea;
            transition: all 0.3s ease;
            cursor: pointer;
        }
        
        .subtopic-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 15px 35px rgba(0,0,0,0.1);
            background: #fff;
        }
        
        .subtopic-title {
            font-size: 1.4em;
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
            padding: 6px 14px;
            border-radius: 25px;
            font-size: 0.85em;
            font-weight: 600;
        }
        
        .papers-list {
            max-height: 350px;
            overflow-y: auto;
            padding-right: 10px;
        }
        
        .paper-item {
            background: white;
            border-radius: 12px;
            padding: 18px;
            margin-bottom: 12px;
            border: 1px solid #e9ecef;
            transition: all 0.3s ease;
        }
        
        .paper-item:hover {
            box-shadow: 0 8px 20px rgba(0,0,0,0.1);
            border-color: #667eea;
            transform: translateX(3px);
        }
        
        .paper-title {
            font-weight: 600;
            color: #2c3e50;
            margin-bottom: 10px;
            font-size: 1.1em;
            line-height: 1.4;
        }
        
        .paper-details {
            display: grid;
            grid-template-columns: auto 1fr auto;
            gap: 12px;
            align-items: center;
            margin-bottom: 12px;
            font-size: 0.95em;
            color: #666;
        }
        
        .paper-journal {
            font-weight: 500;
            color: #495057;
        }
        
        .paper-date {
            background: #e3f2fd;
            color: #1976d2;
            padding: 4px 10px;
            border-radius: 12px;
            font-size: 0.85em;
        }
        
        .paper-summary {
            color: #666;
            font-size: 1em;
            line-height: 1.5;
            margin-bottom: 12px;
        }
        
        .paper-link {
            display: inline-block;
            background: #28a745;
            color: white;
            text-decoration: none;
            padding: 8px 16px;
            border-radius: 20px;
            font-size: 0.9em;
            transition: all 0.3s ease;
        }
        
        .paper-link:hover {
            background: #218838;
            transform: scale(1.05);
        }
        
        .footer {
            text-align: center;
            color: white;
            margin-top: 60px;
            padding: 30px;
            opacity: 0.9;
        }
        
        @media (max-width: 1024px) {
            .chart-grid {
                grid-template-columns: 1fr;
            }
        }
        
        @media (max-width: 768px) {
            .header h1 {
                font-size: 2.5em;
            }
            .subtopics-grid {
                grid-template-columns: 1fr;
            }
            .category-header {
                flex-direction: column;
                align-items: flex-start;
                gap: 20px;
            }
            .category-grid {
                grid-template-columns: 1fr;
            }
        }
        
        /* 스크롤바 스타일링 */
        .papers-list::-webkit-scrollbar {
            width: 8px;
        }
        
        .papers-list::-webkit-scrollbar-track {
            background: #f1f1f1;
            border-radius: 4px;
        }
        
        .papers-list::-webkit-scrollbar-thumb {
            background: #c1c1c1;
            border-radius: 4px;
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
            text("마취학 연구 분류 대시보드 - Anesthesia Research Classification")
        doc.asis(create_enhanced_css())
        doc.asis('<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>')
    
    with tag("body"):
        with tag("div", klass="container"):
            # 향상된 헤더
            with tag("div", klass="header"):
                with tag("h1"):
                    text("🏥 마취학 연구 분류 대시보드")
                with tag("p"):
                    text("Comprehensive Anesthesia Research Classification Dashboard")
                with tag("p", style="font-size: 1em; margin-top: 15px; opacity: 0.8;"):
                    text("23개 주요 마취학 분야별 연구 동향 및 분석")
                with tag("p", style="font-size: 0.9em; margin-top: 10px; opacity: 0.7;"):
                    if metadata.get("analysis_date"):
                        text(f"Last updated: {metadata['analysis_date']}")
                    else:
                        text(f"Last updated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
            
            # 개요 통계 카드
            with tag("div", klass="overview-stats"):
                with tag("div", klass="overview-card"):
                    with tag("div", klass="overview-number papers"):
                        text(str(len(df_papers)))
                    with tag("div", klass="overview-label"):
                        text("분석된 논문 수")

# 차트 섹션 (데이터가 있는 경우에만)
if len(active_categories) > 0:
    with tag("div", klass="charts-section"):
        with tag("h2", klass="charts-title"):
            text("📊 연구 현황 분석")
        
        with tag("div", klass="chart-grid"):
            # 전체 마취학 분야 현황 (활성/비활성)
            with tag("div", klass="chart-container"):
                with tag("div", klass="chart-title"):
                    text("🎯 마취학 분야별 연구 현황")
                with tag("div", id="field-status-chart"):
                    pass  # 차트가 여기에 삽입됨
            
            # 활성 분야의 논문 분포
            with tag("div", klass="chart-container"):
                with tag("div", klass="chart-title"):
                    text("📈 활성 분야별 논문 분포")
                with tag("div", id="active-papers-chart"):
                    pass  # 차트가 여기에 삽입됨
        
        # 세부 주제 분석
        with tag("div", klass="full-width-chart"):
            with tag("div", klass="chart-title"):
                text("🔍 세부 연구 주제 상위 20개")
            with tag("div", id="subtopics-chart"):
                pass  # 차트가 여기에 삽입됨

# 전체 카테고리 개요
with tag("div", klass="category-overview"):
    with tag("h2", klass="category-overview-title"):
        text("🗂️ 마취학 23개 분야 전체 개요")
    
    with tag("div", klass="category-grid"):
        for cat_info in category_structure:
            css_class = "category-item active" if cat_info['status'] == 'active' else "category-item empty"
            with tag("div", klass=css_class):
                with tag("div", klass="category-name"):
                    text(cat_info['category'])
                
                with tag("div", klass="category-stats-inline"):
                    with tag("span", klass="category-stat-badge"):
                        if cat_info['status'] == 'active':
                            text(f"📄 {cat_info['paper_count']}편")
                        else:
                            text("📄 0편")
                    
                    with tag("span", klass="category-stat-badge"):
                        if cat_info['status'] == 'active':
                            text(f"🏷️ {cat_info['subtopic_count']}주제")
                        else:
                            text("🏷️ 0주제")
                    
                    with tag("span", klass="category-stat-badge"):
                        status_text = "✅ 활성" if cat_info['status'] == 'active' else "⏸️ 대기"
                        text(status_text)

# 활성 카테고리 상세 섹션
if len(active_categories) > 0:
    with tag("div", klass="detailed-categories"):
        for cat_info in active_categories:
            category = cat_info['category']
            
            with tag("div", klass="category-section"):
                with tag("div", klass="category-header"):
                    with tag("div", klass="category-title-main"):
                        text(f"📚 {category}")
                    with tag("div", klass="category-stats-header"):
                        with tag("div", klass="category-stat-header"):
                            text(f"📄 {cat_info['paper_count']}편")
                        with tag("div", klass="category-stat-header"):
                            text(f"🏷️ {cat_info['subtopic_count']}개 주제")
                
                with tag("div", klass="subtopics-grid"):
                    for subtopic_info in cat_info['subtopics']:
                        subtopic = subtopic_info['name'] 
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

# 푸터
with tag("div", klass="footer"):
    with tag("p", style="font-size: 1.1em; margin-bottom: 10px;"):
        text("🔬 Generated with Python, Gemini AI & GitHub Pages")
    with tag("p", style="font-size: 1em; margin-bottom: 15px;"):
        text("마취학 연구 동향 분석 및 시각화 대시보드")
    if metadata.get("date_range", {}).get("oldest_formatted"):
        date_range = metadata["date_range"]
        with tag("p", style="font-size: 0.95em; opacity: 0.8;"):
            text(f"📅 논문 발행 기간: {date_range['oldest_formatted']} ~ {date_range['newest_formatted']}")

# JavaScript 및 차트 생성
with tag("script"):
    doc.asis(f"""
    // 차트 데이터 준비
    const categoryData = {df_all_categories.to_json(orient='records')};
    const activeData = {df_active.to_json(orient='records') if len(df_active) > 0 else '[]'};
    const subtopicData = {df_subtopics.to_json(orient='records') if len(df_subtopics) > 0 else '[]'};
    
    // 1. 전체 분야 현황 도넛 차트
    const fieldStatusData = [
        {{
            values: [{len(active_categories)}, {len(empty_categories)}],
            labels: ['활성 연구 분야', '미개발 연구 분야'],
            type: 'pie',
            hole: 0.5,
            marker: {{
                colors: ['#28a745', '#6c757d']
            }},
            textinfo: 'label+percent+value',
            textposition: 'outside',
            hovertemplate: '<b>%{{label}}</b><br>개수: %{{value}}<br>비율: %{{percent}}<extra></extra>'
        }}
    ];
    
    const fieldStatusLayout = {{
        font: {{family: "Arial, sans-serif", size: 12}},
        plot_bgcolor: 'rgba(0,0,0,0)',
        paper_bgcolor: 'rgba(0,0,0,0)',
        margin: {{l: 50, r: 50, t: 20, b: 20}},
        showlegend: true,
        legend: {{
            orientation: 'h',
            x: 0.5,
            xanchor: 'center',
            y: -0.1
        }}
    }};
    
    if (document.getElementById('field-status-chart')) {{
        Plotly.newPlot('field-status-chart', fieldStatusData, fieldStatusLayout, {{responsive: true}});
    }}
    
    // 2. 활성 분야별 논문 분포 (세로 막대 차트)
    if (activeData.length > 0) {{
        const activeSorted = activeData.sort((a, b) => b.paper_count - a.paper_count);
        
        const activePapersData = [{{
            x: activeSorted.map(d => d.category_short),
            y: activeSorted.map(d => d.paper_count),
            type: 'bar',
            marker: {{
                color: activeSorted.map((d, i) => ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8'][i % 5])
            }},
            text: activeSorted.map(d => d.paper_count),
            textposition: 'outside',
            hovertemplate: '<b>%{{x}}</b><br>논문 수: %{{y}}<br>세부주제: %{{customdata}}개<extra></extra>',
            customdata: activeSorted.map(d => d.subtopic_count)
        }}];
        
        const activePapersLayout = {{
            font: {{family: "Arial, sans-serif", size: 11}},
            plot_bgcolor: 'rgba(0,0,0,0)',
            paper_bgcolor: 'rgba(0,0,0,0)',
            margin: {{l: 60, r: 50, t: 20, b: 120}},
            xaxis: {{
                tickangle: -45,
                tickfont: {{size: 10}},
                title: {{text: '연구 분야', font: {{size: 12}}}}
            }},
            yaxis: {{
                title: {{text: '논문 수', font: {{size: 12}}}}
            }},
            showlegend: false
        }};
        
        if (document.getElementById('active-papers-chart')) {{
            Plotly.newPlot('active-papers-chart', activePapersData, activePapersLayout, {{responsive: true}});
        }}
    }}
    
    // 3. 세부 주제 상위 20개 (가로 막대 차트)
    if (subtopicData.length > 0) {{
        const topSubtopics = subtopicData
            .sort((a, b) => b.count - a.count)
            .slice(0, 20)
            .reverse(); // 가로 차트를 위해 역순 정렬
        
        const subtopicsChartData = [{{
            x: topSubtopics.map(d => d.count),
            y: topSubtopics.map(d => d.subtopic.length > 50 ? d.subtopic.substring(0, 47) + '...' : d.subtopic),
            type: 'bar',
            orientation: 'h',
            marker: {{
                color: topSubtopics.map(d => {{
                    const colors = {{'마취전 관리': '#FF6B6B', '마취 약리': '#4ECDC4'}};
                    return colors[d.category_short] || '#45B7D1';
                }})
            }},
            text: topSubtopics.map(d => d.count),
            textposition: 'outside',
            hovertemplate: '<b>%{{y}}</b><br>논문 수: %{{x}}<br>분야: %{{customdata}}<extra></extra>',
            customdata: topSubtopics.map(d => d.category_short)
        }}];
        
        const subtopicsLayout = {{
            height: Math.max(600, topSubtopics.length * 35),
            font: {{family: "Arial, sans-serif", size: 10}},
            plot_bgcolor: 'rgba(0,0,0,0)',
            paper_bgcolor: 'rgba(0,0,0,0)',
            margin: {{l: 300, r: 80, t: 20, b: 60}},
            xaxis: {{
                title: {{text: '논문 수', font: {{size: 12}}}},
                tickfont: {{size: 11}}
            }},
            yaxis: {{
                tickfont: {{size: 9}},
                automargin: true
            }},
            showlegend: false
        }};
        
        if (document.getElementById('subtopics-chart')) {{
            Plotly.newPlot('subtopics-chart', subtopicsChartData, subtopicsLayout, {{responsive: true}});
        }}
    }}
    
    // 인터랙티브 기능
    document.addEventListener('DOMContentLoaded', function() {{
        // 카테고리 아이템 호버 효과
        document.querySelectorAll('.category-item').forEach(item => {{
            item.addEventListener('mouseenter', function() {{
                if (this.classList.contains('active')) {{
                    this.style.borderLeftWidth = '8px';
                    this.style.borderLeftColor = '#28a745';
                }} else {{
                    this.style.borderLeftWidth = '8px';
                    this.style.borderLeftColor = '#6c757d';
                }}
            }});
            
            item.addEventListener('mouseleave', function() {{
                this.style.borderLeftWidth = '2px';
            }});
        }});
        
        // 세부주제 카드 호버 효과
        document.querySelectorAll('.subtopic-card').forEach(card => {{
            card.addEventListener('mouseenter', function() {{
                this.style.borderLeftWidth = '10px';
            }});
            card.addEventListener('mouseleave', function() {{
                this.style.borderLeftWidth = '6px';
            }});
        }});
        
        // 부드러운 등장 애니메이션
        const observerOptions = {{
            threshold: 0.1,
            rootMargin: '0px 0px -50px 0px'
        }};
        
        const observer = new IntersectionObserver(function(entries) {{
            entries.forEach(entry => {{
                if (entry.isIntersecting) {{
                    entry.target.style.opacity = '1';
                    entry.target.style.transform = 'translateY(0)';
                }}
            }});
        }}, observerOptions);
        
        // 애니메이션 대상 요소들
        document.querySelectorAll('.overview-card, .chart-container, .full-width-chart, .category-item, .subtopic-card, .category-section').forEach((el, index) => {{
            el.style.opacity = '0';
            el.style.transform = 'translateY(30px)';
            el.style.transition = 'all 0.8s ease';
            observer.observe(el);
        }});
        
        // 논문 아이템 호버 효과
        document.querySelectorAll('.paper-item').forEach(item => {{
            item.addEventListener('mouseenter', function() {{
                this.style.transform = 'translateX(8px)';
                this.style.boxShadow = '0 12px 25px rgba(0,0,0,0.15)';
            }});
            item.addEventListener('mouseleave', function() {{
                this.style.transform = 'translateX(3px)';
                this.style.boxShadow = '0 8px 20px rgba(0,0,0,0.1)';
            }});
        }});
    }});
    """)

# HTML 저장
output_html = "index.html"  # GitHub Pages를 위해 index.html로 저장
print("💾 향상된 HTML 파일 생성 중...")
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

print("\n🏁 향상된 마취학 연구 분류 대시보드 생성이 완료되었습니다!")
print("📊 새로운 대시보드 주요 기능:")
print("   ✅ 23개 전체 마취학 분야 구조 반영")
print("   ✅ 활성/비활성 분야 구분 시각화")
print("   ✅ 전체 분야 개요 및 현황 파악")
print("   ✅ 개선된 차트 (도넛, 막대, 가로 막대)")
print("   ✅ 향상된 반응형 디자인")
print("   ✅ 부드러운 애니메이션 효과")
print("   ✅ 세부주제별 상세 정보 및 논문 링크")