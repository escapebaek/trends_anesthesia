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
from datetime import datetime, timedelta
import signal
import threading
import numpy as np

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

def extract_first_author(author_string):
    """저자 문자열에서 첫 번째 저자만 추출"""
    if not author_string or author_string in ['N/A', 'Unknown author']:
        return 'N/A'
    
    # 쉼표로 구분된 첫 번째 저자 추출
    first_author = author_string.split(',')[0].strip()
    
    # 'et al.' 제거
    first_author = first_author.replace(' et al.', '').replace(' et al', '')
    
    return first_author if first_author else 'N/A'

def parse_date(date_string):
    """날짜 문자열을 datetime 객체로 변환"""
    if not date_string:
        return None
    
    try:
        # 다양한 날짜 형식 처리
        if '-' in date_string:
            if len(date_string.split('-')) == 2:  # YYYY-MM 형식
                return datetime.strptime(date_string + "-01", "%Y-%m-%d")
            else:  # YYYY-MM-DD 형식
                return datetime.strptime(date_string, "%Y-%m-%d")
        else:
            return datetime.strptime(date_string + "-01-01", "%Y-%m-%d")
    except:
        return None

def get_recent_trend_data(df_papers, months=6):
    """최근 N개월간의 트렌드 데이터 생성"""
    current_date = datetime.now()
    cutoff_date = current_date - timedelta(days=months*30)
    
    # 날짜 파싱
    df_papers['parsed_date'] = df_papers['issue_date'].apply(parse_date)
    recent_papers = df_papers[df_papers['parsed_date'] >= cutoff_date].copy()
    
    if len(recent_papers) == 0:
        return pd.DataFrame()
    
    # 월별 그룹핑
    recent_papers['month_year'] = recent_papers['parsed_date'].dt.to_period('M')
    trend_data = recent_papers.groupby(['month_year', 'category_short']).size().reset_index(name='count')
    trend_data['month_year_str'] = trend_data['month_year'].astype(str)
    
    return trend_data

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

# 2. 개선된 데이터 전처리 
category_stats = []
subtopic_stats = []
all_papers = []

print("🔍 데이터 구조 분석 중...")
print(f"전체 카테고리 수: {len(classified_data)}")

# 데이터 구조 파악을 위한 디버깅
print("\n🐛 JSON 구조 디버깅:")
for i, (category, subtopics) in enumerate(list(classified_data.items())[:2]):
    print(f"  카테고리 {i+1}: {category}")
    print(f"    타입: {type(subtopics)}")
    if isinstance(subtopics, dict):
        print(f"    세부주제 수: {len(subtopics)}")
        for j, (subtopic_name, papers) in enumerate(list(subtopics.items())[:2]):
            print(f"      세부주제 {j+1}: {subtopic_name}")
            print(f"        타입: {type(papers)}")
            if isinstance(papers, list):
                print(f"        논문 수: {len(papers)}")
                if len(papers) > 0:
                    print(f"        첫 번째 논문 키: {list(papers[0].keys()) if isinstance(papers[0], dict) else 'Not a dict'}")
    print()

# 실제 데이터 전처리
for category, subtopics in classified_data.items():
    category_count = 0
    category_subtopics = 0
    
    print(f"\n📂 처리 중: {category}")
    
    if isinstance(subtopics, dict):
        print(f"   세부주제 수: {len(subtopics)}")
        
        for subtopic, papers in subtopics.items():
            if papers and isinstance(papers, list) and len(papers) > 0:
                paper_count = len(papers)
                category_count += paper_count
                category_subtopics += 1
                
                print(f"   - {subtopic}: {paper_count}개 논문")
                
                # 세부주제 통계
                subtopic_stats.append({
                    "category": category,
                    "subtopic": subtopic,
                    "count": paper_count,
                    "category_short": category.split("(")[0].strip()
                })
                
                # 개별 논문 데이터
                for paper in papers:
                    if isinstance(paper, dict):
                        paper_data = paper.copy()
                        paper_data["category"] = category
                        paper_data["subtopic"] = subtopic
                        paper_data["category_short"] = category.split("(")[0].strip()
                        # 저자 이름 개선
                        paper_data["first_author"] = extract_first_author(paper_data.get("author", "N/A"))
                        all_papers.append(paper_data)
    
    print(f"   총 논문 수: {category_count}")
    
    # 논문이 있는 카테고리만 추가 (중요: 0개인 카테고리 제외)
    if category_count > 0:
        category_stats.append({
            "category": category,
            "category_short": category.split("(")[0].strip(),
            "total_papers": category_count,
            "subtopics": category_subtopics
        })

print(f"\n📊 최종 집계:")
print(f"   - 활성 카테고리: {len(category_stats)}개")
print(f"   - 활성 세부주제: {len(subtopic_stats)}개")
print(f"   - 총 논문: {len(all_papers)}개")

# DataFrame 생성 및 데이터 검증
df_categories = pd.DataFrame(category_stats)
df_subtopics = pd.DataFrame(subtopic_stats)
df_papers = pd.DataFrame(all_papers)

# 데이터 타입 확실히 설정 (중요!)
if len(df_categories) > 0:
    df_categories['total_papers'] = df_categories['total_papers'].astype(int)
    df_categories['subtopics'] = df_categories['subtopics'].astype(int)
    
if len(df_subtopics) > 0:
    df_subtopics['count'] = df_subtopics['count'].astype(int)

print(f"\n✅ 데이터프레임 생성 완료:")
print(f"   - 활성 카테고리: {len(df_categories)}개")
print(f"   - 총 세부주제: {len(df_subtopics)}개") 
print(f"   - 총 논문: {len(df_papers)}개")

# 카테고리별 상세 정보 출력
if len(df_categories) > 0:
    print(f"\n📈 카테고리별 논문 수 분포:")
    for _, row in df_categories.sort_values('total_papers', ascending=False).iterrows():
        print(f"   - {row['category_short']}: {row['total_papers']}개 논문, {row['subtopics']}개 세부주제")
    
    print(f"\n📊 총 논문 수 검증: {df_categories['total_papers'].sum()}개")

# 최신 트렌드 데이터 생성 (날짜 정보가 있는 경우에만)
trend_data = pd.DataFrame()
if len(df_papers) > 0 and 'issue_date' in df_papers.columns:
    trend_data = get_recent_trend_data(df_papers, months=12)

# 3. 개선된 색상 팔레트 정의
modern_colors = [
    '#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', 
    '#DDA0DD', '#98D8E8', '#F7DC6F', '#BB8FCE', '#85C1E9',
    '#F8C471', '#82E0AA', '#F1948A', '#AED6F1', '#D7BDE2',
    '#A9DFBF', '#F9E79F', '#D5A6BD', '#AED6F1', '#F4D03F'
]

print("📈 차트 생성 중...")

# 4. HTML 문서 생성
doc, tag, text = Doc().tagtext()

def create_enhanced_css():
    return """
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Inter', sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            color: #333;
            line-height: 1.6;
        }
        
        .container {
            max-width: 1600px;
            margin: 0 auto;
            padding: 20px;
        }
        
        .header {
            text-align: center;
            color: white;
            margin-bottom: 50px;
            padding: 60px 0;
            position: relative;
        }
        
        .header::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
            background: rgba(255,255,255,0.1);
            border-radius: 30px;
            backdrop-filter: blur(10px);
            z-index: -1;
        }
        
        .header h1 {
            font-size: 3.5em;
            font-weight: 800;
            margin-bottom: 15px;
            text-shadow: 2px 2px 8px rgba(0,0,0,0.3);
            background: linear-gradient(45deg, #fff, #f0f0f0);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }
        
        .header p {
            font-size: 1.3em;
            opacity: 0.95;
            font-weight: 300;
        }
        
        .header .subtitle {
            font-size: 1em;
            margin-top: 10px;
            opacity: 0.8;
            font-style: italic;
        }
        
        .dashboard-grid {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 40px;
            margin-bottom: 50px;
        }
        
        .chart-container {
            background: rgba(255,255,255,0.95);
            border-radius: 25px;
            padding: 35px;
            box-shadow: 0 25px 50px rgba(0,0,0,0.15);
            backdrop-filter: blur(20px);
            border: 1px solid rgba(255,255,255,0.3);
            transition: all 0.3s ease;
            position: relative;
            overflow: hidden;
        }
        
        .chart-container::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 4px;
            background: linear-gradient(90deg, #667eea, #764ba2, #667eea);
            background-size: 200% 100%;
            animation: shimmer 3s ease-in-out infinite;
        }
        
        @keyframes shimmer {
            0%, 100% { background-position: 200% 0; }
            50% { background-position: -200% 0; }
        }
        
        .chart-container:hover {
            transform: translateY(-5px);
            box-shadow: 0 35px 70px rgba(0,0,0,0.2);
        }
        
        .full-width {
            grid-column: 1 / -1;
        }
        
        .chart-title {
            font-size: 1.6em;
            font-weight: 700;
            margin-bottom: 25px;
            color: #2c3e50;
            text-align: center;
            position: relative;
            padding-bottom: 15px;
        }
        
        .chart-title::after {
            content: '';
            position: absolute;
            bottom: 0;
            left: 50%;
            transform: translateX(-50%);
            width: 50px;
            height: 3px;
            background: linear-gradient(90deg, #667eea, #764ba2);
            border-radius: 2px;
        }
        
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 25px;
            margin: 40px 0;
        }
        
        .stat-card {
            background: rgba(255,255,255,0.95);
            border-radius: 20px;
            padding: 30px;
            text-align: center;
            box-shadow: 0 15px 35px rgba(0,0,0,0.1);
            border: 1px solid rgba(255,255,255,0.3);
            transition: all 0.3s ease;
            position: relative;
            overflow: hidden;
        }
        
        .stat-card::before {
            content: '';
            position: absolute;
            top: 0;
            left: -100%;
            width: 100%;
            height: 100%;
            background: linear-gradient(90deg, transparent, rgba(255,255,255,0.4), transparent);
            transition: left 0.5s ease;
        }
        
        .stat-card:hover::before {
            left: 100%;
        }
        
        .stat-card:hover {
            transform: translateY(-8px) scale(1.02);
            box-shadow: 0 25px 50px rgba(0,0,0,0.15);
        }
        
        .stat-number {
            font-size: 3em;
            font-weight: 800;
            background: linear-gradient(135deg, #667eea, #764ba2);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            margin-bottom: 15px;
            position: relative;
        }
        
        .stat-label {
            color: #666;
            font-size: 1.2em;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 1px;
        }
        
        .category-section {
            background: rgba(255,255,255,0.95);
            border-radius: 25px;
            margin: 40px 0;
            padding: 40px;
            box-shadow: 0 25px 50px rgba(0,0,0,0.1);
            border: 1px solid rgba(255,255,255,0.2);
        }
        
        .category-header {
            background: linear-gradient(135deg, #667eea, #764ba2);
            color: white;
            padding: 25px 35px;
            border-radius: 20px;
            margin-bottom: 30px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            flex-wrap: wrap;
            position: relative;
            overflow: hidden;
        }
        
        .category-header::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
            background: rgba(255,255,255,0.1);
            opacity: 0;
            transition: opacity 0.3s ease;
        }
        
        .category-header:hover::before {
            opacity: 1;
        }
        
        .category-title {
            font-size: 1.6em;
            font-weight: 700;
            position: relative;
            z-index: 1;
        }
        
        .category-stats {
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
            position: relative;
            z-index: 1;
        }
        
        .category-stat {
            background: rgba(255,255,255,0.25);
            padding: 10px 20px;
            border-radius: 25px;
            font-size: 0.95em;
            font-weight: 600;
            backdrop-filter: blur(10px);
            border: 1px solid rgba(255,255,255,0.3);
        }
        
        .subtopics-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(450px, 1fr));
            gap: 30px;
        }
        
        .subtopic-card {
            background: linear-gradient(145deg, #f8f9fa, #ffffff);
            border-radius: 20px;
            padding: 30px;
            border-left: 6px solid #667eea;
            transition: all 0.4s ease;
            cursor: pointer;
            box-shadow: 0 10px 25px rgba(0,0,0,0.08);
            position: relative;
            overflow: hidden;
        }
        
        .subtopic-card::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
            background: linear-gradient(135deg, rgba(102,126,234,0.05), rgba(118,75,162,0.05));
            opacity: 0;
            transition: opacity 0.3s ease;
        }
        
        .subtopic-card:hover::before {
            opacity: 1;
        }
        
        .subtopic-card:hover {
            transform: translateY(-8px) scale(1.02);
            box-shadow: 0 20px 45px rgba(0,0,0,0.15);
            border-left-width: 8px;
        }
        
        .subtopic-title {
            font-size: 1.4em;
            font-weight: 700;
            color: #2c3e50;
            margin-bottom: 20px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            position: relative;
            z-index: 1;
        }
        
        .paper-count {
            background: linear-gradient(135deg, #667eea, #764ba2);
            color: white;
            padding: 6px 16px;
            border-radius: 25px;
            font-size: 0.85em;
            font-weight: 700;
            box-shadow: 0 4px 15px rgba(102,126,234,0.3);
        }
        
        .papers-list {
            max-height: 400px;
            overflow-y: auto;
            padding-right: 15px;
            position: relative;
            z-index: 1;
        }
        
        .paper-item {
            background: rgba(255,255,255,0.9);
            border-radius: 15px;
            padding: 20px;
            margin-bottom: 15px;
            border: 1px solid rgba(0,0,0,0.05);
            transition: all 0.3s ease;
            backdrop-filter: blur(10px);
        }
        
        .paper-item:hover {
            box-shadow: 0 8px 25px rgba(0,0,0,0.12);
            border-color: #667eea;
            transform: translateX(5px);
            background: rgba(255,255,255,1);
        }
        
        .paper-title {
            font-weight: 700;
            color: #2c3e50;
            margin-bottom: 12px;
            font-size: 1.1em;
            line-height: 1.4;
        }
        
        .paper-details {
            display: grid;
            grid-template-columns: 1fr auto;
            gap: 15px;
            align-items: center;
            margin-bottom: 15px;
            font-size: 0.95em;
            color: #666;
        }
        
        .paper-author-journal {
            display: flex;
            flex-direction: column;
            gap: 5px;
        }
        
        .paper-author {
            font-weight: 600;
            color: #495057;
            font-size: 1.05em;
        }
        
        .paper-journal {
            font-style: italic;
            color: #6c757d;
        }
        
        .paper-date {
            background: linear-gradient(135deg, #e3f2fd, #bbdefb);
            color: #1976d2;
            padding: 6px 12px;
            border-radius: 15px;
            font-size: 0.85em;
            font-weight: 600;
            white-space: nowrap;
        }
        
        .paper-summary {
            color: #666;
            font-size: 1em;
            line-height: 1.5;
            margin-bottom: 15px;
            text-align: justify;
        }
        
        .paper-link {
            display: inline-block;
            background: linear-gradient(135deg, #28a745, #20c997);
            color: white;
            text-decoration: none;
            padding: 8px 16px;
            border-radius: 20px;
            font-size: 0.9em;
            font-weight: 600;
            transition: all 0.3s ease;
            box-shadow: 0 4px 15px rgba(40,167,69,0.3);
        }
        
        .paper-link:hover {
            background: linear-gradient(135deg, #218838, #1ba085);
            transform: translateY(-2px);
            box-shadow: 0 6px 20px rgba(40,167,69,0.4);
        }
        
        .footer {
            text-align: center;
            color: white;
            margin-top: 60px;
            padding: 30px;
            opacity: 0.9;
        }
        
        .footer p {
            margin-bottom: 10px;
        }
        
        @media (max-width: 768px) {
            .dashboard-grid {
                grid-template-columns: 1fr;
                gap: 20px;
            }
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
            .paper-details {
                grid-template-columns: 1fr;
                gap: 10px;
            }
        }
        
        .papers-list::-webkit-scrollbar {
            width: 8px;
        }
        
        .papers-list::-webkit-scrollbar-track {
            background: rgba(0,0,0,0.05);
            border-radius: 4px;
        }
        
        .papers-list::-webkit-scrollbar-thumb {
            background: linear-gradient(135deg, #667eea, #764ba2);
            border-radius: 4px;
        }
        
        .papers-list::-webkit-scrollbar-thumb:hover {
            background: linear-gradient(135deg, #5a6fd8, #6a4190);
        }
        
        .loading-animation {
            opacity: 0;
            transform: translateY(30px);
            animation: fadeInUp 0.6s ease forwards;
        }
        
        @keyframes fadeInUp {
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }
    </style>
    """

# HTML 구조 생성
doc.asis("<!DOCTYPE html>")
with tag("html", lang="ko"):
    with tag("head"):
        doc.asis('<meta charset="UTF-8">')
        doc.asis('<meta name="viewport" content="width=device-width, initial-scale=1.0">')
        with tag("title"):
            text("마취학 연구 분류 대시보드 - Anesthesia Research Trends")
        doc.asis(create_enhanced_css())
        doc.asis('<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>')
    
    with tag("body"):
        with tag("div", klass="container"):
            # 헤더
            with tag("div", klass="header loading-animation"):
                with tag("h1"):
                    text("🏥 마취학 연구 동향 분석")
                with tag("p"):
                    text("Anesthesia Research Trends & Classification Dashboard")
                with tag("p", klass="subtitle"):
                    if metadata.get("analysis_date"):
                        text(f"Last Analysis: {metadata['analysis_date']}")
                    else:
                        text(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
            
            # 통계 카드
            with tag("div", klass="stats-grid loading-animation"):
                with tag("div", klass="stat-card"):
                    with tag("div", klass="stat-number"):
                        text(str(len(df_categories)))
                    with tag("div", klass="stat-label"):
                        text("Research Areas")
                
                with tag("div", klass="stat-card"):
                    with tag("div", klass="stat-number"):
                        text(str(len(df_subtopics)))
                    with tag("div", klass="stat-label"):
                        text("Subtopics")
                
                with tag("div", klass="stat-card"):
                    with tag("div", klass="stat-number"):
                        text(str(len(df_papers)))
                    with tag("div", klass="stat-label"):
                        text("Classified Papers")
                
                with tag("div", klass="stat-card"):
                    with tag("div", klass="stat-number"):
                        text(str(metadata.get("total_papers_analyzed", len(df_papers))))
                    with tag("div", klass="stat-label"):
                        text("Total Analyzed")

# 차트 생성 (데이터가 있는 경우에만)
if len(df_categories) > 0:
    print(f"📊 차트 생성 시작...")
    print(f"   카테고리 데이터: {len(df_categories)}개")
    
    # 1. 수정된 카테고리별 논문 수 바 차트
    chart_data = df_categories.sort_values('total_papers', ascending=True)
    print(f"   바 차트 데이터: {len(chart_data)}개 항목")
    print(f"   데이터 타입 확인: {chart_data['total_papers'].dtype}")
    
    # 바 차트 생성 - text 파라미터 수정
    fig1 = go.Figure()
    fig1.add_trace(go.Bar(
        x=chart_data['total_papers'],
        y=chart_data['category_short'],
        orientation='h',
        text=[f'{val}' for val in chart_data['total_papers']],  # 명시적으로 문자열로 변환
        textposition='outside',
        texttemplate='%{text}',
        marker=dict(
            color=chart_data['total_papers'],
            colorscale='Viridis',
            showscale=False
        ),
        hovertemplate='<b>%{y}</b><br>Papers: %{x}<extra></extra>'
    ))
    
    fig1.update_layout(
        title="📊 Research Distribution by Category",
        xaxis_title="Number of Papers",
        yaxis_title="Research Category",
        height=max(500, len(df_categories) * 50),
        font=dict(family="Inter, Arial, sans-serif", size=12),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        title_font_size=20,
        title_x=0.5,
        margin=dict(l=200, r=100, t=100, b=80),
        yaxis=dict(tickfont=dict(size=11)),
        xaxis=dict(tickfont=dict(size=12), title_font_size=14)
    )

    # 2. 수정된 도넛 차트
    pie_data = df_categories[df_categories['total_papers'] > 0].copy()
    print(f"   파이 차트 데이터: {len(pie_data)}개 항목")
    print(f"   파이 차트 값 확인: {pie_data['total_papers'].tolist()}")
    
    fig2 = go.Figure()
    fig2.add_trace(go.Pie(
        labels=pie_data['category_short'],
        values=pie_data['total_papers'],
        hole=0.5,
        textposition='auto',
        textinfo='percent+label',
        hovertemplate='<b>%{label}</b><br>Papers: %{value}<br>Percentage: %{percent}<extra></extra>',
        textfont_size=11,
        marker=dict(colors=modern_colors[:len(pie_data)])
    ))
    
    fig2.update_layout(
        title="🥧 Research Category Distribution",
        font=dict(family="Inter, Arial, sans-serif", size=12),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        title_font_size=20,
        title_x=0.5,
        showlegend=True,
        legend=dict(orientation="v", yanchor="middle", y=0.5, xanchor="left", x=1.02)
    )

    # 3. 수정된 세부주제 상위 15개 차트
    if len(df_subtopics) > 0:
        top_subtopics = df_subtopics.sort_values('count', ascending=True).tail(15)
        print(f"   세부주제 차트 데이터: {len(top_subtopics)}개 항목")
        print(f"   세부주제 값 확인: {top_subtopics['count'].tolist()}")
        
        # 카테고리별 색상 매핑
        category_color_map = {}
        unique_categories = top_subtopics['category_short'].unique()
        for i, cat in enumerate(unique_categories):
            category_color_map[cat] = modern_colors[i % len(modern_colors)]
        
        colors = [category_color_map[cat] for cat in top_subtopics['category_short']]
        
        fig3 = go.Figure()
        fig3.add_trace(go.Bar(
            x=top_subtopics['count'],
            y=top_subtopics['subtopic'],
            orientation='h',
            text=[f'{val}' for val in top_subtopics['count']],  # 명시적으로 문자열로 변환
            textposition='outside',
            texttemplate='%{text}',
            marker=dict(color=colors),
            hovertemplate='<b>%{y}</b><br>Papers: %{x}<br>Category: %{customdata}<extra></extra>',
            customdata=top_subtopics['category_short']
        ))
        
        fig3.update_layout(
            title="🔍 Top 15 Research Subtopics",
            xaxis_title="Number of Papers",
            yaxis_title="Research Subtopic",
            height=700,
            font=dict(family="Inter, Arial, sans-serif", size=11),
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            title_font_size=20,
            title_x=0.5,
            margin=dict(l=300, r=100, t=100, b=80),
            yaxis=dict(tickfont=dict(size=10)),
            xaxis=dict(tickfont=dict(size=12))
        )
    else:
        print("   ⚠️ 세부주제 데이터가 없습니다.")
        fig3 = None

    # 4. 최신 트렌드 차트 (시간별 논문 발행 동향)
    if len(trend_data) > 0:
        print(f"   트렌드 차트 데이터: {len(trend_data)}개 데이터포인트")
        fig4 = px.line(
            trend_data,
            x='month_year_str',
            y='count',
            color='category_short',
            title="📈 Recent Publication Trends (Last 12 Months)",
            labels={"count": "Number of Papers", "month_year_str": "Month", "category_short": "Category"},
            color_discrete_sequence=modern_colors,
            markers=True
        )
        fig4.update_layout(
            font=dict(family="Inter, Arial, sans-serif", size=12),
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            title_font_size=20,
            title_x=0.5,
            xaxis_title="Publication Month",
            yaxis_title="Number of Papers",
            legend=dict(title="Research Category"),
            hovermode='x unified',
            height=500
        )
        fig4.update_traces(
            mode='lines+markers',
            line=dict(width=3),
            marker=dict(size=8)
        )
    else:
        print("   ⚠️ 트렌드 데이터가 없습니다.")
        fig4 = None

    # HTML에 차트 추가
    with tag("div", klass="dashboard-grid loading-animation"):
        with tag("div", klass="chart-container"):
            with tag("div", klass="chart-title"):
                text("Research Distribution by Category")
            doc.asis(fig1.to_html(full_html=False, include_plotlyjs=False, div_id="category-chart"))
        
        with tag("div", klass="chart-container"):
            with tag("div", klass="chart-title"):
                text("Category Distribution Overview")
            doc.asis(fig2.to_html(full_html=False, include_plotlyjs=False, div_id="category-pie"))
    
    # 세부주제 차트 (데이터가 있는 경우만)
    if fig3 is not None:
        with tag("div", klass="chart-container full-width loading-animation"):
            with tag("div", klass="chart-title"):
                text("Top Research Subtopics")
            doc.asis(fig3.to_html(full_html=False, include_plotlyjs=False, div_id="subtopic-chart"))

    # 트렌드 차트 추가 (데이터가 있는 경우)
    if fig4 is not None:
        with tag("div", klass="chart-container full-width loading-animation"):
            with tag("div", klass="chart-title"):
                text("Publication Trends Over Time")
            doc.asis(fig4.to_html(full_html=False, include_plotlyjs=False, div_id="trend-chart"))

else:
    # 데이터가 없는 경우 오류 메시지 표시
    with tag("div", klass="chart-container full-width loading-animation"):
        with tag("div", style="text-align: center; padding: 50px;"):
            with tag("h2", style="color: #e74c3c;"):
                text("⚠️ 데이터 로딩 오류")
            with tag("p", style="color: #666; font-size: 1.1em;"):
                text("분류된 데이터를 찾을 수 없습니다. JSON 파일 구조를 확인해주세요.")
            with tag("p", style="color: #666; margin-top: 20px;"):
                text("예상 JSON 구조: {'카테고리': {'세부주제': [논문리스트]}}")

print("📊 차트 생성 완료!")

# 카테고리별 상세 섹션
for idx, (_, cat_row) in enumerate(df_categories.sort_values('total_papers', ascending=False).iterrows()):
    category = cat_row['category']
    category_subtopics = df_subtopics[df_subtopics['category'] == category].sort_values('count', ascending=False)
    
    with tag("div", klass="category-section loading-animation", style=f"animation-delay: {idx * 0.1}s"):
        with tag("div", klass="category-header"):
            with tag("div", klass="category-title"):
                text(f"📚 {category}")
            with tag("div", klass="category-stats"):
                with tag("div", klass="category-stat"):
                    text(f"📄 {cat_row['total_papers']} papers")
                with tag("div", klass="category-stat"):
                    text(f"🔍 {cat_row['subtopics']} subtopics")
        
        with tag("div", klass="subtopics-grid"):
            for _, subtopic_row in category_subtopics.iterrows():
                subtopic = subtopic_row['subtopic']
                papers = [p for p in all_papers if p['category'] == category and p['subtopic'] == subtopic]
                
                with tag("div", klass="subtopic-card"):
                    with tag("div", klass="subtopic-title"):
                        with tag("span"):
                            text(subtopic)
                        with tag("span", klass="paper-count"):
                            text(f"{len(papers)}")
                    
                    with tag("div", klass="papers-list"):
                        for paper in papers:
                            with tag("div", klass="paper-item"):
                                with tag("div", klass="paper-title"):
                                    text(paper.get('title', 'No title available'))
                                
                                with tag("div", klass="paper-details"):
                                    with tag("div", klass="paper-author-journal"):
                                        with tag("div", klass="paper-author"):
                                            author = paper.get('first_author', 'N/A')
                                            if author != 'N/A':
                                                text(f"👨‍⚕️ {author}")
                                                if paper.get('author', '').count(',') > 0:
                                                    text(" et al.")
                                            else:
                                                text("👨‍⚕️ Author not available")
                                        
                                        with tag("div", klass="paper-journal"):
                                            text(f"📖 {paper.get('journal', 'Unknown journal')}")
                                    
                                    if paper.get('issue_date'):
                                        with tag("div", klass="paper-date"):
                                            text(f"📅 {paper['issue_date']}")
                                
                                if paper.get('abstract_summary'):
                                    with tag("div", klass="paper-summary"):
                                        summary = paper['abstract_summary']
                                        if len(summary) > 300:
                                            summary = summary[:300] + "..."
                                        text(summary)
                                
                                if paper.get('link'):
                                    with tag("a", href=paper['link'], target="_blank", klass="paper-link"):
                                        text("📄 View on PubMed")

# 푸터 추가
with tag("div", klass="footer"):
    with tag("p", style="font-size: 1.1em; font-weight: 600;"):
        text("🔬 Generated with Python, Gemini AI & Advanced Data Visualization")
    with tag("p", style="font-size: 0.95em;"):
        text("Real-time anesthesia research trends and classification system")
    if metadata.get("date_range", {}).get("oldest_formatted"):
        with tag("p", style="font-size: 0.9em; margin-top: 10px; opacity: 0.8;"):
            date_range = metadata["date_range"]
            text(f"📊 Publication Range: {date_range['oldest_formatted']} ~ {date_range['newest_formatted']}")

# JavaScript 추가
with tag("script"):
    doc.asis("""
    document.addEventListener('DOMContentLoaded', function() {
        const observerOptions = {
            threshold: 0.1,
            rootMargin: '0px 0px -50px 0px'
        };

        const observer = new IntersectionObserver((entries) => {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    entry.target.style.opacity = '1';
                    entry.target.style.transform = 'translateY(0)';
                }
            });
        }, observerOptions);

        document.querySelectorAll('.loading-animation').forEach(el => {
            el.style.opacity = '0';
            el.style.transform = 'translateY(30px)';
            el.style.transition = 'all 0.6s ease';
            observer.observe(el);
        });
    });
    
    document.querySelectorAll('.subtopic-card').forEach(card => {
        card.addEventListener('mouseenter', function() {
            this.style.borderLeftWidth = '10px';
            this.style.borderLeftColor = '#764ba2';
        });
        card.addEventListener('mouseleave', function() {
            this.style.borderLeftWidth = '6px';
            this.style.borderLeftColor = '#667eea';
        });
    });
    
    document.querySelectorAll('.paper-item').forEach(item => {
        item.addEventListener('mouseenter', function() {
            this.style.transform = 'translateX(10px) scale(1.02)';
            this.style.boxShadow = '0 12px 35px rgba(102,126,234,0.2)';
        });
        item.addEventListener('mouseleave', function() {
            this.style.transform = 'translateX(0) scale(1)';
            this.style.boxShadow = '0 8px 25px rgba(0,0,0,0.12)';
        });
    });

    function animateCounters() {
        document.querySelectorAll('.stat-number').forEach(counter => {
            const target = parseInt(counter.textContent);
            const increment = target / 50;
            let current = 0;
            
            const timer = setInterval(() => {
                current += increment;
                if (current >= target) {
                    counter.textContent = target;
                    clearInterval(timer);
                } else {
                    counter.textContent = Math.floor(current);
                }
            }, 30);
        });
    }

    setTimeout(animateCounters, 500);

    console.log('🏥 Anesthesia Research Dashboard loaded successfully!');
    console.log('📊 Total categories:', """ + str(len(df_categories)) + """);
    console.log('📄 Total papers:', """ + str(len(df_papers)) + """);
    """)

# HTML 저장
output_html = "index.html"
print("💾 수정된 HTML 대시보드 생성 중...")
with open(output_html, "w", encoding="utf-8") as f:
    f.write(doc.getvalue())

print(f"✅ 수정된 마취학 분류 대시보드 생성 완료 → {output_html}")
print("\n🔧 주요 수정사항:")
print("   ✓ 바 차트 숫자 표시 문제 해결 (plotly.graph_objects 사용)")
print("   ✓ 파이 차트 카운팅 정확성 개선 (0개 카테고리 제외)")
print("   ✓ 세부주제 차트 정확성 개선")
print("   ✓ 불필요한 저널 차트 제거")
print("   ✓ 데이터 타입 명시적 설정")
print("   ✓ 차트 텍스트 레이블 수정")

# 자동 배포 실행 (기존 코드와 동일)
if AUTO_DEPLOY:
    print("\n🚀 GitHub Pages 자동 배포를 시작합니다...")
    
    if setup_git_repo():
        pages_url = deploy_to_github()
        if pages_url:
            print("🎉 배포가 완료되었습니다!")
            print(f"🌐 대시보드 URL: {pages_url}")
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
    try:
        webbrowser.open("file://" + os.path.abspath(output_html))
        print("🌐 로컬 브라우저에서 대시보드를 열었습니다.")
    except Exception:
        print(f"📁 수동으로 파일을 열어주세요: {os.path.abspath(output_html)}")

print("\n🏁 수정된 마취학 연구 분류 대시보드 생성이 완료되었습니다!")