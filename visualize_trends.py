import json
import pandas as pd
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
import google.generativeai as genai
from typing import Dict, List, Any
import requests
import time
from dotenv import load_dotenv

# .env 파일에서 환경변수 로드
load_dotenv()

# 환경변수에서 API 키 가져오기
API_KEY = os.getenv("GEMINI_API_KEY")

# GitHub 설정 (사용자가 수정해야 할 부분)
GITHUB_REPO_PATH = "."  # 현재 디렉토리가 git 레포지토리라고 가정
GITHUB_REPO_URL = "https://github.com/escapebaek/trends_anesthesia.git"
AUTO_DEPLOY = True      # 자동 배포 여부
AUTO_OPEN_BROWSER = True  # 자동으로 브라우저 열기 여부

# Gemini API 설정
GEMINI_API_KEY = os.getenv('GEMINI_API_KEY')  # 환경변수에서 API 키 읽기
if not GEMINI_API_KEY:
    print("⚠️ GEMINI_API_KEY 환경변수가 설정되지 않았습니다.")
    print("💡 다음 명령어로 설정하세요:")
    print("   export GEMINI_API_KEY='your_api_key_here'")
    GEMINI_API_KEY = input("또는 여기에 API 키를 입력하세요: ").strip()
    if not GEMINI_API_KEY:
        print("❌ API 키가 없으면 차트 생성 기능을 사용할 수 없습니다.")
        USE_GEMINI_CHARTS = False
    else:
        USE_GEMINI_CHARTS = True
else:
    USE_GEMINI_CHARTS = True

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

class GeminiChartGenerator:
    """Gemini API를 사용한 차트 생성 클래스"""
    
    def __init__(self, api_key: str):
        if api_key:
            genai.configure(api_key=api_key)
            self.model = genai.GenerativeModel('gemini-2.0-flash')
            self.enabled = True
        else:
            self.enabled = False
    
    def generate_chart_code(self, data_summary: Dict[str, Any]) -> str:
        """데이터 요약을 바탕으로 차트 생성 코드를 생성"""
        if not self.enabled:
            return ""
        
        prompt = f"""
다음 마취학 연구 데이터를 분석하여 Chart.js를 사용한 인터랙티브 차트 코드를 생성해주세요.

데이터 요약:
- 총 카테고리 수: {data_summary['total_categories']}
- 총 서브토픽 수: {data_summary['total_subtopics']}
- 총 논문 수: {data_summary['total_papers']}
- 카테고리별 논문 수: {json.dumps(data_summary['category_counts'], ensure_ascii=False, indent=2)}
- 연도별 논문 수: {json.dumps(data_summary['yearly_counts'], ensure_ascii=False, indent=2)}
- 상위 저널: {json.dumps(data_summary['top_journals'], ensure_ascii=False, indent=2)}

요구사항:
1. Chart.js를 사용해서 3-4개의 다양한 차트를 생성
2. 카테고리별 분포 (도넛 차트)
3. 연도별 트렌드 (라인 차트)
4. 상위 저널 (바 차트)
5. 반응형 디자인
6. 모던한 색상 팔레트 사용
7. HTML div 요소들과 JavaScript 코드를 모두 포함
8. 차트는 실제 데이터를 사용

응답 형식:
```html
<!-- 차트 컨테이너 HTML -->
<div class="charts-section">
    <div class="chart-container">
        <canvas id="categoryChart"></canvas>
    </div>
    <!-- 추가 차트들... -->
</div>

<script>
// Chart.js 코드
// 실제 데이터를 사용한 차트 생성
</script>
```

한국어 라벨은 그대로 사용하고, 차트가 시각적으로 매력적이고 정보가 명확하게 전달되도록 해주세요.
"""
        
        try:
            print("🤖 Gemini AI에서 차트 코드를 생성 중...")
            response = self.model.generate_content(prompt)
            
            if response and response.text:
                # HTML과 JavaScript 코드 추출
                text = response.text
                
                # ```html 블록에서 코드 추출
                if "```html" in text:
                    code_start = text.find("```html") + 7
                    code_end = text.find("```", code_start)
                    if code_end != -1:
                        return text[code_start:code_end].strip()
                
                # 전체 응답이 코드인 경우
                return text
                
            else:
                print("⚠️ Gemini API 응답이 비어있습니다.")
                return ""
                
        except Exception as e:
            print(f"❌ Gemini API 호출 실패: {e}")
            return ""
    
    def generate_fallback_charts(self, data_summary: Dict[str, Any]) -> str:
        """Gemini API 실패 시 기본 차트 생성"""
        category_data = data_summary['category_counts']
        yearly_data = data_summary['yearly_counts']
        journal_data = data_summary['top_journals']
        
        # 카테고리 데이터를 Chart.js 형식으로 변환
        categories = list(category_data.keys())
        category_values = list(category_data.values())
        
        # 연도 데이터 정렬
        sorted_years = sorted(yearly_data.items())
        years = [str(year) for year, _ in sorted_years]
        yearly_values = [count for _, count in sorted_years]
        
        # 상위 저널 데이터
        top_journals = list(journal_data.items())[:10]  # 상위 10개만
        journal_names = [name for name, _ in top_journals]
        journal_counts = [count for _, count in top_journals]
        
        return f"""
<!-- 차트 섹션 -->
<div class="charts-section">
    <h2 class="section-title">📊 Research Analytics Dashboard</h2>
    
    <div class="charts-grid">
        <div class="chart-container">
            <h3 class="chart-title">카테고리별 논문 분포</h3>
            <canvas id="categoryChart"></canvas>
        </div>
        
        <div class="chart-container">
            <h3 class="chart-title">연도별 논문 트렌드</h3>
            <canvas id="yearlyChart"></canvas>
        </div>
        
        <div class="chart-container full-width">
            <h3 class="chart-title">상위 저널별 논문 수</h3>
            <canvas id="journalChart"></canvas>
        </div>
    </div>
</div>

<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.js"></script>
<script>
// Chart.js가 완전히 로드될 때까지 대기
document.addEventListener('DOMContentLoaded', function() {{
    // Chart.js가 로드되었는지 확인
    if (typeof Chart === 'undefined') {{
        console.error('Chart.js가 로드되지 않았습니다.');
        return;
    }}

    console.log('Chart.js 버전:', Chart.version);

    // 차트 기본 설정
    Chart.defaults.font.family = 'Inter, sans-serif';
    Chart.defaults.color = '#6c757d';

    // 색상 팔레트
    const colors = [
        '#4a90e2', '#50e3c2', '#f39c12', '#e74c3c', '#9b59b6',
        '#1abc9c', '#34495e', '#f1c40f', '#e67e22', '#95a5a6',
        '#3498db', '#2ecc71', '#ff7675', '#a29bfe', '#fd79a8'
    ];

    // 카테고리 도넛 차트
    const categoryCtx = document.getElementById('categoryChart');
    if (categoryCtx) {{
        try {{
            new Chart(categoryCtx, {{{{ // Corrected: Removed extra backslashes
                type: 'doughnut',
                data: {{{{ // Corrected: Removed extra backslashes
                    labels: {json.dumps(categories, ensure_ascii=False)},
                    datasets: [{{{{ // Corrected: Removed extra backslashes
                        data: {category_values},
                        backgroundColor: colors.slice(0, {len(categories)}),
                        borderWidth: 2,
                        borderColor: '#ffffff'
                    }}}}
                }}}},
                options: {{{{ // Corrected: Removed extra backslashes
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {{{{ // Corrected: Removed extra backslashes
                        legend: {{{{ // Corrected: Removed extra backslashes
                            position: 'bottom',
                            labels: {{{{ // Corrected: Removed extra backslashes
                                padding: 20,
                                usePointStyle: true
                            }}}}
                        }}}},
                        tooltip: {{{{ // Corrected: Removed extra backslashes
                            callbacks: {{{{ // Corrected: Removed extra backslashes
                                label: function(context) {{
                                    const total = context.dataset.data.reduce((a, b) => a + b, 0);
                                    const percentage = ((context.parsed * 100) / total).toFixed(1);
                                    return context.label + ': ' + context.parsed + '편 (' + percentage + '%)';
                                }}
                            }}}}
                        }}}}
                    }}}}
                }}}}
            }}}});
            console.log('카테고리 차트 생성 완료');
        }} catch (error) {{
            console.error('카테고리 차트 생성 실패:', error);
        }}
    }}

    // 연도별 라인 차트
    const yearlyCtx = document.getElementById('yearlyChart');
    if (yearlyCtx) {{
        try {{
            new Chart(yearlyCtx, {{{{ // Corrected: Removed extra backslashes
                type: 'line',
                data: {{{{ // Corrected: Removed extra backslashes
                    labels: {json.dumps(years)},
                    datasets: [{{{{ // Corrected: Removed extra backslashes
                        label: '논문 수',
                        data: {yearly_values},
                        borderColor: '#4a90e2',
                        backgroundColor: 'rgba(74, 144, 226, 0.1)',
                        fill: true,
                        tension: 0.4,
                        borderWidth: 3,
                        pointBackgroundColor: '#4a90e2',
                        pointBorderColor: '#ffffff',
                        pointBorderWidth: 2,
                        pointRadius: 6
                    }}}}
                }}}},
                options: {{{{ // Corrected: Removed extra backslashes
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {{{{ // Corrected: Removed extra backslashes
                        y: {{{{ // Corrected: Removed extra backslashes
                            beginAtZero: true,
                            grid: {{{{ // Corrected: Removed extra backslashes
                                color: 'rgba(0,0,0,0.1)'
                            }}}}
                        }}}},
                        x: {{{{ // Corrected: Removed extra backslashes
                            grid: {{{{ // Corrected: Removed extra backslashes
                                color: 'rgba(0,0,0,0.1)'
                            }}}}
                        }}}}
                    }}}},
                    plugins: {{{{ // Corrected: Removed extra backslashes
                        legend: {{{{ // Corrected: Removed extra backslashes
                            display: false
                        }}}},
                        tooltip: {{{{ // Corrected: Removed extra backslashes
                            backgroundColor: 'rgba(0,0,0,0.8)',
                            titleColor: '#ffffff',
                            bodyColor: '#ffffff',
                            borderColor: '#4a90e2',
                            borderWidth: 1
                        }}}}
                    }}}}
                }}}}
            }}}});
            console.log('연도별 차트 생성 완료');
        }} catch (error) {{
            console.error('연도별 차트 생성 실패:', error);
        }}
    }}

    // 저널 바 차트
    const journalCtx = document.getElementById('journalChart');
    if (journalCtx) {{
        try {{
            new Chart(journalCtx, {{{{ // Corrected: Removed extra backslashes
                type: 'bar',
                data: {{{{ // Corrected: Removed extra backslashes
                    labels: {json.dumps(journal_names, ensure_ascii=False)},
                    datasets: [{{{{ // Corrected: Removed extra backslashes
                        label: '논문 수',
                        data: {journal_counts},
                        backgroundColor: colors.slice(0, {len(journal_names)}),
                        borderColor: colors.slice(0, {len(journal_names)}),
                        borderWidth: 1
                    }}}}
                }}}},
                options: {{{{ // Corrected: Removed extra backslashes
                    responsive: true,
                    maintainAspectRatio: false,
                    indexAxis: 'y',
                    scales: {{{{ // Corrected: Removed extra backslashes
                        x: {{{{ // Corrected: Removed extra backslashes
                            beginAtZero: true,
                            grid: {{{{ // Corrected: Removed extra backslashes
                                color: 'rgba(0,0,0,0.1)'
                            }}}}
                        }}}},
                        y: {{{{ // Corrected: Removed extra backslashes
                            grid: {{{{ // Corrected: Removed extra backslashes
                                display: false
                            }}}}
                        }}}}
                    }}}},
                    plugins: {{{{ // Corrected: Removed extra backslashes
                        legend: {{{{ // Corrected: Removed extra backslashes
                            display: false
                        }}}},
                        tooltip: {{{{ // Corrected: Removed extra backslashes
                            backgroundColor: 'rgba(0,0,0,0.8)',
                            titleColor: '#ffffff',
                            bodyColor: '#ffffff'
                        }}}}
                    }}}}
                }}}}
            }}}});
            console.log('저널 차트 생성 완료');
        }} catch (error) {{
            console.error('저널 차트 생성 실패:', error);
        }}
    }}

    console.log('📊 모든 차트가 성공적으로 로드되었습니다!');
}});
</script>
"""

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

def create_data_summary(classified_data: Dict, all_papers: List[Dict]) -> Dict[str, Any]:
    """차트 생성을 위한 데이터 요약 생성"""
    # 카테고리별 논문 수
    category_counts = {}
    for category, subtopics in classified_data.items():
        total = sum(len(papers) for papers in subtopics.values() if isinstance(papers, list))
        if total > 0:
            # 카테고리 이름 단순화
            simple_name = category.split("(")[0].strip()
            category_counts[simple_name] = total
    
    # 연도별 논문 수
    yearly_counts = {}
    for paper in all_papers:
        date_str = paper.get('issue_date', '')
        if date_str:
            try:
                year = int(date_str.split('-')[0])
                yearly_counts[year] = yearly_counts.get(year, 0) + 1
            except:
                pass
    
    # 상위 저널
    journal_counts = {}
    for paper in all_papers:
        journal = paper.get('journal', 'Unknown')
        if journal and journal != 'Unknown':
            journal_counts[journal] = journal_counts.get(journal, 0) + 1
    
    # 상위 10개 저널만 선택
    top_journals = dict(sorted(journal_counts.items(), key=lambda x: x[1], reverse=True)[:10])
    
    return {
        'total_categories': len(category_counts),
        'total_subtopics': sum(len(subtopics) for subtopics in classified_data.values()),
        'total_papers': len(all_papers),
        'category_counts': category_counts,
        'yearly_counts': yearly_counts,
        'top_journals': top_journals
    }

# 1. JSON 로드
json_path = "anesthesia_classified_with_metadata.json"
if not os.path.exists(json_path):
    print(f"❌ {json_path} 파일을 찾을 수 없습니다.")
    print("💡 먼저 analyze_with_gemini.py를 실행하세요.")
    sys.exit(1)

print("📊 분류된 데이터 로드 중...")
with open(json_path, "r", encoding="utf-8") as f:
    full_data = json.load(f)

classified_data = full_data.get("classified_abstracts", {})
metadata = full_data.get("metadata", {})

# 2. 개선된 데이터 전처리 
category_stats = []
subtopic_stats = []
all_papers = []

print("🔍 데이터 구조 분석 중...")
print(f"전체 카테고리 수: {len(classified_data)}")

# 데이터 구조 파악을 위한 디버깅
print("\n🐛 JSON 구조 디버깅:")
for category, subtopics in classified_data.items():
    category_count = sum(len(papers) for papers in subtopics.values() if isinstance(papers, list))
    category_subtopics = len(subtopics)

    if category_count > 0:
        category_stats.append({
            "category": category,
            "category_short": category.split("(")[0].strip(),
            "total_papers": category_count,
            "subtopics": category_subtopics
        })

    for subtopic, papers in subtopics.items():
        if papers and isinstance(papers, list):
            subtopic_stats.append({
                "category": category,
                "subtopic": subtopic,
                "count": len(papers),
                "category_short": category.split("(")[0].strip()
            })
            for paper in papers:
                if isinstance(paper, dict):
                    paper_data = paper.copy()
                    paper_data["category"] = category
                    paper_data["subtopic"] = subtopic
                    paper_data["category_short"] = category.split("(")[0].strip()
                    paper_data["first_author"] = extract_first_author(paper_data.get("author", "N/A"))
                    all_papers.append(paper_data)

print(f"- 총 논문: {len(all_papers)}개")

# DataFrame 생성
df_categories = pd.DataFrame(category_stats)
df_subtopics = pd.DataFrame(subtopic_stats)
df_papers = pd.DataFrame(all_papers)

if not df_categories.empty:
    df_categories['total_papers'] = df_categories['total_papers'].astype(int)
    df_categories['subtopics'] = df_categories['subtopics'].astype(int)

if not df_subtopics.empty:
    df_subtopics['count'] = df_subtopics['count'].astype(int)

print("\n✅ 데이터프레임 생성 완료.")

# 3. 차트 생성을 위한 데이터 요약
data_summary = create_data_summary(classified_data, all_papers)

# 4. Gemini를 사용한 차트 생성
chart_html = ""
if USE_GEMINI_CHARTS:
    print("🤖 Gemini AI로 차트 생성 중...")
    chart_generator = GeminiChartGenerator(GEMINI_API_KEY)
    chart_html = chart_generator.generate_chart_code(data_summary)
    
    if not chart_html:
        print("⚠️ Gemini 차트 생성 실패, 기본 차트 사용")
        chart_html = chart_generator.generate_fallback_charts(data_summary)
else:
    print("📊 기본 차트 생성 중...")
    chart_generator = GeminiChartGenerator(None)  # API 키 없이 초기화
    chart_html = chart_generator.generate_fallback_charts(data_summary)

print("📊 HTML 문서 생성 중...")

# HTML 문서 생성
doc, tag, text = Doc().tagtext()

def create_enhanced_css():
    return """
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700;800&display=swap');
        
        * { margin: 0; padding: 0; box-sizing: border-box; }
        
        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background-color: #f8f9fa;
            color: #343a40;
            line-height: 1.6;
        }
        
        .container {
            max-width: 1800px;
            margin: 0 auto;
            padding: 20px;
        }
        
        .header {
            text-align: center;
            padding: 50px 20px;
            background: linear-gradient(135deg, #4a90e2 0%, #50e3c2 100%);
            border-radius: 24px;
            margin-bottom: 40px;
            color: white;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        }
        
        .header h1 {
            font-size: 3em;
            font-weight: 800;
            margin-bottom: 10px;
            text-shadow: 1px 1px 3px rgba(0,0,0,0.2);
        }
        
        .header p {
            font-size: 1.2em;
            opacity: 0.9;
            font-weight: 300;
        }
        
        .header .subtitle {
            font-size: 0.9em;
            margin-top: 15px;
            opacity: 0.8;
            font-style: italic;
        }
        
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 30px;
            margin-bottom: 40px;
        }
        
        .stat-card {
            background: #ffffff;
            border-radius: 20px;
            padding: 25px;
            text-align: center;
            box-shadow: 0 8px 25px rgba(0,0,0,0.07);
            transition: all 0.3s ease;
        }
        
        .stat-card:hover {
            transform: translateY(-5px) scale(1.03);
            box-shadow: 0 12px 35px rgba(0,0,0,0.1);
        }
        
        .stat-number {
            font-size: 2.8em;
            font-weight: 800;
            color: #4a90e2;
            margin-bottom: 10px;
        }
        
        .stat-label {
            color: #6c757d;
            font-size: 1.1em;
            font-weight: 600;
        }
        
        /* 차트 섹션 스타일 */
        .charts-section {
            background: #ffffff;
            border-radius: 20px;
            margin: 40px 0;
            padding: 40px;
            box-shadow: 0 8px 25px rgba(0,0,0,0.07);
        }
        
        .section-title {
            text-align: center;
            font-size: 2.2em;
            font-weight: 800;
            color: #343a40;
            margin-bottom: 40px;
            background: linear-gradient(135deg, #4a90e2, #50e3c2);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }
        
        .charts-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 30px;
            margin-bottom: 30px;
        }
        
        .chart-container {
            background: #f8f9fa;
            border-radius: 16px;
            padding: 30px;
            box-shadow: 0 5px 15px rgba(0,0,0,0.05);
            transition: all 0.3s ease;
            height: 400px;
        }
        
        .chart-container:hover {
            transform: translateY(-5px);
            box-shadow: 0 12px 30px rgba(0,0,0,0.1);
        }
        
        .chart-container.full-width {
            grid-column: 1 / -1;
            height: 500px;
        }
        
        .chart-title {
            text-align: center;
            font-size: 1.4em;
            font-weight: 700;
            color: #343a40;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid #e9ecef;
        }
        
        .chart-container canvas {
            max-height: 300px;
        }
        
        .chart-container.full-width canvas {
            max-height: 400px;
        }
        
        .category-section {
            background: #ffffff;
            border-radius: 20px;
            margin: 30px 0;
            padding: 30px;
            box-shadow: 0 8px 25px rgba(0,0,0,0.07);
        }
        
        .category-header {
            background: linear-gradient(135deg, #4a90e2, #50e3c2);
            color: white;
            padding: 20px 30px;
            border-radius: 16px;
            margin-bottom: 30px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            flex-wrap: wrap;
        }
        
        .category-title {
            font-size: 1.5em;
            font-weight: 700;
        }
        
        .category-stats {
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
        }
        
        .category-stat {
            background: rgba(255,255,255,0.2);
            padding: 8px 18px;
            border-radius: 20px;
            font-size: 0.9em;
            font-weight: 600;
        }
        
        .subtopics-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(400px, 1fr));
            gap: 25px;
        }
        
        .subtopic-card {
            background: #f8f9fa;
            border-radius: 16px;
            padding: 25px;
            border-left: 5px solid #4a90e2;
            transition: all 0.3s ease;
            cursor: pointer;
        }
        
        .subtopic-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 12px 35px rgba(0,0,0,0.1);
            border-left-color: #50e3c2;
        }
        
        .subtopic-title {
            font-size: 1.3em;
            font-weight: 700;
            color: #343a40;
            margin-bottom: 15px;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        
        .paper-count {
            background: #4a90e2;
            color: white;
            padding: 5px 15px;
            border-radius: 15px;
            font-size: 0.8em;
            font-weight: 700;
        }
        
        .papers-list {
            max-height: 350px;
            overflow-y: auto;
            padding-right: 10px;
        }
        
        .paper-item {
            background: #ffffff;
            border-radius: 12px;
            padding: 15px;
            margin-bottom: 12px;
            border: 1px solid #e9ecef;
            transition: all 0.3s ease;
        }
        
        .paper-item:hover {
            box-shadow: 0 5px 15px rgba(0,0,0,0.08);
            border-color: #4a90e2;
            transform: translateX(5px);
        }
        
        .paper-title {
            font-weight: 700;
            color: #343a40;
            margin-bottom: 8px;
            font-size: 1.05em;
        }
        
        .paper-details {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 10px;
            font-size: 0.9em;
            color: #6c757d;
        }
        
        .paper-author {
            font-weight: 600;
            color: #495057;
        }
        
        .paper-journal {
            font-style: italic;
        }
        
        .paper-date {
            background-color: #e9ecef;
            color: #495057;
            padding: 4px 10px;
            border-radius: 10px;
            font-size: 0.8em;
            font-weight: 600;
        }
        
        .paper-summary {
            color: #495057;
            font-size: 0.95em;
            line-height: 1.5;
            margin-bottom: 12px;
        }
        
        .paper-link {
            display: inline-block;
            background-color: #28a745;
            color: white;
            text-decoration: none;
            padding: 6px 14px;
            border-radius: 12px;
            font-size: 0.85em;
            font-weight: 600;
            transition: all 0.3s ease;
        }
        
        .paper-link:hover {
            background-color: #218838;
            transform: translateY(-2px);
        }
        
        .footer {
            text-align: center;
            color: #6c757d;
            margin-top: 50px;
            padding: 20px;
            font-size: 0.9em;
        }
        
        @media (max-width: 768px) {
            .stats-grid, .subtopics-grid, .charts-grid {
                grid-template-columns: 1fr;
            }
            .header h1 { font-size: 2.2em; }
            .category-header { flex-direction: column; align-items: flex-start; gap: 15px; }
            .chart-container { height: 300px; }
            .chart-container.full-width { height: 350px; }
        }
        
        .papers-list::-webkit-scrollbar { width: 6px; }
        .papers-list::-webkit-scrollbar-track { background: #f1f1f1; border-radius: 3px; }
        .papers-list::-webkit-scrollbar-thumb { background: #ced4da; border-radius: 3px; }
        .papers-list::-webkit-scrollbar-thumb:hover { background: #adb5bd; }
        
        .loading-animation {
            opacity: 0;
            transform: translateY(20px);
            animation: fadeInUp 0.5s ease forwards;
        }
        
        @keyframes fadeInUp {
            to { opacity: 1; transform: translateY(0); }
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
            text("Anesthesia Research Trends")
        doc.asis(create_enhanced_css())
    
    with tag("body"):
        with tag("div", klass="container"):
            # 헤더
            with tag("div", klass="header loading-animation"):
                with tag("h1"):
                    text("🏥 Anesthesia Research Trends")
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

            # 차트 섹션 추가
            if chart_html:
                doc.asis(f'<div class="loading-animation" style="animation-delay: 0.3s">{chart_html}</div>')

# 카테고리별 상세 섹션
for idx, (_, cat_row) in enumerate(df_categories.sort_values('total_papers', ascending=False).iterrows()):
    category = cat_row['category']
    category_subtopics = df_subtopics[df_subtopics['category'] == category].sort_values('count', ascending=False)
    
    with tag("div", klass="category-section loading-animation", style=f"animation-delay: {(idx + 1) * 0.1}s"):
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
        if USE_GEMINI_CHARTS:
            text("🔬 Generated with Python & Gemini AI + Chart.js")
        else:
            text("🔬 Generated with Python & Chart.js")
    with tag("p", style="font-size: 0.95em;"):
        text("Anesthesia research classification system with interactive charts")
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
print("💾 HTML 대시보드 생성 중...")
with open(output_html, "w", encoding="utf-8") as f:
    f.write(doc.getvalue())

print(f"✅ 마취학 분류 대시보드 생성 완료 → {output_html}")
print("\n🔧 새로운 기능:")
if USE_GEMINI_CHARTS:
    print("   ✓ Gemini AI 기반 인터랙티브 차트 생성")
    print("   ✓ Chart.js를 사용한 반응형 차트")
else:
    print("   ✓ Chart.js를 사용한 기본 인터랙티브 차트")
print("   ✓ 카테고리별 도넛 차트")
print("   ✓ 연도별 트렌드 라인 차트")
print("   ✓ 상위 저널 바 차트")
print("   ✓ 반응형 차트 레이아웃")
print("   ✓ 모던한 색상 팔레트")

# 자동 배포 실행
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

print("\n🏁 마취학 연구 분류 대시보드 생성이 완료되었습니다!")
if USE_GEMINI_CHARTS:
    print("🤖 Gemini AI가 생성한 인터랙티브 차트가 포함되었습니다!")
else:
    print("📊 기본 인터랙티브 차트가 포함되었습니다!")
