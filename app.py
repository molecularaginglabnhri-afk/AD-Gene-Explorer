import streamlit as st
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
import mygene
import zipfile
import requests

# ==========================================
# 0. 網頁基本設定
# ==========================================
st.set_page_config(page_title="GSE125050 Gene Expression Explorer", layout="wide")
st.title("🧠 GSE125050 Cell-Type Specific Gene Expression")

st.markdown("""
📌 **Data Source & System Attribution:**
* **Raw Dataset Source:** NCBI GEO ([GSE125050](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125050)) — Cell-type-specific RNA profiles derived from the post-mortem frozen superior frontal gyrus of Alzheimer's disease (AD) patients and non-AD controls.
* **System Development & Data Processing:** Built and processed into this online explorer system by **Wendy Hsu's Lab at NHRI, Taiwan**.
""")

plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# ==========================================
# 1. 自動下載與解壓縮區塊
# ==========================================
@st.cache_resource(show_spinner=False)
def prepare_data():
    if not os.path.exists("extracted_data"):
        # 使用 GitHub Release 網址
        zip_url = "https://github.com/molecularaginglabnhri-afk/AD-Gene-Explorer/releases/download/v1/data.zip"
        
        # 如果壓縮檔不在，就先下載
        if not os.path.exists("data.zip"):
            try:
                r = requests.get(zip_url, stream=True)
                r.raise_for_status()
                with open("data.zip", "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            except Exception as e:
                st.error(f"Download failed: {e}")
                return False
        
        # 下載完後進行解壓縮
        try:
            with zipfile.ZipFile("data.zip", 'r') as zip_ref:
                zip_ref.extractall("extracted_data")
        except Exception as e:
            st.error(f"Extraction failed: {e}")
            return False
            
    return True

with st.spinner("Downloading and extracting dataset (approx. 31MB)... This will only take a moment on the first run."):
    data_ready = prepare_data()

# ==========================================
# 2. 參數與資料集設定 (113個樣本)
# ==========================================
sample_list = [
    ["GSM3561843", "10-26", "Control", "myeloid"], ["GSM3561844", "02-01", "Control", "endothelial"],
    ["GSM3561845", "10-63", "Control", "neuron"], ["GSM3561846", "06-26", "AD", "neuron"],
    ["GSM3561847", "13-40", "Control", "myeloid"], ["GSM3561848", "06-26", "AD", "myeloid"],
    ["GSM3561849", "10-70", "Control", "myeloid"], ["GSM3561850", "07-04", "Control", "endothelial"],
    ["GSM3561851", "04-54", "AD", "endothelial"], ["GSM3561852", "07-28", "Control", "neuron"],
    ["GSM3561853", "04-54", "AD", "neuron"], ["GSM3561854", "13-74", "AD", "neuron"],
    ["GSM3561855", "13-49", "Control", "endothelial"], ["GSM3561856", "02-01", "Control", "astrocyte"],
    ["GSM3561857", "10-70", "Control", "astrocyte"], ["GSM3561858", "08-55", "Control", "neuron"],
    ["GSM3561859", "09-35", "Control", "neuron"], ["GSM3561860", "11-98", "Control", "endothelial"],
    ["GSM3561861", "11-105", "AD", "neuron"], ["GSM3561862", "13-49", "Control", "neuron"],
    ["GSM3561863", "05-58", "AD", "myeloid"], ["GSM3561864", "05-58", "AD", "neuron"],
    ["GSM3561865", "10-80", "AD", "endothelial"], ["GSM3561866", "10-22", "Control", "neuron"],
    ["GSM3561867", "13-66", "AD", "neuron"], ["GSM3561868", "10-63", "Control", "endothelial"],
    ["GSM3561869", "09-20", "Control", "neuron"], ["GSM3561870", "00-54", "AD", "astrocyte"],
    ["GSM3561871", "02-01", "Control", "neuron"], ["GSM3561872", "10-80", "AD", "astrocyte"],
    ["GSM3561873", "03-63", "Control", "neuron"], ["GSM3561874", "03-63", "Control", "endothelial"],
    ["GSM3561875", "02-27", "Control", "endothelial"], ["GSM3561876", "03-59", "AD", "endothelial"],
    ["GSM3561877", "03-59", "AD", "myeloid"], ["GSM3561878", "08-55", "Control", "endothelial"],
    ["GSM3561879", "03-15", "Control", "astrocyte"], ["GSM3561880", "04-54", "AD", "astrocyte"],
    ["GSM3561881", "00-54", "AD", "neuron"], ["GSM3561882", "07-04", "Control", "myeloid"],
    ["GSM3561883", "07-15", "AD", "neuron"], ["GSM3561884", "14-29", "AD", "myeloid"],
    ["GSM3561885", "11-98", "Control", "astrocyte"], ["GSM3561886", "04-54", "AD", "myeloid"],
    ["GSM3561887", "07-04", "Control", "neuron"], ["GSM3561888", "11-105", "AD", "myeloid"],
    ["GSM3561889", "04-44", "AD", "neuron"], ["GSM3561890", "06-26", "AD", "endothelial"],
    ["GSM3561891", "14-29", "AD", "endothelial"], ["GSM3561892", "07-20", "AD", "neuron"],
    ["GSM3561893", "10-70", "Control", "neuron"], ["GSM3561894", "02-27", "Control", "astrocyte"],
    ["GSM3561895", "10-26", "Control", "endothelial"], ["GSM3561896", "03-15", "Control", "neuron"],
    ["GSM3561897", "05-10", "Control", "endothelial"], ["GSM3561898", "07-04", "Control", "astrocyte"],
    ["GSM3561899", "13-47", "AD", "astrocyte"], ["GSM3561900", "13-54", "AD", "neuron"],
    ["GSM3561901", "00-54", "AD", "myeloid"], ["GSM3561902", "13-40", "Control", "neuron"],
    ["GSM3561903", "11-38", "AD", "neuron"], ["GSM3561904", "11-98", "Control", "neuron"],
    ["GSM3561905", "05-10", "Control", "astrocyte"], ["GSM3561906", "05-10", "Control", "neuron"],
    ["GSM3561907", "13-47", "AD", "neuron"], ["GSM3561908", "14-29", "AD", "neuron"],
    ["GSM3561909", "11-105", "AD", "endothelial"], ["GSM3561910", "02-27", "Control", "neuron"],
    ["GSM3561911", "10-22", "Control", "myeloid"], ["GSM3561912", "03-15", "Control", "endothelial"],
    ["GSM3561913", "05-10", "Control", "myeloid"], ["GSM3561914", "10-80", "AD", "myeloid"],
    ["GSM3561915", "03-59", "AD", "neuron"], ["GSM3561916", "10-80", "AD", "neuron"],
    ["GSM3561917", "11-98", "Control", "myeloid"], ["GSM3561918", "02-01", "Control", "myeloid"],
    ["GSM3561919", "10-26", "Control", "astrocyte"], ["GSM3561920", "04-44", "AD", "endothelial"],
    ["GSM3561921", "11-49", "AD", "neuron"], ["GSM3561922", "10-26", "Control", "neuron"],
    ["GSM3561923", "11-49", "AD", "astrocyte"], ["GSM3561924", "07-28", "Control", "myeloid"],
    ["GSM3561925", "02-27", "Control", "myeloid"], ["GSM3561926", "10-63", "Control", "astrocyte"],
    ["GSM3561927", "10-63", "Control", "myeloid"], ["GSM3561928", "09-35", "Control", "endothelial"],
    ["GSM3561929", "03-63", "Control", "myeloid"], ["GSM3561930", "10-22", "Control", "endothelial"],
    ["GSM3561931", "09-35", "Control", "myeloid"], ["GSM3561932", "07-28", "Control", "endothelial"],
    ["GSM3561933", "10-79", "AD", "neuron"], ["GSM3561934", "10-20", "Control", "neuron"],
    ["GSM3561935", "03-21", "AD", "neuron"], ["GSM3561936", "10-20", "Control", "myeloid"],
    ["GSM3561937", "08-41", "AD", "neuron"], ["GSM3561938", "03-21", "AD", "myeloid"],
    ["GSM3561939", "10-39", "Control", "myeloid"], ["GSM3561940", "10-20", "Control", "astrocyte"],
    ["GSM3561941", "08-90", "Control", "endothelial"], ["GSM3561942", "04-33", "AD", "myeloid"],
    ["GSM3561943", "03-21", "AD", "endothelial"], ["GSM3561944", "10-39", "Control", "astrocyte"],
    ["GSM3561945", "09-50", "Control", "neuron"], ["GSM3561946", "10-79", "AD", "endothelial"],
    ["GSM3561947", "10-39", "Control", "neuron"], ["GSM3561948", "09-50", "Control", "endothelial"],
    ["GSM3561949", "08-90", "Control", "astrocyte"], ["GSM3561950", "08-90", "Control", "neuron"],
    ["GSM3561951", "04-33", "AD", "astrocyte"], ["GSM3561952", "04-33", "AD", "neuron"],
    ["GSM3561953", "08-41", "AD", "endothelial"], ["GSM3561954", "10-39", "Control", "endothelial"],
    ["GSM3561955", "03-21", "AD", "astrocyte"]
]

meta_df = pd.DataFrame(sample_list, columns=["GSM_ID", "ID", "Group", "Cell_Type"])
meta_df['Cell_Type'] = meta_df['Cell_Type'].str.capitalize()

# ==========================================
# 3. 核心分析函數 (從解壓縮出來的檔案中搜尋)
# ==========================================
def get_expression(gsm_id, target_id):
    if not data_ready:
        return None
        
    search_path = os.path.join("extracted_data", "**", f"{gsm_id}*.tsv")
    files = glob.glob(search_path, recursive=True)
        
    if not files:
        return None
        
    try:
        df = pd.read_csv(files[0], sep='\t', skiprows=3)
        val = df.loc[df['ID_REF'] == target_id, 'VALUE']
        return val.values[0] if not val.empty else 0.0
    except Exception:
        return None

# ==========================================
# 4. 側邊欄與視覺化執行區
# ==========================================
with st.sidebar:
    st.header("Search Parameters")
    # 將使用者輸入自動轉為大寫並去除空白
    TARGET_GENE_SYMBOL = st.text_input("Enter Gene Symbol:", value="BECN1").strip().upper()
    analyze_btn = st.button("Run Analysis 🚀")

if analyze_btn and TARGET_GENE_SYMBOL:
    if not data_ready:
        st.error("System is still preparing the dataset or failed to download it. Please try again.")
    else:
        with st.spinner(f"Querying Entrez ID for {TARGET_GENE_SYMBOL}..."):
            mg = mygene.MyGeneInfo()
            # 擴大搜尋範圍包含 symbol 與 alias
            query_res = mg.query(TARGET_GENE_SYMBOL, scopes='symbol,alias', fields='entrezgene,symbol,alias', species='human')
            
        TARGET_ID = None
        OFFICIAL_SYMBOL = ""
        
        # 進行嚴格的「完全比對 (Exact Match)」
        if query_res.get('hits'):
            for hit in query_res['hits']:
                if 'entrezgene' not in hit:
                    continue
                    
                hit_symbol = hit.get('symbol', '').upper()
                
                hit_aliases = hit.get('alias', [])
                if isinstance(hit_aliases, str):
                    hit_aliases = [hit_aliases]
                hit_aliases = [a.upper() for a in hit_aliases]
                
                if TARGET_GENE_SYMBOL == hit_symbol or TARGET_GENE_SYMBOL in hit_aliases:
                    TARGET_ID = int(hit['entrezgene'])
                    OFFICIAL_SYMBOL = hit.get('symbol', TARGET_GENE_SYMBOL)
                    break
        
        if TARGET_ID:
            st.success(f"Exact Match Found! **{OFFICIAL_SYMBOL}** (Entrez ID: {TARGET_ID})")
            TARGET_GENE_SYMBOL = OFFICIAL_SYMBOL # 將標題統一替換為官方名稱
            
            with st.spinner("Extracting expression data from 113 files..."):
                meta_df['Expression'] = meta_df['GSM_ID'].apply(lambda x: get_expression(x, TARGET_ID))
                
                df_plot = meta_df.dropna(subset=['Expression']).copy()
                df_plot['Expression'] = pd.to_numeric(df_plot['Expression'])
                
                if df_plot.empty:
                    st.error("No expression data found. Ensure TSV files are inside data.zip.")
                else:
                    fig, ax = plt.subplots(figsize=(12, 7))
                    sns.set_style("ticks")

                    custom_colors = {"Control": "#FFFFFF", "AD": "#FF0000"}
                    
                    sns.barplot(data=df_plot, x="Cell_Type", y="Expression", hue="Group", 
                                palette=custom_colors, edgecolor=".2", capsize=.1, errorbar=('ci', 95), ax=ax)
                    
                    sns.stripplot(data=df_plot, x="Cell_Type", y="Expression", hue="Group", 
                                  dodge=True, color="black", alpha=0.5, size=5, legend=False, ax=ax)

                    cell_types = df_plot['Cell_Type'].unique()
                    for i, cell in enumerate(cell_types):
                        ctrl = df_plot[(df_plot['Cell_Type'] == cell) & (df_plot['Group'] == 'Control')]['Expression']
                        ad = df_plot[(df_plot['Cell_Type'] == cell) & (df_plot['Group'] == 'AD')]['Expression']
                        
                        if len(ctrl) > 1 and len(ad) > 1:
                            t_stat, p_val = ttest_ind(ctrl, ad)
                            if p_val < 0.001: sig = "***"
                            elif p_val < 0.01: sig = "**"
                            elif p_val < 0.05: sig = "*"
                            else: sig = "ns"
                            
                            y_pos = df_plot[df_plot['Cell_Type'] == cell]['Expression'].max() * 0.95
                            ax.text(i, y_pos, sig, ha='center', va='bottom', fontsize=18, fontweight='bold')

                    # 更新標題與 Y 軸標籤
                    ax.set_title(f"{TARGET_GENE_SYMBOL} Expression Levels", fontsize=18, fontweight='bold', pad=15)
                    ax.set_ylabel("Arbitrary units", fontsize=20)
                    ax.set_xlabel("Cell Type", fontsize=20)
                    
                    handles, labels = ax.get_legend_handles_labels()
                    ax.legend(handles[:2], labels[:2], title="Condition", frameon=False, fontsize=20, title_fontsize=20)
                    sns.despine()
                    plt.tight_layout()
                    
                    st.pyplot(fig)
                    
                    with st.expander("Show Raw Data Table"):
                        st.dataframe(df_plot)
        else:
            st.error(f"Cannot find an EXACT match for gene '{TARGET_GENE_SYMBOL}'. Please check the spelling.")
