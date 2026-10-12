import streamlit as st
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
import mygene
import gdown  # 👈 新增：用來下載 Google Drive 檔案

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

# 圖表全域設定
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# ==========================================
# 1. Google Drive 資料下載區塊
# ==========================================
DATA_DIR = "data"
GDRIVE_FOLDER_ID = "1FWAHofVcHj832vjwO7qrTZRv5A0-fuYx" # 從您的網址萃取出的資料夾 ID

# 使用 cache 確保伺服器重啟時只會下載一次
@st.cache_resource(show_spinner=False)
def download_data_from_gdrive():
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)
    
    # 檢查資料夾內是否有足夠的檔案，若少於 113 個則觸發下載
    existing_files = glob.glob(os.path.join(DATA_DIR, "*.tsv"))
    if len(existing_files) < 113:
        try:
            # gdown 會將雲端資料夾內容下載到 DATA_DIR
            url = f"https://drive.google.com/drive/folders/{GDRIVE_FOLDER_ID}?usp=sharing"
            gdown.download_folder(url, output=DATA_DIR, quiet=False, use_cookies=False)
        except Exception as e:
            st.error(f"Failed to download data from Google Drive: {e}")

# 在頁面載入時觸發下載
with st.spinner("Initializing database from Google Drive... This will only take a moment on the first run."):
    download_data_from_gdrive()

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

# ==========================================
# 3. 核心分析函數 (讀取已下載的檔案)
# ==========================================
def get_expression(gsm_id, target_id):
    # 支援尋找下載到 data 裡的檔案
    search_path = os.path.join(DATA_DIR, f"{gsm_id}*.tsv")
    
    # gdown 下載資料夾時，有時會建立多一層子資料夾，所以加個遞迴搜尋保障
    files = glob.glob(search_path)
    if not files:
        files = glob.glob(os.path.join(DATA_DIR, "**", f"{gsm_id}*.tsv"), recursive=True)
        
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
    TARGET_GENE_SYMBOL = st.text_input("Enter Gene Symbol:", value="BECN1").strip().upper()
    analyze_btn = st.button("Run Analysis 🚀")

if analyze_btn and TARGET_GENE_SYMBOL:
    with st.spinner(f"Querying Entrez ID for {TARGET_GENE_SYMBOL}..."):
        mg = mygene.MyGeneInfo()
        query_res = mg.query(TARGET_GENE_SYMBOL, scopes='symbol', fields='entrezgene', species='human')
        
    if query_res['hits'] and 'entrezgene' in query_res['hits'][0]:
        TARGET_ID = int(query_res['hits'][0]['entrezgene'])
        st.success(f"Gene Found! **{TARGET_GENE_SYMBOL}** (Entrez ID: {TARGET_ID})")
        
        with st.spinner("Extracting expression data..."):
            meta_df['Expression'] = meta_df['GSM_ID'].apply(lambda x: get_expression(x, TARGET_ID))
            
            df_plot = meta_df.dropna(subset=['Expression']).copy()
            df_plot['Expression'] = pd.to_numeric(df_plot['Expression'])
            
            if df_plot.empty:
                st.error("No expression data found. Ensure files were successfully downloaded from Google Drive.")
            else:
                # 繪圖區塊
                fig, ax = plt.subplots(figsize=(12, 7))
                sns.set_style("ticks")

                # 自訂顏色：Control 為白色，AD 為紅色
                custom_colors = {"Control": "#FFFFFF", "AD": "#FF0000"}
                
                # 畫柱狀圖 (加上黑色邊框避免白色柱體消失)
                sns.barplot(data=df_plot, x="Cell_Type", y="Expression", hue="Group", 
                            palette=custom_colors, edgecolor=".2", capsize=.1, errorbar=('ci', 95), ax=ax)
                
                # 加上個別樣本點
                sns.stripplot(data=df_plot, x="Cell_Type", y="Expression", hue="Group", 
                              dodge=True, color="black", alpha=0.5, size=5, legend=False, ax=ax)

                # 統計檢定與標記
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

                ax.set_title(f"{TARGET_GENE_SYMBOL} Expression Levels (nRPKM)", fontsize=18, fontweight='bold', pad=15)
                ax.set_ylabel("Normalized Expression (nRPKM)", fontsize=16)
                ax.set_xlabel("Cell Type", fontsize=16)
                
                # 整理圖例
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(handles[:2], labels[:2], title="Condition", frameon=False, fontsize=12, title_fontsize=14)
                sns.despine()
                plt.tight_layout()
                
                st.pyplot(fig)
                
                with st.expander("Show Raw Data Table"):
                    st.dataframe(df_plot)
    else:
        st.error(f"Cannot find gene '{TARGET_GENE_SYMBOL}'. Please check the spelling.")