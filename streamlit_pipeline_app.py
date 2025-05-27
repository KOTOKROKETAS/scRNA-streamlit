import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

# --- Sidebar Navigation ---
st.sidebar.title("Μενού Πλοήγησης")
section = st.sidebar.radio("Μετάβαση σε:", (
    "Αρχική", 
    "Προεπεξεργασία",
    "Clustering & UMAP",
    "Διαφορική Έκφραση",
    "Οπτικοποιήσεις",
    "Στοιχεία Ομάδας"
))

# --- Upload file ---
st.sidebar.markdown("---")
uploaded_file = st.sidebar.file_uploader("Ανέβασμα .h5ad αρχείου", type=["h5ad"])

# --- Load Data ---
@st.cache_data
def load_data(file=None):
    if file is not None:
        adata = sc.read_h5ad(file)
    else:
        adata = sc.read_h5ad("pancreas_data.h5ad")
    return adata

adata = load_data(uploaded_file)
# --- Προσθήκη βασικών QC μεταβλητών αν λείπουν ---
import numpy as np
if 'n_genes' not in adata.obs:
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1
if 'n_counts' not in adata.obs:
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
if 'percent_mt' not in adata.obs and 'MT-' in adata.var_names[0]:
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.obs['percent_mt'] = (adata[:, adata.var['mt']].X.sum(axis=1).A1 / adata.obs['n_counts']) * 100


# --- Αρχική ---
if section == "Αρχική":
    st.title("🔬 scRNA-seq Ανάλυση - Pancreas Dataset")
    st.markdown("""
        Αυτή η εφαρμογή επιτρέπει την ανάλυση δεδομένων μονοκυτταρικής RNA έκφρασης 
        χρησιμοποιώντας pipelines βασισμένα στο Scanpy. Επιλέξτε από το μενού στα αριστερά 
        για να εκτελέσετε διαφορετικά στάδια της ανάλυσης.
    """)

# --- Προεπεξεργασία ---
elif section == "Προεπεξεργασία":
    st.header("🧹 Ποιοτικός Έλεγχος και Φιλτράρισμα")
    st.write("Αριθμός κελιών:", adata.n_obs)
    st.write("Αριθμός γονιδίων:", adata.n_vars)

    # Παράμετροι QC
    min_genes = st.slider("Ελάχιστος αριθμός γονιδίων ανά κύτταρο", 0, 1000, 200)
    max_genes = st.slider("Μέγιστος αριθμός γονιδίων ανά κύτταρο", 1000, 10000, 5000)

    adata_qc = adata[(adata.obs['n_genes'] > min_genes) & (adata.obs['n_genes'] < max_genes), :]
    st.write("Κελιά μετά το φιλτράρισμα:", adata_qc.n_obs)

    fig, ax = plt.subplots()
    sns.histplot(adata_qc.obs['n_genes'], bins=50, ax=ax)
    st.pyplot(fig)

# --- Clustering & UMAP ---
elif section == "Clustering & UMAP":
    st.header("📊 Clustering και UMAP")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.leiden(adata, resolution=1.0)
    sc.tl.umap(adata)

    sc.pl.umap(adata, color=['leiden'], show=False)
    st.pyplot(plt.gcf())

# --- Διαφορική Έκφραση ---
elif section == "Διαφορική Έκφραση":
    st.header("📈 Διαφορική Έκφραση Γονιδίων")
    st.write("Υπολογισμός διαφορικής έκφρασης μεταξύ clusters")
    if 'leiden' not in adata.obs:
        if 'neighbors' not in adata.uns:
            sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.leiden(adata, resolution=1.0)
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, show=False)
    st.pyplot(plt.gcf())

# --- Οπτικοποιήσεις ---
elif section == "Οπτικοποιήσεις":
    st.header("🧬 Οπτικοποιήσεις Γονιδιακής Έκφρασης")
    gene = st.text_input("Εισάγετε όνομα γονιδίου", "INS")
    
    if gene not in adata.var_names:
        st.warning("Το γονίδιο δεν βρέθηκε στα δεδομένα.")
    elif 'leiden' not in adata.obs.columns:
        st.warning("Δεν υπάρχουν διαθέσιμα clusters (leiden). Εκτέλεσε πρώτα clustering.")
    else:
        sc.pl.violin(adata, gene, groupby='leiden', show=False)
        st.pyplot(plt.gcf())
        sc.pl.umap(adata, color=gene, show=False)
        st.pyplot(plt.gcf())

# --- Στοιχεία Ομάδας ---
elif section == "Στοιχεία Ομάδας":
    st.header("👥 Πληροφορίες Ομάδας")
    st.markdown("""
    **Μέλη Ομάδας**
    
    - Παναγιώτης Πεσλής – INF2021185 – Προεπεξεργασία, Pipeline σύνδεση
    - Στυλιανός Ιωαννίδης – INF2021067 – Clustering, UMAP, Οπτικοποιήσεις
    - Χρήστος Χατζηκώστας – Π2020130 – Docker, GitHub, Report LaTeX
    """)