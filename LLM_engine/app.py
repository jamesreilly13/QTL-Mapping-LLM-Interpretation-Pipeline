import streamlit as st
from full_gene_engine import FullGeneEngine
import seaborn as sns
import matplotlib.pyplot as plt

st.set_page_config(layout="wide")
st.title("QTL Full Gene Analysis")

@st.cache_resource
def load_engine():
    return FullGeneEngine()

engine = load_engine()

# ---------------- GROUP SELECTION ----------------
groups = ["All"] + sorted(engine.qtl_data.keys())
selected_group = st.selectbox("Behavioral Group", groups)

# ---------------- PHENOTYPE SELECTION ----------------
if selected_group == "All":
    phenotypes = engine.get_available_phenotypes()
else:
    phenotypes = engine.get_available_phenotypes(selected_group)

selected_phenotype = st.selectbox(
    "Phenotype",
    ["All"] + phenotypes
)

# Convert back to code
if selected_phenotype != "All":
    phenotype_map = engine.phenotype_reverse_map()
    phenotype_code = phenotype_map.get(selected_phenotype)
else:
    phenotype_code = None

# ---------------- QUERY ----------------
query = st.text_input("Ask a biological question")

if st.button("Run") and query:

    group = None if selected_group == "All" else selected_group
    phenotype = None if selected_phenotype == "All" else phenotype_code

    with st.spinner("Analyzing full gene set..."):
        result = engine.answer(query, group, phenotype)

    st.subheader("Biological Interpretation")
    st.write(result["answer"])

    with st.expander("Gene Table", expanded=False):
        if not result["genes"].empty:
            st.dataframe(result["genes"], use_container_width=True)
        else:
            st.write("No gene data.")
    
    with st.expander("GO Pathways", expanded=False):
        if not result["go"].empty:
            st.dataframe(result["go"], use_container_width=True)
        else:
            st.write("No GO data.")

    with st.expander("Variant Summary Table", expanded=True):
        if not result["variants"].empty:
            st.dataframe(result["variants"], use_container_width=True)
        else:
            st.write("No variant summary data.")
    
    with st.expander("Founder Effect Heatmap", expanded=True):

        coef_df = result["coef"]

        if coef_df.empty:
            st.write("No coefficient data.")
        else:

            founders = [
                "A_J","C57BL_6J","129S1_SvImJ",
                "NOD_ShiLtJ","NZO_HlLtJ","CAST_EiJ",
                "PWK_PhJ","WSB_EiJ"
            ]

            available = [c for c in founders if c in coef_df.columns]

            if not available:
                st.write("No founder columns available.")
            else:
                heatmap_data = (
                    coef_df
                    .pivot_table(
                        index="sample",
                        values=available,
                        aggfunc="mean"
                    )
                )

                if heatmap_data.empty:
                    st.write("No heatmap data after filtering.")
                else:
                    fig, ax = plt.subplots(figsize=(14, 6))

                    sns.heatmap(
                        heatmap_data,
                        cmap="coolwarm",
                        center=0,
                        linewidths=0.2,
                        ax=ax
                    )

                    st.pyplot(fig)
    
    with st.expander("PubMed Literature", expanded=False):
        if result["pubmed_text"]:
            st.write(result["pubmed_text"])
        else:
            st.write("No literature found.")

    with st.expander("Gene Literature Evidence Score", expanded=False):
        if "gene_lit_summary" in result:
            st.json(result["gene_lit_summary"])

    