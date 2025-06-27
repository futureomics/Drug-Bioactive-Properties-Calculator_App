import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED

# --- Bioactivity calculation ---
def calculate_properties(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        "SMILES": smiles,
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "RotatableBonds": Descriptors.NumRotatableBonds(mol),
        "QED": QED.qed(mol)
    }

# --- Streamlit app ---
st.title("🧪 Drug Bioactive Properties Calculator")
st.write("Upload a **CSV** file with a column named `smiles` (max 500 molecules).")

uploaded_file = st.file_uploader("Upload CSV File", type=["csv"])

if uploaded_file:
    try:
        df = pd.read_csv(uploaded_file)

        if 'smiles' not in df.columns:
            st.error("CSV must contain a 'smiles' column.")
        else:
            df = df.dropna(subset=['smiles'])
            df = df.head(500)  # Limit to 500 molecules

            with st.spinner("Calculating properties..."):
                results = []
                for smi in df['smiles']:
                    props = calculate_properties(smi)
                    if props:
                        results.append(props)

                if results:
                    result_df = pd.DataFrame(results)
                    st.success(f"Calculated properties for {len(result_df)} molecules.")
                    st.dataframe(result_df)

                    csv = result_df.to_csv(index=False)
                    st.download_button("📥 Download Results as CSV", csv, "drug_properties.csv", "text/csv")

                else:
                    st.warning("No valid SMILES found in the input.")
    except Exception as e:
        st.error(f"Error: {e}")
