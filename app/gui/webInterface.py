import os

import streamlit as st

from app.src.DataHandlers.JSONProteinIO import JSONProteinIO
from app.src.Models.AminoAcidHP import AminoAcidHP
from app.src.Models.Coordinates2D import Coordinates2D
from app.src.Models.Coordinates3D import Coordinates3D
from app.src.Models.Polarity import Polarity
from app.src.Models.ProteinHP import ProteinHP
from app.src.Models.ProteinModel import ProteinModel

DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "data")


def main():
    st.title("Replica Exchange Monte Carlo (REMC) for the AB Initio problem")

    data_loader = JSONProteinIO(os.path.join(DATA_PATH, "proteins.json"))

    try:
        proteins = data_loader.load_proteins(ProteinModel.HYDROPHOBIC_POLAR)
        for protein in proteins:
            st.text(protein)
        st.success("Proteins loaded successfully")

        # data_loader.filename = os.path.join(DATA_PATH, "new_prots.json")
        # data_loader.save_proteins(proteins, ProteinModel.HYDROPHOBIC_POLAR)

    except Exception as e:
        st.error(e)
        st.stop()


if __name__ == "__main__":
    main()
