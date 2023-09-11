import streamlit as st

from app.src.Models.AminoAcidHP import AminoAcidHP
from app.src.Models.Coordinates2D import Coordinates2D
from app.src.Models.Coordinates3D import Coordinates3D
from app.src.Models.Polarity import Polarity
from app.src.Models.ProteinHP import ProteinHP


def main():
    st.title("Replica Exchange Monte Carlo (REMC) for the AB Initio problem")

    aa_list = [
        AminoAcidHP("Alanine", "Ala", Polarity.HYDROPHOBIC),
        AminoAcidHP("Histidine", "His", Polarity.POLAR),
    ]
    protein = ProteinHP("DD-1", aa_list, -20, 2)

    point2D = Coordinates2D((1, 2))
    point3D = Coordinates3D((1, 2, 3))

    st.write(point2D.coordinates)
    st.write(point3D.coordinates)


if __name__ == "__main__":
    main()
