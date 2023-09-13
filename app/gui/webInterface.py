import os

import streamlit as st

from app.gui.ConformationDrawer2D import ConformationDrawer2D
from app.gui.ConformationDrawer3D import ConformationDrawer3D
from app.src.Controllers.ConformationManager import ConformationManager
from app.src.DataHandlers.JSONProteinIO import JSONProteinIO
from app.src.Models.AminoAcidHP import AminoAcidHP
from app.src.Models.Conformation2D import Conformation2D
from app.src.Models.Coordinates2D import Coordinates2D
from app.src.Models.Coordinates3D import Coordinates3D
from app.src.Models.Lattice2D import Lattice2D
from app.src.Models.Lattice3D import Lattice3D
from app.src.Models.Polarity import Polarity
from app.src.Models.ProteinHP import ProteinHP
from app.src.Models.ProteinModel import ProteinModel
from app.src.Optimizers.REMC import REMC

DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "data")


def main():
    # Page Confiug
    st.set_page_config(
        page_title="REMC for AB-Initio",
        page_icon="🧬",
        layout="wide",  # or "wide" pour utiliser l'ecran tout entier .
    )

    hide_streamlit_style = """
                <style>
                #MainMenu {visibility: hidden;}
                footer {visibility: hidden;}
                </style>
                """
    st.markdown(hide_streamlit_style, unsafe_allow_html=True)

    with st.sidebar:
        st.header("Hyperparameters")

        st.text("")
        max_iterations = st.number_input(
            "Max number iterations",
            value=10,
            step=50,
            min_value=10,
            max_value=1000000,
        )

        search_steps = st.number_input(
            "Number of search steps", value=50, step=50, min_value=1, max_value=10000
        )

        st.text("Temperature")
        min_temp = st.number_input(
            "Min", value=100, step=10, min_value=1, max_value=100000
        )
        max_temp = st.number_input(
            "Max", value=200, step=10, min_value=1, max_value=100000
        )

        nb_replicas = st.number_input(
            "Number of replicas", value=5, step=1, min_value=1, max_value=300
        )

        prob_pull_moves = st.number_input(
            "Probability to use pull moves",
            value=0.0,
            step=0.1,
            min_value=0.0,
            max_value=1.0,
        )

    st.title("Replica Exchange Monte Carlo (REMC) for the AB Initio problem")
    st.divider()

    st.subheader("Protein Data Source Selection")

    option = st.selectbox(
        "Would you like to load the proteins used in the article or use custom ones ?",
        ("", "Article proteins", "Custom one"),
    )

    if option == "Article proteins":
        data_loader = JSONProteinIO(os.path.join(DATA_PATH, "proteins.json"))
        proteins = None
        chosen_protein = ""
        chosen_idx = -1
        dims = None

        try:
            proteins = data_loader.load_proteins(ProteinModel.HYDROPHOBIC_POLAR)
            st.success("Proteins loaded successfully")

        except Exception as e:
            st.error(e)
            st.stop()

        st.write("Please choose a protein to run REMC on in the list below")
        protein_names = [""]
        for protein in proteins:
            protein_names.append(protein.name)

        chosen_protein = st.selectbox("Protein", protein_names)

        if chosen_protein != "":
            chosen_idx = protein_names.index(chosen_protein)
            dims = proteins[chosen_idx - 1].recommended_dimension

            st.write("Please type the dimensions of the lattice to use")
            cols = st.columns(3)

            i = 0
            lattice_dims = []

            for col in cols:
                with col:
                    if i < 2:
                        lattice_dims.append(
                            st.number_input(
                                f"Dimension {i}",
                                value=0,
                                step=1,
                                min_value=0,
                                max_value=500,
                            )
                        )
                    if i == 2 and dims == 3:
                        lattice_dims.append(
                            st.number_input(
                                "Dimension {i}",
                                value=0,
                                step=1,
                                min_value=0,
                                max_value=500,
                            )
                        )

                    i += 1

            run_value = st.button("Run REMC")
            run_algo = True
            for dim in lattice_dims:
                if dim == 0 or not run_value:
                    run_algo = False
                    break

            if run_algo:
                lattice = None
                if dims == 2:
                    lattice = Lattice2D((lattice_dims[0], lattice_dims[1]))
                else:
                    lattice = Lattice3D(
                        (lattice_dims[0], lattice_dims[1], lattice_dims[2])
                    )

                conf_manager = ConformationManager(proteins[chosen_idx - 1])

                try:
                    initial_conformation = conf_manager.create_initial_conformation(
                        lattice.dimensions
                    )

                    st.write(
                        "Theoretical protein energy : ", proteins[chosen_idx - 1].e_star
                    )
                    st.write(
                        "Initial Conformation energy: ",
                        initial_conformation.computed_energy,
                    )

                    MC = REMC(
                        search_steps,
                        nb_replicas,
                        min_temp,
                        max_temp,
                        conf_manager,
                        max_iter=max_iterations,
                        rho=prob_pull_moves,
                    )

                    optimal_conformation = initial_conformation
                    with st.spinner(
                        text="Running REMC, this may take a moment depending or the parameters..."
                    ):
                        optimal_conformation = MC.optimize(
                            initial_conformation, proteins[chosen_idx - 1].e_star
                        )

                    st.write(
                        "Optimal energy found by REMC: ",
                        initial_conformation.compute_energy(),
                    )

                    if dims == 2:
                        drawer = ConformationDrawer2D(
                            optimal_conformation, {"H": 1, "P": 0}
                        )
                        drawer.draw()
                    else:
                        drawer = ConformationDrawer3D(
                            optimal_conformation, {"H": "black", "P": "red"}
                        )
                        drawer.draw()

                except Exception as e:
                    st.error(e)
                    print(e)
                    st.stop()

    elif option == "Custom one":
        custom_sequence = st.text_input(
            "Please enter an amino acid sequence. Ex : HPPPHPPPH (uppercase !)",
            value="",
        )

        dimensions = st.selectbox(
            "In which dimension do you want to represent it ?", ["", "2D", "3D"]
        )

        e_star = st.number_input(
            "Please enter the theoretical energy of the protein",
            value=0,
            step=-1,
            min_value=-500,
            max_value=0,
        )

        j = 0
        seq = []
        for amino in custom_sequence:
            if amino != "H" and amino != "P":
                st.error("Invalid amino acid sequence (H or P only))")
                st.stop()

            pol = None
            if amino == "H":
                pol = Polarity.HYDROPHOBIC
            else:
                pol = Polarity.POLAR

            seq.append(
                AminoAcidHP(
                    j,
                    "",
                    "",
                    pol,
                )
            )

            j += 1

        prots = []
        if dimensions == "2D" and e_star != 0:
            prots.append(ProteinHP("Custom-1", seq, e_star, 2))
        elif dimensions == "3D" and e_star != 0:
            prots.append(ProteinHP("Custom-1", seq, e_star, 3))
        else:
            pass

        if len(prots) != 0:
            st.write("Please type the dimensions of the lattice to use")
            cols = st.columns(3)

            i = 0
            lattice_dims = []

            for col in cols:
                with col:
                    if i < 2:
                        lattice_dims.append(
                            st.number_input(
                                f"Dimension {i}",
                                value=0,
                                step=1,
                                min_value=0,
                                max_value=500,
                            )
                        )
                    if i == 2 and dimensions == "3D":
                        lattice_dims.append(
                            st.number_input(
                                "Dimension {i}",
                                value=0,
                                step=1,
                                min_value=0,
                                max_value=500,
                            )
                        )

                    i += 1

            run_value = st.button("Run REMC")
            run_algo = True
            for dim in lattice_dims:
                if dim == 0 or not run_value:
                    run_algo = False
                    break

            lattice = None
            if run_algo:
                if dimensions == "2D":
                    lattice = Lattice2D((lattice_dims[0], lattice_dims[1]))
                else:
                    lattice = Lattice3D(
                        (lattice_dims[0], lattice_dims[1], lattice_dims[2])
                    )
            for cur_prot in prots:
                conf_manager = ConformationManager(cur_prot)

                try:
                    initial_conformation = conf_manager.create_initial_conformation(
                        lattice.dimensions
                    )

                    st.write("Theoretical protein energy : ", cur_prot.e_star)
                    st.write(
                        "Initial Conformation energy: ",
                        initial_conformation.compute_energy(),
                    )

                    MC = REMC(
                        int(search_steps),
                        int(nb_replicas),
                        int(min_temp),
                        int(max_temp),
                        conf_manager,
                        max_iter=max_iterations,
                        rho=prob_pull_moves,
                    )

                    optimal_conformation = initial_conformation
                    with st.spinner(
                        text="Running REMC, this may take a moment depending on the parameters..."
                    ):
                        optimal_conformation = MC.optimize(
                            initial_conformation, cur_prot.e_star
                        )

                    st.write(
                        "Optimal energy found by REMC: ",
                        initial_conformation.computed_energy,
                    )

                    if dimensions == "2D":
                        drawer = ConformationDrawer2D(
                            optimal_conformation, {"H": 1, "P": 0}
                        )
                        drawer.draw()
                    else:
                        drawer = ConformationDrawer3D(
                            optimal_conformation, {"H": "black", "P": "red"}
                        )
                        drawer.draw()

                except Exception as e:
                    # st.error(e)
                    print(e)
                    st.stop()

    else:
        pass


if __name__ == "__main__":
    main()
