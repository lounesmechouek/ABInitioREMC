import copy
import random
from typing import Tuple

from ..Models.Conformation import Conformation
from ..Models.Conformation2D import Conformation2D
from ..Models.Conformation3D import Conformation3D
from ..Models.Coordinates2D import Coordinates2D
from ..Models.Coordinates3D import Coordinates3D
from ..Models.Lattice2D import Lattice2D
from ..Models.Lattice3D import Lattice3D
from ..Models.ProteinHP import ProteinHP


class ConformationManager:
    """Class that manages the conformations of a protein."""

    _protein: ProteinHP  # Protein of the conformations

    _conformations: list[Conformation]  # Conformations of the protein

    def __init__(self, protein: ProteinHP) -> None:
        """Constructor for the ConformationManager class.

        Parameters
        ----------
        protein : ProteinHP
            Protein of the conformations.
        """
        self._protein = protein
        self._conformations = []

    @property
    def protein(self) -> ProteinHP:
        """Getter for the attribute protein of the ConformationManager.

        Returns
        -------
        ProteinHP
            Protein of the conformations.
        """
        return self._protein

    @protein.setter
    def protein(self, protein: ProteinHP) -> None:
        """Setter for the attribute protein of the ConformationManager.

        Parameters
        ----------
        protein : ProteinHP
            Protein of the conformations to be assigned.
        """
        self._protein = protein

    @property
    def conformations(self) -> list[Conformation]:
        """Getter for the attribute conformations of the ConformationManager.

        Returns
        -------
        list[Conformation]
            Conformations of the protein.
        """
        return self._conformations

    @conformations.setter
    def conformations(self, conformations: list[Conformation]) -> None:
        """Setter for the attribute conformations of the ConformationManager.

        Parameters
        ----------
        conformations : list[Conformation]
            Conformations of the protein to be assigned.
        """
        self._conformations = conformations

    def create_initial_conformation(self, lattice_dims=Tuple[int, ...]) -> Conformation:
        """Creates the initial conformation of the protein and adds it to the list of conformations.

        Parameters
        ----------
        lattice_dims : Tuple[int, ...]
            Dimensions of the lattice.

        Returns
        -------
        Conformation
            Initial conformation of the protein.
        """
        if len(lattice_dims) == 2:
            lattice = Lattice2D(lattice_dims)
        elif len(lattice_dims) == 3:
            lattice = Lattice3D(lattice_dims)
        else:
            raise ValueError("The lattice dimensions must be 2 or 3.")

        i = 0
        search_valid_conformation = True
        dict_coords = {}
        last_amino_coords = None

        while search_valid_conformation:
            for amino_acid in self._protein.sequence:
                if i == 0:
                    pos = []
                    # We sample a random position for the first amino acid
                    for axe in range(len(lattice_dims)):
                        pos.append(random.randint(0, lattice_dims[axe] - 1))
                    dict_coords[tuple(pos)] = amino_acid

                    if len(lattice_dims) == 2:
                        lattice.set_cell_value(Coordinates2D(tuple(pos)), True)
                    elif len(lattice_dims) == 3:
                        lattice.set_cell_value(Coordinates3D(tuple(pos)), True)
                    else:
                        raise ValueError("The lattice dimensions must be 2 or 3.")

                    if len(pos) == 2:
                        last_amino_coords = Coordinates2D(tuple(pos))
                    elif len(pos) == 3:
                        last_amino_coords = Coordinates3D(tuple(pos))
                    else:
                        raise ValueError("The lattice dimensions must be 2 or 3.")

                else:
                    position = None
                    # We choose an adjacent position for the next amino acid that is not occupied
                    try:
                        position = lattice.get_random_adjacent_cell(
                            last_amino_coords, list(dict_coords.keys())
                        )

                        dict_coords[position] = amino_acid
                        if len(lattice_dims) == 2:
                            lattice.set_cell_value(Coordinates2D(position), True)
                            last_amino_coords = Coordinates2D(position)
                        elif len(lattice_dims) == 3:
                            lattice.set_cell_value(Coordinates3D(position), True)
                            last_amino_coords = Coordinates3D(position)
                        else:
                            raise ValueError("The lattice dimensions must be 2 or 3.")

                    except Exception:
                        # We could not find any adjacent position for the amino acid
                        dict_coords = {}
                        i = 0
                        last_amino_coords = None
                        lattice.reset_lattice()
                        break
                i += 1

            # We create a conformation and we check if it is valid
            if len(dict_coords) == len(self._protein.sequence):
                if len(lattice_dims) == 2:
                    conformation = Conformation2D(self._protein, lattice, dict_coords)

                elif len(lattice_dims) == 3:
                    conformation = Conformation3D(self._protein, lattice, dict_coords)
                else:
                    raise ValueError("The lattice dimensions must be 2 or 3.")

                if conformation.is_valid():
                    search_valid_conformation = False
                    self._conformations.append(conformation)
                    return conformation
                else:
                    dict_coords = {}
                    i = 0
                    last_amino_coords = None
                    lattice.reset_lattice()

            else:
                dict_coords = {}
                i = 0
                last_amino_coords = None
                lattice.reset_lattice()

    def compute_vhsd_neighbourhood(
        self, conformation: Conformation
    ) -> list[Conformation]:
        """Computes the VHSd neighbourhood of a conformation.

        Parameters
        ----------
        conformation : Conformation
            The conformation to be used to compute the neighbourhood.

        Returns
        -------
        List[Conformation]
            VHSd neighbourhood of the conformation.
        """
        neighbourhood = []

        dim = None
        if len(conformation.lattice.dimensions) == 2:
            dim = 2
        elif len(conformation.lattice.dimensions) == 3:
            dim = 3
        else:
            raise ValueError("The lattice dimensions must be 2 or 3.")

        for coord in conformation.amino_acid_coordinates:
            new_positions = []
            # We first compute the possible end moves
            # End moves are made on the first and last amino acids of the protein

            if coord in [
                conformation.get_amino_acid_coordinates(
                    conformation.protein.sequence[0]
                ).coordinates,
                conformation.get_amino_acid_coordinates(
                    conformation.protein.sequence[
                        len(conformation.protein.sequence) - 1
                    ]
                ).coordinates,
            ]:
                try:
                    if dim == 2:
                        new_end_positions = conformation.lattice.compute_end_moves(
                            Coordinates2D(coord),
                            conformation.protein,
                            conformation.amino_acid_coordinates,
                        )
                    else:
                        new_end_positions = conformation.lattice.compute_end_moves(
                            Coordinates3D(coord),
                            conformation.protein,
                            conformation.amino_acid_coordinates,
                        )

                    for pos in new_end_positions:
                        new_positions.append(pos)
                except Exception as e:
                    # print(e)
                    # print(
                    #    f"No end moves possible for {coord} in {conformation.protein}"
                    # )
                    pass

            # Then we compute corner moves when possible on the rest of the residues
            if coord not in [
                conformation.get_amino_acid_coordinates(
                    conformation.protein.sequence[0]
                ).coordinates,
                conformation.get_amino_acid_coordinates(
                    conformation.protein.sequence[
                        len(conformation.protein.sequence) - 1
                    ]
                ).coordinates,
            ]:
                try:
                    if dim == 2:
                        new_corner_position = conformation.lattice.compute_corner_moves(
                            Coordinates2D(coord),
                            conformation.protein,
                            conformation.amino_acid_coordinates,
                        )
                    else:
                        new_corner_position = conformation.lattice.compute_corner_moves(
                            Coordinates3D(coord),
                            conformation.protein,
                            conformation.amino_acid_coordinates,
                        )
                    for pos in new_corner_position:
                        new_positions.append(pos)
                except Exception as e:
                    # print(e)
                    # print(
                    #    f"No corner moves possible for {coord} in {conformation.protein}"
                    # )
                    pass

            if len(new_positions) > 0:
                for new_pos in new_positions:
                    new_lattice = copy.deepcopy(conformation.lattice)

                    if dim == 2:
                        new_lattice.set_cell_value(Coordinates2D(coord), False)
                        new_lattice.set_cell_value(Coordinates2D(new_pos), True)
                    else:
                        new_lattice.set_cell_value(Coordinates3D(coord), False)
                        new_lattice.set_cell_value(Coordinates3D(new_pos), True)

                    dict_cells = copy.deepcopy(conformation.amino_acid_coordinates)
                    dict_cells[new_pos] = dict_cells.pop(coord)

                    new_conf = None

                    if dim == 2:
                        new_conf = Conformation2D(
                            conformation.protein, new_lattice, dict_cells
                        )
                    else:
                        new_conf = Conformation3D(
                            conformation.protein, new_lattice, dict_cells
                        )

                    self._conformations.append(new_conf)
                    neighbourhood.append(new_conf)
            else:
                # print(f"/!/ No moves possible for {coord} in {conformation.protein}")
                pass

        return neighbourhood
