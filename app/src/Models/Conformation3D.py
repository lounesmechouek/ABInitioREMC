from dataclasses import dataclass
from typing import Tuple

from app.src.Models.Coordinates3D import Coordinates3D

from .AminoAcidHP import AminoAcidHP
from .Conformation import Conformation
from .Lattice3D import Lattice3D
from .Polarity import Polarity
from .ProteinHP import ProteinHP


@dataclass(slots=True)
class Conformation3D(Conformation):
    """Class that represents a conformation in a 3D lattice."""

    def __init__(
        self,
        protein: ProteinHP,
        lattice: Lattice3D,
        amino_acid_coordinates: dict[Tuple[int, int, int], AminoAcidHP],
    ) -> None:
        """Constructor for the Conformation2D class.

        Parameters
        ----------
        lattice : Lattice3D
            Lattice of the conformation.
        """
        self._protein = protein
        self._lattice = lattice
        self._computed_energy = 0

        self._amino_acid_coordinates = amino_acid_coordinates
        try:
            for coord in amino_acid_coordinates.keys():
                self._lattice.set_cell_value(Coordinates3D(coord), True)
        except Exception as e:
            raise e

    def get_amino_acid_coordinates(self, amino_acid: AminoAcidHP) -> Coordinates3D:
        """Gets the coordinates of an amino acid in the conformation.

        Parameters
        ----------
        amino_acid : AminoAcidHP
            Amino acid.

        Returns
        -------
        Coordinates3D
            Coordinates of the amino acid.
        """
        for key, value in self._amino_acid_coordinates.items():
            if value.id == amino_acid.id:
                return Coordinates3D(key)

        raise ValueError("Amino acid not found in the conformation.")

    def is_valid(self) -> bool:
        """Checks if the conformation is valid.

        Returns
        -------
        bool
            True if the conformation is valid, False otherwise.
        """
        valid_conformation = True
        # It is enough to check that each amino acid that is adjacent to the next one in the sequence is adjacent in the lattice space.
        for i in range(len(self._protein.sequence)):
            try:
                cur_coords = self.get_amino_acid_coordinates(self._protein.sequence[i])

                if i == 0:
                    next_coords = self.get_amino_acid_coordinates(
                        self._protein.sequence[i + 1]
                    )
                    # Check if the amino acids are adjacent
                    if not self._lattice.are_adjacent(cur_coords, next_coords):
                        valid_conformation = False
                        break
                elif i == len(self._protein.sequence) - 1:
                    prec_coords = self.get_amino_acid_coordinates(
                        self._protein.sequence[i - 1]
                    )
                    # Check if the amino acids are adjacent
                    if not self._lattice.are_adjacent(cur_coords, prec_coords):
                        valid_conformation = False
                        break
                else:
                    prec_coords = self.get_amino_acid_coordinates(
                        self._protein.sequence[i - 1]
                    )
                    next_coords = self.get_amino_acid_coordinates(
                        self._protein.sequence[i + 1]
                    )
                    # Check if the amino acids are adjacent
                    if (not self._lattice.are_adjacent(cur_coords, next_coords)) or (
                        not self._lattice.are_adjacent(cur_coords, prec_coords)
                    ):
                        valid_conformation = False
                        break
            except ValueError:
                valid_conformation = False
                break

        return valid_conformation

    def compute_energy(self) -> int:
        """Computes the energy of the conformation.

        Returns
        -------
        int
            Energy of the conformation.
        """
        energy = 0
        for i in range(len(self._protein.sequence) - 1):
            for j in range(i + 1, len(self._protein.sequence)):
                if abs(j - i) > 1:  # Not connected in the sequence
                    try:
                        cur_coords = self.get_amino_acid_coordinates(
                            self._protein.sequence[i]
                        )
                        next_coords = self.get_amino_acid_coordinates(
                            self._protein.sequence[j]
                        )
                        if (
                            (self._lattice.are_adjacent(cur_coords, next_coords))
                            and (
                                self._protein.sequence[i].polarity
                                == Polarity.HYDROPHOBIC
                            )
                            and (
                                self._protein.sequence[j].polarity
                                == Polarity.HYDROPHOBIC
                            )
                        ):
                            # Topological H neighbours
                            energy += -1
                    except Exception as e:
                        raise e
        self._computed_energy = energy
        return energy

    def get_topological_neighbours(self, cell: Coordinates3D) -> list[Coordinates3D]:
        """Computes the topological neighbours of a cell.

        Parameters
        ----------
        cell : Coordinates2D
            Cell.

        Returns
        -------
        list[Coordinates2D]
            Topological neighbours of the cell.
        """
        pass
