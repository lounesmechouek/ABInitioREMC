from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Tuple

from .AminoAcidHP import AminoAcidHP
from .Lattice import Lattice
from .ProteinHP import ProteinHP
from .TopoCoordinates import TopoCoordinates


@dataclass(slots=True)
class Conformation(ABC):
    """Abstract class that represents a conformation in a lattice."""

    _protein: ProteinHP  # Protein of the conformation

    _amino_acid_coordinates: dict[
        Tuple[int, ...], AminoAcidHP
    ]  # Coordinates of the amino acids in the conformation

    _lattice: Lattice  # Lattice of the conformation

    _computed_energy: int  # Computed energy of the conformation

    @property
    def protein(self) -> ProteinHP:
        """Getter for the attribute protein of the conformation.

        Returns
        -------
        ProteinHP
            Protein of the conformation.
        """
        return self._protein

    @protein.setter
    def protein(self, protein: ProteinHP) -> None:
        """Setter for the attribute protein of the conformation.

        Parameters
        ----------
        protein : ProteinHP
            Protein of the conformation to be assigned.
        """
        self._protein = protein

    @property
    def amino_acid_coordinates(self) -> dict[Tuple[int, ...], AminoAcidHP]:
        """Getter for the attribute amino_acid_coordinates of the conformation.

        Returns
        -------
        dict[Tuple[int, ...], AminoAcidHP]
            Coordinates of the amino acids in the conformation.
        """
        return self._amino_acid_coordinates

    @amino_acid_coordinates.setter
    def amino_acid_coordinates(
        self, amino_acid_coordinates: dict[Tuple[int, ...], AminoAcidHP]
    ) -> None:
        """Setter for the attribute amino_acid_coordinates of the conformation.

        Parameters
        ----------
        amino_acid_coordinates : dict[Tuple[int, ...], AminoAcidHP]
            Coordinates of the amino acids in the conformation to be assigned.
        """
        self._amino_acid_coordinates = amino_acid_coordinates

    @property
    def lattice(self) -> Lattice:
        """Getter for the attribute lattice of the conformation.

        Returns
        -------
        Lattice
            Lattice of the conformation.
        """
        return self._lattice

    @lattice.setter
    def lattice(self, lattice: Lattice) -> None:
        """Setter for the attribute lattice of the conformation.

        Parameters
        ----------
        lattice : Lattice
            Lattice of the conformation to be assigned.
        """
        self._lattice = lattice

    @property
    def computed_energy(self) -> int:
        """Getter for the attribute computed_energy of the conformation.

        Returns
        -------
        int
            Energy of the conformation.
        """
        return self._computed_energy

    @computed_energy.setter
    def computed_energy(self, computed_energy: int) -> None:
        """Setter for the attribute computed_energy of the conformation.

        Parameters
        ----------
        computed_energy : int
            Energy of the conformation to be assigned.
        """
        self._computed_energy = computed_energy

    @abstractmethod
    def is_valid(self) -> bool:
        """Checks if the conformation is valid.

        Returns
        -------
        bool
            True if the conformation is valid, False otherwise.
        """
        pass

    @abstractmethod
    def compute_energy(self) -> int:
        """Computes the energy of the conformation.

        Returns
        -------
        int
            Energy of the conformation.
        """
        pass

    @abstractmethod
    def get_topological_neighbours(
        self, cell: TopoCoordinates
    ) -> list[TopoCoordinates]:
        """Computes the topological neighbours of a cell.

        Parameters
        ----------
        cell : TopoCoordinates
            Cell.

        Returns
        -------
        list[TopoCoordinates]
            Topological neighbours of the cell.
        """
        pass
