from dataclasses import dataclass

from .AminoAcidHP import AminoAcidHP
from .Protein import Protein


@dataclass(slots=True)
class ProteinHP(Protein):
    """ProteinHP is a class that represents a protein in the HP-Model."""

    _e_star: int  # Optimal (Theoretical) Energy of the protein in the HP-Model.
    _recommended_dimension: int  # Recommended dimension of the protein (2D or 3D) in the HP-Model.

    @property
    def e_star(self) -> int:
        """Getter for the attribute e_star of the protein.

        Returns
        -------
        int
            Optimal (Theoretical) Energy of the protein in the HP-Model.
        """
        return self._e_star

    @e_star.setter
    def e_star(self, e_star: int) -> None:
        """Setter for the attribute e_star of the protein.

        Parameters
        ----------
        e_star : int
            Optimal (Theoretical) Energy of the protein in the HP-Model to be assigned.
        """
        self._e_star = e_star

    @property
    def recommended_dimension(self) -> int:
        """Getter for the attribute recommended_dimension of the protein.

        Returns
        -------
        int
            Recommended dimension of the protein (2D or 3D) in the HP-Model.
        """
        return self._recommended_dimension

    @recommended_dimension.setter
    def recommended_dimension(self, recommended_dimension: int) -> None:
        """Setter for the attribute recommended_dimension of the protein.

        Parameters
        ----------
        recommended_dimension : int
            Recommended dimension of the protein (2D or 3D) in the HP-Model to be assigned.
        """
        self._recommended_dimension = recommended_dimension

    def are_neighbours(
        self, amino_acid_1: AminoAcidHP, amino_acid_2: AminoAcidHP
    ) -> bool:
        """Checks if two amino acids are neighbours in the protein.

        Parameters
        ----------
        amino_acid_1 : AminoAcidHP
            First amino acid.
        amino_acid_2 : AminoAcidHP
            Second amino acid.

        Returns
        -------
        bool
            True if the amino acids are neighbours, False otherwise.
        """
        index_amino1 = None
        index_amino2 = None
        for i, residue in enumerate(self.sequence):
            if residue.id == amino_acid_1.id:
                index_amino1 = i
            elif residue.id == amino_acid_2.id:
                index_amino2 = i

            if index_amino1 is not None and index_amino2 is not None:
                break

        if index_amino1 is None or index_amino2 is None:
            return False

        else:
            neighbours = abs(index_amino1 - index_amino2) == 1
            return neighbours

    def protein_model(self) -> str:
        return "Hydrophobic-Polar"

    def __repr__(self) -> str:
        """Returns a string representation of the protein.

        Returns
        -------
        str
            String representation of the protein.
        """
        return f"{''.join([str(amino_acid) for amino_acid in self.sequence])}"

    def __str__(self) -> str:
        """Returns a string representation of the protein.

        Returns
        -------
        str
            String representation of the protein.
        """
        return f"{''.join([str(amino_acid) for amino_acid in self.sequence])}"
