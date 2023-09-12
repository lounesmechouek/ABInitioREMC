from dataclasses import dataclass

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
