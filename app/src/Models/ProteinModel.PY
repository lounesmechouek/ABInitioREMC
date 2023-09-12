from enum import Enum


class ProteinModel(Enum):
    """Protein models supported by the application."""

    HYDROPHOBIC_POLAR = 0  # Polar amino acid (P)
    OTHER = 1  # Hydrophobic amino acid (H)

    def __str__(self) -> str:
        """Returns a string representation of the model of the protein.

        Returns
        -------
        str
            String representation of the model of the protein.
        """
        return self.name
