from enum import Enum


class Polarity(Enum):
    """Polarity of an amino acid in the HP-Model."""

    POLAR = 0  # Polar amino acid (P)
    HYDROPHOBIC = 1  # Hydrophobic amino acid (H)

    def __str__(self) -> str:
        """Returns a string representation of the polarity of the amino acid.

        Returns
        -------
        str
            String representation of the polarity of the amino acid.
        """
        return self.name
