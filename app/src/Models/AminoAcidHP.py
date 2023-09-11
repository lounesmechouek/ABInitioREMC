from dataclasses import dataclass

from .AminoAcid import AminoAcid
from .Polarity import Polarity


@dataclass(slots=True)
class AminoAcidHP(AminoAcid):
    """AminoAcidHP is a class that represents an amino acid in the HP-Model."""

    _polarity: Polarity  # Polarity of the amino acid (Polar or Hydrophobic)

    @property
    def polarity(self) -> Polarity:
        """Getter for the attribute polarity of the amino acid.

        Returns
        -------
        Polarity
            Polarity of the amino acid.
        """
        return self._polarity

    @polarity.setter
    def polarity(self, polarity: Polarity) -> None:
        """Setter for the attribute polarity of the amino acid.

        Parameters
        ----------
        polarity : Polarity
            Polarity of the amino acid to be assigned.
        """
        self._polarity = polarity

    def acide_model(self) -> str:
        return "Hydrophibic-Polar"
