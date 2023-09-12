from abc import ABC, abstractmethod

from ..Models.Protein import Protein
from ..Models.ProteinModel import ProteinModel


class ProteinIO(ABC):
    """This is an interface used for reading and writing protein data."""

    @abstractmethod
    def load_proteins(self, protein_model: ProteinModel) -> list[Protein]:
        """Loads proteins from a data source.

        Returns
        -------
        list[Protein]
            List of proteins loaded from a data source.
        """
        pass

    @abstractmethod
    def save_proteins(
        self, proteins: list[Protein], protein_model: ProteinModel
    ) -> None:
        """Saves proteins to a destination.

        Parameters
        ----------
        proteins : list[Protein]
            List of proteins to be saved to a destination.
        """
        pass
