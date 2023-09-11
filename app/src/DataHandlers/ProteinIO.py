from abc import ABC, abstractmethod

from ..Models.Protein import Protein


class ProteinIO(ABC):
    """This is an interface used for reading and writing protein data."""

    @abstractmethod
    def load_proteins(self) -> list[Protein]:
        """Loads proteins from a data source.

        Returns
        -------
        list[Protein]
            List of proteins loaded from a data source.
        """
        pass

    @abstractmethod
    def save_proteins(self, proteins: list[Protein]) -> None:
        """Saves proteins to a destination.

        Parameters
        ----------
        proteins : list[Protein]
            List of proteins to be saved to a destination.
        """
        pass
