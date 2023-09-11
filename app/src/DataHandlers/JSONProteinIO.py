import json

from ..Models.Protein import Protein
from ..Models.ProteinHP import ProteinHP
from .ProteinIO import ProteinIO


class JSONProteinIO(ProteinIO):
    """This is an interface used for reading and writing protein data from JSON files."""

    _filename: str

    def __init__(self, filename: str) -> None:
        """Constructor for the JSONProteinIO class.

        Parameters
        ----------
        filename : str
            Name of the JSON file to be read from and written to.
        """
        self._filename = filename

    def load_proteins(self) -> list[Protein]:
        """Loads proteins from a JSON file.

        Returns
        -------
        list[ProteinHP]
            List of proteins loaded from a JSON file.
        """
