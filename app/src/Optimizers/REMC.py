import copy
import math
import random

from ..Controllers.ConformationManager import ConformationManager
from ..Models.Conformation import Conformation
from .SimpleMonteCarlo import SimpleMonteCarlo


class REMC:
    """Class for Replica Exchange Monte Carlo optimization algorithm in the AB-Initio context."""

    _max_iters: int = 100  # Maximum number of iterations
    _phi: int  # Number of search steps.
    _khi: int  # Number of replicas
    _rho: float = 0.0  # Probability to use pull moves
    _conformation_manager: ConformationManager  # Conformation manager
    _tmin: int  # Minimum temperature
    _tmax: int  # Maximum temperature
    _sampled_temperatures: list[float]  # List of replicas

    def __init__(
        self,
        phi: int,
        khi: int,
        tmin: int,
        tmax: int,
        conf_manager: ConformationManager,
        max_iter: int = 100,
        rho: float = 0.0,
    ) -> None:
        """Constructor for the REMC class.

        Parameters
        ----------
        phi : int
            Number of search steps.
        khi : int
            Number of replicas.
        tmin : float
            Minimum temperature.
        tmax : float
            Maximum temperature.
        conf_manager : ConformationManager
            Conformation manager.
        max_iter: int, optional
            Maximum number of iterations, by default 100
        rho : float, optional
            Probability to use pull moves, by default 0.0
        """
        self._max_iters = max_iter
        self.phi = phi
        self._khi = khi
        if tmin > tmax:
            raise ValueError("tmin must be less than tmax.")
        self._tmin = tmin
        self._tmax = tmax
        self._rho = rho
        self._conformation_manager = conf_manager

        self._sampled_temperatures = random.sample(range(tmin, tmax + 1), khi)

    @property
    def phi(self) -> int:
        """Getter for the attribute phi of the REMC class.

        Returns
        -------
        int
            Number of search steps.
        """
        return self._phi

    @phi.setter
    def phi(self, phi: int) -> None:
        """Setter for the attribute phi of the REMC class.

        Parameters
        ----------
        phi : int
            Number of search steps to be assigned.
        """
        self._phi = phi

    @property
    def rho(self) -> float:
        """Getter for the attribute rho of the REMC class.

        Returns
        -------
        float
            Probability to use pull moves.
        """
        return self._rho

    @rho.setter
    def rho(self, rho: float) -> None:
        """Setter for the attribute rho of the REMC class.

        Parameters
        ----------
        rho : float
            Probability to use pull moves to be assigned.
        """
        self._rho = rho

    @property
    def conformation_manager(self) -> ConformationManager:
        """Getter for the attribute conformation_manager of the REMC class.

        Returns
        -------
        ConformationManager
            Conformation manager of the REMC class.
        """
        return self._conformation_manager

    @conformation_manager.setter
    def conformation_manager(self, conformation_manager: ConformationManager) -> None:
        """Setter for the attribute conformation_manager of the REMC class.

        Parameters
        ----------
        conformation_manager : ConformationManager
            Conformation manager to be assigned.
        """
        self._conformation_manager = conformation_manager

    @property
    def khi(self) -> int:
        """Getter for the attribute khi of the REMC class.

        Returns
        -------
        int
            Number of replicas.
        """
        return self._khi

    @khi.setter
    def khi(self, khi: int) -> None:
        """Setter for the attribute khi of the REMC class.

        Parameters
        ----------
        khi : int
            Number of replicas to be assigned.
        """
        self._khi = khi

    @property
    def tmin(self) -> float:
        """Getter for the attribute tmin of the REMC class.

        Returns
        -------
        float
            Minimum temperature.
        """
        return self._tmin

    @tmin.setter
    def tmin(self, tmin: float) -> None:
        """Setter for the attribute tmin of the REMC class.

        Parameters
        ----------
        tmin : float
            Minimum temperature to be assigned.
        """
        self._tmin = tmin

    @property
    def tmax(self) -> float:
        """Getter for the attribute tmax of the REMC class.

        Returns
        -------
        float
            Maximum temperature.
        """
        return self._tmax

    @tmax.setter
    def tmax(self, tmax: float) -> None:
        """Setter for the attribute tmax of the REMC class.

        Parameters
        ----------
        tmax : float
            Maximum temperature to be assigned.
        """
        self._tmax = tmax

    def optimize(self, conformation: Conformation, e_star: int) -> Conformation:
        """Optimizes a conformation using the REMC algorithm.

        Parameters
        ----------
        conformation : Conformation
            Conformation to be optimized.
        e_star : int
            Optimal (Theoretical) Energy of the protein in the HP-Model.

        Returns
        -------
        Conformation
            Optimized conformation.
        """
        MonteCarlo = SimpleMonteCarlo(self.phi)

        optimal_energy = conformation.compute_energy()
        offset = 0
        optimal_replica = None
        previous_replicas = self.khi * [conformation]
        iters = 1

        while (optimal_energy > e_star) and (iters <= self._max_iters):
            print(12 * "--")
            print("offset: ", offset)
            print("optimal_energy: ", optimal_energy)
            replicas = []
            for k in range(self.khi):
                try:
                    replicas.append(
                        MonteCarlo.optimize(
                            previous_replicas[k],
                            self._sampled_temperatures[k],
                            self.conformation_manager,
                        )
                    )
                except Exception as e:
                    raise e

                print(f"replica {k} : {replicas[k].computed_energy}")
                if replicas[k].computed_energy < optimal_energy:
                    optimal_energy = replicas[k].computed_energy
                    optimal_replica = copy.deepcopy(replicas[k])

            i = offset + 1
            while i + 1 < self.khi:
                j = i + 1
                delta = (
                    1 / self._sampled_temperatures[j]
                    - 1 / self._sampled_temperatures[i]
                ) * (replicas[i].computed_energy - replicas[j].computed_energy)
                if delta <= 0:
                    self._sampled_temperatures[i], self._sampled_temperatures[j] = (
                        self._sampled_temperatures[j],
                        self._sampled_temperatures[i],
                    )
                else:
                    q = random.uniform(0, 1)
                    if q <= math.exp(-delta):
                        self._sampled_temperatures[i], self._sampled_temperatures[j] = (
                            self._sampled_temperatures[j],
                            self._sampled_temperatures[i],
                        )
                i += 2

            offset = 1 - offset
            previous_replicas = copy.deepcopy(replicas)
            iters += 1

        return optimal_replica
