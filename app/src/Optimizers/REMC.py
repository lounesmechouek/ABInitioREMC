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
        self._phi = phi
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
    def tmin(self) -> int:
        """Getter for the attribute tmin of the REMC class.

        Returns
        -------
        int
            Minimum temperature.
        """
        return self._tmin

    @tmin.setter
    def tmin(self, tmin: int) -> None:
        """Setter for the attribute tmin of the REMC class.

        Parameters
        ----------
        tmin : int
            Minimum temperature to be assigned.
        """
        self._tmin = tmin

    @property
    def tmax(self) -> int:
        """Getter for the attribute tmax of the REMC class.

        Returns
        -------
        int
            Maximum temperature.
        """
        return self._tmax

    @tmax.setter
    def tmax(self, tmax: int) -> None:
        """Setter for the attribute tmax of the REMC class.

        Parameters
        ----------
        tmax : int
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
        MonteCarlo = SimpleMonteCarlo(self._phi)

        optimal_energy = conformation.compute_energy()
        optimal_replica = copy.deepcopy(conformation)
        replicas = self._khi * [copy.deepcopy(conformation)]

        offset = 0
        iters = 1

        entered = 0

        print(12 * "####")
        print(f"=> Initial energy : {str(optimal_energy)}")
        print(f"=> Initial coords : {str(optimal_replica.amino_acid_coordinates)}")
        print(f"=> Initial temperatures : {str(self._sampled_temperatures)}")

        while (optimal_energy > e_star) and (iters <= self._max_iters):
            print(f"******REMC : ITERATION {iters}/{self._max_iters}*******")
            for k in range(self._khi):
                # We optimise the replicas
                replicas[k] = MonteCarlo.optimize(
                    replicas[k],
                    self._sampled_temperatures[k],
                    self._conformation_manager,
                )

                if replicas[k].computed_energy < optimal_energy:
                    entered += 1
                    optimal_energy = replicas[k].computed_energy
                    optimal_replica = copy.deepcopy(replicas[k])
                    print(f"New optimal energy : {str(optimal_energy)} !")

            i = offset + 1
            while i + 1 < self._khi:
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
                    print(
                        f"=> Temperature swipe made : {str(self._sampled_temperatures)}"
                    )
                else:
                    q = random.uniform(0, 1)
                    if q <= math.exp(-delta):
                        self._sampled_temperatures[i], self._sampled_temperatures[j] = (
                            self._sampled_temperatures[j],
                            self._sampled_temperatures[i],
                        )
                        print(
                            f"=> Temperature swipe made : {str(self._sampled_temperatures)}"
                        )
                i += 2

            iters += 1
            offset = 1 - offset

        print(f"Optimized {entered} times")
        print(f"New Energy : {str(optimal_energy)}")
        print(f"New coords : {str(optimal_replica.amino_acid_coordinates)}")
        return optimal_replica
