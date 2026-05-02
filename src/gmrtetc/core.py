from pathlib import Path
from typing import cast, Literal
from dataclasses import dataclass
from functools import cached_property

import numpy as np
from beartype import beartype
from autoregistry import Registry
from numpy.polynomial import Polynomial
from scipy.interpolate import RegularGridInterpolator as RGI
from astropy.coordinates import SkyCoord, Latitude, Longitude


class ETCError(Exception):
    pass


@beartype
@dataclass
class Observation(Registry, suffix="Observation"):
    band: int
    nant: int
    npol: int
    rastr: str
    decstr: str

    @property
    def info(self):
        try:
            return {
                2: {
                    "bw": 200.0,
                    "fc": 200.0,
                    "fl": 125.0,
                    "fh": 250.0,
                    "bwuse": 50.0,
                    "refgain": 0.33,
                    "confusion": 86.2,
                    "refgbytsys": 0.0013,
                    "sensitivity": Polynomial(
                        [
                            -0.0274432341259116,
                            0.000653316475755705,
                            -5.75221466249264e-06,
                            2.26066919981535e-08,
                            -3.30079139610497e-11,
                        ]
                    ),
                },
                3: {
                    "bw": 200.0,
                    "fc": 400.0,
                    "fl": 260.0,
                    "fh": 500.0,
                    "bwuse": 120.0,
                    "refgain": 0.33,
                    "confusion": 3.8,
                    "refgbytsys": 0.0039,
                    "sensitivity": Polynomial(
                        [
                            -3.94269512191062,
                            0.0609220731323399,
                            -0.000388278427948319,
                            1.30788311514484e-06,
                            -2.45688427095004e-09,
                            2.44226778956594e-12,
                            -1.00460621956708e-15,
                        ]
                    ),
                },
                4: {
                    "bw": 400.0,
                    "fc": 650.0,
                    "fl": 550.0,
                    "fh": 850.0,
                    "bwuse": 200.0,
                    "refgain": 0.33,
                    "confusion": 0.4,
                    "refgbytsys": 0.00379,
                    "sensitivity": Polynomial(
                        [
                            -60.9659547181797,
                            0.52981146165773,
                            -0.00191199335049503,
                            3.66739380599885e-06,
                            -3.94286904604421e-09,
                            2.25267804015136e-12,
                            -5.34321835013018e-16,
                        ]
                    ),
                },
                5: {
                    "bw": 400.0,
                    "fl": 980.0,
                    "fc": 1260.0,
                    "fh": 1500.0,
                    "bwuse": 280.0,
                    "refgain": 0.22,
                    "confusion": 0.04,
                    "refgbytsys": 0.0036,
                    "sensitivity": Polynomial(
                        [
                            -57.79067497317600000,
                            0.283183485351112000,
                            -0.000576846763506177,
                            6.25315209850824e-07,
                            -3.8047941517696e-10,
                            1.23211866985187e-13,
                            -1.65909077237512e-17,
                        ]
                    ),
                },
            }[self.band]
        except KeyError:
            raise ETCError("INVALID BAND. ABORT.")

    @property
    def fl(self) -> float:
        return self.info["fl"]

    @property
    def fh(self) -> float:
        return self.info["fh"]

    @property
    def fc(self) -> float:
        return self.info["fc"]

    @property
    def bw(self) -> float:
        return self.info["bw"]

    @property
    def bwuse(self) -> float:
        return self.info["bwuse"]

    @property
    def refgain(self) -> float:
        return self.info["refgain"]

    @property
    def refgbytsys(self) -> float:
        return self.info["refgbytsys"]

    @property
    def sensitivity(self) -> Polynomial:
        return self.info["sensitivity"]

    @property
    def coords(self) -> SkyCoord:
        return SkyCoord(f"{self.rastr} {self.decstr}")

    @property
    def ra(self) -> Longitude:
        return cast(Longitude, self.coords.ra)

    @property
    def dec(self) -> Latitude:
        return cast(Latitude, self.coords.dec)

    @property
    def galactic(self) -> SkyCoord:
        return cast(SkyCoord, self.coords.galactic)

    @property
    def gl(self) -> Longitude:
        return cast(Longitude, self.galactic.l)

    @property
    def gb(self) -> Latitude:
        return cast(Latitude, self.galactic.b)

    @cached_property
    def skymap(self) -> np.ndarray:
        return np.asarray(
            [
                float(
                    (
                        (Path(__file__).parent / "data" / "t408.dat")
                        .read_text()
                        .replace(" ", "")
                        .replace("\r", "")
                        .replace("\n", "")
                    )[i : i + 5]
                )
                for i in range(0, 181 * 92 * 5, 5)
            ]
        ).reshape(92, 181)

    def tsky(self, f: float) -> float:
        return float(
            -1
            if (
                tsky408 := 10
                ** (
                    RGI(
                        (4 * np.arange(-45, 47), np.arange(-90, 91)),
                        self.skymap,
                        fill_value=7,
                        method="linear",
                        bounds_error=False,
                    )
                )((self.gl, self.gb))
            )
            == 1e7
            else (408.0 / f) ** 2.55 * tsky408
        )

    @property
    def uptime(self) -> float:
        dec = self.dec.rad
        lat = 19.1 / 180.0 * np.pi
        elv = 17.0 / 180.0 * np.pi
        return float(
            np.acos(
                (np.sin(elv) - (np.sin(lat) * np.sin(dec)))
                / (np.cos(lat) * np.cos(dec))
            )
            / (2 * np.pi)
            * 24.0
            * 2.0
        )


@beartype
@dataclass
class PulsarObservation(Observation):
    nf: int
    dm: float
    wint: float
    wscat: float = -1.0
    cdmode: bool = False
    beammode: Literal["IA", "PA", "PC"] = "IA"

    def __post_init__(self):
        if self.wscat < 0.0:
            self.wscat = (
                10
                ** (
                    -6.46
                    + 0.154 * np.log10(self.dm)
                    + 1.07 * np.log10(self.dm) ** 2
                    - 3.86 * np.log10(self.fc / 1000)
                )
                * 1e-3
            )

    @property
    def df(self) -> float:
        return self.bw / self.nf

    @property
    def freqs(self) -> np.ndarray:
        return np.linspace(self.fl, self.fh, self.nf)

    @property
    def usefreqs(self) -> np.ndarray:
        return np.linspace(
            self.fc - self.bwuse / 2.0 + 0.5,
            self.fc + self.bwuse / 2.0,
            self.nf,
        )

    @property
    def wdm(self) -> float:
        return 0.0 if self.cdmode else 8.3e3 * self.dm * self.df / (self.fc**3)

    @property
    def weff(self) -> float:
        return np.sqrt(self.wint**2 + self.wdm**2 + self.wscat**2)

    def sefd(self, f: float) -> float:
        tdef = int(22 * (408 / f) ** 2.55)
        tsky = self.tsky(self.fc) * (self.fc / f) ** 2.55
        refgbytsys = self.refgain / (self.refgain / self.refgbytsys - tdef + tsky)
        gbytsys = self.refgain / (self.refgain / self.sensitivity(f) - tdef + tsky)
        return 1.0 / gbytsys if gbytsys > refgbytsys * 0.5 else 0.0

    @property
    def sumsefd(self) -> float:
        return float(np.sum(np.vectorize(self.sefd)(self.usefreqs) ** 2))


@beartype
@dataclass
class SinglePulseObservation(PulsarObservation):
    def rms(self) -> float:
        try:
            return float(
                np.sqrt(
                    self.sumsefd
                    / self.nf
                    / (self.npol * (self.bwuse * 1e6) * self.weff)
                    / {
                        "IA": self.nant,
                        "PA": self.nant**2,
                        "PC": self.nant * (self.nant - 1),
                    }[self.beammode]
                )
            )
        except KeyError:
            raise ETCError("INVALID BEAM MODE. ABORT.")


@beartype
@dataclass
class FoldedProfileObservation(PulsarObservation):
    tobs: float = 0.0
    period: float = 0.0

    @property
    def duty(self) -> float:
        return self.weff / self.period

    def rms(self) -> float:
        if (self.tobs > 0.0) and (self.period > 0.0):
            try:
                return float(
                    np.sqrt(
                        self.sumsefd
                        / self.nf
                        / (self.npol * (self.bwuse * 1e6) * self.tobs)
                        / {
                            "IA": self.nant,
                            "PA": self.nant**2,
                            "PC": self.nant * (self.nant - 1),
                        }[self.beammode]
                        * (self.duty / (1 - self.duty))
                    )
                )
            except KeyError:
                raise ETCError("INVALID BEAM MODE. ABORT.")
        raise ETCError("PERIOD AND TOBS SHOULD BE > 0. ABORT.")


@beartype
@dataclass
class ContinuumObservation(Observation):
    pass


@beartype
@dataclass
class LineObservation(Observation):
    pass
