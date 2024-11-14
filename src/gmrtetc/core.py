import numpy as np
from typing import Self
from pathlib import Path
from dataclasses import dataclass
from numpy.polynomial import Polynomial
from scipy.interpolate import RegularGridInterpolator as RGI

BANDINFO = {
    2: {
        "bw": 200,
        "fmid": 200,
        "flow": 125,
        "fhigh": 250,
        "bwusable": 50,
        "refgain": 0.33,
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
        "bw": 200,
        "fmid": 400,
        "flow": 260,
        "fhigh": 500,
        "bwusable": 120,
        "refgain": 0.33,
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
        "bw": 400,
        "fmid": 650,
        "flow": 550,
        "fhigh": 850,
        "bwusable": 200,
        "refgain": 0.33,
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
        "bw": 400,
        "flow": 980,
        "fmid": 1260,
        "fhigh": 1500,
        "bwusable": 280,
        "refgain": 0.22,
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
}

DATADIR = (Path(__file__).parent / "data").resolve()


@dataclass
class Observation:
    band: int
    nant: int
    npol: int
    ra: float
    dec: float
    flow: float
    fmid: float
    fhigh: float
    bwusable: float

    @classmethod
    def new(
        cls,
        band: int,
        nant: int,
        npol: int,
        ra: float,
        dec: float,
    ) -> Self:
        if band not in [2, 3, 4, 5]:
            msg = f"Band {band} is invalid. Only Bands 2, 3, 4, and 5 available at the GMRT."
            ValueError(msg)

        if nant > 30:
            msg = f"nant = {nant} > 30. There are only 30 antennas at the GMRT."
            raise ValueError(msg)

        if npol not in [1, 2]:
            msg = "Invalid number of polarisations. Can be only 1 or 2."
            raise ValueError(msg)

        return cls(
            band,
            nant,
            npol,
            ra,
            dec,
            BANDINFO[band]["flow"],
            BANDINFO[band]["fmid"],
            BANDINFO[band]["fhigh"],
            BANDINFO[band]["bwusable"],
        )

    @property
    def galactic(self):
        alpha = np.deg2rad(self.ra)
        delta = np.deg2rad(self.dec)

        gl0 = np.deg2rad(122.9320)
        delta0 = np.deg2rad(27.1284)
        alpha0 = np.deg2rad(192.8595)

        gl = np.rad2deg(
            gl0
            - (
                (np.cos(delta) * np.sin(alpha - alpha0))
                / (
                    (np.sin(delta) * np.cos(delta0))
                    - (np.cos(delta) * np.sin(delta0) * np.cos(alpha - alpha0))
                )
            )
        )

        gb = np.rad2deg(
            np.arcsin(
                np.sin(delta) * np.sin(delta0)
                + np.cos(delta) * np.cos(delta0) * np.cos(alpha - alpha0)
            )
        )

        return gl, gb

    def calctsky(self, freq: float) -> float:
        gl, gb = self.galactic

        with open(f"{DATADIR}/t408.dat") as f:
            text = f.read()
            text = text.replace(" ", "")
            text = text.replace("\r", "")
            text = text.replace("\n", "")

            brange = np.arange(-90, 91)
            lrange = 4 * np.arange(-45, 47)
            skymap = np.array([float(text[i : i + 5]) for i in range(0, len(text), 5)])
            skymap = skymap.reshape((lrange.shape[0], brange.shape[0]))

        t408 = RGI(
            (lrange, brange),
            skymap,
            fill_value=7,
            method="linear",
            bounds_error=False,
        )

        tsky408 = 10 ** t408((gl, gb))
        tsky = -1 if tsky408 == 1e7 else (408.0 / freq) ** 2.55 * tsky408
        return tsky

    def calcgbytsys(self, freq: float) -> float:
        tdef = int(22 * (408 / freq) ** 2.55)
        tsky = self.calctsky(self.fmid) * (self.fmid / freq) ** 2.55

        refgain = BANDINFO[self.band]["refgain"]
        refgbytsys = BANDINFO[self.band]["refgbytsys"]
        gbytsys = BANDINFO[self.band]["sensitivity"](freq)

        gbytsys = refgain / (refgain / gbytsys - tdef + tsky)
        refgbytsys = refgain / (refgain / refgbytsys - tdef + tsky)
        return 1.0 / gbytsys if gbytsys > refgbytsys * 0.5 else 0.0


@dataclass
class ContinuumObservation(Observation):
    def calcsumsefd(self):
        pass

    def calcrms(self):
        pass


@dataclass
class BeamformedObservation(Observation):
    nf: int
    mode: str
    dm: float
    wint: float
    cdmode: bool = False
    wdm: float | None = None
    wscatt: float | None = None

    @property
    def df(self) -> float:
        return self.bwusable / self.nf

    @property
    def freqs(self) -> np.ndarray:
        return np.linspace(self.flow, self.fhigh, self.nf)

    def calcsumsefd(self) -> float:
        sumsefd = 0.0
        for freq in self.freqs:
            sumsefd += self.calcgbytsys(freq)
        return sumsefd

    def calcrms(self) -> float:
        rms = self.calcsumsefd() / (
            self.npol
            * (self.bwusable * 1e6)
            * np.sqrt(
                self.wint**2
                + (0.0 if self.cdmode else 8.3e6 * self.dm * self.df / (self.fmid**3))
                ** 2
                + (
                    10
                    ** (
                        -6.46
                        + 0.154 * np.log10(self.dm)
                        + 1.07 * np.log10(self.dm) ** 2
                        - 3.86 * np.log10(self.fmid / 1000)
                    )
                    if self.wscatt is None
                    else self.wscatt
                )
                ** 2
            )
        )

        if self.mode == "IA":
            rms /= self.nant
        elif self.mode == "PA":
            rms /= self.nant**2
        elif self.mode == "PC":
            rms /= self.nant * (self.nant - 1)
        else:
            raise ValueError("This mode is not available at the GMRT.")
        rms = np.sqrt(rms)

        return rms
