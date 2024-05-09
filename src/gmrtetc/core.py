from dataclasses import dataclass


bandinfo = {
    2: {
        "bw": 200,
        "ubw": 50,
        "fmid": 200,
        "flow": 125,
        "fhigh": 250,
        "refgain": 0.33,
        "refgbytsys": 0.0013,
        "coeffs": [
            -0.0274432341259116,
            0.000653316475755705,
            -5.75221466249264e-06,
            2.26066919981535e-08,
            -3.30079139610497e-11,
        ],
    },
    3: {
        "bw": 200,
        "ubw": 120,
        "fmid": 400,
        "flow": 260,
        "fhigh": 500,
        "refgain": 0.33,
        "refgbytsys": 0.0039,
        "coeffs": [
            -3.94269512191062,
            0.0609220731323399,
            -0.000388278427948319,
            1.30788311514484e-06,
            -2.45688427095004e-09,
            2.44226778956594e-12,
            -1.00460621956708e-15,
        ],
    },
    4: {
        "bw": 400,
        "ubw": 200,
        "fmid": 650,
        "flow": 550,
        "fhigh": 850,
        "refgain": 0.33,
        "refgbytsys": 0.00379,
        "coeffs": [
            -60.9659547181797,
            0.52981146165773,
            -0.00191199335049503,
            3.66739380599885e-06,
            -3.94286904604421e-09,
            2.25267804015136e-12,
            -5.34321835013018e-16,
        ],
    },
    5: {
        "bw": 400,
        "ubw": 280,
        "flow": 980,
        "fmid": 1260,
        "fhigh": 1500,
        "refgain": 0.22,
        "refgbytsys": 0.0036,
        "coeffs": [
            -57.79067497317600000,
            0.283183485351112000,
            -0.000576846763506177,
            6.25315209850824e-07,
            -3.8047941517696e-10,
            1.23211866985187e-13,
            -1.65909077237512e-17,
        ],
    },
}


@dataclass
class Observation:
    kind: str
    mode: str
    band: int
    nant: int
    npol: int
    ra: float
    dec: float
    beammode: str

    @property
    def bw(self):
        pass

    @property
    def bwusable(self):
        pass

    @property
    def chanres(self):
        pass

    @property
    def coords(self):
        pass

    def tsky(self):
        pass

    def gbytsys(self):
        pass

    def sefd(self):
        pass

    def sumsefd(self):
        pass

    def rms(self):
        pass
