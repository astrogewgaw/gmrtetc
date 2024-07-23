from dataclasses import dataclass
import numpy as np
from scipy.interpolate import interp2d
from astropy.coordinates import SkyCoord
import sys 
import re
import math



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
gain =0.33
def readfile(fname='t408.dat'):
	f=open(fname)
	text=f.read()
	text=text.replace('\r','')	
	text=text.replace('\n','')	
	text=text.replace(' ','')	
	skymap=[]	
	for i in range(0, len(text), 5):
		skymap.append(float(text[i:i+5]))
	lrange=4*np.arange(-45,47)
	brange=np.arange(-90,91)
	skymap=np.array(skymap)
	skymap=skymap.reshape((lrange.shape[0],brange.shape[0]))
	return lrange,brange,skymap
def getTsky408(l,b,lrange,brange,skymap):
	t408=interp2d(x=lrange,y=brange,z=skymap.T,
                      fill_value=7,kind='linear')
	return 10**t408(l,b)[0]   


@dataclass
class Observation:
    band: int  #Observation Band
    cr:float   #channel resolution
    nant: int  #number of antennas 
    npol: int  #number of polarizations
    ra: float  #RA of the transient 
    dec: float  #DEC of the transient 
    f1: float    #lower limit of frequency of transient 
    f2: float   #upper limit of frequney of transient 
    dm: float  #Dispersion measure in pc cm^-3
    pw: float  #pulse width in ms 
    w_scatt: str  #scattering width in ms 
    cd: str  # coherent dispersion- takes 'yes' or 'no'
    beammode: str # PA or IA 
    ubw_1: str
    # frequency you want your value of g_t_sys

    @property
    def bw(self):
        
        if self.band == 1: 
            b_w = bandinfo[1]["bw"]
            #fmid = bandinfo[1]["fmid"]
        elif  self.band == 2: 
            b_w = bandinfo[2]["bw"]
            #fmid = bandinfo[2]["fmid"]
        elif  self.band == 3: 
            b_w = bandinfo[3]["bw"]
            #fmid = bandinfo[3]["fmid"]
        elif  self.band == 4: 
            b_w = bandinfo[4]["bw"]
            #fmid = bandinfo[4]["fmid"]
        # elif self.band >= 5: 
        #     print(f"\n Input correct value")
        #     breakpoint

        return b_w
    


    @property
    def bwusable(self):

        if self.band == 1: 
             ubw = bandinfo[1]["ubw"]
        elif self.band == 2: 
             ubw = bandinfo[2]["ubw"]
        elif self.band == 3:
             ubw = bandinfo[3]["ubw"]
        elif self.band == 4:
             ubw = bandinfo[4]["ubw"]
        return ubw
    
    
    def fmid(self):

        if self.band == 1: 
             fmid = bandinfo[1]["fmid"]
        elif self.band == 2: 
             fmid = bandinfo[2]["fmid"]
        elif self.band == 3:
             fmid = bandinfo[3]["fmid"]
        elif self.band == 4:
             fmid = bandinfo[4]["fmid"]
        return fmid
        

    # @property
    # def chanres(self):
    #     pass

    
    def coords(self):
        tt = self.ra.split()
        sec = 0
        if tt[0]:
            sec += int(tt[0][:-1]) * 3600  # Removing 'h' and converting hours to seconds
        if tt[1]:
            sec += int(tt[1][:-1]) * 60  # Removing 'm' and converting minutes to seconds
        if tt[2]:
            sec += float(tt[2][:-1])  # Removing 's' to get seconds
        RA_deg = sec / 240.0  # sec to deg
        RA_rad = RA_deg / 57.29577951308232 # deg to rad
        tt = re.split(r'[d\'"\s]+', self.dec.strip())
        sign = "+"
        if tt[0] == "-":
            sign = "-"
        sec = 0

        if sign == "+":
            if tt[0]:
                sec += int(tt[0]) * 3600
            if tt[1]:
                sec += int(tt[1]) * 60
            if tt[2]:
                sec += float(tt[2])
        elif sign == "-":
            if tt[0]:
                sec += int(tt[0]) * 3600
            if tt[1]:
                sec -= int(tt[1]) * 60
            if tt[2]:
                sec -= float(tt[2])

        DEC_deg = sec / 3600  # arcsec to deg
        DEC_rad = DEC_deg / 57.29577951308232  # deg to rad
        alpha = RA_rad
        delta = DEC_rad
        alpha0 = 192.8595/57.29577951308232
        delta0=  27.1284/57.29577951308232
        longi0= 122.9320/57.29577951308232

        rhs= math.cos(delta) * math.sin(alpha - alpha0)
        rhs= rhs / ( (math.sin(delta) * math.cos(delta0)) - (math.cos(delta) * math.sin(delta0) * math.cos(alpha - alpha0)) )
        longi =  longi0 -  math.atan(rhs)
        longi = longi * 57.29577951308232 #rad to degree

        lhs= math.sin(delta) * math.sin(delta0) + math.cos(delta) * math.cos(delta0) * math.cos(alpha - alpha0)
        lati = math.asin(lhs)
        lati = lati * 57.29577951308232 #rad to degree
        l = longi
        b = lati
        

        return l,b
        

    def tsky(self):
         l,b = Observation.coords(self)
         nu = Observation.fmid(self)
        
         lrange,brange,skymap=readfile()
         Tsky408=getTsky408(l,b,lrange,brange,skymap)
         if(Tsky408==1e7):
            return -1
         Tsky=(408.0/nu)**2.55*Tsky408
         return Tsky #Kelvin

        
    
    def gbytsys(self,freqx):
         fmid = Observation.fmid(self)
         Tsky = Observation.tsky(self)
         band_number = self.band
         t_def = 22*( 408 / freqx) ** (2.55)
         t_def = int(t_def) #convert to integer
         Tsky = Tsky * (fmid/freqx)**(2.55)

         if band_number == 2:
            ref_gby_tsys = 0.0013
            gain = 0.33
            a = -0.0274432341259116
            b = 0.000653316475755705
            c = -5.75221466249264e-06
            d = 2.26066919981535e-08
            e = -3.30079139610497e-11
            x = freqx
            gby_tsys =  a*1 + b*x + c*x*x + d*x*x*x + e*x*x*x*x
            gby_tsys_f = gain/(gain/gby_tsys - t_def + Tsky)
            ref_gby_tsys = gain/(gain/ref_gby_tsys - t_def + Tsky)
            if gby_tsys_f > ref_gby_tsys*0.5:
                return (1.0/gby_tsys_f)
            else:
                return 0.0
         if self.band == 3:
            ref_gby_tsys =0.0039
            gain = 0.33
            a = -3.94269512191062
            b = 0.0609220731323399
            c = -0.000388278427948319
            d = 1.30788311514484e-06
            e = -2.45688427095004e-09
            p = 2.44226778956594e-12
            q = -1.00460621956708e-15
            x = freqx
            gby_tsys = a*1 + b*x + c*x*x + d*x*x*x + e*x*x*x*x + p*x*x*x*x*x + q*x*x*x*x*x*x
            gby_tsys_f = gain/(gain/gby_tsys - t_def + Tsky)
            ref_gby_tsys = gain/(gain/ref_gby_tsys - t_def + Tsky)
            if gby_tsys_f > ref_gby_tsys*0.5:
                return (1.0/gby_tsys_f)
            else:
                return 0.0
         if band_number == 4:
            ref_gby_tsys = 0.00379
            gain = 0.33
            a = -60.9659547181797
            b = 0.52981146165773
            c = -0.00191199335049503
            d = 3.66739380599885e-06
            e = -3.94286904604421e-09
            p = 2.25267804015136e-12
            q = -5.34321835013018e-16
            x = freqx
            gby_tsys= a*1 + b*x + c*x*x + d*x*x*x + e*x*x*x*x + p*x*x*x*x*x + q*x*x*x*x*x*x
            gby_tsys_f = gain/(gain/gby_tsys - t_def + Tsky)
            ref_gby_tsys = gain/(gain/ref_gby_tsys - t_def + Tsky)
            if gby_tsys_f > ref_gby_tsys*0.5:
                return (1.0/gby_tsys_f)
            else:
                return 0.0
        
         return gby_tsys_f
                        

    def sumsefd(self):

        sumsefd = 0.0
        f1 = self.f1
        f2 = self.f2
        ubw = self.bwusable
        fmid = Observation.fmid(self)


        #If you don't know the frquency range of your transient, then put the lower and upper limit both as 0
    

        if f1 == 0: 
            f1 = (fmid - ubw/2) +0.5
        else: 
            f1 = self.f1
        
        if f2 == 0: 
            f2 = (fmid + ubw/2)
        else: 
            f2 = self.f2

        
        for i in np.arange(f1,f2):
            gby_tsys_f = self.gbytsys(i)
            sumsefd = sumsefd + gby_tsys_f*gby_tsys_f
        return sumsefd  #returns the sumsefd for a range 

        

    def rms(self):
        dm = self.dm
        fmid = Observation.fmid(self)
        cr = self.cr
        pw = self.pw
        nant = self.nant
        npol = self.npol
        
        beammode = self.beammode
        w_scatt = self.w_scatt
        w_dm = self.cd

        if  self.ubw_1 == 'nil': 
            ubw_1 = float(self.bwusable)
        else: 
            ubw_1 = self.ubw_1

       
        
        #w_scatt - insert zero if you don't have any value and put value in ms if you have value 
        if w_scatt == 'nil': 
            w_scatt=math.pow(10,(-6.46 + 0.154 * math.log10(dm) + 1.07 * math.pow(math.log10(dm),2) - 3.86 * math.log10(fmid/1000)))
        else: 
            w_scatt = float(self.w_scatt)
        w_scatt_f = w_scatt/1000
    
        #w_dm
        if w_dm == 'yes':
            w_dm = 0
        elif w_dm == 'no': 
            w_dm =  8.3 * math.pow(10,6) * dm * cr/math.pow(fmid,3)
        w_dm_f = w_dm/1000

        pw = pw/1000
        w_eff = math.sqrt( pw*pw + w_dm_f*w_dm_f + w_scatt_f*w_scatt_f )
        sumsefd = Observation.sumsefd(self)
        

        if beammode == "PA": 
            rms = math.sqrt( sumsefd  /( ubw_1 * nant * nant * npol * math.pow(10,6) * w_eff)) 
        if beammode == "IA":
            rms = math.sqrt( sumsefd  /( ubw_1 * nant * npol * math.pow(10,6) * w_eff ))
        if beammode == "PA-IA":
            rms = math.sqrt( sumsefd  /( ubw_1 * (nant-1) * nant * npol * math.pow(10,6) * w_eff ))
        rms = rms*1000
        
        print(f"The rms for the transient is {rms}mJy")
        
     #in mJy (milli-Jansky)



# bandno= Observation(int(input('What is the band_number?')), float(input('channel resolution:')), int(input("number of antennas:")), int(input("number of polarization used:")),  input("what is the RA?"), input("\n what is the dec?"), int(input("Put lower limit value f1 (put 0 if you want the dictionary value of f1)")), int(input("Put upper limit value f2 (put 0 if you want the dictionary value of f2)")), float(input("Dispersion measure:")),float(input("Pulse Width in ms:")), str(input("w_sctt in ms")), str(input("Coherent Dispersion (yes/no)")), str(input("Beammode (IA/PA/PA-IA):")), str(input('usable bw(put nil if you want dictionary value)')))
# h= bandno.rms()
# print(h)

        

        

