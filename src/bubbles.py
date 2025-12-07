import math
import numpy as np
from . import constants as cc
from .SNR_utils import composit_R, RSB_Badry
class GalctoHelio:
    def __init__(self,
                gal_lbd   = None,
                theta_deg = None,
                helio_xyz = None):
        # initialise coordinates
        if gal_lbd is not None:
            self.l_deg, self.b_deg, self.d = gal_lbd
        else:
            self.l_deg, self.b_deg, self.d = [None]*3
        if helio_xyz is not None:
            self.x, self.y, self.z = helio_xyz 
        else:
            self.x, self.y, self.z = self.helio()
        if theta_deg is not None:
            self.theta_deg = theta_deg    
        else:
            self.theta_deg = None
        pass
    def helio(self):
        l_rad = math.radians(self.l_deg)
        b_rad = math.radians(self.b_deg)
        x = self.d * math.cos(b_rad) * math.cos(l_rad)
        y = self.d * math.cos(b_rad) * math.sin(l_rad)
        z = self.d * math.sin(b_rad)
        
        self.x, self.y, self.z = x, y, z
        return self.x, self.y, self.z
    
class bubble(GalctoHelio):
    def __init__(self,gal_lbd=None, helio_xyz=None, theta_deg=None,
                 xi=1.0,name = "No name",known_age = None, DensityProfile=None):
        super().__init__(gal_lbd=gal_lbd, theta_deg=theta_deg,helio_xyz=helio_xyz)
        self.name= name
        self.xi = xi # ism metallicity 
        self.kind = None 
        
        self.known_age = known_age
        self.computed_age = None
        self.density_profile = DensityProfile
        if self.density_profile is None:
            raise ValueError(f"Density profile is not set for this bubble: {self.name}")
        self.local_density()
        pass
    def data_radius(self):
        
        if self.theta_deg is not None:
            theta_rad = math.radians(self.theta_deg)
        else:
            theta_rad = 0.0
        if self.d is None:
            d  =0.0
        else:
            d = self.d
        diameter  = 2.0 * d * math.tan(theta_rad / 2.0)
        return diameter / 2.0
    @staticmethod
    def zero_search(arr,val):
        tmp_val = np.inf
        for idx,elem in enumerate(arr):
            diff = val-elem
            nozeros = True
            if diff < tmp_val: 
                index = idx 
                element = elem
                tmp_val = diff
                if diff <= 0.0:
                    nozeros=False
                    break
        if nozeros:
            print("No zero was found, returning closest number.")
        return index, element 
    def coord(self):
        if self.x is not None:
            return self.x, self.y,self.z
        else: 
            return self.helio()
    def local_density(self, density_profile=None):
        
        #  density_profile is not None:
        #     self.density_profile = density_profile
            
        if self.density_profile is None:
            raise ValueError(f"Density profile is not set for this bubble: {self.name}")
        x, y, z = self.coord()
        self.local_n = self.density_profile.profile(z*cc.pc)/(cc.mH*cc.mu_mol)
        return self.local_n
    def set_kind(self, kind):
        self.kind = kind
    def is_SB(self):
        return self.kind=="SB"
    def is_SNR(self):
        return not self.is_SB()
        
class SNR(bubble):
    def __init__(self, 
                 gal_lbd=None, theta_deg=None, helio_xyz=None,
                 xi=1, name = "No name", 
                 known_age  = None, DensityProfile=None):
    
        super().__init__(gal_lbd=gal_lbd, theta_deg=theta_deg,helio_xyz=helio_xyz, xi=xi, name=name, known_age=known_age, DensityProfile=DensityProfile)
        super().set_kind("SNR")
        
        pass
    def age(self, cs5=10., tmin=1e1, tmax=5e7,npoints=100000):
        if self.known_age is not None: 
            return self.known_age
        else:
            if self.computed_age is not None:
                return self.computed_age
            else:    
                # only for single SNR case
                n=self.local_density()
                radius_obs = self.data_radius()
                time = 10**np.linspace(np.log10(tmin), np.log10(tmax),npoints)
                radii= [composit_R(tt,n = n,xi = self.xi,cs5 = cs5) for tt in time]
                index, rad =  self.zero_search(radii,radius_obs)
                self.computed_age = time[index]
                return time[index]
    def Nsn_analytic(self):
        return 1
    def __repr__(self) -> str:
        return "%20s | n = %8.3g, age = %10.5g  Myr | R = %8.4g, (x,y,z) = (%8.3g,%8.3g,%8.3g)"%(self.name, self.local_n,self.age(),self.data_radius(),self.x,self.y,self.z) #super().__repr__()
    def __str__(self) -> str:
        return "%20s | n = %8.3g, age = %10.5g  Myr | R = %8.4g, (x,y,z) = (%8.3g,%8.3g,%8.3g)"%(self.name, self.local_n,self.age(),self.data_radius(),self.x,self.y,self.z) #super().__repr__()

class SB(bubble):
    def __init__(self, 
                 gal_lbd=None,theta_deg=None,helio_xyz=None,
                 xi=1, name = "No name",
                 known_age  = None, deltat_SNe = 1e6, DensityProfile=None):
        super().__init__(gal_lbd=gal_lbd, theta_deg=theta_deg, helio_xyz=helio_xyz, xi=xi, name=name,known_age= known_age, DensityProfile=DensityProfile)
        super().set_kind("SB")
        
        self.deltat_SNe = deltat_SNe
        self.progenitors = []
        
        pass
    
    def age(self, theta=0.6, tmin=1e1, tmax=5e7,npoints=100000):    
        if self.known_age is not None: 
            return self.known_age
        else:
            if self.computed_age is not None:
                return self.computed_age
            else:
                # only for single SNR case
                n=self.local_density()
                radius_obs = self.data_radius()
                time = 10**np.linspace(np.log10(tmin), np.log10(tmax),npoints)
                radii= [RSB_Badry(tt,n = n,deltat_SNe=self.deltat_SNe, theta=theta) for tt in time]
                index, rad =  self.zero_search(radii,radius_obs)
                self.computed_age = time[index]
                return time[index] #,abs((radius_obs-rad)/radius_obs)
    def Nsn_analytic(self):
        self.Nsn = self.computed_age/(self.deltat_SNe)
        return self.Nsn
    def add_SNe(self,snr):
        self.progenitors+=[snr]
        return
    def add_random_SNe(self,radius=10):
        if self.Nsn is None: self.Nsn_analytic()
        xn = np.random.normal(loc=self.x,scale=radius,size=round(self.Nsn))
        yn = np.random.normal(loc=self.y,scale=radius,size=round(self.Nsn))
        zn = np.random.normal(loc=self.z,scale=radius,size=round(self.Nsn))
        trand = np.random.uniform(-1.,1.)/10.0
        for i in range(0,round(self.Nsn)):
            avg_age = self.age()/self.Nsn*(1+trand)
            #print(self.Nsn)
            self.progenitors += [SNR(None, None,DensityProfile=self.density_profile, helio_xyz=(xn[i],yn[i],zn[i]), name=self.name+"_ SN %02d"%i, known_age=avg_age*(i+1), )]
        return
        
        
    def __repr__(self) -> str:
        MSG = "%20s | n = %8.3g, age = %10.5g  Myr | R = %8.4g, (x,y,z) = (%8.3g,%8.3g,%8.3g), Nsn =%8.3g"%(self.name, self.local_n,self.computed_age,\
            self.data_radius(),self.x,self.y,self.z,self.Nsn )
        if len(self.progenitors)>0:
            for snr in self.progenitors:
                MSG+=" \n       |____ %s"%snr
        return MSG
