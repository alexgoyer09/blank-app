
"""
Nom_du_fichier: fct_etat6.py
Créé le      : 11 juin 2022
Créé par     : Alex Goyer 
Version num   : 6
Modifié le   : 2 juin 2025
"""


import numpy as np
import scipy.optimize
from sympy import *
import pandas as pd
import matplotlib.pyplot as plt
import math



#%% Class qui construit l'identité de la substance
class Substance:
    
    def __init__(self,name,Tc,Pc,w):
        self.name=name
        self.Tc=Tc
        self.Pc=Pc
        self.w=w



#%% Classe mère des équations d'état
class Eos:
    
    R=8.314 # J/mol.K <=> kPa.L/mol/K
    
    def __init__(self):
        
        ...
        

    def fit_data(self,substance):
        self.subs_name=substance.name
        self.Tc=substance.Tc
        self.Pc=substance.Pc
        self.w=substance.w                


    def get_a_b(self):
        
        ...
    
    
    def get_attr(self,P,v,T):
        
        ...


    def get_dynamic_variables(self,**kwargs):
        
        T=symbols("T")
        P=symbols("P")
        v=symbols("v")
        
        for parametre, valeur in kwargs.items():
            if parametre=="P" and valeur!=0.0:
                P=valeur
            elif parametre=="T" and valeur!=0.0:
                T=valeur
            elif parametre=="v" and valeur!=0.0:
                v=valeur
                
        return P,v,T
    

    def get_eos(self,**kwargs):
        
        P,v,T=self.get_dynamic_variables(**kwargs)
        
        a,b=self.get_a_b()
        
        attr=self.get_attr(P,v,T)
        
        eos=Eos.R*T/(v-b)-P-attr
        
        return eos
    

    def get_eos_derivated(self,symbol,**kwargs):

        eos=self.get_eos(**kwargs)
        
        eos_derivated=diff(eos,symbol)
        
        return eos_derivated
    
    
    def get_eos_integrated(self,symbol,limits,**kwargs):
        
        eos=self.get_eos(**kwargs)
        
        eos_integrated=Integral(eos,(symbol,limits[0],limits[-1]))
        
        return eos_integrated.evalf()

    
    def get_roots(self,eos):
        
        roots=solve(eos)
        
        return roots


    def isotherm(self,list_v,T):
       
        list_P=[]
            
        for v in list_v:
            
            eos=self.get_eos(T=T,v=v)
            roots=self.get_roots(eos)        
            
            list_P.append(float(roots[0]))
            
        return list_v,list_P
    
    
    def data_isotherms(self,range_v,n_data_points,list_isotherms):
        
        vmin,vmax=range_v
        list_v=np.linspace(vmin,vmax,n_data_points) 
        
        array_P=np.zeros((len(list_v),len(list_isotherms)))

        for i in range(len(list_isotherms)):
            
            list_v,list_P=self.isotherm(list_v,list_isotherms[i]*self.Tc)
            
            array_P[:,i]=list_P
            
        data_isotherms={"v":list_v,
                        "P":array_P,
                        "T":list_isotherms
                        }
            
        return data_isotherms
      
    
    def get_Psat_iteration_range(self,T):
        
        Psat_range=[]
        
        eos_derivated=self.get_eos_derivated(symbols("v"),T=T)
        roots=solve(eos_derivated)
        
        for r in roots:            
            eos=self.get_eos(T=T,P=symbols("P"),v=re(r))
            P=self.get_roots(eos)

            if re(P[0])>=0:
                Psat_range.append(re(P[0]))
                
        if len(Psat_range)==1:
            Psat_range.append(0)

        return sorted(Psat_range)
    
    
    def find_vl_vg(self,T,P):

        eos=self.get_eos(T=T,P=P)
        roots=self.get_roots(eos)
        
        vl=re(roots[0])
        vg=re(roots[-1])
        
        return vl,vg
        
        
    def calculate_area_under_eos(self,T,P):
        
        vl,vg=self.find_vl_vg(T,P)       
        area=self.get_eos_integrated(symbols("v"),[vl,vg],T=T,P=P)
        
        return area
    
    
    def iteration_Psat(self,T,Psat_range,dP):

        Pmin,Pmax=min(Psat_range),max(Psat_range)
        
        P=Pmax
        Psat=Pmax
        area=self.calculate_area_under_eos(T,P)
        min_area=self.calculate_area_under_eos(T,P)
        
        while area<0 and P>Pmin:
            
            # print(area,min_area,P,Psat)
            area=self.calculate_area_under_eos(T,P)
            
            if abs(area)<abs(min_area) and area<0:
                min_area=area
                Psat=P
                
            P-=dP 
        
        return Psat,min_area
    
    
    def create_list_dP(self,Pmax,precision):
        
        n_Pmax=math.floor(math.log10(abs(Pmax)))
        
        n=int(n_Pmax+precision)
        
        list_dP=np.logspace(-precision,n_Pmax,n+1)
            
        return list_dP[::-1]
    
    
    def raffinate_iteration(self,T,precision):
        
        Psat_range=self.get_Psat_iteration_range(T)
        Pmin,Pmax=min(Psat_range),max(Psat_range)
        
        list_dP=self.create_list_dP(Pmax, precision)
        
        for dP in list_dP:
            
                        
            Psat,min_area=self.iteration_Psat(T,Psat_range,dP)
            Psat_range=[Pmin,Psat]
        
        return round(Psat,precision)
    
    
    def saturation_array(self,list_isotherms):
        
        list_bubble=[]
        list_dew=[]
        list_Psat=[]
        
        for Tr in list_isotherms:
            
            T=Tr*self.Tc
            
            Psat=self.raffinate_iteration(T, 2)
            vl,vg=self.find_vl_vg(T, Psat)
            
            list_bubble.append(vl)
            list_dew.append(vg)
            list_Psat.append(Psat)
            
        data_saturation={"vl":list_bubble,
                        "vg":sorted(list_dew, reverse=True),
                        "Psat":list_Psat
                        }
        
        return data_saturation
        
                
    
    def PvT(self,data_isotherms,data_saturation,limits_v,limits_P):
    # def PvT(self,data_isotherms,limits_v,limits_P):
               
        list_v,array_P,list_isotherms=data_isotherms["v"],data_isotherms["P"],data_isotherms["T"]
        list_vl,list_vg,list_Psat=data_saturation["vl"],data_saturation["vg"],data_saturation["Psat"]
        
        
        fig, ax = plt.subplots(dpi=300)
        
        ax.plot(list_v,array_P[:,],label=["Tr = "+str(round(u,2)) for u in list_isotherms])
        
        ax.plot(list_vl,list_Psat,label="Bubble point")
        ax.plot(list_vg,list_Psat,label="Dew point")
        
        ax.set_title(f"Diagramme PvT de {self.subs_name} ({self.eos_name})")
        ax.set_xlabel("Volume spécifique (L/mol)")
        ax.set_ylabel("Pression (kPa)")

        ax.set_xlim(limits_v)  # Limites de l’axe des X
        ax.set_ylim(limits_P)  # Limites de l’axe des Y
        
        ax.legend()
        
        ax.set_xscale('log')
        ax.grid(True, which='both', color='gray')
        
        return fig
    


#%% Class object de Van der Waals [Vdw]
class Vdw(Eos):
    
    def __init__(self):
        super().__init__()
        self.eos_name="Van der Waals"
    
    
    def get_a_b(self):
        
        a=27*np.power(Eos.R*self.Tc,2)/64/self.Pc   
        b=Eos.R*self.Tc/8/self.Pc                  
        
        return a,b
    
    
    def get_attr(self,P,v,T):
        
        a,b=self.get_a_b()
        
        attr=a/(v**2)                 
        
        return attr    
    
    
    
#%% Class object de Redlich-Kwong [Rk]
class Rk(Eos):
    
    def __init__(self):
        super().__init__()
        self.eos_name="Redlich-Kwong"
    
    
    def get_a_b(self):
        
        a=(1/9/(np.power(2,1/3)-1))*np.power(Eos.R,2)*np.power(self.Tc,2.5)/self.Pc  
        b=(np.power(2,1/3)-1)/3*Eos.R*self.Tc/self.Pc           
        
        return a,b
    
    
    def get_attr(self,P,v,T):
        
        a,b=self.get_a_b()
        
        attr=a/np.power(T,1/2)/v/(v+b) 
        
        return attr
    
    
    
#%% Class object de Peng-Robinson [Pr]
class Pr(Eos):
    
    def __init__(self):
        super().__init__()
        self.eos_name="Peng-Robinson" 
    
    
    def get_a_b(self):
        
        a=0.45724*(Eos.R*self.Tc)**2/self.Pc
        b=0.07780*Eos.R*self.Tc/self.Pc               
        
        return a,b
    
    
    def get_k(self):
        
        k=0.37464+1.54226*self.w-0.26992*(self.w)**2
        
        return k
    
    
    def get_alpha(self,T):
        
        k=self.get_k()
        
        Tr=T/self.Tc
        
        alpha=(1+k*(1-(Tr)**(1/2)))**2
        
        return alpha
    
    
    def get_attr(self,P,v,T):
        
        a,b=self.get_a_b()
        
        alpha=self.get_alpha(T)
        
        attr=a*alpha/(v*(v+b)+b*(v-b))
        
        return attr



#%% Class object de Soave-Redlich-Kwong [Srk]
class Srk(Eos):  
    
    def __init__(self):
        super().__init__()
        self.eos_name="Soave-Redlich-Kwong"   
        
        
    def get_a_b(self):
        
        a=0.42748*(Eos.R*self.Tc)**2/self.Pc
        b=0.08664*Eos.R*self.Tc/self.Pc               
        
        return a,b
    
    
    def get_k(self):
        
        k=0.480+1.574*self.w-0.176*(self.w)**2
        
        return k
    
    
    def get_alpha(self,T):
        
        k=self.get_k()
        
        Tr=T/self.Tc
        
        alpha=(1+k*(1-(Tr)**(1/2)))**2
        
        return alpha
    
    
    def get_attr(self,P,v,T):
        
        a,b=self.get_a_b()
        
        alpha=self.get_alpha(T)
        
        attr=a*alpha/v/(v+b)
        
        return attr   
    
    
    
#%% Class object de Lee-Kesler [Lk]
# class Lk(Eos):
    
    # eos_name="Lee-Kesler"
    
#     def __init__(self):
#         super().__init__()
#         self.eos_title="Lee-Kesler"     
    
#         self.b1=[0.1181193,0.2026579]
#         self.b2=[0.265728,0.331511]
#         self.b3=[0.154790,0.027655]
#         self.b4=[0.030323,0.203488]
#         self.c1=[0.0236744,0.0313385]
#         self.c2=[0.0186984,0.0503618]
#         self.c3=[0,0.016901]
#         self.c4=[0.042724,0.041577]
#         self.d1=[1.55488e-05,4.8736e-05]
#         self.d2=[6.23689e-05,7.40336e-06]
#         self.beta=[0.65392,1.226]
#         self.gamma=[0.060167,0.03754]   


 
    
    
    




compose1=Substance("propane",562,4910,0.152)
# compose2=Substance("ethane",305.4,48.74*100,0.152)
# compose3=Substance("",562,4910,0)

model=Vdw()
model.fit_data(compose1)


range_v=[0.1,5]
n_data_points=500
list_isotherms=[0.5,1]
list_isotherms_Psat=[0.1,0.3,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.87,0.89,0.9,0.92,0.93,0.95,0.97,0.98,1]

# area=model.calculate_area_under_eos(562*0.1,0.01)
# print(area)
# model.get_Psat_iteration_range(560)

# Psat=model.raffinate_iteration(562*0.3, 4)
# print(Psat)

data_isotherms=model.data_isotherms(range_v, n_data_points, list_isotherms)
data_saturation=model.saturation_array(list_isotherms_Psat)
vg,vl,Psat=data_saturation["vl"],data_saturation["vg"],data_saturation["Psat"]

model.PvT(data_isotherms, data_saturation,[0.1,5], [0,1.5*4910])
# model.PvT(data_isotherms, data_saturation,[0.1,5], [-1.5*4910,1.5*4910])


    
    
    
















    
    
    
    