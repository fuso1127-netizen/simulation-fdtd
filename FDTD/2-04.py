import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,show,xlabel,ylabel,title,legend,grid,axis,rcParams,tight_layout
from scipy import real,imag,pi
from RI import WLx, epAg, epAu, RIAu, RIAg

r3 = 25      # コアの半径
s = 1.5      # シェルの半径/コアの半径
r2 = r3*s    # シェルの半径
n1 = 1       # 周辺媒質の屈折率
n2 = RIAu    # 微粒子コアの屈折率
n3 = 1.5     # シェルの屈折率
k = 2 * pi / WLx   # 波数の配列
  
delta = (n2**2) * (2 * (n1**2) * (1 + 2 * (s**3)) + (n3**2) * (2 + s**3)) - 2 * ((n2**2)**2 + (n1**2) * (n3**2)) * (1 - s**3)
b11 = (s**3 / delta) * ((n2**2) * ((n1**2) * (1 + 2*s**3) - (n3**2) * (2 + s**3)) + (2 * (n2**2)**2 - (n1**2) * (n3**2)) * (1 - s**3))

alpha = -4 * pi * r2**3 * (n1**2) * b11    # 分極率
Csca = k**4 / (6 * pi) * abs(alpha)**2     # 散乱断面積
Cabs = k * imag(alpha)          # 吸収断面積
Qsca = Csca / ((r2**2) * pi)    # 散乱効率
Qabs = Cabs / ((r2**2) * pi)   # 吸収効率
 
plt.figure(figsize=(8,6))     
plot(WLx,Qsca, label=r"$Q_{{\rm sca}}$",linewidth = 3.0, color='black')  
plot(WLx,Qabs, label=r"$Q_{{\rm abs}}$",linewidth = 3.0, color='gray')   
xlabel("波長 (nm)",fontsize=22) 
ylabel("散乱・吸収効率",fontsize=22) 

grid(True)    
axis([400,1000,0,20])      
plt.tick_params(labelsize=20) 
#legend(fontsize=20,loc='lower left')
tight_layout()
show()

