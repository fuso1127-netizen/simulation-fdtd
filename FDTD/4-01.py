import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,show,xlabel,ylabel,title,legend,grid,axis,rcParams,tight_layout
from scipy import real,imag,pi,sqrt,log
from RI import WLx, epAg, epAu, RIAu, RIAg

n1 = RIAu    # 回転楕円体の屈折率
n2 = 1       # 周辺媒質の屈折率
a = b = 10   # 非回転軸の長さ〔nm〕
c = 50       # 回転軸の長さ〔nm〕
ee = sqrt(1-(a/c)**2)    # 離心率
lz = (1-ee**2)/ee**2 * (1/(2*ee) * log((1+ee)/(1-ee))-1)   # z 方向 反電場係数
lx = (1-lz)/2            # x 方向 反電場係数
k = 2 * pi / WLx         #真空中の波数

alphax = 4*pi*a*b*c*(n2**2)*((n1**2)-(n2**2))/(3*((n2**2)+lx*((n1**2)-(n2**2))))    # x 方向 分極率
alphaz = 4*pi*a*b*c*(n2**2)*((n1**2)-(n2**2))/(3*((n2**2)+lz*((n1**2)-(n2**2))))    # z 方向 分極率

Csca_x = k**4 / (6 * pi) * abs(alphax)**2     # x 方向 散乱断面積
Cabs_x = k * imag(alphax)      # x 方向 吸収断面積
Qsca_x = Csca_x / (a*a*pi)     # x 方向 散乱効率
Qabs_x = Cabs_x / (a*a*pi)     # x 方向 吸収効率

Csca_z = k**4 / (6 * pi) * abs(alphaz)**2      # z 方向 散乱断面積
Cabs_z = k * imag(alphaz)      # z 方向 吸収断面積
Qsca_z = Csca_z / (a*c*pi)     # z 方向 散乱効率
Qabs_z = Cabs_z / (a*c*pi)     # z 方向 吸収効率

plt.figure(figsize=(8,6)) 
plot(WLx,Qsca_x, label=r"$Q_{{\rm sca},a}$",linewidth = 3.0, color='black') 
plot(WLx,Qabs_x, label=r"$Q_{{\rm abs},a}$",linewidth = 3.0, color='gray') 
xlabel("波長 (nm)",fontsize=22)      # x-axis label
ylabel("散乱・吸収効率",fontsize=22)           # y-axis label
title("散乱・吸収効率 a-軸",fontsize=22)   # Title of the graph 
grid(True)                                 # Show Grid 
axis([300,1000,0,1])                       # Plot Range 
plt.tick_params(labelsize=20)
legend(fontsize=20,loc='upper right')
tight_layout()
show()

plt.figure(figsize=(8,6)) 
plot(WLx,Qsca_z, label=r"$Q_{{\rm sca},c}$",linewidth = 3.0, color='black') 
plot(WLx,Qabs_z, label=r"$Q_{{\rm abs},c}$",linewidth = 3.0, color='gray') 
xlabel("波長 (nm)",fontsize=22)      # x-axis label
ylabel("散乱・吸収効率",fontsize=22)           # y-axis label
title("散乱・吸収効率 c-軸",fontsize=22)   # Title of the graph 
grid(True)                                 # Show Grid 
axis([300,1000,0,50])                     # Plot Range
plt.tick_params(labelsize=20)
legend(fontsize=20,loc='upper left')
tight_layout()
show()
