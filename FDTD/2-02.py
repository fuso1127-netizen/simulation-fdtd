import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,show,xlabel,ylabel,title,legend,grid,axis,rcParams,tight_layout
from scipy import real,imag,pi
from RI import WLx, epAg, epAu, RIAu, RIAg

n1 = 1       # 周辺媒質の屈折率
n2 = RIAg    # 粒子の屈折率
r=25         # 粒子の半径
k = 2 * pi / WLx    # 波数の配列
alpha = 4 * pi * (r**3) * (n1**2) * (n2**2 - n1**2) / (n2**2 + 2 * n1**2)  # 分極率の計算
Csca = k**4 / (6 * pi) * abs(alpha)**2   # 散乱断面積
Cabs = k * imag(alpha)                   # 吸収断面積
Qsca = Csca / ((r**2) * pi)              # 散乱効率
Qabs = Cabs / ((r**2) * pi)              # 吸収効率

plt.figure(figsize=(8,6))        # figure size
plot(WLx,real(RIAg), label="real",linewidth = 3.0, color='black')        # real part of the refractive index of gold
plot(WLx,imag(RIAg), label="imaginary",linewidth = 3.0, color='gray')   # real part of the refractive index of gold
xlabel("波長 (nm)",fontsize=22)       # x-axis label
ylabel("屈折率",fontsize=22)          # y-axis label
title("銀の屈折率",fontsize=22)   # Title of the graph 
grid(True)                                 # Show Grid 
axis([300,700,0,5])                       # Plot Range
plt.tick_params(labelsize=20)  # scale fontsize=18pt
legend(fontsize=20,loc='lower right')
tight_layout()
show()

plt.figure(figsize=(8,6))        # figure size
plot(WLx,real(epAg), label="real",linewidth = 3.0, color='black')        # real part of the refractive index of gold
plot(WLx,imag(epAg), label="imaginary",linewidth = 3.0, color='gray')   # real part of the refractive index of gold
xlabel("波長 (nm)",fontsize=22)       # x-axis label
ylabel("誘電率",fontsize=22)          # y-axis label
title(r"銀の誘電率",fontsize=22)   # Title of the graph 
grid(True)                                 # Show Grid 
axis([300,700,-10,5])                       # Plot Range
plt.tick_params(labelsize=20)  # scale fontsize=18pt
legend(fontsize=20,loc='lower right')
tight_layout()
show()

plt.figure(figsize=(8,6))        # figure size
plot(WLx,Csca/10**4, label=r"$C_{{\rm sca}}$",linewidth = 3.0, color='black')        # real part of the refractive index of gold
plot(WLx,Cabs/10**4, label=r"$C_{{\rm abs}}$",linewidth = 3.0, color='gray')   # real part of the refractive index of gold
xlabel("波長 (nm)",fontsize=22)       # x-axis label
ylabel(r"散乱・吸収断面積(×10000 nm^2)",fontsize=22)          # y-axis label
title(r"散乱・吸収断面積",fontsize=22)   # Title of the graph 
grid(True)                                 # Show Grid 
axis([300,700,0,5])                       # Plot Range
plt.tick_params(labelsize=20)  # scale fontsize=18pt
legend(fontsize=20,loc='lower right')
tight_layout()
show()

plt.figure(figsize=(8,6))        # figure size
plot(WLx,Qsca, label=r"$Q_{{\rm sca}}$",linewidth = 3.0, color='black')        # real part of the refractive index of gold
plot(WLx,Qabs, label=r"$Q_{{\rm abs}}$",linewidth = 3.0, color='gray')   # real part of the refractive index of gold
xlabel("波長 (nm)",fontsize=22)       # x-axis label
ylabel("散乱・吸収効率",fontsize=22)          # y-axis label
title(r"散乱・吸収効率$  (R=25 nm)",fontsize=22)   # Title of the graph 
grid(True)                                 # Show Grid 
axis([300,700,0,20])                       # Plot Range
plt.tick_params(labelsize=20)  # scale fontsize=18pt
legend(fontsize=20,loc='lower right')
tight_layout()
show()
