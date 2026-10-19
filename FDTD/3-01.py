import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,show,xlabel,ylabel,title,legend,grid,axis,tight_layout
from scipy import pi,sqrt,zeros,array,real,imag
from scipy.special import jv,jvp,hankel1,h1vp
from RI import WLx, NumWLx, epAg, epAu, RIAu, RIAg

def h1v(n,x):               # ハンケル関数(短く定義)
    return hankel1(n,x)
def a(n,x,mA,mB):           # 散乱係数 an
    return (mB*jv(n,mB*x)*jvp(n,mA*x)-mA*jv(n,mA*x)*jvp(n,mB*x))/ \
    (mB*jv(n,mB*x)*h1vp(n,mA*x)-mA*h1v(n,mA*x)*jvp(n,mB*x))
def b(n,x,mA,mB):           # 散乱係数 bn
    return (mA*jv(n,mB*x)*jvp(n,mA*x)-mB*jv(n,mA*x)*jvp(n,mB*x))/ \
    (mA*jv(n,mB*x)*h1vp(n,mA*x)-mB*h1v(n,mA*x)*jvp(n,mB*x))
def fn(n,phi):              # 係数 fn
    return (1/k0)*(cos(n*phi)+1j*sin(n*phi))*(pow(-1j,n))

rr = 25      # 円柱の半径
qq = 20      # ベッセル関数の次数

k0 = 2 * pi / WLx            #真空中の波数
m1 = RIAg                    # 円柱の屈折率
m2 = 1.0                     # 周辺媒質の屈折率
x = k0 * m2 * rr             #サイズ パラメータ

Qsca_tm = zeros(NumWLx, dtype=float)
Qext_tm = zeros(NumWLx, dtype=float)
Qsca_te = zeros(NumWLx, dtype=float)
Qext_te = zeros(NumWLx, dtype=float)

for n in range(-qq,qq):
    Qsca_tm = Qsca_tm + (2/x) * abs(b(n,x,m2,m1))**2   # TM 偏光散乱効率
    Qext_tm = Qext_tm + (2/x) * real(b(n,x,m2,m1))     # TM 偏光消光効率
    Qsca_te = Qsca_te + (2/x) * abs(a(n,x,m2,m1))**2   # TE 偏光散乱効率
    Qext_te = Qext_te + (2/x) * real(a(n,x,m2,m1))     # TE 偏光消光効率
        
Qabs_tm = Qext_tm - Qsca_tm
Qabs_te = Qext_te - Qsca_te

plt.figure(figsize=(8,6))            
plot(WLx,Qsca_tm, label=r"$Q_{\rm sca}(\rm TM)$",linewidth = 3.0, color='black')   
plot(WLx,Qsca_te, label=r"$Q_{\rm sca}(\rm TE)$",linewidth = 3.0, color='gray')   
xlabel("波長 (nm)",fontsize=22)  
ylabel("散乱・吸収効率",fontsize=22)    
title("散乱・吸収効率",fontsize=22)     
grid(True)      
axis([300,800,0,5])       
legend(fontsize=20,loc='upper right')     
plt.tick_params(labelsize=20)  
tight_layout()
show() 

