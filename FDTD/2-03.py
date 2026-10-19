import scipy as sp
import scipy.special
import matplotlib as mpl
import matplotlib.pyplot as plt
from RI import WLx, NumWLx, epAu, RIAu
from scipy import pi,arange,zeros,array,real,imag
from matplotlib.pyplot import plot,show,xlabel,ylabel,title,legend,grid, axis
from scipy.special import spherical_jn,spherical_yn

def psi(n,z):                     # 第一種 Riccati-Bessel 関数
    return z*spherical_jn(n,z)
def psiDz(n,z):                   # 第一種 Riccati-Bessel 関数の微分
    return spherical_jn(n,z)+z*spherical_jn(n,z,1)
def xi(n,z):                      # 第三種 Riccati-Bessel 関数
    return z*(spherical_jn(n,z)+1j*spherical_yn(n,z))
def xiDz(n,z):                    # 第三種 Riccati-Bessel 関数の微分
    return (spherical_jn(n,z)+1j*spherical_yn(n,z))+z* \
           (spherical_jn(n,z,1)+1j*spherical_yn(n,z,1))
def a(n,m,x):                   # 散乱係数 an
    return (m*psi(n,m*x)*psiDz(n,x)-psi(n,x)*psiDz(n,m*x))/ \
           (m*psi(n,m*x)*xiDz(n,x)-xi(n,x)*psiDz(n,m*x))
def b(n,m,x):                   # 散乱係数 bn
    return (psi(n,m*x)*psiDz(n,x)-m*psi(n,x)*psiDz(n,m*x))/ \
           (psi(n,m*x)*xiDz(n,x)-m*xi(n,x)*psiDz(n,m*x))

n1 = 1.0    # 周辺媒質の屈折率
n2 = RIAu   # 粒子の屈折率
r = 50      # 粒子の半径
qq = 50     # ベッセル関数の次数

Csca = zeros(NumWLx, dtype=complex)
Cext = zeros(NumWLx, dtype=complex)
Cabs = zeros(NumWLx, dtype=complex)

k0 = 2*pi/WLx                      # 真空中の波数
x = k0*n1*r                        # サイズ パラメータ
m = n2/n1                          # 比屈折率

for n in range(qq):
    Csca = Csca + (2*pi/k0**2)*(2*(n+1)+1)*(abs(a(n+1,m,x)**2)+abs(b(n+1,m,x)**2))
    Cext = Cext + (2*pi/k0**2)*(2*(n+1)+1)*(real(a(n+1,m,x)+b(n+1,m,x)))
    Cabs = Cext - Csca

Qsca = Csca / ((r**2) * pi)   # 散乱効率
Qabs = Cabs / ((r**2) * pi)   # 吸収効率

plot(WLx,Qsca, label=r"$Q_{\rm sca}$",linewidth = 3.0, color='black') 
plot(WLx,Qabs, label=r"$Q_{\rm abs}$",linewidth = 3.0, color='gray')   

xlabel("波長 (nm)",fontsize=22)   
ylabel("散乱・吸収効率",fontsize=22)        
    
grid(True)    
axis([400,800,0,5]) 
legend(fontsize=20,loc='upper right')     
plt.tick_params(labelsize=18)  # Show Legend
show() 
