import scipy as sp
import scipy.special
import matplotlib as mpl
import matplotlib.pyplot as plt
from RI import WLx, NumWLx, epAu, RIAu
from scipy import pi,arange,zeros,array,real,imag
from matplotlib.pyplot import plot,show,xlabel,ylabel,title,legend,grid, axis
from scipy.special import spherical_jn,spherical_yn

def psi(n,z):                     # Riccati-Bessel function of first kind
    return z*spherical_jn(n,z)
def psiDz(n,z):                   # Derivative of Riccati-Bessel function of first kind
    return spherical_jn(n,z)+z*spherical_jn(n,z,1)
def xi(n,z):                      # Riccati-Bessel function of third kind
    return z*(spherical_jn(n,z)+1j*spherical_yn(n,z))
def xiDz(n,z):                    # Derivative of Riccati-Bessel function of third kind
    return (spherical_jn(n,z)+1j*spherical_yn(n,z)) \
            +z*(spherical_jn(n,z,1)+1j*spherical_yn(n,z,1))
def chi(n,z):                   
    return -z*spherical_yn(n,z)
def chiDz(n,z):                   
    return -spherical_yn(n,z)-z*spherical_yn(n,z,1)

def aa(n,m1,m2,x):
    return (m2*psi(n,m2*x)*psiDz(n,m1*x)-m1*psiDz(n,m2*x)*psi(n,m1*x)) \
           /(m2*chi(n,m2*x)*psiDz(n,m1*x)-m1*chiDz(n,m2*x)*psi(n,m1*x))
def bb(n,m1,m2,x):
    return (m2*psi(n,m1*x)*psiDz(n,m2*x)-m1*psi(n,m2*x)*psiDz(n,m1*x)) \
           /(m2*chiDz(n,m2*x)*psi(n,m1*x)-m1*psiDz(n,m1*x)*chi(n,m2*x)) 

def a(n,m1,m2,x,y):
    return (psi(n,y)*(psiDz(n,m2*y)-aa(n,m1,m2,x)*chiDz(n,m2*y))     \
            -m2*psiDz(n,y)*(psi(n,m2*y)-aa(n,m1,m2,x)*chi(n,m2*y)))  \
            /(xi(n,y)*(psiDz(n,m2*y)-aa(n,m1,m2,x)*chiDz(n,m2*y))    \
            -m2*xiDz(n,y)*(psi(n,m2*y)-aa(n,m1,m2,x)*chi(n,m2*y))) 
def b(n,m1,m2,x,y):
    return (m2*psi(n,y)*(psiDz(n,m2*y)-bb(n,m1,m2,x)*chiDz(n,m2*y))  \
            -psiDz(n,y)*(psi(n,m2*y)-bb(n,m1,m2,x)*chi(n,m2*y)))     \
            /(m2*xi(n,y)*(psiDz(n,m2*y)-bb(n,m1,m2,x)*chiDz(n,m2*y)) \
              -xiDz(n,y)*(psi(n,m2*y)-bb(n,m1,m2,x)*chi(n,m2*y))) 

r3 = 100     # コアの半径
s = 1.1      # シェルの半径/コアの半径
r2 = r3*s    # シェルの半径
qq = 20      # ベッセル関数の次数

k0 = 2*pi/WLx       # 真空中の波数
n1 = 1              # 周辺媒質の屈折率
n2 = RIAu           # 微粒子コアの屈折率
n3 = 1.5            # シェルの屈折率

x = k0 * n1 * r3    #サイズ パラメータ(コア)
y = k0 * n1 * r2    # サイズ パラメータ(シェル

m2 = n2 / n1        # 比屈折率(シェル)
m3 = n3 / n1        # 比屈折率(コア)

Csca = zeros(NumWLx, dtype=complex)
Cext = zeros(NumWLx, dtype=complex)
Cabs = zeros(NumWLx, dtype=complex)    

for n in range(qq):
    Csca = Csca + (2*pi / k0**2) * \
           (2 * (n+1)+1) * (abs(a(n+1,m3,m2,x,y))**2 + abs(b(n+1,m3,m2,x,y)**2))
    Cext = Cext + (2*pi / k0**2) * \
           (2 * (n+1)+1) * (real(a(n+1,m3,m2,x,y) + b(n+1,m3,m2,x,y)))
    Cabs = Cext - Csca   

Qsca = Csca / ((r2**2) * pi)    # 散乱効率
Qabs = Cabs / ((r2**2) * pi)    # 吸収効率

plot(WLx,abs(Qsca), label=r"$Q_{\rm sca}$",linewidth = 3.0, color='black')  
plot(WLx,abs(Qabs), label=r"$Q_{\rm abs}$",linewidth = 3.0, color='gray')    

xlabel("波長 (nm)",fontsize=22)     
ylabel("散乱・吸収効率",fontsize=22)         
title(r"散乱・吸収効率",fontsize=22)          
grid(True)  
axis([400,1000,0,10])            
legend(fontsize=20,loc='upper left')        
plt.tick_params(labelsize=18)   
show() 
