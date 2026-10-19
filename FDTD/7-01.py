import matplotlib.pyplot as plt
import numpy as np
from scipy import zeros
from matplotlib.pyplot import plot,show,xlabel,ylabel,title,legend,grid, axis,subplot

f = open("qtable","r")   # "qtable"ファイルをオープン
dat = f.read()           # 文字変数 dat にすべての文字列を読み込む
f.close()                # "qtable"ファイルをクローズ

dat = dat.split("\n")    # 文字変数 dat を分割
datLEN=len(dat)-15       # 分割した dat からヘッダの部分を除いた行数を求める

WLx=zeros(datLEN)
Qext=zeros(datLEN)
Qabs=zeros(datLEN)
Qsca=zeros(datLEN)

DDSversion=dat[0]        # "qtable"の 0 行目  DDSCAT のバージョン
Target=dat[1]            # "qtable"の 1 行目 キーワード(target の種類)
Shape=dat[4]             # "qtable"の 4 行目  target の
NumDipole=dat[5]         # "qtable"の 5 行目 双極子数
aEff=dat[1][1:11]        # "qtable"の 15 行目 有効半径
    
i=14
j=0
while j <= datLEN-1:
    WLx[j]=float(dat[i][12:22])*1000    # 波長の読込み  μm を nm に変換
    Qext[j]=float(dat[i][23:33])        # 消光断面積の読込み  μm を nm に変換
    Qabs[j]=float(dat[i][34:44])        # 吸収断面積の読込み  μm を nm に変換
    Qsca[j]=float(dat[i][45:55])        # 散乱断面積の読込み  μm を nm に変換
    i=i+1
    j=j+1

plot(WLx,Qsca, label=r"$Q_{\rm sca}$")   # Scattering Cross-section
plot(WLx,Qext, label=r"$Q_{\rm ext}$")   # Extinction Cross-section
plot(WLx,Qabs, label=r"$Q_{\rm abs}$")   # Absorption Cross-section

xlabel("波長 (nm)",fontsize=12)      # x-axis label
ylabel("散乱・吸収効率",fontsize=12)        # y-axis label
title("散乱・吸収効率",fontsize=12)        # Title of the graph 
grid(True)                                 # Show Grid 
axis([400,700,0,15])                   # Plot Range 
legend()                                   # Show Legend
show() 
