import h5py
import numpy as np
import math
import scipy
from scipy import signal
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline as spline
from scipy.fftpack import fft, fftshift ,ifft,rfft,fftfreq,rfftfreq
import os
c=2.9979e10
G=6.67408e-8
Msun=1.989e33
Length = G*Msun/c**2
Time = Length/c
Frequency=1/Time

EOS=['ALF2','APR4','GAM2','GNH3','H4','SLy']
MASS=['1200','1225','1250','1275','1300','1325','1350','1375','1400','1425','1450','1475','1500']
Q=['10','08','09']
point=['5000','4500','4000','3500','3000','2500','2000','1500','1000','500']


def fre_do(x,y,mass):
    fd=fft(y)
    N=len(y)
    if (N % 2) == 1:
        N=N+1
    T=x[1]-x[0]
    xf = np.linspace(0.0, 1.0/(2.0*T), int(N/2))
    fq=fftfreq(len(y))
    mask=fq>=0
    fd=2.0*(fd/N)
    fd=fd[mask]
    fd=abs(fd)
    return xf,fd


def analyze(rhM,time,mass):



    peaks,prop=scipy.signal.find_peaks(abs(rhM))
    ampls=rhM[peaks]
    merg=np.amax(abs(ampls))
    merg=np.where(abs(ampls)==merg)
    merg=int(merg[0])
    t1=peaks[merg]

    for mj in range(len(time)):
        if time[mj]>0.0:
            flag='pos'
            t0=mj
            break

    ampl=rhM[t0:]
    tim=time[t0:]

    #ampl=rhM
    #tim=time

    tuk=signal.tukey(len(ampl),0.03)
    dat=ampl*tuk

    fq,fd=fre_do(tim,dat,mass)

    mx=np.where(fd==np.amax(fd))[0][0]
    freq=fq[mx]
    amp=fd[mx]

    return fq,fd,tim,dat


if os.path.exists('results'):
    pass
else:
    os.mkdir('results')

if os.path.exists('results/log'):
    pass
else:
    os.mkdir('results/log')

if os.path.exists('results/linear'):
    pass
else:
    os.mkdir('results/linear')

if os.path.exists('results/3fig'):
    pass
else:
    os.mkdir('results/3fig')

if os.path.exists('results/3fig/linear'):
    pass
else:
    os.mkdir('results/3fig/linear')

if os.path.exists('results/3fig/log'):
    pass
else:
    os.mkdir('results/3fig/log')

for eos in EOS:

    for q in Q:
        for mas in MASS:
            for p in point:
                name=eos+'-q'+q+'-M'+mas+'.h_l2_m2.r500.t-1500_'+p+'.dat'

                try:

                    f=open('data/'+name,'r')

                    lines=f.readlines()[23:]

                    result1=[]
                    result2=[]
                    result3=[]
                    for x in lines:
                        for i in range(len(x.split(' '))):
                            if x.split(' ')[i]!='':
                                result1.append(x.split(' ')[i])
                                for j in range(i+1,len(x.split(' '))):
                                    if x.split(' ')[j]!='':
                                        result2.append(x.split(' ')[j])
                                        for k in range(j+1,len(x.split(' '))):
                                            if x.split(' ')[k]!='':
                                                result3.append(x.split(' ')[k])
                                                break
                                        break
                                break

                    time=[float(i) for i in result1]
                    rh1=[float(i) for i in result2]
                    rh2=[float(i) for i in result3]

                    rh=np.empty(len(rh1))
                    for i in range(len(rh1)):
                        rh[i]=rh1[i]+rh2[i]



                    bn=open('data/BNS/GNH3-q08-M1275.bns')
                    blines=bn.readlines()
                    exec(blines[8])
                    exec(blines[9])
                    mass=mass1+mass2
                    fq,fd,tim,dat=analyze(rh,time,mass)


                    fig=plt.figure()
                    plt.plot(fq*Frequency,fd)
                    plt.xlim(0,5000)
                    plt.xlabel('frequency (Hz)')
                    plt.savefig('results/linear/'+name+'.jpg')
                    plt.close()


                    fig=plt.figure()
                    plt.plot(fq*Frequency,fd)
                    plt.xlim(0,5000)
                    plt.ylim(10**(-6),1)
                    plt.yscale('log')
                    plt.xlabel('frequency (Hz)')
                    plt.savefig('results/log/'+name+'.jpg')
                    plt.close()

                    fig=plt.figure()
                    plt.subplot(212)
                    plt.plot((fq*Frequency),fd)

                    plt.xlim(0,5000)
                    plt.xlabel('Frequency (Hz)')
                    plt.legend(['Postmerger only'])
                    plt.subplot(222)
                    plt.plot(tim,dat)
                    plt.title('Postmerger')
                    plt.subplot(221)
                    plt.plot(time,rh)
                    plt.title('Time Domain')
                    plt.savefig('results/3fig/linear/'+name+'.jpg')
                    plt.close()


                    fig=plt.figure()
                    plt.subplot(212)
                    plt.plot((fq*Frequency),fd)
                    plt.yscale('log')
                    plt.xlim(0,5000)
                    plt.ylim(10**(-6),1)
                    plt.xlabel('Frequency (Hz)')
                    plt.legend(['Postmerger only'])
                    plt.subplot(222)
                    plt.plot(tim,dat)
                    plt.title('Postmerger')
                    plt.subplot(221)
                    plt.plot(time,rh)
                    plt.title('Time Domain')
                    plt.savefig('results/3fig/log/'+name+'.jpg')
                    plt.close()
                except OSError:
                    pass
