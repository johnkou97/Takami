import h5py
import numpy as np
import math
import scipy
from scipy import signal
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline as spline
from scipy.fftpack import fft, fftshift ,ifft,rfft,fftfreq,rfftfreq
import os
import pywt
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


def f20(M,R6):
    return 8.943+4.059*M-1.332*R6-.358*(M**2)-.182*R6*M+.048*(R6**2)

def fspir(M,R8):
    return 6.264+1.929*M-.645*R8+.881*(M**2)-.311*R8*M+.03*(R8**2)

def fpeak(M,R6):
    return 13.822-0.576*M-1.375*R6+.479*(M**2)-.073*R6*M+.044*(R6**2)

def f20_a(M,R6):
    return 9.586+4.09*M-1.427*R6+.048*(M**2)-.261*R6*M+.055*(R6**2)

def fspir_a(M,R8):
    return 5.846+1.75*M-.555*R8+1.002*(M**2)-.316*R8*M+.026*(R8**2)

def fpeak_a(M,R8):
    return 10.942-.369*M-.987*R8+1.095*(M**2)-.201*R8*M+.036*(R8**2)


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


if os.path.exists('results/q1'):
    pass
else:
    os.mkdir('results/q1')

if os.path.exists('results/q1/linear'):
    pass
else:
    os.mkdir('results/q1/linear')

if os.path.exists('results/q1/log'):
    pass
else:
    os.mkdir('results/q1/log')

if os.path.exists('results/spec'):
    pass
else:
    os.mkdir('results/spec')


#q=1
q='10'
for eos in EOS:
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
                        rh[i]=rh1[i]



                    bn=open('data/BNS/'+eos+'-q'+q+'-M'+mas+'.bns')
                    blines=bn.readlines()
                    exec(blines[8])
                    exec(blines[9])
                    mass=mass1+mass2
                    q2=mass1/mass2
                    Mc=pow(q2/pow(1+q2,2),3/5)*mass
                    m_r=np.load('tid_def/'+eos+'.npy')
                    mx=np.amax(m_r[0])
                    idx=np.where(m_r[0]==mx)
                    idx=idx[0][0]
                    cs=spline(m_r[0][1:idx],m_r[1][1:idx])
                    r68=np.zeros((1,2))
                    r68[0,0]=cs(1.6)*Length/1.0e5
                    r68[0,1]=cs(1.8)*Length/1.0e5
                    f_2=f20(Mc,r68[0,0])
                    f_s=fspir(Mc,r68[0,1])
                    f_p=fpeak_a(Mc,r68[0,0])
                    f_0=2*f_p-f_2


                    fq,fd,tim,dat=analyze(rh,time,mass)


                    fig=plt.figure()
                    plt.plot(fq*Frequency,fd)
                    ax=plt.subplot()
                    ax.axvline(x=(f_p*Mc)*1000,color='r',label='peak')
                    ax.axvspan((f_p*Mc)*1000-196, (f_p*Mc)*1000+196, alpha=0.3, color='grey')
                    ax.axvline(x=(f_2*Mc)*1000,color='g',label='2-0')
                    ax.axvspan((f_2*Mc)*1000-229, (f_2*Mc)*1000+229, alpha=0.3, color='yellow')
                    ax.axvline((f_s*Mc)*1000,color='orange',label='spiral')
                    ax.axvspan((f_s*Mc)*1000-286, (f_s*Mc)*1000+286, alpha=0.3, color='cyan')
                    ax.axvline((f_0*Mc)*1000,linestyle="--",color='grey',label='2+0')
                    plt.xlim(0,5000)
                    plt.xlabel('frequency (Hz)')
                    plt.savefig('results/q1/linear/'+name+'.jpg')
                    plt.close()


                    fig=plt.figure()
                    plt.plot(fq*Frequency,fd)
                    ax=plt.subplot()
                    ax.axvline(x=(f_p*Mc)*1000,color='r',label='peak')
                    ax.axvspan((f_p*Mc)*1000-196, (f_p*Mc)*1000+196, alpha=0.3, color='grey')
                    ax.axvline(x=(f_2*Mc)*1000,color='g',label='2-0')
                    ax.axvspan((f_2*Mc)*1000-229, (f_2*Mc)*1000+229, alpha=0.3, color='yellow')
                    ax.axvline((f_s*Mc)*1000,color='orange',label='spiral')
                    ax.axvspan((f_s*Mc)*1000-286, (f_s*Mc)*1000+286, alpha=0.3, color='cyan')
                    ax.axvline((f_0*Mc)*1000,linestyle="--",color='grey',label='2+0')
                    plt.xlim(0,5000)
                    plt.ylim(10**(-6),1)
                    plt.yscale('log')
                    plt.xlabel('frequency (Hz)')
                    plt.savefig('results/q1/log/'+name+'.jpg')
                    plt.close()




                except OSError:
                    pass


#all cases
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
                        rh[i]=rh1[i]



                    bn=open('data/BNS/'+eos+'-q'+q+'-M'+mas+'.bns')
                    blines=bn.readlines()
                    exec(blines[8])
                    exec(blines[9])
                    mass=mass1+mass2
                    q2=mass1/mass2
                    Mc=pow(q2/pow(1+q2,2),3/5)*mass
                    m_r=np.load('tid_def/'+eos+'.npy')
                    mx=np.amax(m_r[0])
                    idx=np.where(m_r[0]==mx)
                    idx=idx[0][0]
                    cs=spline(m_r[0][1:idx],m_r[1][1:idx])
                    r68=np.zeros((1,2))
                    r68[0,0]=cs(1.6)*Length/1.0e5
                    r68[0,1]=cs(1.8)*Length/1.0e5
                    f_2_a=f20_a(Mc,r68[0,0])
                    f_s_a=fspir_a(Mc,r68[0,1])
                    f_p_a=fpeak_a(Mc,r68[0,1])
                    f_0_a=2*f_p_a-f_2_a


                    fq,fd,tim,dat=analyze(rh,time,mass)


                    fig=plt.figure()
                    plt.plot(fq*Frequency,fd)
                    ax=plt.subplot()
                    ax.axvline(x=(f_p_a*Mc)*1000,color='r',label='peak')
                    ax.axvspan((f_p_a*Mc)*1000-196, (f_p_a*Mc)*1000+196, alpha=0.3, color='grey')
                    ax.axvline(x=(f_2_a*Mc)*1000,color='g',label='2-0')
                    ax.axvspan((f_2_a*Mc)*1000-229, (f_2_a*Mc)*1000+229, alpha=0.3, color='yellow')
                    ax.axvline((f_s_a*Mc)*1000,color='orange',label='spiral')
                    ax.axvspan((f_s_a*Mc)*1000-286, (f_s_a*Mc)*1000+286, alpha=0.3, color='cyan')
                    ax.axvline((f_0_a*Mc)*1000,linestyle="--",color='grey',label='2+0')
                    plt.xlim(0,5000)
                    plt.xlabel('frequency (Hz)')
                    plt.savefig('results/linear/'+name+'.jpg')
                    plt.close()


                    fig=plt.figure()
                    plt.plot(fq*Frequency,fd)
                    ax=plt.subplot()
                    ax.axvline(x=(f_p_a*Mc)*1000,color='r',label='peak')
                    ax.axvspan((f_p_a*Mc)*1000-196, (f_p_a*Mc)*1000+196, alpha=0.3, color='grey')
                    ax.axvline(x=(f_2_a*Mc)*1000,color='g',label='2-0')
                    ax.axvspan((f_2_a*Mc)*1000-229, (f_2_a*Mc)*1000+229, alpha=0.3, color='yellow')
                    ax.axvline((f_s_a*Mc)*1000,color='orange',label='spiral')
                    ax.axvspan((f_s_a*Mc)*1000-286, (f_s_a*Mc)*1000+286, alpha=0.3, color='cyan')
                    ax.axvline((f_0_a*Mc)*1000,linestyle="--",color='grey',label='2+0')
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


                    fc = f_p_a
                    dt=(tim[1]-tim[0])*Time*1000
                    band = 2.5
                    wavelet = 'cmor'+str(band)+'-'+str(fc)
                    widths = fc/np.linspace(fc-1.0, fc+1.0, 400)/dt
                    cwtmatr, freqs = pywt.cwt(dat, widths, wavelet, dt)
                    power = abs(cwtmatr)

                    fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
                    ax.pcolormesh(tim, freqs, power,cmap='jet')
                    plt.savefig('results/spec/'+name+'.jpg')
                    plt.close()



                except OSError:
                    pass
