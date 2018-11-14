#!/usr/bin/python

import os
#from matplotlib.pyplot import ion
import matplotlib.pyplot as p
from numpy import *
import numpy as np
from operator import itemgetter
#import os.path
import mpfit
import math
from scipy.stats import f as scipyf
import subprocess
import sys

######### FUNCTIONS FOR THE MAXIMUM #########


def fit(arc,tiempos):
    
    t1 =tiempos[0]#float(raw_input('lim time1: '))
    t2 =tiempos[1]#float(raw_input('lim time2: '))
    
    print('Close the plot')
    p.show(block=True)
    ptos=0

    for i,a in enumerate(dias):
        if a>=t1 and a<=t2:
            ptos+=1
                
    at = np.zeros(ptos)
    att = []
    am = np.zeros(ptos)
    er = np.zeros(ptos)
    j = 0
    
    for i,a in enumerate(dias):
        if a>=t1 and a<=t2:
            att.append((dias[i],magnitud[i],error[i]))

    att.sort(key=itemgetter(0))
    for i in range(len(att)):
        tupla = att[i]
        at[i] = tupla[0]
        am[j] = tupla[1]
        er[j] = tupla[2]
        j+=1
    
    k=0

    for j in range(2,6):
        if ptos>=j+2:
            try:
                aj = polyfit(at,am,j)
                po = poly1d(aj)
                eam = max(error)
                ymaxh = max(am) + eam + 0.05
                yminh = min(am) - eam - 0.05
                x = linspace(at[0]-5,at[len(at)-1]+5,100)
                p.plot(at,am,'.',x,po(x),'-',label='p'+str(j))
                #p.gca().invert_yaxis()
                #p.set_ylim((ymaxh,yminh))
                ax = p.gca()
                ax.set_ylim((ymaxh,yminh))
            except ValueError:
                print('It is not possible to make a polynomial grade'+str(j))
                k = k+1

    if k==3:
        return 0
    else:
        p.legend()
        p.show(block=False)
        maxi = []
        good=0
        des = raw_input('polynomial you choose? ')
        if des == 'p2':
            aj = polyfit(at,am,2) 
            der = polyder(aj)
            r = roots(der)
            for raiz in r:
                val = polyval(aj,raiz)
                maxi.append([raiz,val])
        elif des == 'p3':
            aj = polyfit(at,am,3)
            der = polyder(aj)
            r = roots(der)
            for raiz in r:
                val = polyval(aj,raiz)
                maxi.append([raiz,val])
            good=1
        elif des == 'p4':
            aj = polyfit(at,am,4)
            der = polyder(aj)
            r = roots(der)
            for raiz in r:
                val = polyval(aj,raiz)
                maxi.append([raiz,val])
            good=1
        elif des == 'p5':
            aj = polyfit(at,am,5)
            der = polyder(aj)
            r = roots(der)
            for raiz in r:
                val = polyval(aj,raiz)
                maxi.append([raiz,val])
            good=1
        else:
            print('Invalid')

    print(maxi)
    s=0.
    for err in er:
        s+=err
    erro = s/(len(er))
    
    print('Time, Maximum magnitude, Error:')
    print(np.real(maxi[good][1]),np.real(maxi[good][0]),erro)
    p.show(block=True)
    
    return [np.real(maxi[good][1]),np.real(maxi[good][0]),erro]

########## FUNCTIONS TO FIT ##########
def oneline(x,p2):
    i = 0
    v1 = np.zeros(len(x))
    for i in range(len(x)):
        v1[i] = p2[0]*x[i]+p2[1]
    return v1

def onelinee(x,p1,p2):
    return p1*x+p2

def twoline(x,p1):
    v = np.zeros(len(x))
    i = 0
    for i in range(len(x)):
        if x[i]<p1[3]:
            v[i] = p1[0]*x[i] + p1[2]
        elif x[i]>=p1[3]:
            v[i] = p1[1]*x[i] + p1[3]*(p1[0]-p1[1]) + p1[2]
    return v

def twolinee(x,p1,p2,p3,p4):
    if x<p4:
        return p1*x+p3
    elif x>=p4:
        return p2*x+p4*(p1-p2)+p3

def olivares1(x,pt):
    i = 0
    f = np.zeros(len(x))
    for i in range(len(x)):
        fd = -pt[0]/(1.+np.exp((x[i]-pt[1])/pt[2]))
        l = pt[3]*(x[i]-pt[1])+pt[4]
        f[i] = fd+l
        #f[i] = 10.**(-0.4*y)
    return f

def olivares2(x,pt):
    i = 0
    f = np.zeros(len(x))
    for i in range(len(x)):
        g = -pt[5]*math.exp(-(((x[i]-pt[6])/pt[7])**2))
        fd = -pt[0]/(1.+math.exp((x[i]-pt[1])/pt[2]))
        l = pt[3]*(x[i]-pt[1])+pt[4]
        f[i] = g+fd+l
        #f[i] = 10.**(-0.4*y)
    return f

def olivares1es(x,pt):
    fd = -pt[0]/(1.+math.exp((x-pt[1])/pt[2]))
    l = pt[3]*(x-pt[1])+pt[4]
    f = fd+l
    #f = 10.**(-0.4*y)
    return f

def olivares2es(x,pt):
    g = -pt[5]*math.exp(-(((x-pt[6])/pt[7])**2))
    fd = -pt[0]/(1.+np.exp((x-pt[1])/pt[2]))
    l = pt[3]*(x-pt[1])+pt[4]
    f = fd+l+g
    #f = 10.**(-0.4*y)
    return f


######### FUNCTIONS TO FIT AND USE MPFIT ########
def myfunct(p1, fjac=None, x=None, y=None, err=None):
    model = twoline(x,p1)
    status=0
    return ([status,(y-model)/err])

def myfuncte0(p1, fjac=None, x=None, y=None):
    model = twoline(x,p1)
    status=0
    return ([status,y-model])

def myfunct1(p2, fjac=None, x=None, y=None, err=None):
    model = oneline(x,p2)
    status=0
    return ([status,(y-model)/err])

def myfunct1e0(p2, fjac=None, x=None, y=None):
    model = oneline(x,p2)
    status=0
    return ([status,y-model])

def myfunctoliv1(pt, fjac=None, x=None, y=None, err=None):
    model = olivares1(x,pt)
    status=0
    return([status,(y-model)/err])

def myfunctoliv2(pt, fjac=None, x=None, y=None, err=None):
    model = olivares2(x,pt)
    status=0
    return([status,(y-model)/err])

def chi2(p1,p2):
    s = 0.
    for x,y,yerr in Data:
        teoria = onelinee(x,p1,p2)
        if yerr==0.:
            s += (teoria-y)**2
        else:
            s += (teoria-y)**2/yerr
    return s

def chi22(p1,p2,p3,p4):
    s = 0.
    for x,y,yerr in Data:
        teoria = twolinee(x,p1,p2,p3,p4)
        if yerr==0.:
            s += (teoria-y)**2
        else:
            s += (teoria-y)**2/yerr
    return s

####### EXPLOSION DATES DICTIONARY ########
tablatexp = open('../data/_texp.txt','r')
listatexp = tablatexp.readlines() 
tablatexp.close()

texp = {}
for i in range(1,len(listatexp)):
    lexp2 = listatexp[i].split()
    if lexp2[1]=='nan':
        tiempo = 0.
    else:
        tiempo = float(lexp2[1])
    if lexp2[2]=='nan':
        error = 0.
    else:
        error = float(lexp2[2])
    tipo = lexp2[3]
    texp[lexp2[0]]=[tiempo,error,tipo]

######## FILE TO SAVE #########
tabla = open('tabla.txt','a')
#tabla.write('Sn          Tmax      Mmax   error  ts1    ts2    ts31   ts32   s      error  b      s1     error  s2     error  b1     s3     error  b3     Ttrans    error  Tpt       error  Mend   error  Tend      Mtail  error  Ttail     Comentarios\n')

######## LOOP #########

if len(sys.argv) < 2:
    sna = raw_input('Supernova? ')
else:
    sna = sys.argv[1]
list_files = subprocess.Popen('ls '+sna, stdout=subprocess.PIPE, shell=True)
print(list_files.communicate()[0])
bandaquiere = raw_input('band that you want? ')
print('Supernova: '+sna)
if sna in texp:
      exp = texp[sna]
else:
      exp = [0.,0.,'nan']
print('Texp: '+str(exp[0]))
files = os.listdir(sna)
os.chdir(sna)
#print(files) 
m=0

######## PLOT LIGHT CURVE #########
for b in files:
    n = len(sna)+5
    banda = b[n:]
    datam = []
        
    if banda==bandaquiere:

        f = open(b,'r')
        if os.stat(b).st_size == 0:
            print('Archivo vacio')
            break
        m = 1
        print(f)
        
        lista = f.readlines()
        
        diass = []
        dias = []
        magnitud = []
        error = []
        n = 0
        
        #for j in range(4,len(lista)):
        for j in range(len(lista)):
            lista2 = lista[j].split()
            if lista2[2]=='None':
                diass.append((float(lista2[0]),float(lista2[1]),0.))
                n = 1
            else:
                diass.append((float(lista2[0]),float(lista2[1]),float(lista2[2])))
        
        diass.sort(key=itemgetter(0))
        for j in range(len(diass)):
           print(diass[j])
           dias.append(diass[j][0])
           magnitud.append(diass[j][1])
           if np.isnan(diass[j][2]):
               error.append(0.1)
           else:
               error.append(diass[j][2])
        
        ea = max(error)
        ymax = max(magnitud) + ea + 0.05
        ymin = min(magnitud) - ea - 0.05
        
        if n == 1:
            print('There are no numbers for data errors')
        
        f.close()

        
    else:
        m = m

if np.isfinite(exp[0]) == False: exp=[dias[0]]
print(exp)
if m==0:
    print('Not available band')
    exit()

banda = bandaquiere

if len(sys.argv) > 2:
    print(sys.argv[2])
    select=np.int64(sys.argv[2])
else:
    select=2
#select=0  #for maximum
#select=1  #for s123
#select=2  #for tPT

########## MAXIMUM MAGNITUDES ###########

if select==0:
    magmax = magnitud[0]
    tmagmax = dias[0]
    magmaxe = error[0]
    print(tmagmax, magmax, magmaxe)

    f,(ax1)=p.subplots(1,sharex=True,sharey=False)
    f.subplots_adjust(hspace=0.2)
    ax1.errorbar(np.array(dias)-exp[0],magnitud,yerr = error,fmt='.', label = '%s'%(banda))
    ax1.set_ylabel('Magnitud')
    ax1.set_xlabel('Dias')
    ax1.set_ylim((ymax,ymin))
    ax1.set_xlim((np.min(np.array(dias)-exp[0])-5,np.max(np.array(dias)-exp[0])+5))
    ax1.set_title(sna)
    ax1.legend()
    ax1.scatter(tmagmax-exp[0],magmax,color='red',s=60)#,fmt='.')
    p.show()
    #if raw_input('is it good?( |yes(y)|,no(n)) ')!='n':

    if raw_input('is it good?( |yes(y)|,no(n)) ') =='n':

        f,(ax1)=p.subplots(1,sharex=True,sharey=False)
        f.subplots_adjust(hspace=0.2)
        ax1.errorbar(dias,magnitud,yerr = error,fmt='.', label = '%s'%(banda))
        ax1.set_ylabel('Magnitud')
        ax1.set_xlabel('Dias')
        ax1.set_ylim((ymax,ymin))
        ax1.set_xlim((np.min(dias)-5,np.min(dias)+60))
        ax1.set_title(sna)
        ax1.legend()
        ax1.scatter(tmagmax,magmax,color='red',s=60)#,fmt='.')


        tiempos=[]
        def onclick(event):
            global tiempos
        
            if event.button==1:
                tiempos.append(event.xdata)
                print(str(event.xdata))
            if event.button==2:
                del tiempos[-1]
            if event.button==3:
                tiempos.append(0)
                tiempos.append(0)
        cid = f.canvas.mpl_connect('button_press_event', onclick)
        p.show()
        res = fit(b,tiempos)
        if res == 0:
            print('Cannot be fit with polyfit')
            magmax = nan
            tmagmax = nan
            magmaxe = nan
        else:
            magmax = res[0]
            tmagmax = res[1]
            magmaxe = res[2]
        
        if raw_input('is it good?( |yes(y)|,no(n)) ') =='n':
        
            f,(ax1)=p.subplots(1,sharex=True,sharey=False)
            f.subplots_adjust(hspace=0.2)
            ax1.errorbar(dias,magnitud,yerr = error,fmt='.', label = '%s'%(banda))
            ax1.set_ylabel('Magnitud')
            ax1.set_xlabel('Dias')
            ax1.set_ylim((magmax+0.5,magmax-0.5))
            ax1.set_xlim((np.min(dias)-5,np.min(dias)+160))
            ax1.set_title(sna)
            ax1.legend()
            ax1.scatter(tmagmax,magmax,color='green',s=60)#,fmt='.')
    
            tiempos=[]
            def onclick(event):
                global tiempos
                if event.button==1:
                    tiempos=(event.xdata,event.ydata)
                    print(tiempos)
            cid = f.canvas.mpl_connect('button_press_event', onclick)
            p.show()
            tmagmax = tiempos[0]
            magmax = tiempos[1]
            magmaxe = 0.1
#    if mgmax == 1:
#        if raw_input('fit?(yes,|no|) ')=='yes':
#            rres = fit(b)
#            if rres == 0:
#                print('cannot fit with polyfit')
#            
#            tmax = float(raw_input('Tmax2: '))
#            max2 = float(raw_input('Mmax2: '))
#            errmax2 = float(raw_input('error: '))
#        else:
#            tmax2 = float(raw_input('Tmax2: '))
#            max2 = float(raw_input('Mmax2: '))
#            errmax2 = float(raw_input('error: '))
#    else:
#        tmax2 = nan
#        max2 = nan
#        errmax2 = nan
else:
    magmax = nan
    tmagmax = nan
    magmaxe = nan
tmax2 = nan
max2 = nan
errmax2 = nan
p.close()


if select ==1:
    #replaced with the following piece of code from Thomas:
    f,(ax1)=p.subplots(1,sharex=True,sharey=False)
    f.subplots_adjust(hspace=0.2)
    ax1.errorbar(dias,magnitud,yerr = error,fmt='.', label = '%s'%(banda))
    ax1.set_ylabel('Magnitud')
    ax1.set_xlabel('Dias')
    ax1.set_ylim((ymax,ymin))
    ax1.set_xlim((np.min(dias)-5,np.max(dias)+5))
    ax1.set_title(sna)
    ax1.legend()
    print("Left click put point, right click 0 for each time, click right and left remove last value")
    tiempos=[]
    def onclick(event):
        global tiempos
        if event.button==1:
            tiempos.append(event.xdata)
            print('clicked!'+str(event.xdata))
        if event.button==2:
            del tiempos[-1]
        if event.button==3:
            tiempos.append(0)
            tiempos.append(0)
            tiempos.append(0)
            tiempos.append(0)
    cid = f.canvas.mpl_connect('button_press_event', onclick)
    p.show()
    p.close()
    ts1 = float(tiempos[0])
    ts2 = float(tiempos[1])
    ts31 = float(tiempos[2])
    ts32 = float(tiempos[3])
    print(round(ts1,3),round(ts2,3),round(ts31,3),round(ts32,3))
else:
    ts1 = nan
    ts2 = nan
    ts31 = nan
    ts32 = nan

if select ==1:
    print('Doing s1 and s2')
    #if raw_input('s1,s2? (|yes(y)|,no(n)) ') !='n':

    p.show(block=True)

    Data = []
    xx = []
    nn = 0

    for o,x in enumerate(dias):
        x1 = x-ts1
        xx.append(x1)
        if x>=ts1 and x<=ts2:
            nn += 1

    X = np.zeros(nn)
    Y = np.zeros(nn)
    Ye = np.zeros(nn)
    N = 0
    cero = 0

    for o,x in enumerate(dias):
        if x>=ts1 and x<=ts2:
            X[N] = x-ts1
            Y[N] = magnitud[o]
            if error[o]==0:
                cero = 1
            else:
                Ye[N] = error[o]
            N+=1
            Data.append((x-ts1,magnitud[o],error[o]))

    if cero==0:
        
        p00 = [0.02,16.]
        fa = {'x':X,'y':Y,'err':Ye}
        fab = {'x':X,'y':Y,'err':Ye}

        mm = mpfit.mpfit(myfunct1, p00, functkw=fa)
        print('parameters= ',mm.params)

        st = mm.params[0]
        b = mm.params[1]

        p0 = [st, st, b, (ts2-ts1)*0.2]

        m = mpfit.mpfit(myfunct, p0, functkw=fab)
        print('parameters= ',m.params)
        ttrans = m.params[3]+ts1
        print('ttrans= ',ttrans)

        tiem = linspace(0,ts2-ts1,1000)

        p.errorbar(xx,magnitud,yerr = error,fmt='k.', label = banda)
        p.plot(tiem,twoline(tiem,m.params),'-r',label='Fit dos')
        p.plot(tiem,oneline(tiem,mm.params),'-b',label='Fit una')
        p.ylim((ymax,ymin))
        p.legend()

    else:

        p00 = [0.02,16.]
        fa = {'x':X,'y':Y}
        fab = {'x':X,'y':Y}

        mm = mpfit.mpfit(myfunct1e0, p00, functkw=fa)
        print('parameters= ',mm.params)

        st = mm.params[0]
        b = mm.params[1]

        p0 = [st, st, b, (ts2-ts1)*0.2]

        m = mpfit.mpfit(myfuncte0, p0, functkw=fab)
        print('parameters= ',m.params)

        tiem = linspace(0,ts2-ts1,1000)

        p.errorbar(xx,magnitud,yerr = error,fmt='k.', label = banda)
        p.plot(tiem,twoline(tiem,m.params),'-r',label='Fit dos')
        p.plot(tiem,oneline(tiem,mm.params),'-b',label='Fit una')
        p.ylim((ymax,ymin))
        p.legend()
            
    p.show()

    N = nn # number of data points
    SSR1 = chi2(mm.params[0],mm.params[1]) # chi-square of oneline
    SSR2 = chi22(m.params[0],m.params[1],m.params[2],m.params[3]) # chi-square of twoline

    ftest = ((SSR1-SSR2)/2.)/(SSR2/(N-4.))
    print('oneline s2 = ',mm.params[0]*100.,'err s2= ',mm.perror[0]*100., 'b= ',mm.params[1])
    print('twoline s1= ', m.params[0]*100.,'err s1= ',m.perror[0]*100.,'twoline s2= ',m.params[1]*100.,'err s2= ',m.perror[1]*100., 'b= ', m.params[2])
    print('ttrans= ',m.params[3]+ts1,'err ttrans= ', m.perror[3])
    print('F-test= ',ftest, 'ddf= ',nn-4., 'p-value= ',scipyf.cdf(ftest, 2, nn-4))
    if scipyf.cdf(ftest, 2, nn-4) > 0.95:
        print('2-slopes better')
        good_slope=2
    else:
        print('1-slope better')
        good_slope=1

    if good_slope == 1:
        s1 = nan
        errs1 = nan
        s2 = mm.params[0]*100.
        errs2 = mm.perror[0]*100.
        b1 = mm.params[1]
        s = mm.params[0]*100.
        errs = mm.perror[0]*100.
        b = mm.params[1]
        ttrans = nan
        errttrans = nan
    elif good_slope == 2:
       s1 = m.params[0]*100.
       errs1 = m.perror[0]*100.
       b1 = m.params[2]
       s2 = m.params[1]*100.
       errs2 = m.perror[1]*100.
       s = mm.params[0]*100.
       errs = mm.perror[0]*100.
       b = mm.params[1]
       ttrans = m.params[3]+ts1
       errttrans = m.perror[3]
    else:
       print('Not an option')
else:
   s1 = nan
   errs1 = nan
   s2 = nan
   errs2 = nan
   b1 = nan
   s = nan
   errs = nan
   b = nan
   ttrans = nan
   errttrans = nan

########## s3 ###########
    
print('Doing s3')
if select ==1:
    #if raw_input('s3?(|yes(y)|,no(n)) ')!='n':

    nnn = 0
    cero=0
    xx = []
            
    for o,x in enumerate(dias):
        x1 = x-ts31
        xx.append(x1)
        if x>=ts31 and x<=ts32:
            nnn += 1

    X = np.zeros(nnn)
    Y = np.zeros(nnn)
    Ye = np.zeros(nnn)
    N = 0

    for o,x in enumerate(dias):
        if x>=ts31 and x<=ts32:
            X[N] = x-ts31
            Y[N] = magnitud[o]
            if error[o]==0:
                cero=1
            else:
                Ye[N] = error[o]
            N+=1

    if cero==0:
            
        fb = {'x':X,'y':Y,'err':Ye}
        p01 = [0.01, 13.5]
        m3 = mpfit.mpfit(myfunct1, p01, functkw=fb)

        tiem = linspace(0,ts32-ts31,1000)

        p.errorbar(xx,magnitud,yerr = error,fmt='k.', label = banda)
        p.plot(tiem,oneline(tiem,m3.params),'-b',label='Fit s3')
        p.ylim((ymax,ymin))
        p.legend()
        p.show()


    else:
            
        fb = {'x':X,'y':Y}
        p01 = [0.01, 13.5]
        m3 = mpfit.mpfit(myfunct1e0, p01, functkw=fb)

        tiem = linspace(0,ts32-ts31,1000)

        p.errorbar(xx,magnitud,yerr = error,fmt='k.', label = banda)
        p.plot(tiem,oneline(tiem,m3.params),'-b',label='Fit s3')
        p.ylim((ymax,ymin))
        p.legend()
        p.show()

    print('s3= ',m3.params[0]*100.,'err s3= ',m3.perror[0]*100., 'b3= ', m3.params[1])
    s3 = m3.params[0]*100.
    errs3 = m3.perror[0]*100.
    b3 = m3.params[1]
    
else:
    s3 =nan
    errs3=nan
    b3 = nan



### CALCULAR TPT, MEND, MTAIL
print('Doing tPT')
if select ==2:
    #if raw_input('Tpt?(|yes(y)|,no(n)) ')!='n':
    
    n = len(dias)
    dias1 = np.zeros(n)
    magnitud1 = np.zeros(n)
    error1 = np.zeros(n)
    for k in range(len(dias)):
        dias1[k] = dias[k] - exp[0]
        magnitud1[k] = magnitud[k]
        error1[k] = error[k]
        k+=1
    
    #print(dias, magnitud, error)
    #print(dias1, magnitud1, error1)
    fbb = {'x':dias1,'y':magnitud1,'err':error1}
    p0 = {'value':3., 'fixed':0,'limited':[1,0],'limits':[0.001,10.]} #step FD (a_0)
    p1 = {'value':90., 'fixed':0,'limited':[0,0],'limits':[40.,150.]} #middle of transition phase FD (t_PT)
    p2 = {'value':6., 'fixed':0,'limited':[1,0],'limits':[1.,12.]} #width of the transition phase FD (w_0)
    p3 = {'value':1.000, 'fixed':0,'limited':[1,1],'limits':[0.,3.]} #slope of radioactive decay  
    p4 = {'value':4., 'fixed':0,'limited':[1,0],'limits':[0.,0.]} #zero point at t=t_pt
    p5 = {'value':1.3, 'fixed':0,'limited':[1,1],'limits':[0.001,2.]} #gaussian height
    p6 = {'value':0., 'fixed':0,'limited':[1,1],'limits':[-15,10.]}#center of gaussian
    p7 = {'value':20., 'fixed':0,'limited':[1,0],'limits':[0.00001,3.]} #gaussian width
    pp = [p0,p1,p2,p3,p4,p5,p6,p7]
    vector = linspace(0,dias1[len(dias1)-1],1000)
    vector2 = linspace(min(magnitud1),max(magnitud1),100)
    f = mpfit.mpfit(myfunctoliv2, functkw=fbb, parinfo=pp)
##
    siga = f.perror[0]/(1+math.exp(-30./f.params[2]))
    sigtpt = 0.#f.params[0]*math.exp(-30./f.params[2])*f.perror[1]/(f.params[2]*(math.exp((f.params[1]-30)/f.params[2])+math.exp(-30./f.params[2])))
    sigw =f.params[0]*-30.*math.exp(-30./f.params[2])*f.perror[2]/((f.params[2]*(math.exp(-30./f.params[2])+1))**2)
##
    tpt = f.params[1]+exp[0]
    errtpt = f.perror[1]
    w0 = f.params[2]
    ew0 = f.perror[2]
    a0 = f.params[0]
    ea0 = f.perror[0]
    mtail = olivares2es(f.params[1]+3.0*w0,f.params)
    mtailerr =math.sqrt(((f.params[1]+3.0*w0)*f.perror[3])**2 + f.perror[4]**2)
    mend = olivares2es(f.params[1]-3.0*w0,f.params)
    menderr= math.sqrt(siga**2 + sigtpt**2 + sigw**2)
    vectortpt = np.zeros(len(vector2))
    for m in range(len(vector2)):vectortpt[m]=f.params[1]
##
    print('tpt= ', tpt, 'errtpt= ', errtpt)
    print('mend= ', mend, 'errmend= ', menderr, 'tend= ', tpt-3.0*w0)
    print('mtail= ', mtail, 'errmatil= ', mtailerr, 'ttail= ', tpt+3.0*w0)
    print('w0= ', f.params[2], 'w0err= ', f.perror[2])
    print('a0= ', f.params[0], 'a0err= ', f.perror[0])
    print('s3= ', f.params[3], 's3err= ', f.perror[3], '\n')
    #print('parameters and errors = ', '%e' %f.params[0], '%.2f' %f.perror[0], '%e' %f.params[1], '%.2f' %f.perror[1], '%e' %f.params[2], '%.2f' %f.perror[2], '%e' %f.params[3], '%.2f' %f.perror[3], '%e' %f.params[4],  '%.2f' %f.perror[4],  '%e' %f.params[5], '%.2f' %f.perror[5],  '%e' %f.params[6], '%.2f' %f.perror[6],  '%e' %f.params[7], '%.2f' %f.perror[7])
    #print('tinicio=', exp[0])
    
    p.errorbar(dias1,magnitud1,yerr = error1,fmt='.k',label =banda)
    p.plot(vector,olivares2(vector,f.params),'-b',label= 'Fit Tpt')
    p.plot(vectortpt,vector2,'-r',label='Tpt')
    p.plot(f.params[1]-3.0*w0,mend,'og',label='Mend')
    p.plot(f.params[1]+3.0*w0,mtail,'og',label='Mtail')
    p.ylim((max(magnitud1)+max(error1)+0.5,min(magnitud1)-max(error1)-0.5))
    p.legend()
    p.title(sna)
    p.show()

    if raw_input('is it good?( |yes(y)|,no(n)) ')!='n':
        pass
    else:
        p.errorbar(dias1,magnitud1,yerr = error1,fmt='.k',label =banda)
        p.ylim((max(magnitud1)+max(error1)+0.5,min(magnitud1)-max(error1)-0.5))
        p.legend()
        p.title(sna)
        p.show(block=False)
        tinicio = float(raw_input('initial time? '))
        tfinal = float(raw_input('final time? '))
        p.show(block=True)
        uop = 0
        for k in range(len(dias1)):
            if dias1[k] >= tinicio and dias1[k]<=tfinal:
                uop+=1

        dias2 = np.zeros(uop)
        magnitud2 = np.zeros(uop)
        error2 = np.zeros(uop)

        u=0
        for k in range(len(dias1)):
            if dias1[k] >= tinicio and dias1[k]<=tfinal:
                dias2[u] = dias1[k]
                magnitud2[u] = magnitud1[k]
                error2[u] = error1[k]
                u+=1

        fbb = {'x':dias2,'y':magnitud2,'err':error2}
        p0 = {'value':3., 'fixed':0,'limited':[1,0],'limits':[0.001,10.]} #peldagno fd
        p1 = {'value':90., 'fixed':0,'limited':[1,0],'limits':[40.,150.]} #tiempo de transicion fd (tpt)
        p2 = {'value':6., 'fixed':0,'limited':[1,0],'limits':[1.,12.]} #ancho de la transicion FD
        p3 = {'value':1.0, 'fixed':0,'limited':[1,1],'limits':[0.,5.]} #decaimiento tau
        p4 = {'value':4., 'fixed':0,'limited':[1,0],'limits':[0.,0.]} #pendiente decaimiento
        pp = [p0,p1,p2,p3,p4]
        vector = linspace(0,dias1[len(dias1)-1],1000)
        vector2 = linspace(min(magnitud1),max(magnitud1),100)
        f = mpfit.mpfit(myfunctoliv1, functkw=fbb, parinfo=pp)
        siga = f.perror[0]/(1+math.exp(-30./f.params[2]))
        sigtpt = f.params[0]*math.exp(-30./f.params[2])*f.perror[1]/(f.params[2]*(math.exp((f.params[1]-30)/f.params[2])+math.exp(-30./f.params[2])))
        sigw =f.params[0]*-30.*math.exp(-30./f.params[2])*f.perror[2]/((f.params[2]*(math.exp(-30./f.params[2])+1))**2)

        tpt = f.params[1]+exp[0]
        w0 = f.params[2]
        ew0 = f.perror[2]
        a0 = f.params[0]
        ea0 = f.perror[0]
        errtpt = f.perror[1]
        mtail = olivares1es(f.params[1]+3.0*w0,f.params)
        mtailerr =math.sqrt(((f.params[1]+3.0*w0)*f.perror[3])**2 + f.perror[4]**2)
        mend = olivares1es(f.params[1]-3.0*w0,f.params)
        menderr= math.sqrt(siga**2 + sigtpt**2 + sigw**2)
        vectortpt = np.zeros(len(vector2))
        for m in range(len(vector2)):vectortpt[m]=f.params[1]

        print('tpt= ', tpt, 'errtpt= ', errtpt)
        print('mend= ', mend, 'errmend= ', menderr, 'tend= ', tpt-3.0*w0)
        print('mtail= ', mtail, 'errmatil= ', mtailerr, 'ttail= ', tpt+3.*w0)
        print('w0= ', f.params[2], 'w0err= ', f.perror[2])
        print('a0= ', f.params[0], 'a0err= ', f.perror[0])
        print('s3= ', f.params[3]*100., 's3err= ', f.perror[3]*100.,'\n')
        #print('parameters and errors = ', '%e' %f.params[0], '%.2f' %f.perror[0], '%e' %f.params[1], '%.2f' %f.perror[1], '%e' %f.params[2], '%.2f' %f.perror[2], '%e' %f.params[3], '%.2f' %f.perror[3], '%e' %f.params[4],  '%.2f' %f.perror[4])

        p.errorbar(dias1,magnitud1,yerr = error1,fmt='.k',label =banda)
        p.plot(vector,olivares1(vector,f.params),'-b',label= 'Fit Tpt')
        p.plot(vectortpt,vector2,'-r',label='Tpt')
        p.plot(f.params[1]-3.0*w0,mend,'og',label='Mend')
        p.plot(f.params[1]+3.0*w0,mtail,'og',label='Mtail')
        p.ylim((max(magnitud1)+max(error1)+0.5,min(magnitud1)-max(error1)-0.5))
        p.legend()
        p.title(sna)
        p.show()
        
        if raw_input('is it good?( |yes(y)|,no(n)) ')!='n':
            pass
        else:
            this,(ax1)=p.subplots(1,sharex=True,sharey=False)
            this.subplots_adjust(hspace=0.2)
            ax1.errorbar(dias1,magnitud1,yerr = error1,fmt='.k', label = '%s'%(banda))
            ax1.set_ylabel('Magnitud')
            p.ylim((max(magnitud1)+max(error1)+0.5,min(magnitud1)-max(error1)-0.5))
            ax1.set_xlabel('Dias')
            ax1.set_title(sna)
            ax1.legend()
            print("Left click put point, right click 0 for each time, click right and left remove last value")
            equits=[]
            def onclick(event):
                global tiempos
                equits.append(event.xdata)
                print(event.xdata)
            cid = this.canvas.mpl_connect('button_press_event', onclick)
            p.show()
            p.close()

            print(equits)
            fbb = {'x':dias2,'y':magnitud2,'err':error2}
            p0 = {'value':3., 'fixed':0,'limited':[1,0],'limits':[0.001,10.]} #peldagno fd
            p1 = {'value':equits[0]+(equits[1]-equits[0])/2., 'fixed':0,'limited':[1,0],'limits':[40.,150.]} #tiempo de transicion fd (tpt)
            p2 = {'value':(equits[1]-equits[0])/6., 'fixed':1,'limited':[1,0],'limits':[1.,12.]} #ancho de la transicion FD
            p3 = {'value':1.0, 'fixed':0,'limited':[1,1],'limits':[0.,5.]} #decaimiento tau
            p4 = {'value':4., 'fixed':0,'limited':[1,0],'limits':[0.,0.]} #pendiente decaimiento
            pp = [p0,p1,p2,p3,p4]
            vector = linspace(0,dias1[len(dias1)-1],1000)
            vector2 = linspace(min(magnitud1),max(magnitud1),100)
            f = mpfit.mpfit(myfunctoliv1, functkw=fbb, parinfo=pp)
            siga = f.perror[0]/(1+math.exp(-30./f.params[2]))
            sigtpt = f.params[0]*math.exp(-30./f.params[2])*f.perror[1]/(f.params[2]*(math.exp((f.params[1]-30)/f.params[2])+math.exp(-30./f.params[2])))
            sigw =f.params[0]*-30.*math.exp(-30./f.params[2])*f.perror[2]/((f.params[2]*(math.exp(-30./f.params[2])+1))**2)
            
            tpt = f.params[1]+exp[0]
            w0 = f.params[2]
            ew0 = f.perror[2]
            a0 = f.params[0]
            ea0 = f.perror[0]
            errtpt = f.perror[1]
            mtail = olivares1es(f.params[1]+3.0*w0,f.params)
            mtailerr =math.sqrt(((f.params[1]+3.0*w0)*f.perror[3])**2 + f.perror[4]**2)
            mend = olivares1es(f.params[1]-3.0*w0,f.params)
            menderr= math.sqrt(siga**2 + sigtpt**2 + sigw**2)
            vectortpt = np.zeros(len(vector2))
            for m in range(len(vector2)):vectortpt[m]=f.params[1]
            
            print('tpt= ', tpt, 'errtpt= ', errtpt)
            print('mend= ', mend, 'errmend= ', menderr, 'tend= ', tpt-3.0*w0)
            print('mtail= ', mtail, 'errmatil= ', mtailerr, 'ttail= ', tpt+3.*w0)
            print('w0= ', f.params[2], 'w0err= ', f.perror[2])
            print('a0= ', f.params[0], 'a0err= ', f.perror[0])
            print('s3= ', f.params[3], 's3err= ', f.perror[3],'\n')
            #print('parameters and errors = ', '%e' %f.params[0], '%.2f' %f.perror[0], '%e' %f.params[1], '%.2f' %f.perror[1], '%e' %f.params[2], '%.2f' %f.perror[2], '%e' %f.params[3], '%.2f' %f.perror[3], '%e' %f.params[4],  '%.2f' %f.perror[4])
            
            p.errorbar(dias1,magnitud1,yerr = error1,fmt='.k',label =banda)
            p.plot(vector,olivares1(vector,f.params),'-b',label= 'Fit Tpt')
            p.plot(vectortpt,vector2,'-r',label='Tpt')
            p.plot(f.params[1]-3.0*w0,mend,'og',label='Mend')
            p.plot(f.params[1]+3.0*w0,mtail,'og',label='Mtail')
            p.ylim((max(magnitud1)+max(error1)+0.5,min(magnitud1)-max(error1)-0.5))
            p.legend()
            p.title(sna)
            p.show()

else:
    tpt = nan
    errtpt =nan
    mtail = nan
    mtailerr = nan
    mend = nan
    menderr = nan                   
    w0 = nan
    ew0 = nan
    a0 = nan
    ea0 = nan
########## SAVE #############

ptmagmax = '%.2f' %tmagmax
pmagmax = '%.2f' %magmax
pmagmaxe = '%.2f' %magmaxe
ps1 = '%.2f' %s1
perrs1 = '%.2f' %errs1
ps2 = '%.2f' %s2
perrs2 = '%.2f' %errs2
pb1 = '%.2f'%b1
ps = '%.2f' %s
perrs = '%.2f' %errs
pb = '%.2f' %b
pttrans = '%.2f' %ttrans
perrttrans = '%.2f' %errttrans
ps3 = '%.2f' %s3
perrs3 = '%.2f' %errs3
pb3 = '%0.2f' %b3
ptpt = '%.2f' %tpt
perrtpt = '%.2f' %errtpt
pmtail = '%.2f' %mtail
perrmtail = '%.2f' %mtailerr
pmend = '%.2f' %mend
perrmend = '%.2f' %menderr
w0= '%.2f' %w0
ew0= '%.2f' %ew0
a0= '%.2f' %a0
ea0= '%.2f' %ea0
pmax2 = '%.2f' %max2
pmax2err = '%.2f' %errmax2
ptmax2 = '%.2f' %tmax2
comen = raw_input('comments? ')
coment = tipo+','+comen 

#if raw_input('Save?(|yes|,no) ') !='n':
tabla.write(sna.ljust(12,' ')+ptmagmax.ljust(10,' ')+pmagmax.ljust(7,' ')+pmagmaxe.ljust(7,' ')+str(ts1).ljust(15,' ')+str(ts2).ljust(15,' ')+str(ts31).ljust(15,' ')+str(ts32).ljust(15,' ')+ps.ljust(7,' ')+perrs.ljust(7,' ')+pb.ljust(7,' ')+ps1.ljust(7,' ')+perrs1.ljust(7,' ')+ps2.ljust(7,' ')+perrs2.ljust(7,' ')+pb1.ljust(7,' ')+ps3.ljust(7,' ')+perrs3.ljust(7,' ')+pb3.ljust(7,' ')+pttrans.ljust(10,' ')+perrttrans.ljust(7,' ')+ptpt.ljust(10,' ')+perrtpt.ljust(10,' ')+pmend.ljust(10,' ')+perrmend.ljust(10,' ')+pmtail.ljust(10,' ')+perrmtail.ljust(10,' ')+w0.ljust(10,' ')+ew0.ljust(10,' ')+a0.ljust(10,' ')+ea0.ljust(10,' ')+pmax2.ljust(7,' ')+pmax2err.ljust(7,' ')+ptmax2.ljust(10,' ')+coment+'\n')
tabla.close()
