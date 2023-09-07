# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 23:24:34 2022

@author: diego
"""

import numpy as np
import pylab as pl
import time
from datetime import date # ordenar los tiempos
from datetime import datetime, timedelta 
import math

def toYearFraction(date): # estas funciones pasan la fecha a decimal
    def sinceEpoch(date):  # returns seconds since epoch
        return time.mktime(date.timetuple())
    
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year + 1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed / yearDuration

    return date.year + fraction

def extract_data(ts_file, t_cut = False): # esta funcion extrae los datos
    name = ts_file[-8:-4]
    data = np.loadtxt(ts_file, usecols=(1, 2, 3, 4, 5, 6, 7, 8), dtype=float)
    year = data[:, 0]
    days = data[:, 1]
    dates = [datetime(year=int(yr), month=1, day=1) + timedelta(days=int(days[i])-1)
             for i, yr in enumerate(year)]
    dates = np.array(dates)
    t = np.array([toYearFraction(date_i) for date_i in dates])
    if t_cut != False:
        ind = np.where(t >= t_cut)
    else:
        ind = np.where(t >= t[0])
    ew = data[ind, 2].T.flatten()
    ns = data[ind, 3].T.flatten()
    up = data[ind, 4].T.flatten()
    sigm_ew = data[ind, 5].T.flatten() 
    sigm_ns = data[ind, 6].T.flatten() 
    sigm_up = data[ind, 7].T.flatten() 
    datas = np.vstack((t,ew, ns, up, sigm_ew, sigm_ns, sigm_up)).T

    return {'station' : name , 'data': datas}

datos =extract_data("EMAT.txt", False) # toma todos los datos como un array
print(datos['data'][:,1]) #primera fila de matriz
fig, axs = pl.subplots(3)
fig.suptitle('BN05')
axs[0].plot(datos['data'][:,0], datos['data'][:,1],color="cornflowerblue", linewidth=1.0, linestyle="dotted") 
axs[0].set_ylabel('E-W [mm]')
axs[0].set_xlabel('Fecha [años]')
# en axs[0] ploteamos primera columna de datos con la fecha que son los años
axs[1].plot(datos['data'][:,0], datos['data'][:,2],color="cornflowerblue", linewidth=1.0, linestyle="dotted")
axs[1].set_ylabel('N-S [mm]')
axs[2].plot(datos['data'][:,0], datos['data'][:,3],color="cornflowerblue", linewidth=1.0, linestyle="dotted")
axs[2].set_ylabel('Up [mm]')
axs[2].set_xlabel('Fecha [años]')
# para data[:,4 o 5 se ven cosas con picks como delta]
axs[0].grid(True)
#axs[0].axs(figsize=(8,6),facecolor='white')


## Primera función  LINEAL
c1=datos['data'][:,0] #entra mi vector de tiempo y devuelve matriz de diseño
ala=np.zeros((c1.shape[0],1)) # vector de zeros del largo de mis datos fecha
d=np.ones([c1.shape[0],1]) # vector de unos que van al comienzo de la matriz por def
for i in range(c1.shape[0]):
        ala[i,0]=c1[i]-c1[0]   
    #ala vector de delta t
    
G_trend=np.concatenate((d,ala),axis=1) # junto vector de unos y delta t en una sola matriz que es mi G_trend v cte (linea)

#G_trend será siempre el mismo por como se define la funcion de fechas mas arriba (el mismo para la misma base de datos)    

# ejemplo uso de funcion G_linea=gtrend(datos['data'][:,0])

# Luego debemos sacar el m_est (por minimos cuadrados con peso) para asi
# obtener el d_est que son los datos del modelo

def d_est(sigma,g,datos):
    # g mi matriz de diseño
    # datos-datos de la componente 
    # sigma el error asociado a esa componente #Ej sigma_ew=datos['data'][:,4] 
#sig
    sigma_1= np.diag(sigma**2)
    sigma1=np.linalg.inv(sigma_1) #calcula el inverso multiplicativo de una matriz
    first=np.dot(np.transpose(g),sigma1)
    parent=np.linalg.inv(np.dot(first,g))
    seg=np.dot(parent,np.transpose(g))
    ter=np.dot(seg,sigma1)
    #ew_date=datos['data'][:,1]
    m=np.dot(ter,datos) # este es el m_est
    e=np.dot(g,m) # multiplocacion G m_est y nos da d_station
    return e # nos retorna el d_est

#sigma_ew=datos['data'][:,4]
#dat=datos['data'][:,1] #E-W
# d_est(sigma_ew,G_trend,dat)

dat_ew=datos['data'][:,1]
dat_ns=datos['data'][:,2]
dat_up=datos['data'][:,3]

sigma_ew=datos['data'][:,4]
sigma_ns=datos['data'][:,5]
sigma_up=datos['data'][:,6]


ew=d_est(sigma_ew,G_trend,dat_ew)
ns=d_est(sigma_ns,G_trend,dat_ns)
up=d_est(sigma_up,G_trend,dat_up)



t=datos['data'][:,0] # es lo mismo que c1 de mas arriba pero lo cambio a t -_-
G_ss=np.zeros((t.shape[0],4))

for i in range(t.shape[0]):
    trigo=np.array([[math.sin(2*math.pi*t[i])],[math.cos(2*math.pi*t[i])],[math.sin(4*math.pi*t[i])],[math.cos(4*math.pi*t[i])]])          
    for j in range(4):
        G_ss[i,j]=trigo[j]

#Estos son los estacionales

ew1=d_est(sigma_ew,G_ss,dat_ew)
ns1=d_est(sigma_ns,G_ss,dat_ns)
up1=d_est(sigma_up,G_ss,dat_up)



# print(datos['data'][:,1]) #primera fila de matriz
# fig, axs = pl.subplots(3)
# fig.suptitle('LVI1')

# axs[0].plot(datos['data'][:,0], datos['data'][:,1],color="cornflowerblue", linewidth=1.0, linestyle="dotted") 
# axs[0].plot(datos['data'][:,0], ew,color="red", linewidth=1.0) 
# axs[0].plot(datos['data'][:,0], ew1,color="m", linewidth=1.0)
# axs[0].set_ylabel('E-W [mm]')
# axs[0].set_xlabel('Fecha [años]')
# # en axs[0] ploteamos primera columna de datos con la fecha que son los años
# axs[1].plot(datos['data'][:,0], datos['data'][:,2],color="cornflowerblue", linewidth=1.0, linestyle="dotted")
# axs[1].plot(datos['data'][:,0], ns,color="red", linewidth=1.0)
# axs[1].plot(datos['data'][:,0], ns1,color="m", linewidth=1.0)
# axs[1].set_ylabel('N-S [mm]')

# axs[2].plot(datos['data'][:,0], datos['data'][:,3],color="cornflowerblue", linewidth=1.0, linestyle="dotted")
# axs[2].plot(datos['data'][:,0], up,color="red", linewidth=1.0)
# axs[2].plot(datos['data'][:,0], up1,color="m", linewidth=1.0)

# axs[2].set_ylabel('Up [mm]')
# axs[2].set_xlabel('Fecha [años]')

#Ahora SALTOS

#SALTOS y saltos por antena
#aqui debemos saber la fecha exacta del sismo o cambio de antena
t_illapel=toYearFraction(datetime(2015,9,16))
x=np.where(t == t_illapel) #aqui esto no sirve xq no tiene esa fecha
t_27f=toYearFraction(datetime(2010,2,27))
x2=np.where(t == t_27f)

h0=np.zeros((t.shape[0],1))
for i in range(t.shape[0]):
    if i>=829:   # el uno es la posicion de mi salto 
        h0[i]=1
        #primer salto posicion 1695
    # segundo salto 2261

h1=np.zeros((t.shape[0],1))
for i in range(t.shape[0]):
    if i>=1768:   # el uno es la posicion de mi salto 
        h1[i]=1
        #primer salto posicion 1695
    # segundo salto 2261

# h2=np.zeros((t.shape[0],1))
# for i in range(t.shape[0]):
#     if i>=2262:   # el uno es la posicion de mi salto 
#         h2[i]=1
#         #primer salto posicion 1695
#     #segundo salto 2261
# salto=np.concatenate((h1,h2),axis=1)
G_salto=np.concatenate((h0,h1),axis=1)


ew2=d_est(sigma_ew,G_salto,dat_ew)
ns2=d_est(sigma_ns,G_salto,dat_ns)
up2=d_est(sigma_up,G_salto,dat_up)



# G_total0=np.concatenate((G_trend,G_ss,G_salto),axis=1)
#PARA IR PROBANDO pero falta G_post

# ew_total=d_est(sigma_ew,G_total0,dat_ew)
# ns_total=d_est(sigma_ns,G_total0,dat_ns)
# up_total=d_est(sigma_up,G_total0,dat_up)



## Ahora le agregamos POST-SISMICO

t=datos['data'][:,0]
l0=np.zeros((t.shape[0],1))
for i in range(t.shape[0]):
    if i>=829:   # el uno es la posicion de mi salto
        tau=5/365.25 #con 10 funciona perfect, con 1 funciona mjr, mientras mas pequeño el numerador crece pick de N-S
        l0[i]=np.log(1+(t[i]-t[829])/tau)

l1=np.zeros((t.shape[0],1))
for i in range(t.shape[0]):
    if i>=1768:   # el uno es la posicion de mi salto
        tau=1000/365.25 #con 10 funciona perfect, con 1 funciona mjr, mientras mas pequeño el numerador crece pick de N-S
        l1[i]=np.log(1+(t[i]-t[1727])/tau)


G_post=np.concatenate((l0,l1),axis=1)


ew3=d_est(sigma_ew,G_post,dat_ew)
ns3=d_est(sigma_ns,G_post,dat_ns)
up3=d_est(sigma_up,G_post,dat_up)


fig, axs = pl.subplots(3)
fig.suptitle('Estación EMAT',fontsize=14)

#axs[0].plot(datos['data'][:,0], datos['data'][:,1],color="cornflowerblue", linewidth=1.0, linestyle="dotted") 
axs[0].plot(datos['data'][:,0], ew,color="red", linewidth=1.0,label='Lineal') 
axs[0].plot(datos['data'][:,0], ew1,color="m", linewidth=1.0,label='Estacional')
axs[0].plot(datos['data'][:,0], ew2,color="k", linewidth=1.0,label='Saltos')
axs[0].plot(datos['data'][:,0], ew3,color="c", linewidth=1.0,label='Post')
axs[0].plot(t_illapel,0,marker='*', ms=6, mfc='yellow',mec='k',label='Terremoto')
axs[0].plot(t_27f,0,marker='*', ms=6, mfc='yellow',mec='k')
#ms-tamaño simbolo, mfc-color simbolo,mec-colorancho #buscar verde claro lightgreen
axs[0].set_ylabel('E-W [mm]',fontsize=9)
axs[0].set_xlabel('Fecha [años]',fontsize=9)
axs[0].set_title('Componentes Modelo de Trayectoria [E-W]',fontsize=10)
#axs[0].legend(fontsize=5)

# en axs[0] ploteamos primera columna de datos con la fecha que son los años
#axs[1].plot(datos['data'][:,0], datos['data'][:,2],color="cornflowerblue", linewidth=1.0, linestyle="dotted")
axs[1].plot(datos['data'][:,0], ns,color="red", linewidth=1.0)
axs[1].plot(datos['data'][:,0], ns1,color="m", linewidth=1.0)
axs[1].plot(datos['data'][:,0], ns2,color="k", linewidth=1.0)
axs[1].plot(datos['data'][:,0], ns3,color="c", linewidth=1.0)
axs[1].plot(t_illapel,0,marker='*', ms=6, mfc='yellow',mec='k') #plot en la fecha de illapel,0
axs[1].plot(t_27f,0,marker='*', ms=6, mfc='yellow',mec='k') #plot en la fecha de illapel,0

axs[1].set_ylabel('N-S [mm]',fontsize=9)
axs[1].set_xlabel('Fecha [años]',fontsize=9)
axs[1].set_title('Componentes Modelo de Trayectoria [N-S]',fontsize=10)

#axs[2].plot(datos['data'][:,0], datos['data'][:,3],color="cornflowerblue", linewidth=1.0, linestyle="dotted")
axs[2].plot(datos['data'][:,0], up,color="red", linewidth=1.0)
axs[2].plot(datos['data'][:,0], up1,color="m", linewidth=1.0)
axs[2].plot(datos['data'][:,0], up2,color="k", linewidth=1.0)
axs[2].plot(datos['data'][:,0], up3,color="c", linewidth=1.0)
axs[2].plot(t_illapel,0,marker='*', ms=6, mfc='yellow',mec='k')
axs[2].plot(t_27f,0,marker='*', ms=6, mfc='yellow',mec='k')

axs[2].set_xlabel('Fecha [años]',fontsize=9)
axs[2].set_ylabel('Up [mm]',fontsize=9)
axs[2].set_title('Componentes Modelo de Trayectoria [Up]',fontsize=10)

#pl.axis('tight')
pl.tight_layout()
fig.legend(fontsize=7)
#pl.show()
#pl.grid() # se le puede poner true dentro y cambia la grilla
#pl.legend( fontsize=15)
#plt.figure(figsize=(8,6),facecolor='white')















#Ahora vamos a probar todos juntos (MODELO COMPLETO)
G_total=np.concatenate((G_trend,G_ss,G_salto,G_post),axis=1)

ew_total=d_est(sigma_ew,G_total,dat_ew)
ns_total=d_est(sigma_ns,G_total,dat_ns)
up_total=d_est(sigma_up,G_total,dat_up)

fig, axs = pl.subplots(3)
fig.suptitle('Estación EMAT',fontsize=14)

axs[0].plot(datos['data'][:,0], datos['data'][:,1],'b.',markersize=2,label='Datos') 
axs[0].plot(datos['data'][:,0], ew_total,color="red", linewidth=1,label='Modelo') 
axs[0].plot(t_illapel,0,marker='*', ms=6, mfc='c',mec='c',label='Terremoto')
axs[0].plot(t_27f,0,marker='*', ms=6, mfc='c',mec='c')

#ms-tamaño simbolo, mfc-color simbolo,mec-colorancho #buscar verde claro lightgreen
axs[0].set_ylabel('E-W [mm]',fontsize=9)
#axs[0].set_xlabel('Fecha [años]',fontsize=9)
axs[0].set_title('Ajuste Modelo de Trayectoria [E-W]',fontsize=10)
#axs[0].legend(loc='upper right',fontsize=5)

# en axs[0] ploteamos primera columna de datos con la fecha que son los años
axs[1].plot(datos['data'][:,0], datos['data'][:,2],'b.',markersize=2)
axs[1].plot(datos['data'][:,0],ns_total ,color="red", linewidth=1)
axs[1].set_ylabel('N-S [mm]')
axs[1].plot(t_illapel,0,marker='*', ms=6, mfc='c',mec='c')
axs[1].plot(t_27f,0,marker='*', ms=6, mfc='c',mec='c')
#ms-tamaño simbolo, mfc-color simbolo,mec-colorancho #buscar verde claro lightgreen
#axs[1].set_xlabel('Fecha [años]',fontsize=9)
axs[1].set_title('Ajuste Modelo de Trayectoria [N-S]',fontsize=10)
#axs[1].legend(fontsize=5)

axs[2].plot(datos['data'][:,0], datos['data'][:,3],'b.',markersize=2)
axs[2].plot(datos['data'][:,0], up_total,color="red", linewidth=1.0)

axs[2].set_ylabel('Up [mm]')
axs[2].set_xlabel('Fecha [años]')
axs[2].plot(t_illapel,155,marker='*', ms=6, mfc='c',mec='c')
axs[2].plot(t_27f,155,marker='*', ms=6, mfc='c',mec='c')
axs[2].set_title('Ajuste Modelo de Trayectoria [Up]',fontsize=10)
#axs[2].legend(fontsize=5)
pl.tight_layout()
fig.legend(fontsize=8)
pl.figure(figsize=(6, 12), dpi=80)

#Grafico RESIDUALES--

res_ew=dat_ew-ew_total
res_ns=dat_ns-ns_total
res_up=dat_up-up_total


fig, axs = pl.subplots(3,2)
fig.suptitle('Estación EMAT',fontsize=14)

axs[0,0].plot(datos['data'][:,0], datos['data'][:,1],'c.',markersize=2,label='Datos') 
axs[0,0].plot(datos['data'][:,0], ew_total,color="brown", linewidth=1,label='Modelo') 
axs[0,0].plot(t_illapel,0,marker='*', ms=6, mfc='c',mec='c',label='Terremoto')
#ms-tamaño simbolo, mfc-color simbolo,mec-colorancho #buscar verde claro lightgreen
axs[0,0].set_ylabel('E-W [mm]',fontsize=9)
#axs[0].set_xlabel('Fecha [años]',fontsize=9)
axs[0,0].set_title('Modelo de Trayectoria',fontsize=10)

axs[0,0].legend(loc='upper right',fontsize=6)

# en axs[0] ploteamos primera columna de datos con la fecha que son los años
axs[1,0].plot(datos['data'][:,0], datos['data'][:,2],'c.',markersize=2)
axs[1,0].plot(datos['data'][:,0],ns_total ,color="brown", linewidth=1)
axs[1,0].set_ylabel('N-S [mm]')
axs[1,0].plot(t_illapel,0,marker='*', ms=6, mfc='c',mec='c')
#ms-tamaño simbolo, mfc-color simbolo,mec-colorancho #buscar verde claro lightgreen
#axs[1].set_xlabel('Fecha [años]',fontsize=9)
#axs[1].set_title('Ajuste Modelo de Trayectoria [N-S]',fontsize=10)
#axs[1].legend(fontsize=5)

axs[2,0].plot(datos['data'][:,0], datos['data'][:,3],'c.',markersize=2)
axs[2,0].plot(datos['data'][:,0], up_total,color="brown", linewidth=1.0)

axs[2,0].set_ylabel('Up [mm]')
axs[2,0].set_xlabel('Fecha [años]')
axs[2,0].plot(t_illapel,0,marker='*', ms=6, mfc='c',mec='c')
#axs[2].set_title('Ajuste Modelo de Trayectoria [Up]',fontsize=10)

#Columna 2
axs[0,1].plot(datos['data'][:,0], res_ew,'k.',markersize=2,label='Residuales') 
#axs[0,1].plot(datos['data'][:,0], ew_total,color="red", linewidth=1,label='Modelo') 
axs[0,1].plot(t_illapel,0,marker='*', ms=6, mfc='c',mec='c',label='Terremoto')
#ms-tamaño simbolo, mfc-color simbolo,mec-colorancho #buscar verde claro lightgreen
axs[0,1].set_ylabel('E-W [mm]',fontsize=9)
#axs[0].set_xlabel('Fecha [años]',fontsize=9)
axs[0,1].set_title('Residuales',fontsize=10)
axs[0,1].legend(loc='upper right',fontsize=6)

#axs[0].legend(loc='upper right',fontsize=5)

# en axs[0] ploteamos primera columna de datos con la fecha que son los años
axs[1,1].plot(datos['data'][:,0], res_ns,'k.',markersize=2)
#axs[1,1].plot(datos['data'][:,0],ns_total ,color="red", linewidth=1)
axs[1,1].set_ylabel('N-S [mm]')
axs[1,1].plot(t_illapel,0,marker='*', ms=6, mfc='c',mec='c')
#ms-tamaño simbolo, mfc-color simbolo,mec-colorancho #buscar verde claro lightgreen
#axs[1].set_xlabel('Fecha [años]',fontsize=9)
#axs[1].set_title('Ajuste Modelo de Trayectoria [N-S]',fontsize=10)
#axs[1].legend(fontsize=5)

axs[2,1].plot(datos['data'][:,0], res_up,'k.',markersize=2)
#axs[2,1].plot(datos['data'][:,0], up_total,color="red", linewidth=1.0)

axs[2,1].set_ylabel('Up [mm]')
axs[2,1].set_xlabel('Fecha [años]')
axs[2,1].plot(t_illapel,0,marker='*', ms=6, mfc='c',mec='c')


#axs[2].set_title('Ajuste Modelo de Trayectoria [Up]',fontsize=10)
pl.tight_layout()
pl.figure(figsize=(15, 12), dpi=100)


