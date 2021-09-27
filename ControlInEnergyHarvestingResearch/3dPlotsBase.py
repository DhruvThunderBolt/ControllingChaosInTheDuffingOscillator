import numpy as np
from scipy.integrate import odeint
'''import scipy.integrate as integrate'''
import matplotlib.pyplot as plt
import matplotlib
import math
from mpl_toolkits import mplot3d
import statistics
import sys
import operator
import collections
import time


font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 15}

matplotlib.rc('font', **font)

#%% Setting parameters
# Andre's regime with varying g
gam = 0.02
om = 1.0
kappa = 0.5
xi = 0.05
lam = 0.05

C = xi/kappa
R = kappa/lam/xi


# om = 1.3# these values are for happy whale
# xi = 0.2
# C = 1
# R = 0.98
# g=0.57
# gam = 0.125

#%% Define function integrated


def Harvest(x,t):

    # Assigning symbols to elements of x for ease
    q = x[0] #q is first state
    p = x[1]
    v = x[2]
    Edrive = x[3]
    Ediss = x[4]

    # Renaming some parameters -- but this seems like a re-write, see 25 onwards
    # lam = 1/(R*C)
    # xi = C*(xi/C)
    #Symbolic equations of motion
    dqdt = p
    dpdt = -q**3+q + g*math.cos(om*t)-2*gam*p+xi*v
    dvdt = -v/(R*C)-(xi/C)*p
    #Some other time integrals we want to compute
    dEdrive_dt = g*math.cos(om*t)*p #energy absorbed from drive
    dEdiss_dt = 2*gam*p*p # energy disspated through Gamma
    #    dEC = C*v*(-lam*v-(xi/C)*p) #energy in the capacitor
    dER = v*v/R #Energy in the resistor
    #    dOsc = p*(-q * q * q +  q + g * math.cos(om * t) - 2 * gam * p + xi * v) + (
    #                - q * q * q + q) * p # Energy IN the oscillator
    return [dqdt,dpdt,dvdt,dEdrive_dt,dEdiss_dt, dER]

#%% Creating empty lists to save values
dissList = []
driveList = []
effList = []
harvestList = []
fig = plt.figure(figsize=plt.figaspect(0.5))
testList = [26017, 27000, 28570]
loop = 1
for i in testList:
    g = i/10000.0
    # R=i/100

    t0 = time.time() #Calling computer clock

    x0 = [1,0,0,0,0,0] #Initial values. Change here.

    totPoints = 1000000
    periods = 1000
    pointPerPeroid = totPoints/periods

    t = np.linspace(0, periods*(2*math.pi)/om, totPoints)
    #Going from 0 to 50000, with 100000 points of data.

    x = odeint(Harvest, x0, t) #(Function, initial condition, time)
    # np.save('BaseGPeriod19', x)
    # print("saved")
    #Going to produce five columns of results, first colum is q, then p , then v
    numOfPoints = 980000 #This is the transient number of points to be rejected right here

    q = x[:,0][numOfPoints:] #Starting FROM numOfPoints
    p = x[:,1][numOfPoints:]
    v = x[:,2][numOfPoints:]
    Edrive = x[:,3][numOfPoints:]
    Ediss = x[:,4][numOfPoints:]
    # Ecap = x[:,5][600000:]
    ER = x[:,5][numOfPoints:]
    # EOsc = x[:,7][600000:]

    #Utility function being defined on the fly for averaging energy throughput
    #   def Average(lst):
    #        return sum(lst) / len(lst)
    #Where do we use this?

    HEnergyinonedrive = (ER[-1]-ER[-(totPoints-(numOfPoints+1))])/((totPoints-numOfPoints)/pointPerPeroid)
    #160 because 200-40, Harvested energy in one drive (takes the last value subtracts before transient,
    #then divides by number of periods)
    Energyinonedrive = (Edrive[-1]-Edrive[-(totPoints-(numOfPoints+1))])/((totPoints-numOfPoints)/pointPerPeroid) #Driven energy
    DissEnergyinonedrive = (Ediss[-1]-Ediss[-(totPoints-(numOfPoints+1))])/((totPoints-numOfPoints)/pointPerPeroid)
    enEffNum = HEnergyinonedrive/Energyinonedrive


    dissList.append(DissEnergyinonedrive)
    driveList.append(Energyinonedrive)
    harvestList.append(HEnergyinonedrive)
    effList.append(enEffNum)

    # Data saved
    # Nice plotting set up

    ax = fig.add_subplot(1, 3, loop, projection='3d')

    ax.scatter3D(q, p, v, c=v)
    ax.set_xlabel('q')
    ax.set_ylabel('p')
    ax.set_zlabel('v')
    ax.view_init(15, 45)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.subplots_adjust(hspace=0.4, wspace=0.4)
    t1 = time.time()
    print(t1-t0)
    if loop == 1:
        plt.title("a)", loc='left')
    elif loop == 2:
        plt.title("b)", loc='left')
    else:
        plt.title("c)", loc='left')
    loop = loop + 1

plt.savefig('Fig.9.png',bbox_inches='tight', dpi = 100)
plt.close('all')