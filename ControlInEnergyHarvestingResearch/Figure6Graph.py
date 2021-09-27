import numpy as np
from scipy.integrate import odeint
'''import scipy.integrate as integrate'''
import matplotlib.pyplot as plt
import matplotlib
import math
import sympy
import statistics
import sys
import operator
import collections
import time
import random

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 15}

matplotlib.rc('font', **font)

f = np.load("10000PerBaseGPeriod19.npy")

#%% Setting parameters
# Andre's regime with varying g
gam = 0.02
om = 1.0
kappa = 0.5
xi = 0.05
lam = 0.05
# gWant = 2.6014
factor = 10
C = 0.05/0.5
R = 0.5/0.05/0.05
# print(R)

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


# def Correct(x,t):
#     # Assigning symbols to elements of x for ease
#     q = x[0] #q is first state
#     p = x[1]
#     v = x[2]
#     Edrive = x[3]
#     Ediss = x[4]
#
#
#     # Renaming some parameters -- but this seems like a re-write, see 25 onwards
#     # lam = 1/(R*C)
#     # xi = C*(xi/C)
#     #Symbolic equations of motion
#     dqdt = p
#     dpdt = -q**3+q + g*math.cos(om*t)-2*gam*p+xi*v
#     dvdt = -v/(R*C)-(xi/C)*p
#     #Some other time integrals we want to compute
#     dEdrive_dt = g*math.cos(om*t)*p #energy absorbed from drive
#     dEdiss_dt = 2*gam*p*p # energy disspated through Gamma
#     #    dEC = C*v*(-lam*v-(xi/C)*p) #energy in the capacitor
#     dER = v*v/R #Energy in the resistor
#     #    dOsc = p*(-q * q * q +  q + g * math.cos(om * t) - 2 * gam * p + xi * v) + (
#     #                - q * q * q + q) * p # Energy IN the oscillator
#     return [dqdt,dpdt,dvdt,dEdrive_dt,dEdiss_dt, dER]
# def Correct(x,t):
#     t1 = t[0]
#     dt = t[1]-t[0]
#     # Assigning symbols to elements of x for ease
#     q = x[0] #q is first state
#     p = x[1]
#     v = x[2]
#     Edrive = x[3]
#     Ediss = x[4]
#
#     qActNext = f[periods*pointPerPeroid-1][0]
#     pActNext = f[periods*pointPerPeroid-1][1]
#     vActNext = f[periods*pointPerPeroid-1][3]
#     # Renaming some parameters -- but this seems like a re-write, see 25 onwards
#     # lam = 1/(R*C)
#     # xi = C*(xi/C)
#     #Symbolic equations of motion
#     dqdt = p
#     dpdt = -q**3+q + g*math.cos(om*t)-2*gam*p+xi*v
#     dvdt = -v/(R*C)-(xi/C)*p
#     #Some other time integrals we want to compute
#     dEdrive_dt = g*math.cos(om*t)*p #energy absorbed from drive
#     dEdiss_dt = 2*gam*p*p # energy disspated through Gamma
#     #    dEC = C*v*(-lam*v-(xi/C)*p) #energy in the capacitor
#     dER = v*v/R #Energy in the resistor
#     #    dOsc = p*(-q * q * q +  q + g * math.cos(om * t) - 2 * gam * p + xi * v) + (
#     #                - q * q * q + q) * p # Energy IN the oscillator
#     return [dqdt,dpdt,dvdt,dEdrive_dt,dEdiss_dt, dER]

#%% Creating empty lists to save values
dissList = []
driveList = []
effList = []
harvestList = []
#averagePoincare = []
# forceTot = 0

# testList = [26017]
loop = 0
testList = [25504, 25534, 26023, 26050, 26340, 27000, 26150, 26220, 26350, 26407, 26510, 28570, 26630, 26680, 26750, 26830, 26900]
# testList = [26017]
averagesF = []
for j in testList:

    g = j/10000.0
    gOr = j/10000.0
    if j == 26104:
        t0 = time.time()

        totPoints = 1000000*factor
        periods = 1000*factor
        pointPerPeroid = int(totPoints/periods)
        stepSize = (2*math.pi)/om/pointPerPeroid
        force = []
        forceAbsolute = 0
        checks = 20
        x_lis = []
        previousStart = 0
        forceAv = 0
        forceAbsAv = 0
        x = f
        force = np.zeros(periods)
    else:
        # R=i/100

        t0 = time.time() #Calling computer clock

        x0 = [1,0,0,0,0,0] #Initial values. Change here.

        totPoints = 1000000*factor
        periods = 1000*factor
        pointPerPeroid = int(totPoints/periods)
        stepSize = (2*math.pi)/om/pointPerPeroid
        force = []
        forceAbsolute = 0
        checks = 20
        x_lis = []
        previousStart = 0
        # for i in range(checks):
        #     k = i+1
        for i in range(periods):
            # g = g + random.uniform(-0.0005, 0.0005)
            # gam = gam + random.uniform(-0.0005, 0.0005)
            start = i*(2*math.pi)/om
            end = (i+1)*(2*math.pi)/om
            # print((end-previousStart)/pointPerPeroid)
            t = np.linspace(start, end, pointPerPeroid)

            z = odeint(Harvest, x0, t) #(Function, initial condition, time)
            z_list = z.tolist()
            # print(z_list[0],z_list[1],z_list[3])
            for j in range(z.__len__()-1):
                x_lis.append(z_list[j])

            # print(z)`

            t = np.linspace(end-stepSize, end+stepSize, 2)

            y = odeint(Harvest, z[pointPerPeroid-1], t)
            # print(x_lis[x_lis.__len__()-1])
            # print(y)
            y_list = y.tolist()
            qY = y_list[0][0]
            pY = y_list[0][1]
            vY = y_list[0][2]
            if i+1 < periods:
                qActNext = f[(i+1)*pointPerPeroid][0]
                pActNext = f[(i+1)*pointPerPeroid][1]
                vActNext = f[(i+1)*pointPerPeroid][3]
            else:
                qActNext = f[(i+1)*pointPerPeroid-1][0]
                pActNext = f[(i+1)*pointPerPeroid-1][1]
                vActNext = f[(i+1)*pointPerPeroid-1][3]
            dq = qActNext-qY
            dp = pActNext-pY
            dv = vActNext-vY
            x_lis.append([qActNext, pActNext, vActNext, y[0][3], y[0][4], y[0][5]])
            # print(y)
            # print(z)
            # print(y[1])
            # print(dq, dp, stepSize)
            force.append(dp/stepSize)
            # forceTot = forceTot+(dp/stepSize+dq/(stepSize**2))**2
            forceAbsolute = forceAbsolute + abs(dp/stepSize)
            x0 = x_lis[x_lis.__len__()-2]
            # print(force)
        # print(x_lis[i*pointPerPeroid-1], f[i*pointPerPeroid-1])
        x = np.array(x_lis)
        # forceAv = forceTot/periods
        forceAbsAv = forceAbsolute/periods
    # np.save('BaseGPeriod19', x)
    # print("saved")
    #Going to produce five columns of results, first colum is q, then p , then v
    numOfPoints = 980000*factor #This is the transient number of points to be rejected right here

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

    averagesF.append(forceAbsAv)
    loop = loop+1
    print(loop)

fig, axes = plt.subplots(1, 1, figsize=(20, 15), tight_layout=True)
# fig, axes = plt.subplots(1,1,figsize=(20, 15), tight_layout=True)
col = 4
row = 3

interP = plt.subplot(col,row,1)
values = []
for i in testList:
    values.append(float(i/10000))
plt.scatter(values, averagesF)
interP.set_title("Inter-Period Push", loc = 'center')
plt.xlabel("Driving Force Amplitude")
plt.ylabel("Absolute Value of Average Force")

# What is AveragePoincare ?
pointslist = [25504, 25534, 26023, 26050, 26340, 27000, 26150, 26220, 26350, 26407, 26510, 28570, 26630, 26680, 26750, 26830, 26900]
# pointslist = [26017]
#%% What is this module -- are we scanning over g here ?
averageForces = []
for k in pointslist:
    g = k/10000.0
    gOr = k/10000.0
    t0 = time.time() #Calling computer clock

    x0 = [1,0,0,0,0,0] #Initial values. Change here.

    totPoints = 1000000*factor
    periods = 1000*factor
    pointPerPeroid = int(totPoints/periods)
    stepSize = (2*math.pi)/om/pointPerPeroid
    force = []
    forceAbsolute = 0
    checks = 20
    x_lis = []
    previousStart = 0
    # for i in range(checks):
    #     k = i+1
    for i in range(periods):
        # g = g + random.uniform(-0.0005, 0.0005)
        # gam = gam + random.uniform(-0.0005, 0.0005)
        start = i*(2*math.pi)/om
        end = (i+1)*(2*math.pi)/om
        # print((end-previousStart)/pointPerPeroid)
        t = np.linspace(start, end, pointPerPeroid)

        z = odeint(Harvest, x0, t) #(Function, initial condition, time)
        z_list = z.tolist()
        # print(z_list[0],z_list[1],z_list[3])
        for j in range(z.__len__()-1):
            x_lis.append(z_list[j])

        # print(z)`

        t = np.linspace(end-stepSize, end+stepSize, 2)

        y = odeint(Harvest, z[pointPerPeroid-1], t)
        # print(x_lis[x_lis.__len__()-1])
        # print(y)
        y_list = y.tolist()
        qY = y_list[0][0]
        pY = y_list[0][1]
        vY = y_list[0][2]
        if i+1 < periods:
            qActNext = f[(i+1)*pointPerPeroid][0]
            pActNext = f[(i+1)*pointPerPeroid][1]
            vActNext = f[(i+1)*pointPerPeroid][3]
        else:
            qActNext = f[(i+1)*pointPerPeroid-1][0]
            pActNext = f[(i+1)*pointPerPeroid-1][1]
            vActNext = f[(i+1)*pointPerPeroid-1][3]
        dq = qActNext-qY
        dp = pActNext-pY
        dv = vActNext-vY
        x_lis.append([qActNext, pActNext, vActNext, y[0][3], y[0][4], y[0][5]])
        # print(y)
        # print(z)
        # print(y[1])
        # print(dq, dp, stepSize)
        force.append(dp/stepSize)
        # forceTot = forceTot+(dp/stepSize+dq/(stepSize**2))**2
        forceAbsolute = forceAbsolute + abs(dp/stepSize)
        x0 = x_lis[x_lis.__len__()-1]
        # print(force)
    # print(x_lis[i*pointPerPeroid-1], f[i*pointPerPeroid-1])
    x = np.array(x_lis)
    # forceAv = forceTot/periods
    forceAbsAv = forceAbsolute/periods
    averageForces.append(forceAbsAv)
    loop = loop + 1
    print(loop)

terminal = plt.subplot(col,row,2)
gs = []
for i in pointslist:
    gs.append(float(i/10000))
plt.scatter(gs, averageForces)
terminal.set_title("Terminal Push", loc = 'center')
plt.xlabel("Driving Force Amplitude")
plt.ylabel("Absolute Value of Average Force")


plt.savefig('Average Gs Parametric Drift',bbox_inches='tight', dpi = 100)

plt.close('all')