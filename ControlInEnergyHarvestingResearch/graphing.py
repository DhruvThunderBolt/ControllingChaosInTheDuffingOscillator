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

# Control method for terminal push
# Tested amplitudes
testList = [28570]
loop = 1
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

        totPoints = 1000000*1
        periods = 1000*1
        pointPerPeroid = int(totPoints/periods)
        stepSize = (2*math.pi)/om/pointPerPeroid
        force = []
        forceAbsolute = 0
        checks = 20
        x_lis = []
        previousStart = 0

        for i in range(periods):
            # Amplitude drift implemented for parametric drift
            g = g + random.uniform(-0.0005, 0.0005)

            # start and endpoints of one oscillation
            start = i*(2*math.pi)/om
            end = (i+1)*(2*math.pi)/om

            t = np.linspace(start, end, pointPerPeroid)

            z = odeint(Harvest, x0, t)
            z_list = z.tolist()

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
    print(forceAbsAv)
    # np.save('BaseGPeriod19', x)
    # print("saved")
    #Going to produce five columns of results, first colum is q, then p , then v
    numOfPoints = 980000*1 #This is the transient number of points to be rejected right here

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

    fig, axes = plt.subplots(1, 1, figsize=(20, 15), tight_layout=True)
    # fig, axes = plt.subplots(1,1,figsize=(20, 15), tight_layout=True)
    col = 4
    row = 3

    phase = plt.subplot(col,row,1)
    # phase = plt.subplot(col,row,loop + 1)
    plt.scatter(q,p, s=1)

    #Translating to Poincare
    xs = []
    ys = []
    xs = [x[int(totPoints/periods)*i,0] for i in range(periods)]
    ys = [x[int(totPoints/periods)*i,1] for i in range(periods)]
    plt.scatter(xs[990*1:],ys[990*1:], color="red")
    #Poincare

    #   averagePoincare.append([statistics.mean(xs[9900:]), statistics.mean(ys[9900:])])
    #   What the heck is this ? Oh something we were doing earlier, not to worry
    # print(qmin, qmax, pmin, pmax)
    plt.xlabel('q')
    # plt.xlim(qmin-5, qmax+5)
    plt.ylabel('p')
    # plt.ylim(pmin-5, pmax+5)
    plt.axis([np.amin(q)-1, np.amax(q)+1, np.amin(p)-1, np.amax(p)+1])
    # plt.axis([-1.75,1.75,-1, 1])
    # if(loop == 1):
    phase.set_title('a) Phase Plot for 1,000 Oscillations', loc = 'left')
    # else:
    #     phase.set_title('c)                     g = {:.5f}'.format(gOr), loc = 'left')
    # # phase.set_title('Phase Space')
    #
    # LineGraphF = plt.subplot(col,row,loop+4)
    # plt.plot(np.linspace(0,periods*2*math.pi/om, periods),force)
    # plt.xlabel('Time')
    # plt.ylabel('Force')
    # # plt.xlim(0,periods*2*math.pi/om)
    # # plt.hist(v, bins=500, density = True)
    # # plt.xlabel('v')
    # # plt.ylabel(r'$P_v$')
    # if loop == 1:
    #     LineGraphF.set_title('e)            Force v. Time for g = {:.5f}'.format(gOr), loc = 'left')
    # else:
    #     LineGraphF.set_title('f)            Force v. Time for g = {:.5f}'.format(gOr), loc = 'left')
    # # HistogramV.set_title("Histogram of v")
    t1 = time.time()
    print(t1-t0)



    plt.subplots_adjust(hspace=0.4, wspace=0.4)
    # plt.suptitle("Om = {}, xi = {}, C = {}, R = {}, g = {}, gam = {}".format(om, xi, C, R, g, gam), fontsize = 25)
    # plt.savefig('10000PeriodsHistPlotxi{:.2f}g{:.4f}r{:.2f}om{:.3f}param.png'.format(xi,g,R,om),bbox_inches='tight', dpi = 100)
    # plt.show()
    # plt.close('all')
    loop = loop + 1

# What is AveragePoincare ?
# pointslist = [25504, 25534, 26023, 26050, 26340, 27000, 26150, 26220, 26350, 26407, 26510, 28570, 26630, 26680, 26750, 26830, 26900]
pointslist = [28570]
#%% What is this module -- are we scanning over g here ?
averageForces = []
for k in pointslist:
    g = k/10000.0
    gOr = k/10000.0
    if k == 26104:
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
            g = g + random.uniform(-0.0005, 0.0005)
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
        # print(x)
        # print(x.__len__())
    print(forceAbsAv)
    # np.save('GCorrectiveDATAxi{:.2f}g{:.4f}r{:.2f}om{:.3f}param.png'.format(xi,g,R,om), x)
    #Going from 0 to 50000, with 100000 points of data.
    #print(x)

    #Going to produce five columns of results, first colum is q, then p , then v
    numOfPoints = 980000*factor #This is the transient number of points to be rejected right here

    q = x[:,0][numOfPoints:] #Starting FROM numOfPoints
    p = x[:,1][numOfPoints:]
    v = x[:,2][numOfPoints:]
    # print(v)

    Edrive = x[:,3][numOfPoints:]
    Ediss = x[:,4][numOfPoints:]
    # Ecap = x[:,5][600000:]
    ER = x[:,5][numOfPoints:]
    # EOsc = x[:,7][600000:]

    # for i in range(periods):
    #     start = i*(2*math.pi)/om
    #     end = (i+1)*(2*math.pi)/om
    #     # print((end-previousStart)/pointPerPeroid)
    #     t = np.linspace(start, end, pointPerPeroid)
    #
    #     z = odeint(Harvest, x0, t) #(Function, initial condition, time)
    #     z_list = z.tolist()
    #     # print(z_list[0],z_list[1],z_list[3])
    #     if i%checkEveryXPeriods==0:
    #         for j in range(z.__len__()-1):
    #             x_lis.append(z_list[j])
    #         # print(z)
    #
    #         t = np.linspace(end-stepSize, end+stepSize, 2)
    #
    #         y = Correct(z, pointPerPeroid, i, end, stepSize)
    #         # print(type(y))
    #         x_lis.append(y)
    #         # print(y)
    #         # print(z)
    #         # print(y[1])
    #         x0 = y
    #     else:
    #         for j in range(z.__len__()):
    #             x_lis.append(z_list[j])
    #     print(i)

    # for i in range(periods):
    #     periodsAcheck = int(periods/checks)
    #     start = i*(2*math.pi)/om
    #     end = (i+1)*(2*math.pi)/om
    #     # print((end-previousStart)/pointPerPeroid)
    #     t = np.linspace(start, end, pointPerPeroid)
    #
    #     z = odeint(Harvest, x0, t) #(Function, initial condition, time)
    #     z_list = z.tolist()
    #     # print(z_list[0],z_list[1],z_list[3])
    #     for j in range(z.__len__()-1):
    #         x_lis.append(z_list[j])
    #     # print(z)
    #
    #     t = np.linspace(end-stepSize, end+stepSize, 2)
    #
    #     y = Correct(z, pointPerPeroid, i, end, stepSize)
    #     # print(type(y))
    #     x_lis.append(y)
    #     # print(y)
    #     # print(z)
    #     # print(y[1])
    #     x0 = y
    #     print(i)
    #     # print(y)
    # for i in range(checks):
    #
    #     start = i*periodsAcheck*(2*math.pi)/om
    #     end = (i+1)*periodsAcheck*(2*math.pi)/om
    #     # print((end-previousStart)/pointPerPeroid)
    #     t = np.linspace(start, end, pointPerPeroid*periodsAcheck)
    #
    #     z = odeint(Harvest, x0, t) #(Function, initial condition, time)
    #     # print(f[i*periodsAcheck+10][0])
    #     z_list = z.tolist()
    #     # print(z_list[0],z_list[1],z_list[3])
    #     for j in range(z.__len__()-1):
    #         x_lis.append(z_list[j])
    #     # print(z)
    #     # print(x_lis[i*periodsAcheck+10][0])
    #     t = np.linspace(end-stepSize, end+stepSize, 2)
    #
    #     y = Correct(z, pointPerPeroid, (i+1)*periodsAcheck-1, end, stepSize) # need to adjust periodsize x value for periods a check
    #     # print(type(y))
    #     x_lis.append(y)
    #     # print(y)
    #     # print(z)
    #     # print(y[1])
    #     # print(x_lis[i*periodsAcheck+10][0])
    #     x0 = y
    #     print(i)
    #     # print(y)

    #Utility function being defined on the fly for averaging energy throughput
    #   def Average(lst):
    #        return sum(lst) / len(lst)
    #Where do we use this?

    # HEnergyinonedrive = (ER[-1]-ER[-(totPoints-(numOfPoints+1))])/((totPoints-numOfPoints)/pointPerPeroid)
    # #160 because 200-40, Harvested energy in one drive (takes the last value subtracts before transient,
    # #then divides by number of periods)
    # Energyinonedrive = (Edrive[-1]-Edrive[-(totPoints-(numOfPoints+1))])/((totPoints-numOfPoints)/pointPerPeroid) #Driven energy
    # DissEnergyinonedrive = (Ediss[-1]-Ediss[-(totPoints-(numOfPoints+1))])/((totPoints-numOfPoints)/pointPerPeroid)
    # enEffNum = HEnergyinonedrive/Energyinonedrive
    #
    #
    # dissList.append(DissEnergyinonedrive)
    # driveList.append(Energyinonedrive)
    # harvestList.append(HEnergyinonedrive)
    # effList.append(enEffNum)

    # Data saved
    # Nice plotting set up

    col = 4
    row = 3

    phase = plt.subplot(col,row,2)
    plt.scatter(q,p, s=1)

    #Translating to Poincare
    xs = []
    ys = []
    xs = [x[int(totPoints/periods)*i,0] for i in range(periods)]
    ys = [x[int(totPoints/periods)*i,1] for i in range(periods)]
    plt.scatter(xs[990*factor:], ys[990*factor:], color="red")
    #Poincare

    #   averagePoincare.append([statistics.mean(xs[9900:]), statistics.mean(ys[9900:])])
    #   What the heck is this ? Oh something we were doing earlier, not to worry
    # print(qmin, qmax, pmin, pmax)
    plt.xlabel('q')
    # plt.xlim(qmin-5, qmax+5)
    plt.ylabel('p')
    # plt.ylim(pmin-5, pmax+5)
    plt.axis([np.amin(q)-1, np.amax(q)+1, np.amin(p)-1, np.amax(p)+1])
    # plt.axis([-1.75,1.75,-1, 1])
    phase.set_title('b) Phase Plot for 10,000 Oscillations'.format(gOr), loc = 'left')
    # phase.set_title('Phase Space')



    # EHistplt = plt.subplot(col,row,2)
    # plt.hexbin(q,p, extent=[np.amin(q)-1, np.amax(q)+1, np.amin(p)-1, np.amax(p)+1])
    # # plt.hexbin(q,p, extent=[-1.75,1.75, -1, 1])
    # plt.xlabel('q')
    # plt.ylabel('p')
    # EHistplt.set_title('b)', loc = 'left')
    # # EHistplt.set_title("Histogram of Phase Space")
    #
    #
    # # Histogram = plt.subplot(col,row,3)
    # # plt.hist(p, bins=500, density = True)
    # # # fig = plt.figure()
    # # # ax = fig.add_axes([0,0,1,1])
    # # # langs = ['C', 'C++', 'Java', 'Python', 'PHP']
    # # # students = [23,17,35,29,12]
    # # Histogram.bar(["Average Force Squared"],[forceAv])
    # # # plt.show()
    # # # plt.xlabel('p')
    # # # plt.ylabel(r'$P_p$')
    # # Histogram.set_title('c)', loc = 'left')
    # # # Histogram.set_title("Histogram of p")
    #
    # Histogram = plt.subplot(col,row,3)
    # plt.hist(p, bins=500, density = True)
    # # fig = plt.figure()
    # # ax = fig.add_axes([0,0,1,1])
    # # langs = ['C', 'C++', 'Java', 'Python', 'PHP']
    # # students = [23,17,35,29,12]
    # Histogram.bar(["Magnitude of Average Force "], [forceAbsAv])
    # # plt.show()
    # # plt.xlabel('p')
    # # plt.ylabel(r'$P_p$')
    # Histogram.set_title('c)', loc = 'left')
    # # Histogram.set_title("Histogram of p")
    #
    # capacitor = [x * xi for x in v]
    #
    # CapSpace = plt.subplot(col, row, 4)
    # plt.scatter(q, v,s=0.5, )
    # plt.xlabel('q')
    # plt.ylabel('v')
    # plt.axis([np.amin(q)-1, np.amax(q)+1, np.amin(v)-1, np.amax(v)+1])
    # # plt.axis([-1.5,1.5, -0.2,0.2])
    # CapSpace.set_title('d)', loc = 'left')
    # # CapSpace.set_title("Capacitor Space")
    #
    # HistCapacitor = plt.subplot(col,row,5)
    # # plt.hexbin(q,capacitor, extent=[-3.5,3.5, -0.1,0.1])
    # # plt.hexbin(p,v, extent=[-1,1, -0.2,0.2])
    # plt.hexbin(v, p, extent=[np.amin(v)-1, np.amax(v)+1, np.amin(p)-1, np.amax(p)+1])
    # plt.xlabel('v')
    # plt.xlim(np.amin(v)-1, np.amax(v)+1)
    # plt.ylabel('p')
    # # plt.axis([-1.8,1.8, -0.4,0.4])
    # HistCapacitor.set_title('e)', loc = 'left')
    # # HistCapacitor.set_title("Histogram of Capacitor Space")
    #
    LineGraphF = plt.subplot(col,row,3)
    plt.plot(np.linspace(0,periods*2*math.pi/om, periods),force)
    plt.xlabel('Time')
    plt.ylabel('Force')
    # plt.xlim(0,periods*2*math.pi/om)
    # plt.hist(v, bins=500, density = True)
    # plt.xlabel('v')
    # plt.ylabel(r'$P_v$')
    LineGraphF.set_title('c) Force v. Time for 10,000 Oscillations'.format(gOr), loc = 'left')
    # HistogramV.set_title("Histogram of v")

    # damping = [x * -2 * gam for x in p]
    # # What the heck is this ?
    # # drove = [g * math.cos(om * t) for t in np.linspace(0, 100*(2*math.pi)/om, (totPoints-numOfPoints))]
    # drove = [g * math.cos(om * t) for t in np.linspace(0, (totPoints-numOfPoints)/pointPerPeroid*(2*math.pi)/om, (totPoints-numOfPoints))]

    # driving = plt.subplot(col, row, 7)
    # driving = plt.axes(projection='3d')
    # driving.scatter3D(p,q,v)
    # driving.set_xlabel('p')
    # driving.set_ylabel('q')
    # driving.set_zlabel('v')
    # # plt.zlabel("z")
    # # driving.xlim(np.amin(p)-1,np.amax(p)+1)
    # # driving.ylim(np.amin(q)-1,np.amax(q)+1)
    # # driving.zlim(np.amin(v)-1,np.amax(v)+1)
    # # plt.axis([-1.75,1.75, -1,1])
    # driving.set_title('g)', loc = 'left')
    # driving.set_title('3d Projection of Phase Space')


    # HistDrive = plt.subplot(col,row,8)
    # plt.hexbin(p,drove, extent=[np.amin(p)-1,np.amax(p)+1, np.amin(drove)-1, np.amax(drove)+1])
    # # plt.hexbin(p,drove, extent=[-1.75,1.75, -1,1])
    # plt.xlabel('p')
    # plt.ylabel(r'$g\mathrm{cos}(\omega t)$')
    # HistDrive.set_title('h)', loc = 'left')
    # HistDrive.set_title("Histogram of Driving Space")



    # Histogramdrive = plt.subplot(col,row,9)
    # labels = [r'$E_R$',r'$E_{Drive}$',r'$E_{Diss}$']
    # barNum = [round(HEnergyinonedrive,3),round(Energyinonedrive,3),round(DissEnergyinonedrive,3)]
    # x = np.arange(len(labels))
    # width = 0.35
    # Histogramdrive.bar(x, barNum, width)
    # Histogramdrive.set_xticks(x)
    # Histogramdrive.set_xticklabels(labels)
    # plt.ylim(top=barNum[0]+2)
    # plt.ylabel('Average Energy per Period')
    # Histogramdrive.set_title('i)', loc = 'left')

    t1 = time.time()
    print(t1-t0)



    plt.subplots_adjust(hspace=0.4, wspace=0.4)
    # plt.suptitle("Om = {}, xi = {}, C = {}, R = {}, g = {}, gam = {}".format(om, xi, C, R, g, gam), fontsize = 25)
    # plt.savefig('GCorrectionxi{:.2f}g{:.4f}r{:.2f}om{:.3f}param.png'.format(xi,g,R,om),bbox_inches='tight', dpi = 100)
    # plt.show()
    # plt.close('all')
    averageForces.append(forceAbsAv)


plt.savefig('Fig.8.png',bbox_inches='tight', dpi = 100)
plt.close('all')

# plt.title("Force against g values")
# plt.xlabel("g")
# plt.ylabel("Absolute Value of Average Force")
# plt.scatter(pointslist, averageForces)
# plt.savefig("Average Gs")
# plt.close('all')