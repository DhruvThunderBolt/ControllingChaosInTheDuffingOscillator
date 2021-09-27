import numpy as np
from scipy.integrate import odeint
'''import scipy.integrate as integrate'''
import matplotlib.pyplot as plt
import matplotlib
import math
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
#averagePoincare = []

# What is AveragePoincare ?

#%% What is this module -- are we scanning over g here ?
pointsList = [26014]
# pointsList = [98,99]
for i in pointsList:
    # print(i)
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

    fig, axes = plt.subplots(1,1,figsize=(20, 15), tight_layout=True)
    col = 4
    row = 3

    phase = plt.subplot(col,row,1)
    plt.scatter(q,p, s=1)

    #Translating to Poincare
    xs = []
    ys = []
    xs = [x[int(totPoints/periods)*i,0] for i in range(periods)]
    ys = [x[int(totPoints/periods)*i,1] for i in range(periods)]
    plt.scatter(xs[990:],ys[990:], color="red")
    #Poincare

    #   averagePoincare.append([statistics.mean(xs[9900:]), statistics.mean(ys[9900:])])
    #   What the heck is this ? Oh something we were doing earlier, not to worry

    plt.xlabel('q')
    plt.ylabel('p')
    plt.axis([-3.5, 3.5, -5,5])
    # plt.axis([-1.75,1.75,-1, 1])
    phase.set_title('a)', loc = 'left')
    # phase.set_title('Phase Space')



    EHistplt = plt.subplot(col,row,2)
    plt.hexbin(q,p, extent=[-3.5, 3.5, -5,5])
    # plt.hexbin(q,p, extent=[-1.75,1.75, -1, 1])
    plt.xlabel('q')
    plt.ylabel('p')
    EHistplt.set_title('b)', loc = 'left')
    # EHistplt.set_title("Histogram of Phase Space")


    Histogram = plt.subplot(col,row,3)
    plt.hist(p, bins=500, density = True)
    plt.xlabel('p')
    plt.ylabel(r'$P_p$')
    Histogram.set_title('c)', loc = 'left')
    # Histogram.set_title("Histogram of p")

    capacitor = [x * xi for x in v]

    CapSpace = plt.subplot(col, row, 4)
    plt.scatter(q, v,s=0.5, )
    plt.xlabel('q')
    plt.ylabel('v')
    plt.axis([-3.5,3.5, -3,3])
    # plt.axis([-1.5,1.5, -0.2,0.2])
    CapSpace.set_title('d)', loc = 'left')
    # CapSpace.set_title("Capacitor Space")

    HistCapacitor = plt.subplot(col,row,5)
    # plt.hexbin(q,capacitor, extent=[-3.5,3.5, -0.1,0.1])
    # plt.hexbin(p,v, extent=[-1,1, -0.2,0.2])
    plt.hexbin(v,p, extent=[-3.5, 3.5, -5,5])
    plt.xlabel('v')
    plt.ylabel('p')
    # plt.axis([-1.8,1.8, -0.4,0.4])
    HistCapacitor.set_title('e)', loc = 'left')
    # HistCapacitor.set_title("Histogram of Capacitor Space")

    HistogramV = plt.subplot(col,row,6)
    plt.hist(v, bins=500, density = True)
    plt.xlabel('v')
    plt.ylabel(r'$P_v$')
    HistogramV.set_title('f)', loc = 'left')
    # HistogramV.set_title("Histogram of v")

    damping = [x * -2 * gam for x in p]
    # What the heck is this ?
    # drove = [g * math.cos(om * t) for t in np.linspace(0, 100*(2*math.pi)/om, (totPoints-numOfPoints))]
    drove = [g * math.cos(om * t) for t in np.linspace(0, (totPoints-numOfPoints)/pointPerPeroid*(2*math.pi)/om, (totPoints-numOfPoints))]

    driving = plt.subplot(col, row, 7)
    plt.scatter(q, drove, c='orange',s=1)
    plt.xlabel("q")
    plt.ylabel(r'$g\mathrm{cos}(\omega t)$')
    plt.xlim(-3,3)
    plt.ylim(-3,3)
    # plt.axis([-1.75,1.75, -1,1])
    driving.set_title('g)', loc = 'left')
    # driving.set_title('Driving Force in Phase Space')


    HistDrive = plt.subplot(col,row,8)
    plt.hexbin(p,drove, extent=[-4.5,4.5, -3,3])
    # plt.hexbin(p,drove, extent=[-1.75,1.75, -1,1])
    plt.xlabel('p')
    plt.ylabel(r'$g\mathrm{cos}(\omega t)$')
    HistDrive.set_title('h)', loc = 'left')
    # HistDrive.set_title("Histogram of Driving Space")



    Histogramdrive = plt.subplot(col,row,9)
    labels = [r'$E_R$',r'$E_{Drive}$',r'$E_{Diss}$']
    barNum = [round(HEnergyinonedrive,3),round(Energyinonedrive,3),round(DissEnergyinonedrive,3)]
    x = np.arange(len(labels))
    width = 0.35
    Histogramdrive.bar(x, barNum, width)
    Histogramdrive.set_xticks(x)
    Histogramdrive.set_xticklabels(labels)
    plt.ylim(top=1.25)
    plt.ylabel('Average Energy per Period')
    Histogramdrive.set_title('i)', loc = 'left')

    t1 = time.time()
    print(t1-t0)



    plt.subplots_adjust(hspace=0.4, wspace=0.4)
    # plt.suptitle("Om = {}, xi = {}, C = {}, R = {}, g = {}, gam = {}".format(om, xi, C, R, g, gam), fontsize = 25)
    # plt.savefig('HistPlotxi{:.2f}g{:.4f}r{:.2f}om{:.3f}param.png'.format(xi,g,R,om),bbox_inches='tight', dpi = 100)
    plt.show()
    # plt.close('all')

# np.savetxt('paperOrbits/averages.dat', averagePoincare)



# plt.subplots(1, 1, figsize=(20, 15))
# col = 4
# row = 1



# dissEnergy = plt.subplot(col, row, 1)
# plt.plot(np.linspace(2.55, 2.61, len(dissList)), dissList)
# plt.axis([2.55,2.61,8.5,12.5])
# plt.xticks(fontsize= 15)
# plt.yticks(fontsize= 15)
# plt.xlabel('R', fontsize = 18)
# plt.ylabel('Ediss',fontsize = 18)
# plt.grid(True)
# plt.xticks(np.arange(2.55, 2.61, 0.01), rotation = 'vertical')


# driveEnergy = plt.subplot(col, row, 2)
# plt.plot(np.linspace(2.55, 2.61, len(driveList)), driveList)
# plt.axis([2.55,2.61,8.5,13])
# plt.xticks(fontsize= 15)
# plt.yticks(fontsize= 15)
# plt.xlabel('R',fontsize = 18)
# plt.ylabel('Edrive',fontsize = 18)
# plt.grid(True)
# plt.xticks(np.arange(2.55, 2.61, 0.01), rotation = 'vertical')


# harvestEnergy = plt.subplot(col, row, 3)
# plt.plot(np.linspace(2.55, 2.61, len(harvestList)), harvestList)
# plt.axis([2.55,2.61,0.154,0.166])
# plt.xticks(fontsize= 15)
# plt.yticks(fontsize= 15)
# plt.xlabel('R',fontsize = 18)
# plt.ylabel('EV',fontsize = 18)
# plt.grid(True)
# plt.xticks(np.arange(2.55, 2.61, 0.01), rotation = 'vertical')



# plt.subplots_adjust(hspace=0.5, wspace=0.3)
# # plt.suptitle("Om = {}, xi = {}, C = {}, g = {}, gam = {}".format(om, xi, C, g, gam), fontsize=25)
# plt.savefig('finalcountdown/AndreScanFinal4xi{:.2f}g{:.2f}om{:.2f}.png'.format(xi, g, om), bbox_inches='tight', dpi=100)
# plt.close('all')

# plt.plot(np.linspace(2.55, 2.61, len(effList)), effList)
# plt.axis([2.55,2.61, 0.013, 0.018])
# plt.xticks(fontsize= 15)
# plt.yticks(fontsize= 15)
# plt.xlabel('R',fontsize = 18)
# plt.ylabel('Energy Efficiency',fontsize = 18)
# plt.grid(True)
# plt.xticks(np.arange(2.55, 2.61, 0.01), rotation = 'vertical')
# plt.savefig('finalcountdown/AndreScanFinal4withEffxi{:.2f}g{:.2f}om{:.2f}.png'.format(xi, g, om), bbox_inches='tight', dpi=100)

# np.savetxt('finalcountdown/energyListFlat.dat', [DissEnergyinonedrive,Energyinonedrive,HEnergyinonedrive,enEffNum])