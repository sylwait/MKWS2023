import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

def combust(T: float, P: float, spcies_dict: dict, reactor: str, dt: float= 0.0005, till_ignition: bool=True):
    gas = ct.Solution('gri30.yaml')
    gas.TPX = T, P, spcies_dict
    # Typ reaktora wg symulacji
    if reactor == 'constant temperature':
        r = ct.IdealGasReactor(gas)
    
    elif reactor == 'constant pressure':
        r = ct.IdealGasConstPressureReactor(gas)
    
    else:
        raise TypeError(f"Improper reactor type '{reactor}'")
    
    sim = ct.ReactorNet([r])
    # Stany
    states = ct.SolutionArray(gas, extra=['time_in_ms', 'time_in_sec'])

    # Symulacja dla 10 sec = dt*n_steps
    simulation_time = 10
    # dt rozmiar kroku
    n_steps = int(simulation_time/dt)
    ignition_delay = 0.0
    
    # Warunek sprawdzający, czy doszło do samozapłonu 

    ignited = False
    
    time = 0.0
    for n in range(n_steps):
        time += dt
        sim.advance(time)
        states.append(r.thermo.state, time_in_sec= time, time_in_ms= time*1e3)
        
        if ignited == False:
            if states.T[n] >= (T + 400):
                ignition_delay = time
                ignited = True
                if till_ignition == True:
                    break

    return gas, states, ignition_delay

# PART I-const temperature
fig1, axs = plt.subplots(2,1, figsize=(10, 10), layout='tight')

igd_lst = []
P = np.linspace(1,5,9)
for t in P:
    _, states, igd = combust(1250, t*ct.one_atm, {'CH4':1, 'O2':2, 'N2': 7.52}, 'constant temperature', dt= 2e-4)
    axs[0].plot(states.time_in_ms, states.T, '.-', label= f'{t: .3} atm, I.D. = {igd*1000: .5} ms')
    igd_lst.append(igd)

axs[0].set_xlabel('Time [ms]$\longrightarrow$')
axs[0].set_ylabel("Flame Temperature [K]$\longrightarrow$")
axs[0].grid(linestyle='-.')
axs[0].legend()

axs[1].plot(np.array(igd_lst)*1000, P, '.-')
axs[1].set_xlabel('Auto-ignition time [ms]$\longrightarrow$')
axs[1].set_ylabel("Pressure [atm]$\longrightarrow$")
axs[1].grid(linestyle='-.')

fig1.suptitle(f"Auto-ignition of methane at 1250 K and pressure variation")

# PART I-const pressure
fig2, axs = plt.subplots(2,1, figsize=(10, 10), layout='tight')

igd_lst = []
T = np.linspace(950, 1450, 9)
for t in T:
    _, states, igd = combust(t, 5*ct.one_atm, {'CH4':1, 'O2':2, 'N2': 7.52}, 'constant pressure')
    axs[0].plot(states.time_in_ms, states.T, '.-', label= f'{t} K, I.D. = {igd*1000: .5} ms')
    igd_lst.append(igd)

axs[0].set_xlabel('Time [ms]$\longrightarrow$')
axs[0].set_ylabel("Flame Temperature [K]$\longrightarrow$")
axs[0].grid(linestyle='-.')
axs[0].legend()

axs[1].plot(np.array(igd_lst)*1000, T, '.-')
axs[1].set_xlabel('Auto-ignition time [ms]$\longrightarrow$')
axs[1].set_ylabel("Initial temperature [K]$\longrightarrow$")
axs[1].grid(linestyle='-.')

fig2.suptitle(f"Auto-ignition of methane at 5 atm and initial temperature variation")

plt.show()

#Porównanie czasu zapłonu wg składu mieszanki (moli metanu)

fig4, axs = plt.subplots(1,1, figsize=(10, 10), layout='tight')

sim_iters = 10000
times = np.zeros(sim_iters)
gas = ct.Solution('gri30.yaml')
data = np.zeros([sim_iters,4])

def check_time(tab):
    delta = 400.                   #sprawdzenie kiedy doszlo do zaplonu
    for i in range(sim_iters-1):
        if (tab[i+1,0] - tab[i,0]) > delta:
            break
    return i

fuel_mole = np.array([0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
P = np.linspace(1,5,9)
ignition_time = np.zeros([fuel_mole.size,P.size])
t_counter = 0
temperature = 1250
for t in P:
    j = 0
    for m in fuel_mole:
        gas.TPX = temperature, t*ct.one_atm, 'CH4:%f, O2:1, N2:3.76' % m
        r = ct.IdealGasReactor(gas)
        sim = ct.ReactorNet([r])
        time = 0
        for n in range(sim_iters):
            time += 5.e-4
            sim.advance(time)
            data[n,0] = r.T
            times[n] = time
        ignition_time[j, t_counter] = times[check_time(data)]
        j+=1
    t_counter += 1

print('%15s    %15s' % ('CH4 [mole]', 'Ignition time [s]'))
       
for i in range(fuel_mole.size):
    print('%15.2f %15.3f' % (fuel_mole[i], ignition_time[i,0]))
        
t_counter = 0
for t in P:
    axs.plot(fuel_mole, ignition_time[:,t_counter], '.-', label= f'{t: .3} atm')
    t_counter += 1

axs.set_xlabel('Fuel quantity [mole]$\longrightarrow$')
axs.set_ylabel("Ignition time [s]$\longrightarrow$")
axs.grid(linestyle='-.')
axs.legend()
fig4.suptitle(f"Ignition time at %dK with varying fuel quantity and pressure" % temperature)
    
plt.show()
