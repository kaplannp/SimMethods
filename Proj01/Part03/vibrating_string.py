# coding: utf-8
"""Solving the vibrating string in Python with FEniCS

Author: Juan Luis Cano Rodríguez <juanlu@pybonacci.org>

Adapted by John Stratton <strattja@whitman.edu>

References
----------
* Zachmanoglou, E. C., and Dale W. Thoe. "Introduction to Partial Differential
  Equations with Applications". New York: Dover Publications, 1986.

"""
import numpy as np
from numpy import cos, sin, pi
from scipy.optimize import root_scalar
from scipy.optimize import minimize
from scipy.optimize import differential_evolution


def euler(t, dt, c, gamma, pos, vel):
    '''
    This function executes a series of approximations using eulers method.
    It updates the state t/dt times at increments of dt.
    params:
      t: total time to pass
      dt: steps in which to evaluate system over time to pass
      c: constant. The bigger the c, the bigger the change in velocity over the
        same sized timestep, dt
      gamma: damping constant. Should be negative. If positive, you're exploding
      pos: array of positions of each unit on the string
      vel: array of velocities at each unit on the string
    Postconditions:
      pos and vel are updated

    '''
    for timestep in range(int(t/dt)):

        #calculate first order derivative force contribution
        pos += vel*dt # update position based on the change in velocity
        pos_totalforce = np.zeros(pos.shape)
        #offset by 1, so summing (n) = (n) - (n+1)
        pos_totalforce[:-1] += (pos[1:]-pos[:-1])
        pos_totalforce[-1] += -pos[-1]
        #summing (n+1) = (n+1) - (n)
        pos_totalforce[1:] += (pos[:-1]-pos[1:])
        pos_totalforce[0] += -pos[0]
        #at the end of the day n = [(n) - (n+1)] + [(n) - (n-1)]
        #only at this point is total force, actually total force
        pos_totalforce = c*+pos_totalforce 

        #calculate second order derivative force contribution

        vel_totalforce = vel*gamma # update the force from velocity 
        # print((vel - abs(vel - vel_totalforce) > 0).any()) false as expected
        # contribution

        totalforce = pos_totalforce+vel_totalforce 
        vel += totalforce*dt 

def getPeriod(dt, c, pos, vel):
    '''
    params:
      dt: steps in which to evaluate system over time to pass
      c: constant. The bigger the c, the bigger the change in velocity over the
        same sized timestep, dt
      pos: array of positions of each unit on the string
      vel: array of velocities at each unit on the string
    returns: the period of a cycle
    Postconditions:
      pos and vel have been modified to reflect state at next transition
    '''
    period = 0;
    vel0 = vel[vel.shape[0]//2]
    euler(dt+.0000001, dt, c, pos, vel)
    vel1 = vel[vel.shape[0]//2]
    while(vel0*vel1 >= 0):
        #while the velocity has not changed sign
        #first round, vel0 == 0, so vel0*vel1 is greater than 0
        euler(dt+.0000001, dt, c, pos, vel)
        vel0 = vel1
        vel1 = vel[vel.shape[0]//2]
        #update period with another dt
        period += dt
    return period*2 #period is only one half cycle

if __name__ == '__main__':# s
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib import animation
    #THE TWIDDLES
    L = .3 # m the length of the string (violin g)
    t0 = 0.0 # s start time
    simdt = 1e-4 # s granularity of timestep
    c = 1.5e5 # spring constant. bigger is more force
    gamma = -1e0 # the force developed from velocity
    granularity_line = 101 # this is the number of segments to divide string
    #target for optimization
    targetFreq = 2

    #Plotting Twiddles
    graph_ylim = 0.05  # m 
    framedt = simdt*500 # s granularity of video

    #initialization
    x = np.linspace(0, L, num=granularity_line)
    pos = (1/2 - np.abs(x-L/2)/L)/50 # who knows what this is.
    vel = np.zeros(granularity_line)

    euler(framedt, simdt, c, gamma, pos, vel)
    def anim(i):
        global pos, vel
        euler(framedt, simdt, c, gamma, pos, vel)
        line.set_data(x, pos)
        return line,

    fig = plt.gcf()
    line, = plt.plot(x, pos)
    plt.ylim(-graph_ylim * 2, graph_ylim * 2)
    animate = animation.FuncAnimation(fig, anim, range(40), repeat=True, blit=False)
    plt.show()
    #writervideo = animation.FFMpegWriter(fps=60) 
    #animate.save('string.mp4', fps=90, writer=writervideo)
