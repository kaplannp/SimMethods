import numpy as np
import scipy
import scipy.integrate
import matplotlib.pyplot as plt

#initial conditions
#position = np.asarray([0.0, 1.0])
#velocity = np.asarray([1.0, 1.0])
#px, py, vx, vy
state = np.asarray([0.0, 1.0, 1.0, 2.0])

def derivative(t, y, g):
  velocity = y[2:]
  windresist = (-0.1)*np.linalg.norm(velocity)**2
  # unit vector in direction of wind resistance
  windresist = windresist * ((-1*velocity)/np.linalg.norm(velocity))

  #  windresist = 0.1*velocity**2
  dvdt = np.asarray([0.0, g]) + windresist

  #compute dxdt
  result = np.empty(4)
  result[:2] = velocity
  result[2:] = dvdt
  return result

def plotBallY(states, t_eval):
    plt.plot(t_eval, states[:,1])
    plt.ylabel("height (m)")
    plt.xlabel("time (s)")
def plotBallXY(states, t_eval):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    py = states[:,1]
    px = states[:,0]
    ax.plot3D(px, py, t_eval)
    ax.set_xlabel("xpos (m)")
    ax.set_ylabel("ypos (m)")
    ax.set_zlabel("time (t)")


endt = 0.6
t_eval = np.linspace(0.0, endt, 20);
result = scipy.integrate.solve_ivp(derivative, (0.0, endt), state, t_eval=t_eval,
        args=(-9.8,))
states = result['y'].T
plotBallY(states, t_eval)
plt.show()
plotBallXY(states, t_eval)
plt.show()
#print("Symbolic answer: pos=", pos_t, "speed=", vel_t)
final_result = states[-1,:]
print(final_result)

