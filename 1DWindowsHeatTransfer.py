import numpy as np
import matplotlib.pyplot as plt


def startup():
    ini_temp = float(input('Enter the initial temperature: '))
    nodes = int(input('Enter the number of nodes: '))
    # Function
    temp = np.ones(nodes) * ini_temp
    temp_in = int(input('Enter the inside temperature: '))
    temp_out = int(input('Enter the outside temperature: '))

    return temp, nodes, temp_out, temp_in


def simulation(temp, nodes, temp_out, temp_in):
    # Variables
    k = 0.8
    rho = 2500
    cp = 800
    h_in = 5
    h_out = 20
    length = 0.2
    dx = length / (nodes - 1)
    dt = 20
    iterations = 10000

    # Calculations
    kappa = (k * dt) / (rho * (dx**2) * cp)
    alpha = (2 * dt) / (rho * dx * cp)
    kdx = k / dx
    temp_plot = np.zeros(iterations)

    for x in range(0, iterations):
        temp_next = temp.copy()

        for i in range(1, (nodes - 1)):
            temp_next[i] = temp[i] + ((temp[i-1] - (2 * temp[i]) + temp[i+1]) * kappa)
        temp_next[0] = ((h_in * (temp_in - temp[0]) + kdx * (temp[0+1] - temp[0])) * alpha) + temp[0]
        temp_next[nodes-1] = ((h_out * (temp_out - temp[nodes-1]) + kdx * (temp[nodes - 2] - temp[nodes-1])) * alpha) + temp[nodes-1]
        temp = temp_next
        temp_plot[x] = temp[2]


    temp = list(np.around(np.array(temp), 2))
    print(temp)
    return temp, temp_plot, iterations, dt


def plot(temp_plot, iterations, dt):
    time = range(iterations)
    time = np.array(time) * (dt / 3600)
    plt.plot(time, temp_plot)
    plt.ylabel('Temperature [C]')
    plt.xlabel('Time [h]')
    plt.show()


if __name__ == '__main__':
    temp, nodes, temp_out, temp_in = startup()
    temp_list, temp_plot, iterations, dt = simulation(temp, nodes, temp_out, temp_in)
    plot(temp_plot, iterations, dt)