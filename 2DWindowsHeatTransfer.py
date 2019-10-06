import numpy as np
import matplotlib.pyplot as plt


def startup():
    print('Using the FEM-method it is possible to simulate the heat distribution through a glass window')
    length_x = float(input('Enter the thickness of the glass pane [m]: '))
    t_in = float(input('Enter the inside temperature of the room [C]: '))
    t_out = float(input('Enter the outside temperature of the surrounding [C]: '))
    ini_t = float(input('Give a initial temperature to help the simulation [C]: '))
    length_y = length_x * 5
    dx = 0.005
    dy = 0.005
    nodes_x = int((length_x / dx) + 1)
    nodes_y = int((length_y / dy) + 1)
    # Function
    t_list = np.ones((nodes_y, nodes_x)) * ini_t


    return t_list, nodes_x, nodes_y, t_out, t_in, length_x


def simulation(t_list, nodes_x, nodes_y, t_out, t_in, length_x):
    # Variables
    k = 0.8
    rho = 2500
    cp = 800
    h_in = 5
    h_out = 20
    #length = 1
    dx = length_x / (nodes_x - 1)
    dy = length_x / (nodes_y - 1)
    dt = 10
    iterations = 100

    # Calculations
    kappa = dt / (rho * dx * dy * cp)
    alpha_2dt = (2 * dt) / (rho * dx * dy * cp)
    kdx_dy = (k * dx) / dy
    kdy_dx = (k * dy) / dx
    kdx_2dy = (k * dx) / (2 * dy)
    kdy_2dx = (k * dy) / (2 * dx)

    t_list[0:(nodes_y), 0] = t_in
    t_list[0:(nodes_y), (nodes_x - 1)] = t_out

    # Startup conditions
    print('Your initial conditions where the numbers represent temperature at each node [C]')
    print(t_list)
    t_plot = list(np.zeros(iterations))
    for x in range(0, iterations):
        t_next_list = t_list.copy()
        for n in range(1, (nodes_x - 1)):
            for m in range(1, (nodes_y - 1)):
                # Inside interior boundary calculations
                # t_next_list[m, 0] = (alpha_2dt * ((kdx_2dy * (t_list[m, n + 1] - 2 * t_list[m, n] + t_list[m, n - 1])) + (kdy_dx * (t_list[m + 1, n] - t_list[m, n])) + (h_in * dy) * (t_in - t_list[m, n]))) + t_list[m, n]
                # Outside interior boundary calculations
                # t_next_list[m, nodes_x - 1] = (alpha_2dt * ((kdx_2dy * (t_list[m, n + 1] - 2 * t_list[m, n] + t_list[m, n - 1])) + (kdy_dx * (t_list[m + 1, n] - t_list[m, n])) + (h_out * dy) * (t_out - t_list[m, n]))) + t_list[m, n]
                # Interior nodes
                t_next_list[m, n] = (((kdy_dx * (t_list[m + 1, n] - (2 * t_list[m, n]) + t_list[m - 1, n])) + (kdx_dy * (t_list[m, n + 1] - (2 * t_list[m, n]) + t_list[m, n - 1]))) * kappa) + t_list[m, n]

        t_list = t_next_list


        for y in range(0, nodes_y):
            t_list[y, 0:nodes_x] = list(np.around(np.array(t_list[y, 0:nodes_x]), 1))
        t_plot[x] = t_list[int(nodes_y / 2), 1]
    print('Simulation results')
    print(t_list)

    return t_list, iterations, dt, t_plot


def plot_result(t_plot, iterations, dt):
    time = range(iterations)
    time = np.array(time) * (dt / 3600)
    plt.plot(time, t_plot)
    plt.ylabel('Temperature [C]')
    plt.xlabel('Time [h]')
    plt.show()


if __name__ == '__main__':
    t_list, nodes_x, nodes_y, t_out, t_in, length_x = startup()
    t_list, iterations, dt, t_plot = simulation(t_list, nodes_x, nodes_y, t_out, t_in, length_x)
    plot_result(t_plot, iterations, dt)




    #temp_list, temp_plot, iterations, dt = simulation(temp, nodes, temp_out, temp_in)
    #plot(temp_plot, iterations,