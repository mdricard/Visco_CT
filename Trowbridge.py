import numpy as np
import matplotlib.pyplot as plt
from BiomechTools import low_pass, zero_crossing, max_min, simpson_nonuniform, critically_damped
import math

class Trowbridge:
    n_reps = 0
    energy_absorbed = np.zeros(20)
    energy_returned = np.zeros(20)
    peak_torque = np.zeros(20)
    stiffness = np.zeros(20)

    def __init__(self, fn):
        with open(fn) as infile:
            temp = infile.readline()
            header = temp.split(',')
            self.n = int(header[0])
            self.hor_limb_wt = float(header[1])
            self.sampling_rate = 1000

        data = np.loadtxt(fn, delimiter=',', skiprows=12)
        self.pt = data[:, 0]
        self.emg1 = data[:, 1]
        self.emg2 = data[:, 2]
        self.emg3 = data[:, 3]
        self.tor = data[:, 4]
        self.pos = data[:, 5]#/57.2957795131 # convert to radians
        self.angle = self.pos / 57.2957795131
        self.vel = data[:, 6]
        self.smooth_tor = []
        self.smooth_pos = []
        self.smooth_vel = []
        self.rep_start = np.zeros(20, dtype=np.int32)
        self.max_loc = np.zeros(20, dtype=np.int32)
        self.rep_end = np.zeros(20, dtype=np.int32)


    def plot_all(self):
        plt.plot(self.pt, self.smooth_pos)
        plt.show()

    def filter_data(self):
        self.smooth_tor = low_pass(self.tor, self.sampling_rate, 5)
        self.smooth_pos = low_pass(self.pos, self.sampling_rate, 5)
        self.smooth_vel = critically_damped(self.vel, self.sampling_rate, 5)

    def get_min(self, first_pt, last_pt):
        min_location = first_pt
        min_val = self.smooth_pos[first_pt]
        for i in range(first_pt, last_pt):
            if self.smooth_pos[i] < min_val:
                min_location = i
                min_val = self.smooth_pos[i]
        return min_location

    def get_max(self, first_pt, last_pt):
        max_location = first_pt
        max_val = self.smooth_pos[first_pt]
        for i in range(first_pt, last_pt):
            if self.smooth_pos[i] > max_val:
                max_location = i
                max_val = self.smooth_pos[i]
        return max_location
    
    def plot_angle_deg_rad(self):
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.plot(self.pt, self.angle, 'g-')
        ax2.plot(self.pt, self.pos, 'b-')
        ax1.set_xlabel('time (s)')
        ax1.set_ylabel('Angle (rad)', color='g')
        ax2.set_ylabel('Angle (deg)', color='b')   
        plt.legend(['Angle (rad)', 'Angle (deg)']) 
        plt.show()

    def find_rep(self, show_graph):
        p_18, r_or_f = zero_crossing(self.smooth_pos, 18, 10, self.n - 2)
        rep = 0
        cnt = len(p_18)
        for i in range(0, cnt-3, 2):
            self.rep_start[rep] = self.get_min(p_18[i], p_18[i + 1])
            self.max_loc[rep] = self.get_max(p_18[i + 1], p_18[i + 2])
            self.rep_end[rep] = self.get_min(p_18[i + 2], p_18[i + 3])
            rep = rep + 1
        self.n_reps = rep
        if show_graph:
            plt.plot(self.pt, self.smooth_pos)
            for rep in range(self.n_reps):
                plt.vlines(self.rep_start[rep], 10, 40, linestyles='dotted', colors='black')
                plt.vlines(self.max_loc[rep], 10, 40, linestyles='dashed', colors='red')
                plt.vlines(self.rep_end[rep], 10, 40, linestyles='dotted', colors='black')
            for i in range(len(p_18)):
                plt.scatter(p_18[i], self.smooth_pos[p_18[i]])
            plt.show()

    def analyze_reps(self):
        for rep in range(self.n_reps):
            load_area = simpson_nonuniform(self.smooth_pos[self.rep_start[rep] : self.max_loc[rep]], self.smooth_tor[self.rep_start[rep] : self.max_loc[rep]])
            self.energy_returned[rep] = -1.0 * simpson_nonuniform(self.smooth_pos[self.max_loc[rep] : self.rep_end[rep]], self.smooth_tor[self.max_loc[rep] : self.rep_end[rep]])
            self.energy_absorbed[rep] = load_area - self.energy_returned[rep]
            self.peak_torque[rep] = self.smooth_tor[self.max_loc[rep]]
            self.stiffness[rep] = (self.smooth_tor[self.max_loc[rep]] - self.smooth_tor[self.rep_start[rep]]) / (
                np.radians(self.smooth_pos[self.max_loc[rep]] - self.smooth_pos[self.rep_start[rep]]))


    def graph_rep(self, rep):
        x = self.smooth_pos[self.rep_start[rep] : self.rep_end[rep]]
        y = self.smooth_tor[self.rep_start[rep] : self.rep_end[rep]]
        plt.plot(x, y)
        plt.scatter(x[0], y[0])     # identify the start of the repetition
        plt.show()

    def save_rep(self, rep, rep_name):
        p = self.smooth_pos[self.rep_start[rep] - 200: self.rep_end[rep] + 200] 
        t = self.smooth_tor[self.rep_start[rep] - 200: self.rep_end[rep] + 200]
        t = self.smooth_tor[self.rep_start[rep] - 200: self.rep_end[rep] + 200]
        v = self.smooth_vel[self.rep_start[rep] - 200: self.rep_end[rep] + 200]
        pnts = np.arange(self.rep_start[rep]- 200, self.rep_end[rep] + 200, 1)
        np.savetxt(rep_name, np.column_stack((pnts, p, t, v)),  fmt='%.6f', delimiter=',', newline='\n', header="point, angle, torque, velocity")

    def graph_all_reps(self):
        for rep in range(self.n_reps):
            x = self.smooth_pos[self.rep_start[rep] : self.rep_end[rep]]
            y = self.smooth_tor[self.rep_start[rep] : self.rep_end[rep]]
            plt.plot(x, y)
            plt.ylabel('Torque (Nm)')
            plt.xlabel('Position (degrees)')
            plt.grid(True)
            plt.legend(['Rep 1', 'Rep 2', 'Rep 3', 'Rep 4', 'Rep 5', 'Rep 6', 'Rep 7', 'Rep 8', 'Rep 9', 'Rep 10'])
        plt.show()

    def print_results(self):
        for rep in range(self.n_reps):
            print('Peak Torque (Nm): ' + format(self.peak_torque[rep], '.2f'))
            print('Stiffness (Nm/rad): ' + format(self.stiffness[rep], '.2f'))
            print('Energy Absorbed (Nm): ' + format(self.energy_absorbed[rep], '.2f'))
            print('Energy Returned (Nm): ' + format(self.energy_returned[rep], '.2f'))
