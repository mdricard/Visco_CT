import numpy as np
import matplotlib.pyplot as plt
from BiomechTools import low_pass, zero_crossing, max_min, simpson_nonuniform, critically_damped, residual_analysis
import os.path
from os import path


class Visco:
    n_reps = 0
    energy_absorbed = np.zeros(20)
    energy_returned = np.zeros(20)
    peak_torque = np.zeros(20)
    stiffness = np.zeros(20)

    def __init__(self, fn): #, subject, cond, rom, trial
        with open(fn) as infile:
            temp = infile.readline()
            header = temp.split(',')
            self.n = int(header[0])
            self.hor_limb_wt = float(header[1])
            self.sampling_rate = 1000

        data = np.genfromtxt(fn, delimiter=',', skip_header=12)
        self.pt = data[:, 0]
        self.ga = data[:, 1]
        self.so = data[:, 2]
        self.ta = data[:, 3]
        self.tor = data[:, 4]
        self.pos = data[:, 5] #/57.2957795131 # convert to radians
        #self.angle = self.pos / 57.2957795131
        self.vel = data[:, 6]


        self.smooth_tor = []
        self.smooth_pos = []
        self.smooth_vel = []
        self.rep_start = np.zeros(20, dtype=np.int32)
        self.max_loc = np.zeros(20, dtype=np.int32)
        self.rep_end = np.zeros(20, dtype=np.int32)
        self.peak_torque_pt = np.zeros(20, dtype=np.int32)
        self.stiffness_start_pt = np.zeros(20, dtype=np.int32)

    def plot_all(self):
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.plot(self.pt, self.smooth_tor, 'b-')
        ax2.plot(self.pt, self.smooth_pos, 'r-')
        ax1.set_ylabel('Torque')
        ax2.set_ylabel('Angle') 
        fig.tight_layout()
        plt.grid(True)
        plt.legend(['Torque', 'Angle'])
        plt.show()

    def filter_data(self):
        self.smooth_tor = low_pass(self.tor, self.sampling_rate, 3)
        self.smooth_pos = low_pass(self.pos, self.sampling_rate, 3)
        self.smooth_vel = critically_damped(self.vel, self.sampling_rate, 10)

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

    def get_max_tor(self, rep):
        max_location = self.rep_start[rep]
        max_val = self.smooth_tor[self.rep_start[rep]]
        self.peak_torque_pt[rep] = self.rep_start[rep]     # sometimes the peak torque occurs at first pt, why?
        for i in range(self.rep_start[rep], self.max_loc[rep]):
            if self.smooth_tor[i] > max_val:
                max_val = self.smooth_tor[i]
                self.peak_torque_pt[rep] = i
                self.peak_torque[rep] = self.smooth_tor[i]
        print('Peak Torque  = ',self.peak_torque[rep])
        print('Pk torque Point = ', self.peak_torque_pt[rep])
        print('Rep Start Pt ', self.rep_start[rep])
        print('Max Angle Point ', self.max_loc[rep])
        print('Frank Burns eats worms')

    def find_rep(self, show_graph):
        p_18, r_or_f = zero_crossing(self.smooth_pos, reference_value=20, start=10, stop=self.n - 2)
        if r_or_f[0] == 'falling':      # make sure first point is falling as it pos passes reference (135 deg)
            start_pt = 0
            cnt = len(p_18)
        else:
            start_pt = 1
            cnt = len(p_18) - 1
        rep = 0
        cntr = start_pt       # cntr is used to check for another rep after exiting for loop
        #print(p_18)
        #print('count: ', cnt)
        for i in range(start_pt, cnt-3, 2):
            self.rep_start[rep] = self.get_min(p_18[i], p_18[i + 1])
            self.max_loc[rep] = self.get_max(p_18[i + 1], p_18[i + 2])
            #self.peak_torque_pt[rep] = self.get_max(p_18[i + 1], p_18[i + 2])
            #self.peak_torque[rep] = self.smooth_tor[self.peak_torque_pt[rep]]
            self.rep_end[rep] = self.get_min(p_18[i + 2], p_18[i + 3])
            cntr = cntr + 2
            rep = rep + 1
        # now use cntr to check for one last rep
        #if cntr <= (len(p_18) - 3):
        #    self.rep_start[rep] = self.rep_end[rep - 1]
        #    self.max_loc[rep] = self.get_max(p_18[cntr + 1], p_18[cntr + 2])
        #    #self.peak_torque_pt[rep] = self.get_max(p_18[i + 1], p_18[i + 2])
        #    #self.peak_torque[rep] = self.smooth_tor[self.peak_torque_pt[rep]]
        #    self.rep_end[rep] = self.get_min(p_18[cntr + 2], p_18[cntr + 3])
        #    rep = rep + 1
        self.n_reps = rep
        if show_graph:
            plt.plot(self.pt, self.smooth_pos)
            for rep in range(self.n_reps):
                plt.vlines(self.rep_start[rep], 10, 40, linestyles='solid', colors='green')
                plt.vlines(self.max_loc[rep], 10, 40, linestyles='dashed', colors='red')
                plt.vlines(self.rep_end[rep], 10, 40, linestyles='dotted', colors='black')
            for i in range(len(p_18)):
                plt.scatter(p_18[i], self.smooth_pos[p_18[i]])
            plt.show()

    def find_stiffness_start_pt(self, rep, n_degrees):
        # n_degrees is the number of degrees below the endpoint (typically 5 deg)
        i = self.max_loc[rep] - 10
        while self.smooth_pos[i] > (self.smooth_pos[self.max_loc[rep]] - n_degrees):
            i = i - 1
        self.stiffness_start_pt[rep] = i

    def analyze_reps(self):
        for rep in range(self.n_reps):
            self.find_stiffness_start_pt(rep, 5.0)
            load_area = simpson_nonuniform(self.smooth_pos[self.rep_start[rep] : self.max_loc[rep]], self.smooth_tor[self.rep_start[rep] : self.max_loc[rep]])
            self.energy_returned[rep] = -1.0 * simpson_nonuniform(self.smooth_pos[self.max_loc[rep] : self.rep_end[rep]], self.smooth_tor[self.max_loc[rep] : self.rep_end[rep]])
            self.energy_absorbed[rep] = load_area - self.energy_returned[rep]
            self.get_max_tor(rep)
            self.find_stiffness_start_pt(rep, 10.00)
            self.stiffness[rep] = (self.peak_torque[rep] - self.smooth_tor[self.stiffness_start_pt[rep]]) / (
                np.radians(self.smooth_pos[self.peak_torque_pt[rep]] - self.smooth_pos[self.stiffness_start_pt[rep]]))

    def graph_residual(self):
        residual = residual_analysis(self.tor, self.sampling_rate, 0.5, 20.0)
        freq = np.arange(0.5, 20, 0.5)
        plt.plot(freq, residual)
        plt.show()

    def graph_rep(self, rep):
        x = self.smooth_pos[self.rep_start[rep]: self.rep_end[rep]]
        y = self.smooth_tor[self.rep_start[rep]: self.rep_end[rep]]
        max_tor_pt = self.peak_torque_pt[rep] - self.rep_start[rep]
        pos_minus_5_pt = self.stiffness_start_pt[rep] - self.rep_start[rep]
        plt.plot(x, y)
        plt.scatter(x[max_tor_pt], y[max_tor_pt])
        plt.scatter(x[pos_minus_5_pt], y[pos_minus_5_pt])
        plt.ylabel('Torque (Nm)')
        plt.xlabel('Angle (deg)')
        plt.grid(True)
        plt.show()

    def plot_torque_mmg(self, rep):
        fig, ax1 = plt.subplots()
        ax1.plot(self.smooth_tor[self.rep_start[rep]: self.rep_end[rep]])
        ax2 = ax1.twinx()
        ax2.plot(self.mmg[self.rep_start[rep]: self.rep_end[rep]])
        ax1.set_ylabel('Torque')
        ax2.set_ylabel('mmg')
        fig.tight_layout()
        plt.grid(True)
        plt.show()

    def plot_six(self):
        for rep in range(6):
            fig, ax1 = plt.subplots(2, 3, rep + 1)
            ax1.plot(self.smooth_tor[self.rep_start[rep]: self.rep_end[rep]])
            ax2 = ax1.twinx()
            ax2.plot(self.mmg[self.rep_start[rep]: self.rep_end[rep]])
            ax1.set_ylabel('Torque')
            ax2.set_ylabel('mmg')
            plt.grid(True)
            fig.tight_layout()
        plt.grid(True)
        plt.ylabel('Torque (Nm)')
        plt.xlabel('Angle (r)')
        plt.show()

    def save_rep(self, rep, rep_name):
        p = self.smooth_pos[self.rep_start[rep] : self.rep_end[rep]]
        t = self.smooth_tor[self.rep_start[rep] : self.rep_end[rep]]
        v = self.smooth_vel[self.rep_start[rep] : self.rep_end[rep]]
        mh = self.MHemg[self.rep_start[rep] : self.rep_end[rep]]
        vl = self.VLemg[self.rep_start[rep] : self.rep_end[rep]]
        mmg = self.mmg[self.rep_start[rep] : self.rep_end[rep]]

        np.savetxt(rep_name, np.column_stack((p, t, v, mh, vl, mmg)),  fmt='%.6f', delimiter=',', newline='\n', header="angle, torque, velocity, MH, VL, MMG")

    def graph_all_reps(self):
        for rep in range(self.n_reps):
            x = self.smooth_pos[self.rep_start[rep] : self.rep_end[rep]]
            y = self.smooth_tor[self.rep_start[rep] : self.rep_end[rep]]
            max_tor_pt = self.peak_torque_pt[rep] - self.rep_start[rep]
            pos_minus_5_pt = self.stiffness_start_pt[rep] - self.rep_start[rep]
            plt.plot(x, y)
            plt.scatter(x[max_tor_pt], y[max_tor_pt])
            plt.scatter(x[pos_minus_5_pt], y[pos_minus_5_pt])
            plt.grid
            plt.ylabel('Torque (Nm)')
            plt.xlabel('Angle (r)') 
        plt.show()

    def print_results(self):
        for rep in range(self.n_reps):
            print('Peak Torque (Nm): ' + format(self.peak_torque[rep], '.2f'))
            print('Stiffness (Nm/rad): ' + format(self.stiffness[rep], '.2f'))
            print('Energy Absorbed (Nm): ' + format(self.energy_absorbed[rep], '.2f'))
            print('Energy Returned (Nm): ' + format(self.energy_returned[rep], '.2f'))

    def save_stats_long(self, stat_file_path):
        fn = stat_file_path + 'Loran Stats LONG.csv'
        with open(fn, 'a') as stat_file:
            for rep in range(self.n_reps):
                stat_file.write(
                    self.subject + ',' + self.cond + ',' + self.rom + ',' + self.trial + ',' + str(rep) + ',' + str(self.peak_torque[rep]) + ',' + str(self.stiffness[rep]) + ',' + str(self.energy_absorbed[rep]) + ',' + str(self.energy_returned[rep]) + '\n')
        stat_file.close()
