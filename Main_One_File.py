#from Trowbridge import *
from Visco import *
from os import listdir
from os.path import isfile, join

mypath = 'D:/CT REP Data/'
fn = 'S1PS10WDPre.txt'
file_name = mypath + fn
p1 = Visco(file_name)
#p1 = Trowbridge(file_name) #, subject, cond, rom, trial)
#p1.graph_residual()
#p1.plot_angle_deg_rad()
p1.filter_data()
p1.plot_all()
p1.find_rep(True)  # change to True to view start, peak torque, end for each rep

p1.analyze_reps()

p1.graph_rep(5)
#p1.save_rep(8, 'D:/Rep 8.csv')
p1.print_results()
#p1.save_stats_long(stat_file_path)
p1.graph_all_reps()