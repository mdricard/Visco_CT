from Visco import *
from os import listdir
from os.path import isfile, join

mypath = 'D:/CT REP Data/'

stat_file_path = mypath
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

csv_string = '.csv'
csvFiles = [i for i in onlyfiles if csv_string in i]

subject_string = 'S10'
subject_files = [i for i in csvFiles if subject_string in i]
print(subject_files)

cond_string = 'C2'
subj_cond_files = [i for i in subject_files if cond_string in i]

rom_string = '100%'
subj_cond_rom_files = [i for i in subj_cond_files if rom_string in i]

for i in subj_cond_rom_files:
    line = i.split()
    subject = line[0]
    subject = subject[1:]
    cond = line[1]
    cond = cond[1:]
    rom = line[2]
    trial = line[3]
    trial = trial[1]
    file_name = mypath + i
    p1 = Visco(file_name, subject, cond, rom, trial)
    #p1.graph_residual()
    #p1.plot_all()
    p1.filter_data()
    p1.find_rep(False)  # change to True to view start, peak torque, end for each rep
    p1.analyze_reps()

    #p1.graph_rep(5)
    #p1.save_rep(8, 'D:/Rep 8.csv')
    p1.print_results()
    p1.save_stats_long(stat_file_path)
    p1.graph_all_reps()
    #subject = subj_cond_files[i].split()

#if __name__ == '__main__':
 #   file_name = 'D:/Loran Stiffness/S3 C0 100%  T2.csv'
  #  p1 = Visco(file_name, subject, cond, rom, trial)
   # #p1.graph_residual()
    #p1.plot_all()
    #p1.filter_data()
    #p1.find_rep(False)  # change to True to view start, peak torque, end for each rep
    #p1.analyze_reps()

    #p1.graph_rep(5)
    #p1.save_rep(8, 'D:/Rep 8.csv')
    #p1.print_results()
    #p1.save_stats_long(stat_file_path)
    #p1.graph_all_reps()
