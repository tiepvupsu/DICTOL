import sys 
import os, csv, subprocess
from time import strftime 
"""
submit a script, function:
1. script (no input) 
python submitjob.py 'test1'
2. function 
python submitjob.py 'LRSDL_top' "\\'mySynthetic\\'" 100 5 5 0.1 0.1 0.1 
"""

mem = '8'
walltime = '23:00:00' 

args = sys.argv 


function_name = args[1]
# if function_name is a file name (remove `.m`)
if function_name.find('.') > 0:
    function_name = function_name.split('.')[0]
print function_name, type(function_name)

# print 'Number of arguments:', len(args), 'arguments.'
# print 'Argument List:', str(args)
def getTime():
    return strftime("%m%d_%H%M")
##
jobid_fn = 'jobid.csv'
with open(jobid_fn, "rb") as csvfile:
    reader = csv.reader(csvfile)
    for job in reader:
        joibid_count =  int(job[0])
        # print joibid_count + 1
joibid_count += 1 
pbs_file = str(joibid_count) + '_'+ function_name+'_' + getTime()

nargin = len(args)

######### Write pbs file 
pbs_target = open(pbs_file, 'w+')
pbs_target.write('#PBS -l nodes=1:ppn=1\n')
pbs_target.write('#PBS -l walltime='+walltime+'\n')
pbs_target.write('#PBS -l pmem='+mem+'gb'+'\n')
pbs_target.write('#PBS -j oe\n')
pbs_target.write('#PBS -o pbsout/\n')

pbs_target.write('echo " "\n')
pbs_target.write('echo " "\n')
pbs_target.write('echo "Job started on `hostname` at `date`"\n')
pbs_target.write('module load matlab\n')
pbs_target.write('cd '+ os.getcwd() +'\n')

matlab_str = 'matlab -nosplash -nodisplay -singleCompThread -r ' + \
    '"'+function_name +'('
if nargin == 2:
    matlab_str += ')"\n'
else: 
    for i in range(2, nargin - 1):
        matlab_str += '${arg'+str(i)+'},'
    matlab_str += '${arg'+str(nargin-1)+'})"\n'

pbs_target.write(matlab_str)
pbs_target.write('echo " "\n')
pbs_target.write('echo "Job Ended at `date`"\n')
pbs_target.write('echo " "\n')
pbs_target.close()
########## submit the pbs file 
## prepare qsub_command
qsub_command = 'qsub '
if nargin > 2:
    qsub_command += '-v '
    for i in range(2, nargin - 1):
        qsub_command += 'arg'+str(i)+'='+args[i] +','
    qsub_command += 'arg'+str(nargin-1)+'='+ args[nargin-1]+' '


qsub_command += pbs_file 
## submit qsub_command 
print qsub_command 
# if run_flag:
exit_status = subprocess.call(qsub_command, shell=True)
if exit_status is 1:  # Check to make sure the job submitted
    print "Job {0} failed to submit".format(qsub_command) 

print "Done submitting jobs!"
print "Remove pbs file...",
os.remove(pbs_file)


jobid_target = open(jobid_fn, 'w+')
jobid_target.write(str(joibid_count))
