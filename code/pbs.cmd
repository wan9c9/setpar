#!/bin/sh
###
### This is a sample PBS job script using MPICH-2
### To use MPICH-2
### you must add the following line to the end of your .bashrc file
### . /share1/bin/mpich2.sh
### you also have to make a .mpd.conf file containing MPD secret word 
### > cd $HOME
### > touch .mpd.conf
### > chmod 600 .mpd.conf
### > echo "MPD_SECRETWORD=any-secret-word" >> .mpd.conf

### Job name
#PBS -N tpar
### Declare job non-rerunable
#PBS -r n


### Queue name (parallel or oneday)

### parallel is the queue for running production jobs.
### Each job in this queue can use 1-3 nodes.
### Parallel jobs will be favoured by the system.

#####PBS -q oneday
#PBS -q stat-wkli

### Wall time required. This example is 30 min
#PBS -l walltime=8:00:00

### Number of nodes 

### The following means 1 node and 1 core. 
### Clearly, this is for serial job
###PBS -l nodes=1:ppn=1

### The following means 1 nodes required. Processor Per Node=8, 
### i.e., total 8 cores to be allocated. 
### ppn (Processor Per Node) can be either 1 or 2 or 4 or 8.
#PBS -l nodes=1:ppn=1

### Another example 
### 2 nodes required. Processor per node=8, total 16 cores 
### need to be allocated.  
###PBS -l nodes=2:ppn=8


#The following stuff will be executed in the first allocated node.
#Please don't modify it

echo $PBS_JOBID : `wc -l < $PBS_NODEFILE` CPUs allocated: `cat $PBS_NODEFILE`
cd $PBS_O_WORKDIR    
PATH=$PBS_O_PATH
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
NNODES=`uniq $PBS_NODEFILE | wc -l`


# Create nodefile (with duplicates removed) for use by mpdboot

JID=`echo ${PBS_JOBID} | sed "s/.gridpoint.gridpoint//" `
NODEFILE=/tmp/node.job.$JID
MACHFILE=/tmp/machine.job.$JID
cat $PBS_NODEFILE | /usr/bin/uniq | sed -e "s/.hku.hk//g"  > $NODEFILE
cat $PBS_NODEFILE | /usr/bin/uniq | sed -e "s/.hku.hk//g" | /bin/awk '{print $1":8"}' > $MACHFILE

# Define MPICH path (Uncomment if neccessary)
##MPIPATH=/share1/mpich2-1.2.1/bin/
##MPIPATH=/share1/pgi/linux86-64/2010/mpi2/mpich/bin/

### Boot the MPI2 engine.
echo "mpdboot...."
${MPIPATH}mpdboot -v -n ${NNODES} --rsh=/usr/bin/rsh --file=${NODEFILE} 


echo "mpdtrace...."
${MPIPATH}mpdtrace -l
echo ===========================================================
echo "Job Start  Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

# Run the parallel MPI executable "a.out"
${MPIPATH}mpiexec -n $NPROCS ./main.exe >> ${PBS_JOBNAME}.`echo ${PBS_JOBID} | sed "s/.gridpoint.gridpoint//" `

echo "Job Finish Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

# Shut down the MPI2 engine and exit the PBS script.
echo "mpdallexit...."
${MPIPATH}mpdallexit

rm -f  ${NODEFILE}
rm -f  ${MACHFILE}

exit 0

