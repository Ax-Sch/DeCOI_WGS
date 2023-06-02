echo "Setting up Spark Cluster."

# ADJUST THESE PATHS:
spark_master_path="/home/aschmidt/spark-3.1.1-bin-hadoop3.2/sbin/start-master.sh"
spark_worker_path="/home/aschmidt/spark-3.1.1-bin-hadoop3.2/sbin/start-worker.sh"
slurm_id=""

export LC_ALL="en_US.UTF-8"
#export PATH=$(echo "$PATH" | sed -e "s>/gpfs01/share/conda/condabin:/gpfs01/share/conda/bin:>>")

# ADJUST THIS PATH:
export PYTHONPATH="/home/aschmidt/miniconda3/lib/python3.9/site-packages"
HAIL_HOME=$(pip3 show hail | grep Location | awk -F' ' '{print $2 "/hail"}')
export SPARK_CLASSPATH=$HAIL_HOME"/backend/hail-all-spark.jar"		
#SPARK_LOCAL_DIRS=$TMPDIR

# If variables are not set, assign default values; potentially MODIFY to the queue name / constraints of your slurm settings:
queue=${queue:-medium}
hours_to_run=${hours_to_run:-48}
num_workers=${num_workers:-5}
DRIVE_MEM=${DRIVE_MEM:-24}
EXEC_MEM=${EXEC_MEM:-23}
DRIVE_MEM_SUB=$((DRIVE_MEM - 2))
EXEC_MEM_SUB=$((EXEC_MEM - 2))
CPUS_PER_TASK=${CPUS_PER_TASK:-4}

echo "Env variables set."
echo "Start spark master in this job"

#echo $SPARK_LOCAL_DIRS
$spark_master_path

echo "Add addtitional workers, save job-ids and add a bash trap to cancel them, when this process stops"

# you somehow need to figure out the network-address of your master node (this job here). The network address is needed for the communication between master and slaves.
# The part "spark://${SLURMD_NODENAME}.hpc.uni-bonn:7077" has to be set to your [master-network-address], e.g.: spark://[master-network-address]:7077
for i in $(seq 1 $num_workers)
do
   slurm_id="$slurm_id "$(sbatch --time=${hours_to_run}:00:00 -J spark_worker --ntasks=1 --partition="$queue" \
   --mem=${EXEC_MEM}G --cpus-per-task=${CPUS_PER_TASK} \
   --chdir=$(pwd) \
   -x $worker_nodes_excluded \
   --account=ag_ihg_ludwig \
   --wrap "source $(pwd)/workflow/scripts/LOCAL_DIR.sh; $spark_worker_path spark://${SLURMD_NODENAME}.hpc.uni-bonn:7077; sleep ${hours_to_run}h" | awk '{print $NF}')
done


echo "started workers as jobs:"
echo $slurm_id
trap "scancel $slurm_id" EXIT

first_slurm_id=$(echo $slurm_id | cut -d " " -f1)
echo "waiting for first worker to start, job-id: $first_slurm_id"
first_job_active=0

while [ $first_job_active -eq 0 ]
do
echo "first worker not yet active, job-id: $first_slurm_id"
sleep 30s
first_job_active=$(sacct | \
{ grep "$first_slurm_id" || true; } | \
{ grep "RUNNING" || true; } | wc -l)
done

# MODIFY the last line as above (see line 36).
spark_submit_command="spark-submit \
--jars $HAIL_HOME/backend/hail-all-spark.jar \
--conf spark.driver.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
--conf spark.executor.extraClassPath=./backend/hail-all-spark.jar \
--conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
--conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator \
--conf spark.driver.memory=${DRIVE_MEM_SUB}G \
--conf spark.executor.memory=${EXEC_MEM_SUB}G \
--master spark://${SLURMD_NODENAME}.hpc.uni-bonn:7077 "
