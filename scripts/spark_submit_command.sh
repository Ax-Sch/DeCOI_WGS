echo "Setting up Spark Cluster."

spark_master_path="/ceph01/scratch/aschmidt/spark/spark-3.1.1-bin-hadoop3.2/sbin/start-master.sh"
spark_worker_path="/ceph01/scratch/aschmidt/spark/spark-3.1.1-bin-hadoop3.2/sbin/start-worker.sh"
slurm_id=""

export LC_ALL="en_US.UTF-8"
export PATH=$(echo "$PATH" | sed -e "s>/gpfs01/share/conda/condabin:/gpfs01/share/conda/bin:>>")
export PYTHONPATH="/home/aschmidt/scratch/conda/lib/python3.9/site-packages"
HAIL_HOME=$(pip3 show hail | grep Location | awk -F' ' '{print $2 "/hail"}')
export SPARK_CLASSPATH=$HAIL_HOME"/backend/hail-all-spark.jar"		

# If variables are not set, assign default values
queue=${queue:-batch}
hours_to_run=${hours_to_run:-12}
num_workers=${num_workers:-5}

echo "Env variables set."
echo "Start spark master in this job"

$spark_master_path

echo "Add addtitional workers, save job-ids and add a bash trap to cancel them, when this process stops"

for i in $(seq 1 $num_workers)
do
   slurm_id="$slurm_id "$(sbatch --time=${hours_to_run}:00:00 -J spark_worker --ntasks=1 --partition="$queue" \
   --mem=24G --cpus-per-task=4 \
   --chdir=$(pwd) \
   --gres localtmp:140G \
   -x $worker_nodes_excluded \
   --wrap "source $(pwd)/scripts/LOCAL_DIR.sh; $spark_worker_path spark://${SLURMD_NODENAME}.cluster.imbie:7077; sleep ${hours_to_run}h" | awk '{print $NF}')
done

echo "started workers as jobs:"
echo $slurm_id
trap "scancel $slurm_id" EXIT
sleep 20s

spark_submit_command="spark-submit \
--jars $HAIL_HOME/backend/hail-all-spark.jar \
--conf spark.driver.extraClassPath=$HAIL_HOME/backend/hail-all-spark.jar \
--conf spark.executor.extraClassPath=./backend/hail-all-spark.jar \
--conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
--conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator \
--conf spark.driver.memory=40G \
--conf spark.executor.memory=40G \
--master spark://${SLURMD_NODENAME}.cluster.imbie:7077 "
