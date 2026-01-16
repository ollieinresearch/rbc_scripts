cp /home/ollie/scratch/rbc_scripts/run_sim.sh .
cp /home/ollie/scratch/rbc_scripts/post_process.sh .
cp /home/ollie/scratch/rbc_scripts/plotting.sh .

if [ ! -f "$PWD/params.conf" ]; then
    cp /home/ollie/scratch/rbc_scripts/default.conf .
fi