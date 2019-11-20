# .bashrc
# set locale
export LC_CTYPE="en_US.utf-8" 

# add miniconda/bin to the PATH when in execute node
case "${HOSTNAME}" in
    med-login*)
        ;;
    *)
        export PATH=$HOME/work/miniconda/bin:$HOME/.local/bin:$PATH:
        export TMPDIR=$HOME/scratch/tmp
        # execute code at start up ::
        conda activate base
        echo "You are now in an execute node"
        ;;
esac

# Function qstat2() -----------------------------------------------------------
#
# Print output similar to qstat but do not truncate the "name" column.

qstat2()
{
    {
        echo -e "job-ID\tprior\tname\tuser\tstate\tsubmit/start at\tqueue\tslots\tja-task-ID"
        qstat -xml "$@" \
        | tr '\n' ' ' \
        | sed 's#<job_list[^>]*>#\n#g' \
        | sed 's#<[^>]*>##g' \
        | grep " "  \
        | sed 's/^ \+//g' \
        | sed -e 's/ \+/\t/g'
    } \
    | column -t -s $'\t'
}

export -f qstat2
