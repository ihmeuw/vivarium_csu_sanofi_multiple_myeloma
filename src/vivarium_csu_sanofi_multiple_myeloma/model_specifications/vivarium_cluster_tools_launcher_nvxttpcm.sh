
    export VIVARIUM_LOGGING_DIRECTORY=/ihme/costeffectiveness/results/vivarium_csu_sanofi_multiple_myeloma/v6.0_treatment_effect/united_states_of_america/2021_06_13_21_56_32/logs/2021_06_13_21_56_32_run/worker_logs
    export PYTHONPATH=/ihme/costeffectiveness/results/vivarium_csu_sanofi_multiple_myeloma/v6.0_treatment_effect/united_states_of_america/2021_06_13_21_56_32:$PYTHONPATH

    /share/code/collijk/miniconda3/envs/vivarium-mm/bin/rq worker -c settings --name ${JOB_ID}.${SGE_TASK_ID} --burst         -w "vivarium_cluster_tools.psimulate.distributed_worker.ResilientWorker"         --exception-handler "vivarium_cluster_tools.psimulate.distributed_worker.retry_handler" vivarium

    