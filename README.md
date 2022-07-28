# HPCCF-HPC2

Configuration code for HPC2. For now, this only contains the `modulefiles`, which, in production, 
live at `/software/modules/4.6.1/ucdhpc-20.04/modulefiles`. Short term, we can set up Actions to
spin up a Ubuntu VM with `modules` installed that runs `module display` on changed modules to
check that they are at least syntatically correct; more interesting CI will require automating
the actual software builds as well.
