source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setups
setup gcc v6_4_0
setup boost v1_66_0 -q e15:prof
setup cmake v3_13_1

swdir_fnal=/dune/app/users/rodriges/sw
swdir_ox=/home/rodrigues/dune/daq/primsim/external
swdir=""
[ -e "$swdir_fnal" ] && swdir=$swdir_fnal
[ -e "$swdir_ox" ] && swdir=$swdir_ox
[ -z "$swdir" ] && echo Found no external software directory
prepend LD_LIBRARY_PATH $swdir/lib
prepend PATH $swdir/bin
prepend MANPATH $swdir/share/man
unset swdir swdir_fnal swdir_ox
