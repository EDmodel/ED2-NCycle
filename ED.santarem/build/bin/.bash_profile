# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

PATH=$PATH:$HOME/bin

export PATH

NCARG_ROOT=

export NCARG_ROOT

PATH=$PATH:$HOME/bin:$NCARG_ROOT/bin

export PATH

MANPATH=$MANPATH:$NCARG_ROOT/man

module load openmpi

module load hdf5

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/dmedvigy/hdf5-1.8.10-patch1/hdf5/lib

export LD_LIBRARY_PATH
