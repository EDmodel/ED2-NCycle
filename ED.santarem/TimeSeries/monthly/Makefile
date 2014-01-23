F_COMP=mpif90
F_OPTS=-O3 -static -xHost -traceback
#F_OPTS=-check bounds -traceback -check uninit -ftrapuv -static
PAR_INCS=
OUTILS=
HDF5_DIR=/home/dmedvigy/hdf5-1.8.10-patch1/hdf5
HDF5_INCS=-I$(HDF5_DIR)/include -I$(HDF5_DIR)/lib
HDF5_LIBS=-L$(HDF5_DIR)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -lm
INCLUDES=$(HDF5_INCS)
F_COMMAND=$(F_COMP) -c $(F_OPTS) $(INCLUDES)

mts:
	$(F_COMMAND) -DOLAM_MPI hdf5_f2f.F90
	$(F_COMMAND) -DOLAM_MPI hdf5_utils.F90
	$(F_COMMAND) daily_ts_template.f90
	$(F_COMP) -o daily_ts_template hdf5_utils.o hdf5_f2f.o daily_ts_template.o $(HDF5_LIBS)

eval:
	$(F_COMMAND) -DOLAM_MPI hdf5_f2f.F90
	$(F_COMMAND) -DOLAM_MPI hdf5_utils.F90
	$(F_COMMAND) monthly_ts_eval.f90
	$(F_COMP) -o monthly_ts_eval hdf5_utils.o hdf5_f2f.o monthly_ts_eval.o $(HDF5_LIBS)