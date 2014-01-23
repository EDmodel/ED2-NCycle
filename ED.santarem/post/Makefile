F_DOPTS=-g -free -prec-div -fno-omit-frame-pointer -check all -inline-debug-info -traceback -debug extended -check uninit -fp-stack-check -ftrapuv -auto
F_SOPTS=-xT -O3 -fno-alias -free -traceback

F_OPTS=$(F_SOPTS)

HDF5_DIR=/tigress-hsm/dmedvigy/hdf5-1.8.7/hdf5
HDF5_INCS=-I$(HDF5_DIR)/include -I$(HDF5_DIR)/lib
HDF5_LIBS=-L$(HDF5_DIR)/lib -lhdf5 -lm -lhdf5_fortran -lhdf5 -lhdf5_hl -lz


F_COMP=mpif90

tsy:
	$(F_COMP) -c $(F_OPTS) hdf5_utils.f90
	gcc -c -O3 -DUNDERSCORE -DLITTLE -DPC_LINUX1 hdf5_f2c.c $(HDF5_INCS)
	$(F_COMP) -c $(F_OPTS) timeseries-Y.f90
	$(F_COMP) -o timeseries-Y hdf5_utils.o hdf5_f2c.o timeseries-Y.o -i-static $(F_OPTS) $(HDF5_INCS) $(HDF5_LIBS)


tsm:
	$(F_COMP) -c $(F_OPTS) hdf5_utils.f90
	gcc -c -O3 -DUNDERSCORE -DLITTLE -DPC_LINUX1 hdf5_f2c.c $(HDF5_INCS)
	$(F_COMP) -c $(F_OPTS) timeseries-M.f90
	$(F_COMP) -o timeseries-M hdf5_utils.o hdf5_f2c.o timeseries-M.o -i-static $(F_OPTS) $(HDF5_INCS) $(HDF5_LIBS)

asites:
	$(F_COMP) -c $(F_OPTS) hdf5_utils.f90
	gcc -c -O3 -DUNDERSCORE -DLITTLE -DPC_LINUX1 hdf5_f2c.c $(HDF5_INCS)
	$(F_COMP) -c $(F_OPTS) agg_sites.f90
	$(F_COMP) -o agg_sites hdf5_utils.o hdf5_f2c.o agg_sites.o -i-static $(F_OPTS) $(HDF5_INCS) $(HDF5_LIBS)


