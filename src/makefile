#
# Makefile for FIRTEZ-dz
#-------------------------------------------
# A. Pastor Yabar & J.M. Borrero August 2019
#
#
#===========================================
# To be edited by the USER:
#
LIBS = \
        -L /Users/apy/libraries/lib -llapack -lblas \
        -L /opt/local/include/ -lfftw3
#
INCS = -I/opt/local/include/
#
LDFLAGS = 
#
#===========================================
# From here on, you should not need to edit anything...
#
#
# Compiler dependent:
# 
ifeq ($(COMP),intel)
# Intel compiler:
  #FC = mpiifort
  FC = mpif90
  #CC = mpiicc
#
  FPP = -fpp
  # Debugging: FOPT = -fopenmp -traceback -debug all -check all -g -O0
#
  FR8 = -freal-4-real-8
#
else
# Gfortran compiler:
  FC = mpif90
  #CC = mpicc
#
  FPP = -cpp
  # Debugging: FOPT = -fopenmp -fbacktrace -Wall -fcheck=all -g -finit-real=nan -finit-integer=0 -O0
#
  FR8 = -fdefault-real-8
#
endif

#
# Common to both compiler
# 
FOPT=-O3 -fopenmp
F90F = -c
F77F = -c
#
#
#

OBJ = \
	user_mpi.o \
	lm.o \
	checknan.o \
	voigt.o \
        cons_param.o \
	time_param.o \
	invert_param.o \
	height_handler.o \
	phys_param.o \
	grid_param.o \
	coupled_param.o \
        forward_param.o \
	user_fftw3.o \
        atm_param.o \
	code_modes.o \
	derivvar.o \
	allocate_utils.o \
	misc.o \
	splines.o \
	log.o \
        atom_database.o \
	model_spatial_extension.o \
	join_model.o \
        damping.o \
        lines_database.o \
	magsplit.o \
	absorption_matrix.o \
	chemical_equilibrium.o \
	line_opacity.o \
	cont_opacity.o \
	normalization.o \
	rtesolver.o \
	pre_post_duties.o \
	input_file.o \
        hydrostatic_tools.o \
	forward.o \
	coupled_inversion.o \
	coupled_matrix_inversion.o \
	get_dmodel.o \
	invert3d.o \
	main.o

user_mpi.o: user_mpi.f90
	$(FC) $(FPP) $(F90F) $(FOPT) user_mpi.f90

lm.o: lm.f90
	$(FC) $(F90F) $(FOPT) lm.f90

checknan.o: checknan.f90
	$(FC) $(F90F) $(FOPT) checknan.f90

voigt.o: voigt.f
	$(FC) $(F77F) $(FOPT) $(FR8) voigt.f

cons_param.o: cons_param.f90
	$(FC) $(F90F) $(FOPT) cons_param.f90

time_param.o: time_param.f90 cons_param.o
	$(FC) $(F90F) $(FOPT) time_param.f90

invert_param.o: invert_param.f90 cons_param.o
	$(FC) $(F90F) $(FOPT) invert_param.f90

height_handler.o: height_handler.f90 cons_param.o
	$(FC) $(F90F) $(FOPT) height_handler.f90

phys_param.o: phys_param.f90 cons_param.o
	$(FC) $(F90F) $(FOPT) phys_param.f90

grid_param.o: grid_param.f90 cons_param.o
	$(FC) $(F90F) $(FOPT) grid_param.f90

coupled_param.o: coupled_param.f90 cons_param.o
	$(FC) $(F90F) $(FOPT) coupled_param.f90

forward_param.o: forward_param.f90 cons_param.o
	$(FC) $(F90F) $(FOPT) forward_param.f90

user_fftw3.o: user_fftw3.f90 cons_param.o
	$(FC) $(F90F) $(FOPT) $(INCS) user_fftw3.f90

atm_param.o: atm_param.f90 cons_param.o grid_param.o \
        forward_param.o
	$(FC) $(F90F) $(FOPT) atm_param.f90

code_modes.o: code_modes.f90 cons_param.o
	$(FC) $(F90F) $(FOPT) code_modes.f90

derivvar.o: derivvar.f90 cons_param.o
	$(FC) $(F90F) $(FOPT) derivvar.f90

allocate_utils.o: allocate_utils.f90 cons_param.o
	$(FC) $(F90F) $(FOPT) allocate_utils.f90

misc.o: misc.f90 cons_param.o allocate_utils.o
	$(FC) $(F90F) $(FOPT) misc.f90

splines.o: splines.f90 cons_param.o
	$(FC) $(F90F) $(FOPT) splines.f90

log.o: log.f90 cons_param.o
	$(FC) $(F90F) $(FOPT) log.f90

atom_database.o: atom_database.f90 cons_param.o
	$(FC) $(F90F) $(FOPT) atom_database.f90

model_spatial_extension.o: model_spatial_extension.f90 \
        cons_param.o grid_param.o forward_param.o \
        phys_param.o allocate_utils.o user_mpi.o \
        invert_param.o
	$(FC) $(F90F) $(FOPT) model_spatial_extension.f90

join_model.o: join_model.f90 cons_param.o grid_param.o \
        forward_param.o phys_param.o
	$(FC) $(F90F) $(FOPT) join_model.f90

damping.o: damping.f90 cons_param.o splines.o\
        code_modes.o log.o derivvar.o atom_database.o
	$(FC) $(F90F) $(FOPT) damping.f90

lines_database.o: lines_database.f90 cons_param.o \
        forward_param.o allocate_utils.o log.o damping.o
	$(FC) $(F90F) $(FOPT) lines_database.f90

magsplit.o: magsplit.f90 cons_param.o \
        forward_param.o allocate_utils.o \
        code_modes.o
	$(FC) $(F90F) $(FOPT) magsplit.f90

absorption_matrix.o: absorption_matrix.f90 cons_param.o \
        forward_param.o allocate_utils.o code_modes.o \
        atm_param.o derivvar.o magsplit.o atom_database.o \
        lines_database.o
	$(FC) $(F90F) $(FOPT) absorption_matrix.f90

chemical_equilibrium.o: chemical_equilibrium.f90 cons_param.o \
        code_modes.o invert_param.o derivvar.o atom_database.o \
        lines_database.o
	$(FC) $(F90F) $(FOPT) chemical_equilibrium.f90

cont_opacity.o: cont_opacity.f90 cons_param.o \
        code_modes.o derivvar.o atom_database.o \
        lines_database.o chemical_equilibrium.o
	$(FC) $(F90F) $(FOPT) cont_opacity.f90

normalization.o: normalization.f90 cons_param.o \
        forward_param.o allocate_utils.o \
        chemical_equilibrium.o cont_opacity.o
	$(FC) $(F90F) $(FOPT) normalization.f90

line_opacity.o: line_opacity.f90 cons_param.o \
        code_modes.o derivvar.o atom_database.o \
        lines_database.o chemical_equilibrium.o
	$(FC) $(F90F) $(FOPT) line_opacity.f90

rtesolver.o: rtesolver.f90 cons_param.o \
        code_modes.o atm_param.o absorption_matrix.o \
        phys_param.o
	$(FC) $(F90F) $(FOPT) rtesolver.f90

coupled_inversion.o: coupled_inversion.f90 coupled_param.o \
	invert_param.o grid_param.o forward_param.o \
	allocate_utils.o height_handler.o
	$(FC) $(F90F) $(FOPT) coupled_inversion.f90

coupled_matrix_inversion.o: coupled_matrix_inversion.f90 \
	allocate_utils.o misc.o cons_param.o
	$(FC) $(F90F) $(FOPT) coupled_matrix_inversion.f90
#Ready
input_file.o: input_file.f90 cons_param.o grid_param.o \
        forward_param.o allocate_utils.o \
        code_modes.o user_mpi.o log.o invert_param.o \
        lines_database.o coupled_param.o
	$(FC) $(F90F) $(FOPT) input_file.f90
#Ready
pre_post_duties.o: pre_post_duties.f90 grid_param.o \
        forward_param.o phys_param.o allocate_utils.o \
        code_modes.o user_mpi.o log.o invert_param.o \
        atom_database.o lines_database.o absorption_matrix.o \
        phys_param.o time_param.o misc.o \
        normalization.o coupled_param.o coupled_inversion.o \
	user_fftw3.o input_file.o cons_param.o
	$(FC) $(F90F) $(FOPT) pre_post_duties.f90

hydrostatic_tools.o: hydrostatic_tools.f90 cons_param.o \
        chemical_equilibrium.o
	$(FC) $(F90F) $(FOPT) hydrostatic_tools.f90

forward.o: forward.f90 cons_param.o grid_param.o \
        forward_param.o phys_param.o \
        code_modes.o user_mpi.o damping.o invert_param.o \
        atm_param.o chemical_equilibrium.o cont_opacity.o \
        absorption_matrix.o phys_param.o rtesolver.o \
        time_param.o hydrostatic_tools.o line_opacity.o
	$(FC) $(F90F) $(FOPT) forward.f90

get_dmodel.o: get_dmodel.f90 cons_param.o grid_param.o \
        forward_param.o phys_param.o \
        code_modes.o lm.o user_mpi.o invert_param.o \
        atm_param.o phys_param.o forward.o \
	coupled_inversion.o lm.o height_handler.o
	$(FC) $(F90F) $(FOPT) get_dmodel.f90

invert3d.o: invert3d.f90 phys_param.o \
        code_modes.o user_mpi.o invert_param.o \
        forward.o misc.o get_dmodel.o \
        model_spatial_extension.o coupled_inversion.o \
	get_dmodel.o
	$(FC) $(F90F) $(FOPT) invert3d.f90

main.o: main.f90 user_mpi.o invert3d.o
	$(FC) $(F90F) $(FOPT) main.f90

firtez-dz: $(OBJ)
	$(FC) $(OBJ) $(LDFLAGS) $(LIBS) -o firtez-dz.x -fopenmp

clean:
	rm -f *.o *.a *.x *~ *.mod
