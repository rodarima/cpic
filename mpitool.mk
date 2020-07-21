MPITOOL?=mpitool

$(if $(MPI_COMPILER),,$(error Required variable MPI_COMPILER not defined))
$(if $(MPI_LINKER),,$(error Required variable MPI_LINKER not defined))

MPI_INCLUDE:=$(shell $(SHELL) $(MPITOOL) --compile $(MPI_COMPILER) ^-I)
$(if $(MPI_INCLUDE),,$(error $(MPITOOL) failed))

MPI_LDFLAGS:=$(shell $(SHELL) $(MPITOOL) --link $(MPI_LINKER) ^-L)
$(if $(MPI_LDFLAGS),,$(error $(MPITOOL) failed))

MPI_LDLIBS:=$(shell $(SHELL) $(MPITOOL) --link $(MPI_LINKER) ^-l)
$(if $(MPI_LDLIBS),,$(error $(MPITOOL) failed))
