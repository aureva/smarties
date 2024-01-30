.POSIX:
.SUFFIXES: .o .cpp .f90
.PHONY: clean all

PREFIX = $(SMARTIES_ROOT)
CXXFLAGS = -O3 -g
FCFLAGS = -O3 -L/com/Modules/modulefiles/fftw -lfftw3 -fallow-argument-mismatch -fbacktrace #-fwhole-file -fcheck=all -fbacktrace -Wall  -Wline-truncation  -Wcharacter-truncation -Wsurprising  -Waliasing #-fimplicit-none  -fwhole-file -fcheck=all -fbacktrace -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing -Wimplicit-interface  -Wunused-parameter   -pedantic 

FC = mpifort
CXX = mpicxx
LINK = $(FC)

M_CXXFLAGS = -I$(PREFIX)/include -L/com/Modules/modulefiles/fftw -lfftw3
M_FCFLAGS = -I$(PREFIX)/include -I.
MPI_EXTRA_LIB = -lstdc++
M_LDFLAGS = -L$(PREFIX)/lib -lsmarties -L$(PREFIX)/extern/lib $(MPI_EXTRA_LIB)
M_CHFLAGS = -L/com/Modules/modulefiles/fftw -lfftw3

C = \
	smarties.o\
        smarties_stat.o\
	scalar_slice.o \
	avgslice.o calcbeta.o \
        convec.o derivwall2.o ddx.o ddy.o \
        ddz_uv_p.o ddz_uv.o ddz_w.o \
        dealias1.o dealias2.o \
        divstress_uv.o divstress_w.o filt_da.o \
        filter_2laa.o filter_2dsl.o \
        openfiles.o output_tauw.o output_actions.o\
        optim_lag.o optim_scl_lag.o \
        scalar_rhs.o rmsdiv.o \
        root.o stepbl_w.o root8.o newroots.o \
        sgs_stag.o spectrum.o step_s.o stepbl_uv.o  \
        tridag.o wallstress2.o surf_flux_gabls3.o \
        s_zeroslice.o zeroslice.o output_average.o\
        interp3D.o lagrng_sd.o lagrng_sd_s.o \
        update1.o update3.o update9.o surf_flux.o \
        optim_pln.o optim_scl_pln.o \
        optim_lag_dyn.o optim_scl_lag_dyn.o filter_la.o \
        lagrng_dyn.o lagrng_dyn_s.o \
        caldrag.o surf_flux2.o \
        inflow.o  inflowS.o mgm.o mgm_scalar.o \
        press_stag_parallel.o amd.o amd_scalar.o

O = \
app_main.o\
main.o\
smarties_extern.o\

M = \
app_main\


all: $C $M


pre-build:
ifeq (,$(wildcard ./smarties.f90))
	ln -s $(SMARTIES_ROOT)/include/smarties.f90
endif

pre-build2:
ifeq (,$(wildcard ./smarties_extern.cpp))
	ln -s $(SMARTIES_ROOT)/source/smarties/smarties_extern.cpp
endif


$M: $O

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(M_CXXFLAGS) $<

.f.o:
	$(FC) -c $(FCFLAGS) $(M_FCFLAGS) $(M_CHFLAGS) $< -lmpifort

.f90.o:
	$(FC) -c $(FCFLAGS) $(M_FCFLAGS) $(M_CHFLAGS) $< -lmpifort 

.o:
	$(LINK) $C $O $(LDFLAGS) $(M_LDFLAGS) $(M_CHFLAGS) -o $@


clean:; rm -f -- $M $O $C *.mod


