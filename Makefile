# Universal machine independent makefile for ROMS model 
#========== ======= =========== ======== === ==== =====
# Set machine dependent definitions and rules.

#include Makedefs/Makedefs.Intel
#include Makedefs/Makedefs.IntelEM64T
#include Makedefs/Makedefs.Trestles_MPI
#include Makedefs/Makedefs.Intel_MPI
#include Makedefs/Makedefs.hoffman2
include Makedefs.Gordon_MPI
#include  Makedefs/Makedefs.CMPI
#include  Makedefs/Makedefs.Tungsten
#include Makedefs/Makedefs.Altix
#include Makedefs/Makedefs.p690
#include Makedefs/Makedefs.PGI
#include Makedefs/Makedefs.GCC
#include Makedefs/Makedefs.g95
#include Makedefs/Makedefs.Fujitsu




#                                      All source code files are 
# Configuration for SCRUM/ROMS model:  sorted into eight  groups,
#============== === ========== ======  separated by blanc lines:
#
#  1) main driving part;          5) vertical mixing schemes; 
#  2) 2D time stepping engine;    6) on-fly model diagnostic routines;
#  3) 3D time stepping engine;    7) netCDF I/O routines;
#  4) sea-water EOS routines;     8) model forcing routines; 
#

# Few modifications for the SUNTANS part, JK, Oct/30/2009
#         MSRCS  : mod_transfer.F
#         KRNSRC : extract.F distribute.F
#         IOSRC  : read_transfer.F wrt_transfer.F

KRNSRC = \
	main.F		step2D_FB.F	read_inp.F	sta_par.F\
	timers.F	setup_kwds.F	check_switches1.F\
	check_srcs.F	check_switches2.F\
	set_weights.F	set_scoord.F	init_scalars.F	init_arrays.F\
	setup_grid1.F	setup_grid2.F	visc_rescale.F	set_nudgcof.F\
	analytical.F	ana_grid.F	ana_initial.F\
\
	zetabc.F	u2dbc_im.F	v2dbc_im.F	obc_volcons.F\
\
	pre_step3d.F	step3d_t.F	step3d_uv1.F	step3d_uv2.F\
	prsgrd.F	rhs3d.F	        set_depth.F	omega.F\
	visc3d_GP.F	t3dmix.F	t3dmix_ISO.F	ana_vmix.F\
	u3dbc_im.F	v3dbc_im.F	t3dbc_im.F\
\
	rho_eos.F	ab_ratio.F	alfabeta.F\
\
	lmd_vmix.F      lmd_kpp.F	lmd_swr_frac.F\
	bvf_mix.F\
\
	diag.F		wvlcty.F	grid_stiffness.F        slopes.F\
\
	exchange.F	mpi_setup.F	mpi_exchange4SR.F   mpi_exchange8WA1.F\
\
	lenstr.F	get_date.F


IOSRC = nf_fread.F	nf_read_bry.F	set_cycle.F	checkdims.F\
	insert_node.F	closecdf.F	put_global_atts.F\
	get_grid.F	get_initial.F	def_grid.F	wrt_grid.F\
	def_rst.F	wrt_rst.F	def_his.F	wrt_his.F\
	set_avg.F	wrt_avg.F	def_station.F	wrt_station.F\
\
	get_forces.F	get_all_bry.F	get_tclima.F	get_uclima.F\
	get_ssh.F	get_smflux.F	get_stflux.F	get_srflux.F\
	get_sst.F	get_dqdt.F	get_sss.F\
\
	def_diags_uv.F	wrt_diags_uv.F\
	set_diags_uv_avg.F	wrt_diags_uv_avg.F\
\
	def_diags_vrt.F	wrt_diags_vrt.F\
	set_diags_vrt.F	set_diags_vrt_avg.F	wrt_diags_vrt_avg.F\
\
	def_diags_ts.F	wrt_diags_ts.F\
	set_diags_ts_avg.F	wrt_diags_ts_avg.F\
\
	def_diags_pv.F	wrt_diags_pv.F\
	set_diags_pv_avg.F	wrt_diags_pv_avg.F\
\
	def_diags_eddy.F\
	set_diags_eddy_avg.F	wrt_diags_eddy_avg.F

FLTSRC = floats/init_floats.F	floats/init_arrays_floats.F	floats/get_initial_floats.F\
	floats/interp_rho.F	floats/rhs_floats.F	floats/step_floats.F\
	floats/def_floats.F	floats/wrt_floats.F
#
  SRCS = $(KRNSRC) $(IOSRC) 
#FLOATS ----  SRCS = $(KRNSRC) $(IOSRC) $(FLTSRC) 
   RCS = $(SRCS:.F=.f)
  OBJS = $(RCS:.f=.o)
  SBIN = roms
 LROMS = libroms.a
 LROMS2 = 2/$(LROMS)
 LROMS3 = 2/$(LROMS) 3/$(LROMS)
 LROMS4 = 2/$(LROMS) 3/$(LROMS) 4/$(LROMS)
#
# Eecutable file.
#========== =====
#
$(SBIN): $(OBJS)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o a.out $(OBJS) $(LCDF) $(LMPI)
	mv a.out $(SBIN)
2level: $(OBJS) $(LROMS2)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o a.out $(OBJS) $(LROMS2) $(LCDF) $(LMPI)
	mv a.out $(SBIN)
3level: $(OBJS) $(LROMS3)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o a.out $(OBJS) $(LROMS3) $(LCDF) $(LMPI)
	mv a.out $(SBIN)
4level: $(OBJS) $(LROMS4)
	$(LDR) $(FFLAGS) $(LDFLAGS) -o a.out $(OBJS) $(LROMS4) $(LCDF) $(LMPI)
	mv a.out $(SBIN)
#
# Multilevel libraries
#
$(LROMS): $(OBJS)
	/bin/rm -f $(LROMS)
	ar r $(LROMS) $(OBJS)
2/$(LROMS): $(SRCS)
	cd 2; make -f ./Makefile $(LROMS); cd ..
#
# Everything
#===========
all: tools depend $(SBIN)
#
# Auxiliary utility programs and List of Dependecies:
# ========= ======= ======== === ==== == ============
# 
#  SUNTANS : TOOLS have only mpc now. Add other utilities like partit, ncjoin
#             etc later.
#
 TOOLS = mpc cross_matrix cppcheck srcscheck checkkwds\
	 partit ncjoin ncrename nccp nccpstruct ncextract redefs

tools: $(TOOLS)

    TMP = .

mpc: mpc.F Makefile
	$(CPP) $(CPPFLAGS) -P mpc.F > $(TMP)/mpc.f
	$(LDR) $(FFLAGS) $(LDFLAGS) -o mpc $(TMP)/mpc.f
cross_matrix: TOOLS/cross_matrix.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o cross_matrix TOOLS/cross_matrix.o
partit: TOOLS/partit.o insert_node.o lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o partit TOOLS/partit.o insert_node.o lenstr.o $(LCDF)
ncjoin: TOOLS/ncjoin.o lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o ncjoin TOOLS/ncjoin.o lenstr.o $(LCDF)
ncrename: TOOLS/ncrename.o lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o ncrename TOOLS/ncrename.o lenstr.o $(LCDF)
nccp: TOOLS/nccp.o lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o nccp TOOLS/nccp.o lenstr.o $(LCDF)
nccpstruct: TOOLS/nccpstruct.o lenstr.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o nccpstruct TOOLS/nccpstruct.o lenstr.o $(LCDF)
ncextract: TOOLS/ncextract.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o ncextract TOOLS/ncextract.o
cppcheck: TOOLS/cppcheck.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o cppcheck TOOLS/cppcheck.o
srcscheck: TOOLS/srcscheck.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o srcscheck TOOLS/srcscheck.o
checkkwds: TOOLS/checkkwds.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o checkkwds TOOLS/checkkwds.o
redefs: TOOLS/redefs.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o redefs TOOLS/redefs.o


checkdefs: check_switches1.F setup_kwds.F
check_switches1.F: cppcheck cppdefs.h
	./cppcheck
check_srcs.F: srcscheck Makefile
	./srcscheck
setup_kwds.F: checkkwds read_inp.F
	./checkkwds

depend: checkdefs cross_matrix
	./cross_matrix *.F

#
# Target to create tar file.
# ====== == ====== === =====
#
tarfile: clean
	tar cvf roms.tar Make* *.F *.h *.inc etalon_data.* README.*  *.in *.in.*
#
# Cleaning targets
# ======== =======
#
rmtools:
	/bin/rm -f $(TOOLS)
clean:
	/bin/rm -rf *.f *.o floats/*.f floats/*.o TOOLS/*.f TOOLS/*.o *.a *.s *.trace *.mod *~ $(COMP_FILES) roms 
allclean: clean
	cd 2; make -f ./Makefile clean; cd ..
	cd 3; make -f ./Makefile clean; cd ..
	cd 4; make -f ./Makefile clean; cd ..
clobber: clean
	/bin/rm -f check_switches1.F setup_kwds.F check_srcs.F
	/bin/rm -f $(SBIN) $(TOOLS)
	/bin/rm -f core core.*
#
# Include automatically generated dependency list:
#
include Make.depend

