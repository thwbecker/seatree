include makefile.system

USE_FLAGS = -O2 $(F_EXT_SOURCE_FLAG)  -fdollar-ok

all:	dirs progs

dirs: 
	if [ ! -s $(ARCH)/ ];then \
		mkdir $(ARCH)/; \
	fi


progs: $(ARCH)/voxint $(ARCH)/data_a2b $(ARCH)/joint_lsqr_vx_isodamp \
	$(ARCH)/vox2roughness $(ARCH)/mapview_3d

# voxint
VOX_OBJS = $(ARCH)/voxint.m1.o $(ARCH)/delazs.m1.o $(ARCH)/gceq.m1.o \
	$(ARCH)/param.m1.o $(ARCH)/path_vx.m1.o $(ARCH)/span.m1.o \
	$(ARCH)/isqre.m1.o $(ARCH)/linint.m1.o \
	$(ARCH)/rsoinc.m1.o 


LSQR_OBJS =  $(ARCH)/joint_lsqr_vx_isodamp.m2.o $(ARCH)/main.m2.o \
	$(ARCH)/readmatrix.m2.o $(ARCH)/gradamp.m2.o 


MV3D_OBJS = $(ARCH)/mapview_3d.m3.o $(ARCH)/param.m3.o 

$(ARCH)/voxint: $(VOX_OBJS)
	$(F77) $(FLAGSUSE) -o $(ARCH)/voxint $(VOX_OBJS)

$(ARCH)/data_a2b: data_a2b.f
	$(F77) $(USE_FLAGS) -o $(ARCH)/data_a2b data_a2b.f

$(ARCH)/mapview_3d: $(MV3D_OBJS)
	$(F77) $(USE_FLAGS) -o $(ARCH)/mapview_3d $(MV3D_OBJS)

$(ARCH)/vox2roughness: $(ARCH)/vox2roughness.m2.o 
	$(F77) $(USE_FLAGS) -o $(ARCH)/vox2roughness $(ARCH)/vox2roughness.m2.o 


$(ARCH)/joint_lsqr_vx_isodamp: $(LSQR_OBJS)
	$(F77) $(USE_FLAGS)  $(LSQR_OBJS) -o $(ARCH)/joint_lsqr_vx_isodamp



$(ARCH)/%.o: %.f
	$(F77) -c $(USE_FLAGS) $< -o $(ARCH)/$*.o

$(ARCH)/%.m1.o: src/Mat/%.f
	$(F77) -c $(USE_FLAGS) $< -o $(ARCH)/$*.m1.o

$(ARCH)/%.m2.o: src/Inv/%.f
	$(F77) -c $(USE_FLAGS) $< -o $(ARCH)/$*.m2.o

$(ARCH)/%.m3.o: src/Vis/%.f
	$(F77) -c $(USE_FLAGS) $< -o $(ARCH)/$*.m3.o

$(ARCH)/%.m2.o: src/Inv/%.f90
	$(F90) -c $(USE_FLAGS) $< -o $(ARCH)/$*.m2.o

$(ARCH)/%.o: %.f90
	$(F90) -c $(USE_FLAGS) $< -o $(ARCH)/$*.o


clean:
	rm $(ARCH)/*
