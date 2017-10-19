include ../makefile.system


OBJS = $(ARCH)/matrix_bw_vx.o $(ARCH)/delazs.o $(ARCH)/param.o \
	$(ARCH)/path_vx.o $(ARCH)/voxels.o $(ARCH)/aniso_cor.o \
	$(ARCH)/summodel.o $(ARCH)/parametrize.o

all: dirs progs


dirs: 
	if [ ! -s $(ARCH)/ ];then \
		mkdir $(ARCH)/; \
	fi;\
	if [ ! -s matrices/ ];then \
		mkdir matrices/; \
	fi

clean:
	rm -f $(ARCH)/*.o $(ARCH)/matrix_bw_vx $(ARCH)/read_mat


progs: $(ARCH)/matrix_bw_vx $(ARCH)/read_mat


$(ARCH)/matrix_bw_vx: $(OBJS)
	$(F77) $(ADDFLAGS)  $(FFLAGS) $(OBJS) -o $(ARCH)/matrix_bw_vx 

$(ARCH)/read_mat: $(ARCH)/read_mat.o
	$(F77) $(ADDFLAGS)  $(FFLAGS) $(ARCH)/read_mat.o -o $(ARCH)/read_mat 


$(ARCH)/%.o: %.f
	$(F77) -c $(FFLAGS) $(ADDFLAGS) $< -o $(ARCH)/$*.o

