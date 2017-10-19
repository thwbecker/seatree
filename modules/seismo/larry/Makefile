include makefile.system

all:	dirs progs

dirs: 
	if [ ! -s $(ARCH)/ ];then \
		mkdir $(ARCH)/; \
	fi



progs: $(ARCH)/blk_matrix_ata $(ARCH)/blk_cholesky \
	$(ARCH)/blk_lsqr $(ARCH)/blk2gmt $(ARCH)/blk2gmtz

$(ARCH)/%: $(ARCH)/%.o
	$(F77)  $(ARCH)/$*.o -o $(ARCH)/$*


$(ARCH)/%.o: %.f
	$(F77) -c $(FFLAGS) $(F_EXT_SOURCE_FLAG) $< -o $(ARCH)/$*.o

clean:
	rm $(ARCH)/*.o $(ARCH)/*.a $(ARCH)/blk_matrix_ata $(ARCH)/blk_cholesky \
		$(ARCH)/blk_lsqr $(ARCH)/blk2gmt	
