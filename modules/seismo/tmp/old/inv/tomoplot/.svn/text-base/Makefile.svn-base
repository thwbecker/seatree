include ../../makefile.system


all: ../$(ARCH)/mapview_3d ../$(ARCH)/sol2xyz


clean: 
	rm ../$(ARCH)/mapview_3d


../$(ARCH)/mapview_3d: ../$(ARCH)/mapview_3d.o ../$(ARCH)/param.o
	$(F77) $(FFLAGS) $(ADDFLAGS)  -o ../$(ARCH)/mapview_3d \
		../$(ARCH)/mapview_3d.o ../$(ARCH)/param.o

../$(ARCH)/mapview_3d.o: mapview_3d.f
	$(F77) $(FFLAGS) $(ADDFLAGS)  -c -o ../$(ARCH)/mapview_3d.o mapview_3d.f

../$(ARCH)/sol2xyz: ../$(ARCH)/sol2xyz.o ../$(ARCH)/param.o
	$(F77) $(FFLAGS) $(ADDFLAGS)  -o ../$(ARCH)/sol2xyz \
		../$(ARCH)/sol2xyz.o ../$(ARCH)/param.o

../$(ARCH)/sol2xyz.o: sol2xyz.f
	$(F77) $(FFLAGS) $(ADDFLAGS)  -c -o ../$(ARCH)/sol2xyz.o sol2xyz.f

../$(ARCH)/param.o: ../../mat/param.f
	$(F77) $(FFLAGS) $(ADDFLAGS)  -c -o ../$(ARCH)/param.o ../../mat/param.f
