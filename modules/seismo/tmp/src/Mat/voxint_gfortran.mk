FLAGSSUSE= -ffixed-line-length-none -fdollar-ok -O3 
FC=gfortran
SRCDIR=.

#EXECUTABLE FILE
voxint: voxint_$(FC).mk param.o path_vx.o linint.o \
	voxint.o span.o isqre.o delazs.o gceq.o rsoinc.o
	$(FC) $(FLAGSSUSE) -o voxint voxint.o delazs.o gceq.o \
	param.o path_vx.o span.o isqre.o linint.o rsoinc.o
	rm *.o

#OBJECT FILES
voxint.o: voxint.f
	$(FC) $(FLAGSSUSE) -c -o voxint.o voxint.f

param.o: $(SRCDIR)/param.f
	$(FC) $(FLAGSSUSE) -c -o param.o $(SRCDIR)/param.f

path_vx.o: $(SRCDIR)/path_vx.f
	$(FC) $(FLAGSSUSE) -c -o path_vx.o $(SRCDIR)/path_vx.f

delazs.o: $(SRCDIR)/delazs.f
	$(FC) $(FLAGSSUSE) -c -o delazs.o $(SRCDIR)/delazs.f
	
span.o: $(SRCDIR)/span.f
	$(FC) $(FLAGSSUSE) -c -o span.o $(SRCDIR)/span.f
	
isqre.o: $(SRCDIR)/isqre.f
	$(FC) $(FLAGSSUSE) -c -o isqre.o $(SRCDIR)/isqre.f
	
gceq.o: $(SRCDIR)/gceq.f
	$(FC) $(FLAGSSUSE) -c -o gceq.o $(SRCDIR)/gceq.f

linint.o: $(SRCDIR)/linint.f
	$(FC) $(FLAGSSUSE) -c -o linint.o $(SRCDIR)/linint.f

rsoinc.o: $(SRCDIR)/rsoinc.f
	$(FC) $(FLAGSSUSE) -c -o rsoinc.o $(SRCDIR)/rsoinc.f

clean:
	rm *.o
