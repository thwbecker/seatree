include makefile.include

#
# all modules are controlled by the global parmater file
# global_para.dat
#
all: model data inv

tar: clean_all tarfile

clean: clean_all

#
# just build the binaries
#
build_only:
	cd makemodel;\
	make tools; \
	cd ../makedata; \
	make tools; \
	cd ../inversion; \
	make tools; \
	cd ..

#
# generate the checkerboard input model
#
model:
	cd makemodel;\
	make ; \
	cd ..
#
# generate data
#
data:
	cd makedata;\
	make ; \
	cd ..
#
# generate inversion
#
inv:
	cd inversion;\
	make ; \
	cd ..


clean_all:
	rm *~;\
	cd makemodel; make clean_all; cd ..;\
	cd makedata; make clean_all; cd ..;\
	cd inversion; make clean_all; cd ..

clean:
	rm *~;\
	cd makemodel; make clean; cd ..;\
	cd makedata; make clean; cd ..;\
	cd inversion; make clean; cd ..


tarfile:
	rm syn2d.tgz;\
	cd ..;	\
	tar --create --verbose --gzip --exclude '.svn' --exclude 'old\/' \
		--exclude '.gmt' --exclude 'RCS' --file syn2d.tgz syn2d/ ; \
	cd syn2d;\
	mv ../syn2d.tgz . ;\

