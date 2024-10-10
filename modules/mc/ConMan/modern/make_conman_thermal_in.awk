#
# scripts to generate input for ConMan 2020 CIG version
#
# simple thermal convection example
#
BEGIN{
    #
    # default parameters
    #
    #
    if(nel=="")			# number of elements in z
	nel=50;
    if(aspect=="")		# aspect ratio
	aspect=1.0;
    if(rayleigh=="")		# rayleigh number
	rayleigh = 1e6;
    if(heating=="")		# internal heating
	heating=0;
    if(viscosityE=="")		# non-dimensional activation energy
	viscosityE=0.0;
    if(ntimestep=="")
	ntimestep = 1000;		# total timestep
    if(nprintstep=="")
	nprintstep=   100;		# output every nprintstep steps

    
    viscosityP=0.0;		# depth dependece

    if(print_geom=="")		# print
	print_geom=0;		# general input or geometry
    


    tmax = 10.0;			# maxium time (in diffusion time)
    dt_out = 0.01;		# for time series
    tsave=1.0;			# time controlled output
    
    nelz = nel;
    nelx = nel*aspect;
    
    nodex=nelx+1;
    nodez=nelz+1;

    numnp=nodex*nodez;

    iflow=1;			# 1: execute

    necho=0;			# 0: terse 1: verbose

    inrstr=0;			# 0: conductive start 1: restart fro unit 16
    iorstr=0;			# 0: no restart 1: restart

    Tpert=0.1			# temp pert
    Tbot=1.0;
    Ttop=0.0;
    
    # number of edge nodes for nusselt smoother (top and bottom)
    nodebn=nodex*2;

    # rheology
    ntimvs=6;			# 0: const visc 3: Arrheinius 4: diff/disl
                                # 5: Stein & Hansen 6: Frank Kaminetskii

    nwrap=0;			# for periodic BC

    itype=4;			# 1: ALA 2: TALA 3: EBA 4: Boussinesq approx

    isolve=1;			# 1: explicit 2: implicit 3: Picard
    #isolve=2;			# 1: explicit 2: implicit 3: Picard


    numat=1;
    visc0[1]=1;penalty[1]=1e7;diff[1]=1.0;ra[1]=rayleigh;
    dmhu[1]=heating;viscE[1]=viscosityE;viscV[1]=viscosityV;
    ecut[1]=1e6;

    if(!print_geom){
	#
	# print general input file
	# 
	printf("%i by %i element thermal convection problem\n",nelx,nelz);
	print("#Nds   X   Z  ck echo rrst wrst nus tdvf  wr");
	print(numnp,nelx,nelz,iflow,necho,inrstr,iorstr,nodebn,ntimvs,nwrap,itype,isolve);
	print("time step information")
	print(ntimestep,1.0);
	print("output information")
	print(nprintstep,tmax, dt_out, tsave,1.0);
	print("velocity boundary condition flags: IFCMT,DELNXTLN")
	print("bnode   enode   incr bcf1 bcf2")
	# free slip
	print(1,             nelx*nodez+1,nodez,0,1);	# bottom
	print(nelx*nodez+1,(nelx+1)*nodez,1,    1,0);	# right
	print(1,                    nodez,1,    1,0)	# left
	print(nodez,       (nelx+1)*nodez,nodez,0,1);	# top
	# pin the corners
	print(1,1,            1,1,1); # BL
	print(nodez,nodez,    1,1,1); # TL
	print(nelx*nodez+1,nelx*nodez+1, 1,1,1); # BR
	print((nelx+1)*nodez,(nelx+1)*nodez,1,1,1); # TR
	print(0,0,0,0,0);		   # end VBC

	print("temperature boundary condition flags")
	print(1,nelx*nodez+1,nodez,1);	# bottom
	print(nodez,(nelx+1)*nodez,nodez,1);	# top
	print(0,0,0,0,0);		   # end TBC

	print("bndy info (top - bottom rows)");
	print(1,nelx*nodez+1,nodez);	# bottom
	print(nodez,(nelx+1)*nodez,nodez);	# top
	print(0,0,0);

	print("bndy info (2nd from top - 2nd from bottom rows)");
	print(2,nelx*nodez+2,nodez);	# bottom+1
	print(nodez-1,(nelx+1)*nodez-1,nodez);	# top-1
	print(0,0,0);
	
	print("initial condition information");
	print(Tpert,1.0,aspect);

	print("equation of state information")
	print(0.0, 273.0,  1000.0, 1.0, 1.0);

	print("element information")
	print(numat, 0);
	print("viscosity")
	for(i=1;i<=numat;i++)
	    printf("%e\n",visc0[i]);
	print("penalty number")
	for(i=1;i<=numat;i++)
	    printf("%e\n",penalty[i]);
	print("diffusivity (always one)")
	for(i=1;i<=numat;i++)
	    printf("%e\n",diff[i]);
	print("Rayleigh number")
	for(i=1;i<=numat;i++)
	    printf("%e\n",ra[i]);
	print("internal heating parameter")
	for(i=1;i<=numat;i++)
	    printf("%g\n",dmhu[i]);
	print("activation energy")
	for(i=1;i<=numat;i++)
	    printf("%g\n",viscE[i]);

	print("activation volume")
	for(i=1;i<=numat;i++)
	    printf("%g\n",viscV[i]);

	
    }else{
	# print geometry
        # geom
	print("coordinates");
	printf("%i\t%i\t%.1f %.1f\n",1,4,0,0)	# BL
	printf("%i\t%i\t%.1f %.1f\n",nelx*nodez+1,1,  aspect,0)	# BR
	printf("%i\t%i\t%.1f %.1f\n",(nelx+1)*nodez,1,aspect,1)	# TR
	printf("%i\t%i\t%.1f %.1f\n",nodez,1,0,1)	# TL
	# generation
	print(nodex-1,nodez,nodez-1,1) # 
	print(0,0,0,0);
	print("velocity boundary conditions (non-zero)");
	print(0,0,0,0);
	print("temperature boundary conditions (non-zero)");
	printf("%i\t%i\t%.1f\n",1,    2,  Tbot);	# bottom temperature
	printf("%i\t%i\t%.1f\n",nelx*nodez+1,0,Tbot);
	print(nodex-1,nodez);
	print(0,    0, 0.0);	# should be 0,0 to end group?
	print("element connectivity and material groups");
	print(1,1,1, 1,nodez+1,nodez+2,2); # lower left element, sorted CCW
	print(nelx,nelz,nodez,nelz,1,1);
	print(0,      0,      0,      0,      0,      0,      0); # end

    }


}
