c .. Begin mem_structure.inc: .mem file database structures
	structure /datetime/
	  integer*4		yr
	  integer*4		mo
	  integer*4		dy
	  integer*4		hr
	  integer*4		mn
	  real*4		sec			! SECONDS
	end structure

c .. hid: EVENT ID SUMMARY DATA STRUCTURE
	structure /hid/
	  character*8		typ			! EVENT TYPE
	  integer*4		id			! MEM ID
	  integer*4		gid			! GRM ID
	  integer*4		aid			! ARKIVE ID
	  integer*4		atap			! ARKIVE TAPE
	  integer*4		afil			! ARKIVE FILE
	  integer*4		nst			! COUNT OF HST STRUCTURES
	  integer*4		npn			! COUNT OF HPN STRUCTURES
	  integer*4		npx			! COUNT OF HPX STRUCTURES
	  integer*4		ncd			! COUNT OF HCD STRUCTURES
	  integer*4		naf			! COUNT OF HAF STRUCTURES
	  integer*4		nam			! COUNT OF HAM STRUCTURES
	end structure

c .. hfix : EVENT HYPOCENTER FIXED STRUCTURE
	structure /hfix/
	  real*8		t0			! ORIGIN TIME CUSP SECS
	  record /datetime/	t			! ORIGIN TIME
	  integer*4		msk			! EVENT TYPE; FIX
	  real*4		lat			! LATITUDE
	  real*4		lon			! LONGITUDE
	  real*4		z			! DEPTH (- DOWN)
	end structure

c .. hhyp : EVENT HYPOCENTER SUMMARY STRUCTURE
	structure /hhyp/
	  real*8		t0			! ORIGIN TIME CUSP SECS
	  record /datetime/	t			! ORIGIN TIME
	  integer*4		msk			! EVENT TYPE; FIX
	  real*4		lat			! LATITUDE
	  real*4		lon			! LONGITUDE
	  real*4		z			! DEPTH (- DOWN)
	  real*4		rms			! RMS OF SOLUTION
	  integer*4		np			! NUMBER OF PHASES IN SOLUTION
	  real*4		gap			! AZIMUTHAL GAP
	  real*4		dmn			! DISTANCE TO CLOSEST STN
	  real*4		elt			! ERROR LAT
	  real*4		eln			! ERROR LON
	  real*4		ez			! ERROR DEPTH
	  real*4		et			! ERROR TIME
	end structure

c .. hmag: EVENT MAGNITUDE SUMMARY STRUCTURE
	structure /mag/
	  real*4		m			! MAGNITUDE
	  integer*4		np			! NUMBER OF PHASES
	  real*4		rms			! RMS OF CALCULATION
	end structure

	structure /hmag/
	  record /mag/		md			! CODA DURATION MAGNITUDE
	  record /mag/		mc			! CODA AMPLITUDE MAGNITUDE
	  record /mag/		ml			! WOOD-ANDERSON MAGNITUDE
	  record /mag/		mh			! HELICORDER MAGNITUDE
	  record /mag/		ms			! SURFACEWAVE MAGNITUDE
	  record /mag/		mb			! BODYWAVE TELESEISMIC
	  record /mag/		mw			! ENERGY MAGNITUDE
	  record /mag/		mz			! CODA DURATION LOGAIN
	end structure

c .. hrmk: REMARKS DESCRIPTIVE OF EVENT
	structure /hrmk/
	  character*120		qry
	  character*120		rmk
	  character*120		flt
	end structure

c .. hset: EVENT SET SPECIFIC DATA STRUCTURE
	structure /hset/
	  integer*4		is			! HST STRUCTURE #
	  integer*4		set			! SET NUMBER
	  real*8		t0			! SET START TIME CUSP SECS
	  record /datetime/	t		! SET START DATE-TIME (REFERENCE TIME)
	  character*4 		net			! NETWORK NAME
	  character*4 		dev			! DEVICE NAME
	  real*4		dt			! SECS/SAMPLE
	  integer*4		inc			! BYTES/SAMPLE
	  integer*4		b			! DIGITIZER BITS
	  integer*4		mc			! MAX COUNTS
	  real*4		vm			! MAX VOLTS
	  integer*4		syn			! TIME CODE SYNCH IDENTIFIER
	end structure

c .. hpin: EVENT PIN SPECIFIC DATA STRUCTURE
	structure /hpin/
	  integer*4		is			! HST STRUCTURE #
	  integer*4		ip			! HPN STRUCTURE #
	  integer*4		set			! SET NUMBER
	  integer*4		pin			! PIN NUMBER
	  character*4		net			! NET NAME 3 BYTES
	  character*8		nam			! SITE NAME 6 BYTES
	  character*4		typ			! COMPONENT 3 BYTES
	  integer*4		rtc			! DIGITIZER TIME CNTS TO 1ST SAMPLE
	  integer*4		key		! OFFSET TO START OF SEISMOGRAM IN BYTES
	  integer*4		n			! LENGTH OF TRACE (BYTES)
	  integer*4		bias			! DC OFFSET OF TRACE
	end structure

c .. hpix: pIN PHASE ARRIVAL TIME DATA STRUCTURE
	structure /hpix/
	  integer*4		ip			! HPN STRUCTURE #
	  integer*4		set			! SET NUMBER
	  integer*4		pin			! PIN NUMBER
	  integer*4		rtc			! SAMPLE COUNT OF ARRIVAL
	  real*8		t0			! ARRIVAL TIME CUSP SECS
	  record/datetime/	t			! TIME OF ARRIVAL
	  character*8		phz			! PHASE DESCRIPTOR
	  integer*4		az			! AZIMUTH TO STATION ( 0 = n)
	  integer*4		ia			! TAKE-OFF ANGLE ( 0 = up)
	  real*4		x			! DISTANCE TO STN  (km)
	  real*4		ttr			! TRAVELTIME RESIDUAL (s)
	  real*4		tc			! DELAY TIME CORRECTION (s)
	  real*4		err			! TIMING ERROR ESTIMATE (s)
	end structure

c .. hcod: pIN CODA DURATION DATA STRUCTURE
	structure /hcod/
	  integer*4		ip			! HPN STRUCTURE #
	  integer*4		set			! SET NUMBER
	  integer*4		pin			! PIN NUMBER
	  integer*4		rtc			! SAMPLE COUNT AT START OF CODA
	  integer*4		neq			! NUMBER OF CODA WINDOWS
	  real*4		amp			! AMPLITUDE OF s IN DIGITAL COUNTS
	  real*4		afx			! NOMINAL MINIMUN AMPLITUDE
	  real*4		qfx			! FIXED CODA DECAY CONSTANT
	  real*4		afr			! FREE AMPLITUDE
	  real*4		qfr			! FREE CODA DECAY CONSTANT
	  real*4		rms			! RESIDUAL OF FIT
	  real*4		tau			! LENGTH OF CODA
	  real*4		rbb			! AMPLITUDE OF FINAL SAMPLE
	  character*8		phz			! PHASE DESCRIPTOR
	end structure

c .. hamf: pIN AMPLITUDE/PERIOD DATA STRUCTURE (MENLO, CIT)
	structure /hamf/
	  integer*4		ip			! HPN STRUCTURE #
	  integer*4		set			! SET NUMBER
	  integer*4		pin			! PIN NUMBER
	  integer*4		rtc1			! SAMPLE COUNT AT 1ST Z.C.
	  real*4		wt			! analysis weight
	  real*4		amp1			! PEAK-TO-PEAK AMP
	  real*4		per1			! PERIOD (S)
	  character*8		phz1			! PHASE DESCRIPTOR
	end structure

c .. hamp: WOOD-ANDERSON AMPLITUDE/PERIOD DATA STRUCTURE (CALTECH)
	structure /hamp/
	  integer*4		ip			! HPN STRUCTURE #
	  integer*4		set			! SET NUMBER
	  integer*4		pin			! PIN NUMBER
	  real*4		amp1			! PEAK-TO-PEAK AMP
	  real*4		per1			! PERIOD (S)
	  character*8		phz1			! PHASE DESCRIPTOR
	end structure

c .. End of mem_structure.inc

