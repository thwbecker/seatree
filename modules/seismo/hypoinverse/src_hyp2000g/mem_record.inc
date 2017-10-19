c .. begin mem_record.inc
c: caxp: record contains one or more misaligned fields.
	integer*4 maxpin
	parameter (maxpin = 1536)
	integer*4 maxset
	parameter (maxset = 64)
	record /hid/ hid		!EVENT (SUMMMARY) DATA
	record /hhyp/ hhy	! caxp:
	record /hfix/ hfx	! caxp:
	record /hmag/ hmg
	record /hrmk/ hrm
	record /hset/ hst(maxset)	!SET DATA
	record /hpin/ hpn(maxpin)	!PIN DATA INCLUDING STA NAME
	record /hpix/ hpx(maxpin*2)
	record /hcod/ hcd(maxpin)
	record /hamf/ haf(maxpin)
	record /hamp/ ham(maxpin)
	common /tuprec/  hid, hhy, hfx, hmg, hrm, hst, hpn, hpx,
	1	hcd, haf, ham
c .. end of mem_record.inc
