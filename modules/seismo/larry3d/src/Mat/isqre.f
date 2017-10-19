	function isqre(lat,lon,nsqrs,nsqtot,nlatzones,n,eq_incr)
c----finds the index of the square where (lat,lon) is
	implicit real*8 (a-h,o-z)
	real*8 lat,lon,loc_incr
	dimension nsqrs(nlatzones),nsqtot(nlatzones+1)
	lazone=(90.-lat)/eq_incr+1
	if((90.-lat).gt.180.)lazone=nlatzones
	if((90.-lat).gt.181.)stop "problems in function isqre"
	if(lazone.gt.nlatzones)then
	   print*,"problems in function isqre, latitude",lazone,lat
	   stop
	endif
	if(lon.lt.0.)lon=360.+lon
	loc_incr=360./float(nsqrs(lazone))
	isqre=(lon/loc_incr)+1
	isqre=isqre+nsqtot(lazone)
	if(isqre.gt.n)then
	   print*,"problems in function isqre, longitude",n
	   stop
	endif
	return
	end
	
