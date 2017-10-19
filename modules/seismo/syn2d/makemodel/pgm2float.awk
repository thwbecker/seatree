#
# conver ASCII PGM format to a list of floats
#
# if cstring=1, will print GMT grd conversion string
BEGIN{
    n=0;
    if(cstring == "")
	cstring =0;
}
{
  if(NR==1){
	type=$1;
  } else if(NR==3){
    nx=$1;ny=$2;
    if(cstring == 1)
      printf("xyz2grd -R1/%i/1/%i -ZTL -I1/1\n",nx,ny);
    else if(cstring == 2)
      print(nx,ny);
  }else if(NR==4){
    bscale=$1;
    bscaleh = bscale /2;
  }else if((!cstring)&&(NR>4)){
    n++;
    printf("%11g\n",($1-bscaleh)/bscaleh);
  }
}
