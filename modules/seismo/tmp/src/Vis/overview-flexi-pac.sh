#!/bin/bash
layxcol=5
layinc=2
# arguments: 1 model, 2 cpt file, 3 horiz size of voxels, 4 remove average y/n
model=$1
color=$2
size=$3
ave=$4
rtop=$5
rbot=$6
nlay=$7
awk '{print $1, $2*100}' $model > tmp.txt
dh=`echo "($rtop - $rbot)/$nlay"|bc`
mode=1 #i.e. type of plot (see fortran code)
#binfold=/home/boschil/SUWA/THREEDEE/VOXELS/dy_inv_twb_version/bin/mapview_3d
#binfold=/home1/boschil/BOWA/TOMO/inversion
binfold=.
$binfold/mapview_3d << !
$size
0
$nlay
tmp.txt
"$color"
$ave
$rtop, $rbot
$mode
!
l1=1
l2=`echo "$l1 + ($layxcol * $layinc)"|bc`
l3=`echo "$l2 + ($layxcol * $layinc)"|bc`
ydw=3.1
yup=`echo "$ydw * ($layxcol - 1)"|bc`
gmtset BASEMAP_TYPE plain LABEL_FONT_SIZE 14 LABEL_FONT 6 ANNOT_FONT_SIZE 12 ANNOT_FONT 6
layern=`echo $l1 |awk '{printf ("%02d",$1)}'`
echo plotting layer $layern
psxy layer_$layern -P -A1000 -X1 -Y26.5 -V -JH180/6 -Bwesn -R0/360/-90/90 -K -m > $model.ps 
pscoast -V -R -O -K -JH -W1 -Dc >> $model.ps
psxy plate_boundaries.codes -W3/150/150/150 -m -O -K -V -R -JH >> $model.ps
depth=`echo "($layern * $dh) - ($dh / 2)"|bc`
pstext -W255/255/255Owhite -JQ0/6 -R -O -K -V -N << ! >> $model.ps
5 -85 12 0 6 BL $depth km
!
layern=`expr $l1 + $layinc`
laymax=`echo "($l1 + ($layxcol * $layinc)) - 1"|bc`
while [ "$layern" -le "$laymax" ];do
layern=`echo $layern |awk '{printf ("%02d",$1)}'`
echo plotting layer $layern
psxy layer_$layern -P -A -Y-$ydw -V -JH -Bwesn -R -O -K -m >> $model.ps 
pscoast -V -R -O -K -JH -W1 -Dc >> $model.ps
psxy plate_boundaries.codes -W3/150/150/150 -m -O -K -V -R -JH >> $model.ps
depth=`echo "($layern * $dh) - ($dh / 2)"|bc`
pstext -W255/255/255Owhite -JQ0/6 -R -O -K -V -N << ! >> $model.ps
5 -85 12 0 6 BL $depth km
!
layern=`expr $layern + $layinc`
done
#
layern=$l2
layern=`echo $layern |awk '{printf ("%02d",$1)}'`
echo plotting layer $layern
psxy layer_$layern -P -A -X6.1 -Y$yup -V -JH -Bwesn -R -O -K -m >> $model.ps 
pscoast -V -R -O -K -JH -W1 -Dc >> $model.ps
psxy plate_boundaries.codes -W3/150/150/150 -m -O -K -V -R -JH >> $model.ps
depth=`echo "($layern * $dh) - ($dh / 2)"|bc`
pstext -W255/255/255Owhite -JQ0/6 -R -O -K -V -N << ! >> $model.ps
5 -85 12 0 6 BL $depth km
!
#for layern in 05 06; do
layern=`expr $l2 + $layinc`
laymax=`echo "($l2 + ($layxcol * $layinc)) - 1"|bc`
while [ "$layern" -le "$laymax" ];do
layern=`echo $layern |awk '{printf ("%02d",$1)}'`
echo plotting layer $layern
psxy layer_$layern -P -A -Y-$ydw -V -JH -Bwesn -R -O -K -m >> $model.ps 
pscoast -V -R -O -K -JH -W1 -Dc >> $model.ps
psxy plate_boundaries.codes -W3/150/150/150 -m -O -K -V -R -JH >> $model.ps
depth=`echo "($layern * $dh) - ($dh / 2)"|bc`
pstext -W255/255/255Owhite -JQ0/6 -R -O -K -V -N << ! >> $model.ps
5 -85 12 0 6 BL $depth km
!
layern=`expr $layern + $layinc`
done
layern=$l3
layern=`echo $layern |awk '{printf ("%02d",$1)}'`
echo plotting layer $layern
psxy layer_$layern -P -A -X6.1 -Y$yup -V -JH -Bwesn -R -O -K -m >> $model.ps 
pscoast -V -R -O -K -JH -W1 -Dc >> $model.ps
psxy plate_boundaries.codes -W3/150/150/150 -m -O -K -V -R -JH >> $model.ps
depth=`echo "($layern * $dh) - ($dh / 2)"|bc`
pstext -W255/255/255Owhite -JQ0/6 -R -O -K -V -N << ! >> $model.ps
5 -85 12 0 6 BL $depth km
!
#for layern in 08; do
layern=`expr $l3 + $layinc`
laymax=`echo "($l3 + ($layxcol * $layinc)) - 1"|bc`
while [ "$layern" -le "$laymax" ];do
layern=`echo $layern |awk '{printf ("%02d",$1)}'`
echo plotting layer $layern
psxy layer_$layern -P -A -Y-$ydw -V -JH -Bwesn -R -O -K -m >> $model.ps 
pscoast -V -R -O -K -JH -W1 -Dc >> $model.ps
psxy plate_boundaries.codes -W3/150/150/150 -m -O -K -V -R -JH >> $model.ps
depth=`echo "($layern * $dh) - ($dh / 2)"|bc`
pstext -W255/255/255Owhite -JQ0/6 -R -O -K -V -N << ! >> $model.ps
5 -85 12 0 6 BL $depth km
!
layern=`expr $layern + $layinc`
done
psscale -C$color -O -K -V -D-3.05/-1/10/.5h -S -B0.5/:"@~d@~v/v (%)":>>$model.ps
echo "now plot profile of mean heterogeneity vs. depth"
awk '{if($1 != ">")print $1,$2}' vprofile.txt |\
psxy -W2 -X-10.5 -Y-6 -R0/2900/-1/1 -JX16/4 -O -Ba500:"depth (km)":/a0.2:"mean @~d@~v/v (%)":WS >> $model.ps
for file in `ls layer*`; do
rm $file
done
./twb_modifybb $model.ps
convert $model.ps $model.jpg
