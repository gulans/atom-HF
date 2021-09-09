i=36
sed 's/_Rmin_/1d-5/g' <$i.template > $i.template.tmp;
sed -i 's/_Rmax_/18/g' $i.template.tmp;

sed -i 's/_f1_/1/g' $i.template.tmp;
sed -i 's/_f1w_/1d0/g' $i.template.tmp;

sed -i 's/_f2_/12/g' $i.template.tmp;
sed -i 's/_f2w_/1d0/g' $i.template.tmp;

sed -i 's/_f3_/0/g' $i.template.tmp;
sed -i 's/_f3w_/0d0/g' $i.template.tmp;

sed -i 's/_hyboverride_/.false./g' $i.template.tmp;
sed -i 's/_HFw_/0d0/g' $i.template.tmp;
sed -i 's/_HFsrw_/0d0/g' $i.template.tmp;
sed -i 's/_HFsrp_/0d0/g' $i.template.tmp;

sed -i 's/_grid_/10/g' $i.template.tmp ;



for n in {0..90};
do
Ngrid=$((100+n*10))	
sed "s/_Ngrid_/$Ngrid/g" < $i.template.tmp > $i-$n.input ;

done


