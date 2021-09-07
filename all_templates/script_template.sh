
for i in {1..2};
do
sed -i 's/_Rmin_/1d-5/g' $i.template;
sed -i 's/_Rmax_/18/g' $i.template;
sed -i 's/_Ngrid_/400/g' $i.template;
done

for i in {3..10};
do
sed -i 's/_Rmin_/1d-5/g' $i.template;
sed -i 's/_Rmax_/24/g' $i.template;
sed -i 's/_Ngrid_/500/g' $i.template;
done

for i in {11..18};
do
sed -i 's/_Rmin_/1d-5/g' $i.template;
sed -i 's/_Rmax_/24/g' $i.template;
sed -i 's/_Ngrid_/800/g' $i.template;
done

for i in {19..36};
do
sed -i 's/_Rmin_/1d-6/g' $i.template;
sed -i 's/_Rmax_/24/g' $i.template;
sed -i 's/_Ngrid_/1000/g' $i.template;
done

for i in {37..54};
do
sed -i 's/_Rmin_/1d-6/g' $i.template;
sed -i 's/_Rmax_/24/g' $i.template;
sed -i 's/_Ngrid_/1200/g' $i.template;
done

for i in {55..86};
do
sed -i 's/_Rmin_/1d-6/g' $i.template;
sed -i 's/_Rmax_/24/g' $i.template;
sed -i 's/_Ngrid_/1200/g' $i.template;
done

for i in {87..92};
do
sed -i 's/_Rmin_/1d-6/g' $i.template;
sed -i 's/_Rmax_/24/g' $i.template;
sed -i 's/_Ngrid_/1200/g' $i.template;
done


for i in {1..92};
do
sed -i 's/_f1_/402/g' $i.template;
sed -i 's/_f1w_/1d0/g' $i.template;

sed -i 's/_f2_/0/g' $i.template;
sed -i 's/_f2w_/0d0/g' $i.template;

sed -i 's/_f3_/0/g' $i.template;
sed -i 's/_f3w_/0d0/g' $i.template;

sed -i 's/_hyboverride_/.false./g' $i.template;
sed -i 's/_HFw_/0d0/g' $i.template;
sed -i 's/_HFsrw_/0d0/g' $i.template;
sed -i 's/_HFsrp_/0d0/g' $i.template;

sed -i 's/_grid_/3/g' $i.template;

done

