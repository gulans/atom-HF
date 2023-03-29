import os
import numpy as np

def get_sigma(el):
    mass={2:4,10:20,18:40,20:40,30:65,36:84,54:131,86:222}
    try:
        m=mass[el]
    except:
        m=0
        print("Nav zināma elementa kodola massa.", el)
    #visscher1997
    rn=(0.836*m**(1/3)+0.570)*1e-15 #in meters
    rn=(0.863*m**(1/3)+0.571)*1e-15 #in meters
    a0=0.529177249*1e-10 #in meters
    rn=rn/a0 #in a.u.
    ksi=3/(2*rn**2)
    sigma=1/((2*ksi)**(1/2))
    return str(sigma)


def get_sigma2(el):
    sigmas={2:2.07007858208354E-05,
       10:3.10511442005306E-05,
       18:3.73988997033959E-05,
       36:4.61329490496609E-05,
       54:5.25765779034729E-05,
       4:2.51999840416E-05,
       12:3.26395164729E-05,
       20:3.74326850917E-05,
       30:4.29653071204E-05,
       38:4.67304500420E-05,
       48:5.02386312822E-05,
       56:5.32764990728E-05,
       80:5.96114776527E-05,
       86:6.14472567899E-05}



    return(sigmas[el])

list=[2,10,18,36,54]

list=[4,12,20,30,38,48,56,80,86]


list=[2,4,10,12,18,20,30,36,38,48,54,56,80,86]
list=[80]
nf='run'    #new folder
tf='all_templates/clasic_config'   #template folder
tf='all_templates/occ'

generate=True
run=True


for rm in np.arange(1e-4,0,-1e-5):
  if generate:


    for el in list:
        os.system("mkdir -p {}".format(nf))
        os.system("cat {}/head.template {}/{}.template > {}/{}.run ".format(tf,tf,el,nf,el))
        #os.system("cp {}/{}.template {}/{}.run".format(tf,el,nf,el))

    for el in list:
        Rmin="{}".format(rm)
        Rmax='30'
        Ngrid='3000'
        sigma=get_sigma2(el)
        sigma="0d0"

        f1="101"
        f1w="1d0"
        f1par='2 \\n0.8040d0\\n0.2195164512209d0'
        f2="130"
        f2w="1d0"

        f3="0"
        f3w="0d0"

        hyboverride=".false."
        HFw="0"
        HFsrw="0"
        HFsrp="0"

        grid="5"
        rel="2"
        os.system("sed -i 's/_Z_/{}/g' {}/{}.run;".format(str(el)+"d0",nf,el))
        os.system("sed -i 's/_Rmin_/{}/g' {}/{}.run;".format(Rmin,nf,el))
        os.system("sed -i 's/_Rmax_/{}/g' {}/{}.run;".format(Rmax,nf,el))
        os.system("sed -i 's/_Ngrid_/{}/g' {}/{}.run;".format(Ngrid,nf,el))
        os.system("sed -i 's/_sigma_/{}/g' {}/{}.run;".format(sigma,nf,el))
        os.system("sed -i 's/_f1_/{}/g' {}/{}.run;".format(f1,nf,el))
        os.system("sed -i 's/_f1w_/{}/g' {}/{}.run;".format(f1w,nf,el))
        os.system("sed -i 's/_f1par_/{}/g' {}/{}.run;".format(f1par,nf,el))
        os.system("sed -i 's/_f2_/{}/g' {}/{}.run;".format(f2,nf,el))
        os.system("sed -i 's/_f2w_/{}/g' {}/{}.run;".format(f2w,nf,el))
        os.system("sed -i 's/_f3_/{}/g' {}/{}.run;".format(f3,nf,el))
        os.system("sed -i 's/_f3w_/{}/g' {}/{}.run;".format(f3w,nf,el))
        os.system("sed -i 's/_hyboverride_/{}/g' {}/{}.run;".format(hyboverride,nf,el))
        os.system("sed -i 's/_HFw_/{}/g' {}/{}.run;".format(HFw,nf,el))
        os.system("sed -i 's/_HFsrw_/{}/g' {}/{}.run;".format(HFsrw,nf,el))
        os.system("sed -i 's/_HFsrp_/{}/g' {}/{}.run;".format(HFsrp,nf,el))
        os.system("sed -i 's/_grid_/{}/g' {}/{}.run;".format(grid,nf,el))
        os.system("sed -i 's/_rel_/{}/g' {}/{}.run;".format(rel,nf,el))

  if run:
    for el in list:
        os.system("./atomHF <{}/{}.run".format(nf,el))

