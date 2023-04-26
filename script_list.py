import os
import numpy as np


list_of_names=['he','be','ne','mg','ar','ca','zn','kr','sr','cd','xe','ba','hg','rn','og']

list_of_names=['xe']


def get_rn(el):
    mass={2:4,10:20,18:40,20:40,30:65,36:84,54:131,86:222}
    try:
        m=mass[el]
    except:
        m=0
        print("Nav zināma elementa kodola massa.", el)
    #visscher1997
    rn=(0.836*m**(1/3)+0.570)*1e-15 #in meters
    a0=0.529177249*1e-10 #in meters
    rn=rn/a0 #in a.u.
    ksi=3/(2*rn**2)
    return str(rn)

def get_rn2(el):
    rn_list={2:3.5854812559E-05,
        4:4.3647652710E-05,
        10:5.3782158977E-05,
        12:5.6533300866E-05,
        18:6.4776793917E-05,
        20:6.4835312443E-05,
        30:7.4418094895E-05,
        36:7.9904610996E-05,
        38:8.0939513733E-05,
        48:8.7015861884E-05,
        54:9.1065303450E-05,
        56:9.2277603243E-05,
        80:1.0325010801E-04,
        86:1.0642977075E-04,
        118:1.1581939044E-04}

      # {2:2.07007858208354E-05,
      # 10:3.10511442005306E-05,
      # 18:3.73988997033959E-05,
      # 36:4.61329490496609E-05,
      # 54:5.25765779034729E-05,
      # 4:2.51999840416E-05,
      # 12:3.26395164729E-05,
      # 20:3.74326850917E-05,
      # 30:4.29653071204E-05,
      # 38:4.67304500420E-05,
      # 48:5.02386312822E-05,
      # 56:5.32764990728E-05,
      # 80:5.96114776527E-05,
      # 86:6.14472567899E-05}
    return(rn_list[el])


def name_to_number(el_name):
    num_list={'he':2,
        'be':4,
        'ne':10,
        'mg':12,
        'ar':18,
        'ca':20,
        'zn':30,
        'kr':36,
        'sr':38,
        'cd':48,
        'xe':54,
        'ba':56,
        'hg':80,
        'rn':86,
        'og':118}
    return(num_list[el_name])

def list_to_list(names):
    numbers=[]
    for n in names:
        i=name_to_number(n)
        numbers.append(i)
    return (numbers)




list=list_to_list(list_of_names)

nf='run'    #new folder
tf='all_templates/occ'

generate=True
run=True

if generate:

    for el in list:
        os.system("mkdir -p {}".format(nf))
        os.system("cat {}/head.template {}/{}.template > {}/{}.run ".format(tf,tf,el,nf,el))
        #os.system("cp {}/{}.template {}/{}.run".format(tf,el,nf,el))


    for el in list:
        Rmin='5d-8'
        Rmax='30'
        Ngrid='2200'
        sigma=get_rn2(el)
        #sigma="0d0"
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

        grid="3"
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

##create a copy of input file with a name
if(True):
    for i in range(len(list)):
        os.system("cp {}/{}.run {}/{}.run".format(nf,list[i],nf,list_of_names[i]))

if run:
    for el in list:
        os.system("./atomHF <{}/{}.run".format(nf,el))

