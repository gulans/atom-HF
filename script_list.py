import os
import numpy as np
# This script creates input files by combining all_templates/occ/head.template with
# all_templates/occ/Z.template (where Z is atomic umber of the element) and replaces _X_
# with all the input parameters (described in pdf). The input files are being stored in
# "run" folder.

#list of elements to create input files:
list_of_names=['he','be','ne','mg','ar','ca','zn','kr','sr','cd','xe','ba','hg','rn','og']

list_of_names=['ne']

#calcucates <r_nucleus> from the mass of the atom (#visscher1997)
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
    return(rn_list[el])


def name_to_number(el_name):
    num_list={'h':1,'he':2,'li':3,'be':4,'b':5,'c':6,'n':7,'o':8,
            'f':9,'ne':10,'na':11,'mg':12,'al':13,'si':14,'p':15,'s':16,
            'cl':17,'ar':18,'k':19,'ca':20,'sc':21,'ti':22,'v':23,'cr':24,
            'mn':25,'fe':26,'co':27,'ni':28,'cu':29,'zn':30,'ga':31,'ge':32,
            'as':33,'se':34,'br':35,'kr':36,'rb':37,'sr':38,'y':39,'zr':40,
            'nb':41,'mo':42,'tc':43,'ru':44,'rh':45,'pd':46,'ag':47,'cd':48,
            'in':49,'sn':50,'sb':51,'te':52,'i':53,'xe':54,'cs':55,'ba':56,
            'la':57,'ce':58,'pr':59,'nd':60,'pm':61,'sm':62,'eu':63,'gd':64,
            'tb':65,'dy':66,'ho':67,'er':68,'tm':69,'yb':70,'lu':71,'hf':72,
            'ta':73,'w':74,'re':75,'os':76,'ir':77,'pt':78,'au':79,'hg':80,
            'tl':81,'pb':82,'bi':83,'po':84,'at':85,'rn':86,'fr':87,'ra':88,
            'ac':89,'th':90,'pa':91,'u':92,'np':93,'pu':94,'am':95,'cm':96,
            'bk':97,'cf':98,'es':99,'fm':100,'md':101,'no':102,'lr':103,'rf':104,
            'db':105,'sg':106,'bh':107,'hs':108,'mt':109,'og':118}
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
run=False

if generate:

    for el in list:
        os.system("mkdir -p {}".format(nf))
        os.system("cat {}/head.template {}/{}.template > {}/{}.run ".format(tf,tf,el,nf,el))
        #os.system("cp {}/{}.template {}/{}.run".format(tf,el,nf,el))


    for el in list:
        Rmin='1d-8'
        Rmax='30'
        Ngrid='2200'
        sigma=get_rn2(el)
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
        rel="0"
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

