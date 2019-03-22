    RR=read_sl_fit('locals/'+gal[i]+'_'+snname[i]+'_sp1.C11.gm.CCM.BN',l,sp,spfit,wei,x_j,mcor_j,z_j,age_j)
  
    readcol,'locals/'+gal[i]+'_'+snname[i]+'_sp1.dat',l1,f1,e1,/silent
     f1n=f1
    im=dblarr(1,n_elements(f1))
    eim=dblarr(1,n_elements(f1))
    
    tabinv,l,l1,in1
    f3 = sp - spfit
    f5 = sp / spfit
    f1=bspline_interpol(f3,4,in1)
    f1n=bspline_interpol(f5,4,in1)
   
    im[0,*]=f1
    imn=im
    imn[0,*]=f1n
    eim[0,*]=e1
    
    mass=sl_stellar_mass('locals/'+gal[i]+'_'+snname[i]+'_sp1.C11.gm.CCM.BN',z[i])
    
    OH_0= 4.9E-4 ; Metallicite solaire - Allende, Prieto 2001
    wx=x_j/total(x_j)
    wm=Mcor_j/total(Mcor_j)
    mx=mssm(alog10(age_j),err=wx,/iswei)
    mm=mssm(alog10(age_j),err=wm,/iswei)
    cube_age_l=mx[0]
    cube_age_m=mm[0]
    cube_age_l_e=mx[1]
    cube_age_m_e=mm[1]
    mx=mssm(z_j,err=wx,/iswei)
    mm=mssm(z_j,err=wm,/iswei)
    cube_Z_l=mx[0]
    cube_Z_m=mm[0]
    cube_Z_l_e=mx[1]
    cube_Z_m_e=mm[1]
    cube_av=RR.av
    
    in1=where(age_j lt 0.3e9)
    in2=where(age_j gt 0.3e9 and age_j lt 2.4e9)
    in3=where(age_j ge 2.4e9)
    
    cube_age_l_y=TOTAL((x_j[in1]/total(x_j[in1]))*alog10(age_j[in1]))
    cube_age_m_y=TOTAL((Mcor_j[in1]/total(Mcor_j[in1]))*alog10(age_j[in1]))
    cube_Z_m_y=TOTAL(Mcor_j[in1]*Z_j[in1]/total(Mcor_j[in1]))
    cube_Z_l_y=TOTAL(x_j[in1]*Z_j[in1]/total(x_j[in1]))
    
    cube_age_l_i=TOTAL((x_j[in2]/total(x_j[in2]))*alog10(age_j[in2]))
    cube_age_m_i=TOTAL((Mcor_j[in2]/total(Mcor_j[in2]))*alog10(age_j[in2]))
    cube_Z_m_i=TOTAL(Mcor_j[in2]*Z_j[in2]/total(Mcor_j[in2]))
    cube_Z_l_i=TOTAL(x_j[in2]*Z_j[in2]/total(x_j[in2]))
    
    cube_age_l_o=TOTAL((x_j[in3]/total(x_j[in3]))*alog10(age_j[in3]))
    cube_age_m_o=TOTAL((Mcor_j[in3]/total(Mcor_j[in3]))*alog10(age_j[in3]))
    cube_Z_m_o=TOTAL(Mcor_j[in3]*Z_j[in3]/total(Mcor_j[in3]))
    cube_Z_l_o=TOTAL(x_j[in3]*Z_j[in3]/total(x_j[in3]))