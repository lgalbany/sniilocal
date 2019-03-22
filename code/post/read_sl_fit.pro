function read_sl_fit,file,l,sp,fit,wei,xl,xm,z,age,mini

  on_error,2

  t1=''
  t2=''
tt3=''
  ;t3=dblarr(4)
  t4=dblarr(9)

  get_lun,lun
  openr,lun,file

  REPEAT BEGIN

    readf,lun,t1
    tt=strsplit(t1,count=cnt,/extract)

    if cnt eq 2 then begin

      if tt[1] eq '[N_base]' then nb=long(tt[0])
      if tt[1] eq '[chi2/Nl_eff]' then chi2=double(tt[0])
     
      if tt[1] eq '[Nl_obs]' then begin
           n=long(tt[0])
	   l=dblarr(n)
	   sp=dblarr(n)
	   fit=dblarr(n)
	   wei=dblarr(n)
	   for i=0,n-1 do begin
	      readf,lun,tt3
              t3=strsplit(tt3,' ',count=ccnntt,/extract) 
              if string(t3[1]) eq '**********' then print,t3
              if string(t3[1]) eq '**********' then print,'changed'
              if string(t3[1]) eq '**********' then t3[1]=0.0
	      l[i]=double(t3[0])
	      sp[i]=double(t3[1])
	      fit[i]=double(t3[2])
	      wei[i]=double(t3[3])
	   endfor
      endif

    endif

    if cnt gt 2 then begin      
      if tt[1] eq '[fobs_norm' then scale=double(tt[0])
      if tt[1] eq '[Mini_tot' then m_initot=double(tt[0])
      if tt[1] eq '[Mcor_tot' then m_cortot=double(tt[0])
      if tt[1] eq '[v0_min' then v0=double(tt[0])
      if tt[1] eq '[vd_min' then vd=double(tt[0])
      if tt[1] eq '[AV_min' then AV=double(tt[0])
      if tt[1] eq '[YAV_min' then YAV=double(tt[0])

      
      if tt[2] eq 'x_j(%)' then begin
	   xl=dblarr(nb)
	   xm=dblarr(nb)
	   z=dblarr(nb)
	   age=dblarr(nb)
	   mini=dblarr(nb)
	   l2m=dblarr(nb)
	   varext=dblarr(nb)
	   mstars=dblarr(nb)
	   
	   for i=0,nb-1 do begin
	      readf,lun,t4
	       xl[i]=t4[1]
	     mini[i]=t4[2]
    	   xm[i]=t4[3]
	      age[i]=t4[4]
	        z[i]=t4[5]
	      l2m[i]=t4[6]
	   varext[i]=t4[7]
	   mstars[i]=t4[8]
	      
	   endfor
      endif
      
      if cnt eq 5 then begin
        if tt[1] eq '[S/N' and tt[3] eq 'S/N' then s2n=double(tt[0]) 
      endif
      
      
    endif
    

  ENDREP UNTIL eof(lun)
  
  sp *= scale
  fit *= scale
  
  close,lun
  free_lun,lun

  pars={lambda:l,spec:sp,fit:fit,weight:wei,xl:xl,xm:xm,z:z,age:age,mini:mini,l2m:l2m,varext:varext, $
        mstars:mstars,AV:AV,YAV:YAV,v0:v0,vd:vd,m_cortot:m_cortot,m_initot:m_initot,nbase:nb,chi2:chi2,s2n:s2n}

return,pars
end


 
