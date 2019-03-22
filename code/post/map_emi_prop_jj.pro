PRO map_emi_prop_jj,name,path,suff_read,redshift,sn_tresh,scale,sn_o3_rat,invor,suff1_write,make_ps=make_ps,verify=verify,voronoi=voronoi,_extra=e

astrolib

suffix='.fits'

ids=['[N II] 6548',$
'Halpha 6562.8',$
'[N II] 6583',$
'[O II] 3727',$
'[Ne III] 3868',$  
'Hbeta 4861',$
'Ne I 4931',$     
'[O III] 4959',$   
'[O III] 5007',$   
'[O I] 5577',$
'He I 5875',$
'[O I] 6300',$
'[S III] 6312',$   
'[O I] 6364',$     
'Hdelt 4101',$
'Heps 3970',$   
'Hgamma 4340',$    
'[O III] 4363',$      
'[S II] 6716',$    
'[S II] 6731']    

lam=[6548.05,$
 6562.80,$
 6583.45,$
 3727.43,$
 3868.76,$ 
 4861.33 ,$
 4930.944,$
 4958.911,$
 5006.843,$
 5577.34 ,$
 5875.62, $
 6300.304, $
 6312.06 ,$
 6363.776,$
 4101.7,$
 3970.1,$
 4340.47 ,$
 4363.209,$
 6716.440,$
 6730.815]

 lam1=[3727.43,4958.911,5006.843,4340.47,4861.33,6562.80,6548.05,6583.45,6300.304,6716.440,6730.815,5875., 4101.7]
 label=['OII' ,'OIIIa' ,'OIIIb' ,'Hg'   ,'Hb'   ,'Ha'   ,'NIIa' ,'NIIb' ,'OI'    ,'SIIa'  ,'SIIb','HeI','Hd']

 one=dblarr(13)
 one[*]=1.0

 load_map,path+name+suff_read+'_OII'+suffix,sn_tresh,f_OII,e_OII,_extra=e
 load_map,path+name+suff_read+'_OIIIb'+suffix,sn_tresh,f_OIIIb,e_OIIIb,_extra=e
 load_map,path+name+suff_read+'_OIIIa'+suffix,sn_tresh,f_OIIIa,e_OIIIa,_extra=e
 load_map,path+name+suff_read+'_Hg'+suffix,sn_tresh,f_Hg,e_Hg,_extra=e
 load_map,path+name+suff_read+'_Hb'+suffix,sn_tresh,f_Hb,e_Hb,_extra=e
 load_map,path+name+suff_read+'_Ha'+suffix,sn_tresh,f_Ha,e_Ha,_extra=e
 load_map,path+name+suff_read+'_NIIa'+suffix,sn_tresh,f_NIIa,e_NIIa,_extra=e
 load_map,path+name+suff_read+'_NIIb'+suffix,sn_tresh,f_NIIb,e_NIIb,_extra=e
 load_map,path+name+suff_read+'_OI'+suffix,sn_tresh,f_OI,e_OI,_extra=e
 load_map,path+name+suff_read+'_SIIa'+suffix,sn_tresh,f_SIIa,e_SIIa,_extra=e
 load_map,path+name+suff_read+'_SIIb'+suffix,sn_tresh,f_SIIb,e_SIIb,_extra=e
 
 load_map,path+name+suff_read+'_NeIII'+suffix,sn_tresh,f_NeIII,e_NeIII,_extra=e
 load_map,path+name+suff_read+'_OIII4363'+suffix,sn_tresh,f_OIII4363,e_OIII4363,_extra=e
 load_map,path+name+suff_read+'_HeI'+suffix,sn_tresh,f_HeI,e_HeI,_extra=e
 load_map,path+name+suff_read+'_NeI'+suffix,sn_tresh,f_NeI,e_NeI,_extra=e
 load_map,path+name+suff_read+'_OI5577'+suffix,sn_tresh,f_OI5577,e_OI5577,_extra=e
 load_map,path+name+suff_read+'_SIII'+suffix,sn_tresh,f_SIII,e_SIII,_extra=e
 load_map,path+name+suff_read+'_OI6364'+suffix,sn_tresh,f_OI6364,e_OI6364,_extra=e
 load_map,path+name+suff_read+'_Hd'+suffix,sn_tresh,f_Hd,e_Hd,_extra=e
 load_map,path+name+suff_read+'_He'+suffix,sn_tresh,f_He,e_He,_extra=e

;;;;;; END VORONOI
s=size(f_OII,/dim)

cube=dblarr(13)
e_cube=dblarr(13)

   cube[0]=f_OII
 e_cube[0]=e_OII
   cube[1]=f_OIIIa
 e_cube[1]=e_OIIIa
   cube[2]=f_OIIIb
 e_cube[2]=e_OIIIb
   cube[3]=f_hg
 e_cube[3]=e_hg
   cube[4]=f_hb
 e_cube[4]=e_hb
   cube[5]=f_ha
 e_cube[5]=e_ha
   cube[6]=f_NIIa
 e_cube[6]=e_NIIa
   cube[7]=f_NIIb
 e_cube[7]=e_NIIb
   cube[8]=f_OI
 e_cube[8]=e_OI
   cube[9]=f_SIIa
 e_cube[9]=e_SIIa
  cube[10]=f_SIIb
e_cube[10]=e_SIIb
  cube[11]=f_HeI
e_cube[11]=e_HeI
  cube[12]=f_Hd
e_cube[12]=e_Hd

cube_results=fltarr(2)
ebvim=fltarr(1)
e_ebvim1=fltarr(1)
e_ebvim2=fltarr(1)


ohim=fltarr(1)
e_ohim=fltarr(1)
qim=ohim
qim[*]=!values.f_nan
errq=qim

t_r32=ohim

;Extinction

iha=5
ihb=4

ebvim[*]=!values.f_nan
e_ebvim1[*]=!values.f_nan
e_ebvim2[*]=!values.f_nan
mask=fltarr(1)
mask[*]=!values.f_nan

m=400

ebv=0.01*(dindgen(m+1)-100)
rat=dblarr(m+1)

for i=0,400 do begin 
    fm_unred,[4861.33,6562.8],[1.,1.],-ebv[i],f 
    rat[i]=f[0]/f[1]/2.86
endfor

;forprint,ebv,rat

imrat=reform((cube[ihb]/cube[iha]))
dimc=size(cube,/dim)
imrat=dblarr(1)
imrat[0]=cube[ihb]/cube[iha]

errHb_Ha= reform(sqrt( (e_cube[ihb]/cube[iha])^2+ (cube[ihb]*e_cube[iha]/cube[iha]^2)^2))
errHb_Ha=dblarr(1)
errHb_Ha[0]= reform(sqrt( (e_cube[ihb]/cube[iha])^2+ (cube[ihb]*e_cube[iha]/cube[iha]^2)^2))
in=where(reform((cube[ihb]/cube[iha])/errHb_Ha) ge 2.)
mask[in]=1

writefits,path+name+suff1_write+'_HaHb_mask.fits',mask

if keyword_set(verify) then print,minmax(errHb_Ha,/nan)

cube_results[0]=imrat
cube_results[1]=errHb_Ha
writefits,path+name+suff1_write+'_HaHb.fits',cube_results;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_HaHb.fits',max=0.5,min=0

if keyword_set(verify) then print,'Hb/Ha at center and error: ',minmax(imrat[15:25,15:25])/2.86,minmax(errHb_Ha[15:25,15:25]);,minmax(rat)

in=sort(rat)
rat=rat[in]
ebv=ebv[in]

    if finite(imrat) then begin
      ebvim=spline(rat,ebv,imrat,/double)     
      e_ebvim1=spline(rat,ebv,imrat-errHb_Ha,/double)
      e_ebvim2=spline(rat,ebv,imrat+errHb_Ha,/double)
                 
      fm_unred,lam1,one,ebvim,f 
      fm_unred,lam1,one,e_ebvim1,f1
      fm_unred,lam1,one,e_ebvim2,f2
      e_f=0.5*abs(f1-f2)

      ;;Extinction correct
      for k=0,12 do begin
          e_cube[k] = sqrt(((f[k]>1.0)*e_cube[k])^2+(e_f[k]*cube[k])^2)
          cube[k] *= (f[k]>1.0)
      endfor
     endif

f_OII	 =    reform(cube[0])
e_OII	 =  reform(e_cube[0])
f_OIIIa  =    reform(cube[1])
e_OIIIa  =  reform(e_cube[1])
f_OIIIb  =   reform( cube[2])
e_OIIIb  =  reform(e_cube[2])
f_hg	 =   reform( cube[3])
e_hg	 =  reform(e_cube[3])
f_hb	 =    reform(cube[4])
e_hb	 =  reform(e_cube[4])
f_ha	 =    reform(cube[5])
e_ha	 =  reform(e_cube[5])
f_NIIa   =   reform( cube[6])
e_NIIa   =  reform(e_cube[6])
f_NIIb   =    reform(cube[7])
e_NIIb   =  reform(e_cube[7])
f_OI	 =    reform(cube[8])
e_OI	 =  reform(e_cube[8])
f_SIIa   =    reform(cube[9])
e_SIIa   =  reform(e_cube[9])
f_SIIb   =   reform(cube[10])
e_SIIb   = reform(e_cube[10])
f_HeI    =   reform(cube[11])
e_HeI    = reform(e_cube[11])
f_Hd     =   reform(cube[12])
e_Hd     = reform(e_cube[12])

mkhdr,hjj,cube,/image
sxaddpar,hjj,'CONTENT','EXT-CORR LINE FLUX: GAUSSIAN FIT','erg/s/A/cm^2'
sxaddpar,hjj,'POS01','OII','erg/s/A/cm^2'
sxaddpar,hjj,'POS02','OIII4959','erg/s/A/cm^2'
sxaddpar,hjj,'POS03','OIII5007','erg/s/A/cm^2'
sxaddpar,hjj,'POS04','Hg','erg/s/A/cm^2'
sxaddpar,hjj,'POS05','Hb','erg/s/A/cm^2'
sxaddpar,hjj,'POS06','Ha','erg/s/A/cm^2'
sxaddpar,hjj,'POS07','NII6548','erg/s/A/cm^2'
sxaddpar,hjj,'POS08','NII6583','erg/s/A/cm^2'
sxaddpar,hjj,'POS09','OI6300','erg/s/A/cm^2'
sxaddpar,hjj,'POS10','SII6716','erg/s/A/cm^2'
sxaddpar,hjj,'POS11','SII6731','erg/s/A/cm^2'
sxaddpar,hjj,'POS12','HeI5875','erg/s/A/cm^2'
sxaddpar,hjj,'POS13','Hd','erg/s/A/cm^2'
writefits,path+name+suff1_write+'_JJ.fits',cube,hjj
writefits,path+name+suff1_write+'_JJ.fits',e_cube,hjj,/app

if keyword_set(verify) then print,'Hb/Ha at center and error: ',minmax(ebvim[15:25,15:25]),minmax(e_ebvim[15:25,15:25])

av1=3.1*ebvim
errav1=abs(3.1*0.5*(e_ebvim1-e_ebvim2))

if keyword_set(verify) then print,minmax(errav1,/nan)

in=where(av1 lt 0, cn)

if cn gt 0 then begin 
   Av1[in] = 0.
endif   

cube_results[0]=Av1
cube_results[1]=errAv1
writefits,path+name+suff1_write+'_AV.fits',cube_results
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_AV.fits'

if not keyword_set(metalonly) then begin
; Electronic density
N_e=f_SIIa/f_SIIb
err_Ne= sqrt((e_SIIa/f_SIIb)^2+(f_SIIb*e_SIIa/(f_SIIb^2))^2)

cube_results[0]=N_e
cube_results[1]=err_Ne
writefits,path+name+suff1_write+'_SaSb.fits',cube_results
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_SaSb.fits'
endif
;Diagram diagnostic

OIII_Hb=alog10(f_OIIIb/f_Hb)
NII_Ha=alog10(f_NIIb/f_Ha)
e_OIII_Hb=sqrt((e_OIIIb/f_Hb)^2+(e_Hb*f_OIIIb/(f_Hb^2))^2)
e_NII_Ha=sqrt((e_NIIb/f_Ha)^2+(e_Ha*f_NIIb/(f_Ha^2))^2)

err_NII_Ha=abs(e_NII_Ha/(10.0^NII_Ha*alog(10.)))
err_OIII_Hb=abs(e_OIII_Hb/(10.0^OIII_Hb*alog(10.)))

if keyword_set(verify) then print,'OIII/Hb at center and error: ',minmax((f_OIIIb/f_Hb)[15:25,15:25]),minmax(e_OIII_Hb[15:25,15:25]);,minmax(rat)
if keyword_set(verify) then print,'NII/Hb at center and error: ',minmax((f_NIIb/f_Ha)[15:25,15:25]),minmax(e_NII_Ha[15:25,15:25]);,minmax(rat)

cube_results[0]=OIII_Hb
cube_results[1]=err_OIII_Hb
writefits,path+name+suff1_write+'_OIII_Hb.fits',cube_results;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_OIII_Hb.fits'

cube_results[0]=NII_Ha
cube_results[1]=err_NII_Ha
writefits,path+name+suff1_write+'_NII_Ha.fits',cube_results;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_NII_Ha.fits'

if keyword_set(make_ps) then begin
          ;hardcopy,/open,file=Path+name+'_BPT.eps',/color
  begin_ps,path+name+suff1_write+'_BPT.eps'
  loadct,12
  plotsym,0,1,/fill
  plot,NII_Ha,OIII_Hb,psym=8,xrange=[-1.5,0.5],yrange=[-1.5,1],xsty=1,ysty=1,ythick=5,xthick=5,thick=4,xtitle=textoidl('Log([NII]6583/H\alpha)',font=0),ytitle=textoidl('log([OIII]5007/H\beta)',font=0),$
       font=0,charsize=1.5,charthick=3

                       ;oploterror,NII_Ha,OIII_Hb,err_NII_Ha,err_OIII_Hb,psym=8
  oplot,(indgen(24)*0.1-2),0.61/((indgen(24)*0.1-2)-0.47)+1.19,linestyle=1,thick=7,COLOR=fsc_color('red') ;Kewley 
  oplot,(indgen(20)*0.1-2),0.61/((indgen(20)*0.1-2)-0.05)+1.3,linestyle=2,thick=7,COLOR=fsc_color('red') ;Kauffman
  !p.thick=3
 end_ps
;hardcopy,/close
endif 

mask_AGN=fltarr(1)

in=where(finite(OIII_Hb+NII_Ha),cn)
if cn gt 0 then mask_AGN[in]=1.

;if n_elements(WHERE(OIII_Hb GE 0.61/(NII_Ha-0.47)+1.19)) then begin
;   print,'No Agn' 
;endif else begin
   in=where(OIII_Hb GE 0.61/(NII_Ha-0.05)+1.3,cn)
   if cn gt 0 then mask_AGN[in]=2.
;mask_AGN(WHERE(OIII_Hb GE 0.61/(NII_Ha-0.05)+1.3))=(WHERE(OIII_Hb GT 0.61/(NII_Ha-0.05)+1.3))*0.+'Nan';kauffman
;mask_AGN(WHERE(OIII_Hb GE 0.61/(NII_Ha-0.47)+1.19))=(WHERE(OIII_Hb GT 0.61/(NII_Ha-0.47)+1.19))*0.+'Nan';kewley

   in=where(OIII_Hb GE 0.61/(NII_Ha-0.47)+1.19,cn)
   if cn gt 0 then mask_AGN[in]=3.
;endelse
writefits,path+name+suff1_write+'_AGN.fits',mask_AGN;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_AGN.fits',min=0.5



;Map ionization O2O3
O2O3=alog10(f_OII/f_OIIIb)
e_O2O3=sqrt((e_OII/f_OIIIb)^2+(e_OIIIb*f_OII/(f_OIIIb^2))^2)
err_O2O3=abs(e_O2O3/(10.0^O2O3*alog(10.)))

cube_results[0]=O2O3
cube_results[1]=err_O2O3

writefits,path+name+suff1_write+'_O2O3.fits',cube_results;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_O2O3.fits'

;Map ionization N2O2
N2O2=alog10(f_NIIb/f_OII)
e_N2O2=sqrt((e_NIIb/f_OII)^2+(e_OII*f_NIIb/(f_OII^2))^2)
err_N2O2=abs(e_N2O2/(10.0^N2O2*alog(10.)))

cube_results[0]=N2O2
cube_results[1]=err_N2O2

writefits,path+name+suff1_write+'_N2O2.fits',cube_results;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_N2O2.fits'

;Map N/O - Perez-Montero & Contini 2009

cube_results[0]=0.93*N2O2-0.2
cube_results[1]=0.93*err_N2O2

writefits,path+name+suff1_write+'_NO.fits',cube_results;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_NO.fits'

;;;;;; metallicities

;Metallicity N2 Raffa 2013

cube_results[0]=8.75+0.46*NII_Ha;*mask_AGN
cube_results[1]=0.46*err_NII_Ha

writefits,path+name+suff1_write+'_OH_R13_N2.fits',cube_results;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_OH_R13_N2.fits'

;Metallicity O3N2 Raffa 2013
O3N2=alog10((f_OIIIb/f_Hb)/(f_NIIb/f_Ha))
e_O3N2=sqrt( (e_OIII_Hb/(f_NIIb/f_Ha) )^2+(e_NII_Ha*(f_OIIIb/f_Hb)/(f_NIIb/f_Ha)^2)^2)
err_O3N2=abs(e_O3N2/(10.0^O3N2*alog(10.)))

cube_results[0]= 8.53-0.22*O3N2
cube_results[1]= 0.22*err_O3N2

writefits,path+name+suff1_write+'_OH_R13_O3N2.fits',cube_results;*mask_AGN;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_OH_R13_O3N2.fits'

;Metallicity N2 Pettini 2004

cube_results[0]=(9.37+2.03*NII_Ha+1.26*NII_Ha^2+0.32*NII_Ha^3);*mask_AGN
cube_results[1]=err_NII_Ha*abs(2.03+2.52*NII_Ha+0.96*NII_Ha^2);*mask_AGN

writefits,path+name+suff1_write+'_OH_N2.fits',cube_results;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_OH_N2.fits'


;Metallicity O3N2 Pettini 2004
O3N2=alog10((f_OIIIb/f_Hb)/(f_NIIb/f_Ha))
e_O3N2=sqrt( (e_OIII_Hb/(f_NIIb/f_Ha) )^2+(e_NII_Ha*(f_OIIIb/f_Hb)/(f_NIIb/f_Ha)^2)^2)
err_O3N2=abs(e_O3N2/(10.0^O3N2*alog(10.)))

cube_results[0]= 8.73-0.32*O3N2
cube_results[1]= 0.32*err_O3N2

if keyword_set(verify) then print,'Met_O3N2 at center and error: ',minmax(8.73-0.32*o3n2[15:25,15:25]),minmax(0.32*err_o3n2[15:25,15:25]);,minmax(rat)


writefits,path+name+suff1_write+'_OH_O3N2.fits',cube_results;*mask_AGN;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_OH_O3N2.fits'


;Metallicity D16 Dopita 16
f_SII=f_SIIa+f_SIIb
e_SII=sqrt(e_SIIa^2+e_SIIb^2)
NII_SII=alog10(f_NIIb/f_SII)
e_NII_SII=sqrt((e_NIIb/f_SII)^2+(e_SII*f_NIIb/(f_SII^2))^2)
err_NII_SII=abs(e_NII_SII/(10.0^NII_SII*alog(10.)))

;D16=8.77 + NII_SII + 0.264 * NII_Ha
D16=8.77 + NII_SII + 0.264 * NII_Ha + 0.45 * (NII_SII + 0.264 * NII_Ha + 0.3)^5
err_D16=sqrt( e_NII_SII^2 + (0.264*e_NII_Ha)^2  ) ;+ todo el show del ^5

cube_results[0]=D16
cube_results[1]= err_D16
writefits,path+name+suff1_write+'_OH_D16.fits',cube_results;*mask_AGN;,header

;Metallicity O3N2  Perez-Montero & Contini 2009

cube_results[0]= 8.33-0.31*O3N2-0.35*(0.93*N2O2-0.2)
cube_results[1]= sqrt((0.31*err_O3N2)^2+(0.35*0.93*err_N2O2)^2)

writefits,path+name+suff1_write+'_OH_O3N2_NO.fits',cube_results;*mask_AGN;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_OH_O3N2_NO.fits'


; Metallicity R23 Kobulnicky

ii=-1

repeat begin
  ii++
  load_map,path+name+suff_read+'_OIIIb'+suffix,sn_o3_rat-ii,f_OIIIb1,e_OIIIb1,_extra=e
  load_map,path+name+suff_read+'_OIIIa'+suffix,sn_o3_rat-ii,f_OIIIa1,e_OIIIa1,_extra=e

  o3_rat=f_OIIIa1/f_OIIIb1
  in=where(finite(o3_rat),cn)
endrep until (cn gt 5 or sn_o3_rat-ii lt 2) or ii ge 0

if cn gt 5 then begin 
   m_o3_rat=median(o3_rat[in])
   e_m_o3_rat=median(abs(o3_rat[in]-m_o3_rat))/0.6745
   if keyword_set(verify) then print,m_o3_rat,e_m_o3_rat
endif else begin
   m_o3_rat=0.3
   if cn eq 0 then e_m_o3_rat=0.3 else e_m_o3_rat=median(abs(o3_rat[in]-m_o3_rat))/0.6745       ;  if cn eq 0 then e_m_o3_rat=0.3 else 
endelse

in=where(~ finite(f_OIIIa))
f_OIIIa[in]=m_o3_rat*f_OIIIb[in]

e_OIIIa[in]=sqrt((m_o3_rat*e_OIIIb[in])^2+(e_m_o3_rat*f_OIIIb[in])^2)

R23=alog10((f_OII+f_OIIIa+f_OIIIb)/f_Hb)	
err_ratioR23=sqrt((e_OII^2+e_OIIIa^2+e_OIIIb^2)/f_Hb^2 + ((10.0^R23/f_Hb)*e_Hb)^2 ) 
err_R23=err_ratioR23/(10.0^R23*alog(10.))

cube_results[0]=R23
cube_results[1]=err_R23
writefits,path+name+suff1_write+'_R32.fits',cube_results;*mask_AGN;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_R32.fits'

O32=alog10((f_OIIIa+f_OIIIb)/f_OII)
e_O3O2=sqrt((e_OIIIb/f_OII)^2+(e_OII*f_OIIIb/(f_OII^2))^2+(e_OIIIa/f_OII)^2+(e_OII*f_OIIIa/(f_OII^2))^2)
err_O32=abs(e_O3O2/(10.0^O32*alog(10.)))


cube_results[0]=0.8*o32-3.02
cube_results[1]=0.8*o32-3.02

writefits,path+name+suff1_write+'_U.fits',cube_results;*mask_AGN;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_U.fits'

if keyword_set(verify) then print,'O3O2 at center and error ',minmax(O32[15:25,15:25],/nan),minmax(err_O32[15:25,15:25],/nan)

N2O2=alog10(f_NIIb/f_OII)

nan=WHERE(~finite(n2o2),nnan)
in=where(finite(n2o2))
up=WHERE(N2O2[in] GE -1.2,nup)
low=WHERE(N2O2[in] lt -1.2,nlow)


if keyword_set(verify) then print,nup,nlow,nnan

if nup gt 0 then begin
   ohim[in[up]]=8.7
   e_ohim[in[up]]=0.1
   t_r32=R23[in[up]]
   e_t_r32=err_R23[in[up]]
   for i=0,2 do begin 
      qim[in[up]]=ionization_parameter(O32[in[up]],ohim[in[up]],err_O32[in[up]],e_ohim[in[up]],e_qim)
      errq[in[up]]=e_qim
; if keyword_set(verify) then print,'KK04 at center and error, Q error: ',minmax(ohim[15:25,15:25],/nan),minmax(e_ohim[15:25,15:25],/nan),minmax(qim[15:25,15:25],/nan),minmax(errq[15:25,15:25],/nan)
      
      
      ohim[in[up]]=9.72-0.777*t_r32-0.951*t_r32^2-0.072*t_r32^3-0.811*t_r32^4 $
          -qim[in[up]]*(0.0737-0.0713*t_r32-0.141*t_r32^2+0.0373*t_r32^3-0.058*t_r32^4)
 
      tt1=(-0.777-2*0.951*t_r32-3*0.072*t_r32^2-4*0.811*t_r32^3)
      tt2=(-0.0713-2*0.141*t_r32+3*0.0373*t_r32^2-4*0.058*t_r32^3)
      e_ohim[in[up]]=sqrt( (tt1-qim[in[up]]*tt2)^2*e_t_r32^2+errq[in[up]]^2*(0.0737-0.0713*t_r32-0.141*t_r32^2+0.0373*t_r32^3-0.058*t_r32^4)^2) 
     
   endfor   
endif

if nlow gt 0 then begin
   ohim[in[low]]=8.2
   e_ohim[in[low]]=0.1
   t_r32=R23[in[low]]   
   e_t_r32=err_R23[in[low]]
   for i=0,2 do begin 
      qim[in[low]]=ionization_parameter(O32[in[low]],ohim[in[low]],err_O32[in[low]],e_ohim[in[low]],e_qim)      
      errq[in[low]]=e_qim

      ohim[in[low]]=9.40+4.65*t_r32-3.17*t_r32^2-qim[in[low]]*(0.272+0.547*t_r32-0.513*t_r32^2)
      
      tt1=(4.65-2*3.17*t_r32)
      tt2=(0.547-2*0.513*t_r32)
      e_ohim[in[low]]=sqrt( (tt1-qim[in[low]]*tt2)^2*e_t_r32^2+errq[in[low]]^2*(0.272+0.547*t_r32-0.513*t_r32^2)^2)
     
   endfor
endif

in=where(finite(qim))
qim[in]-=alog10(3e10)
cube_results[0]=qim
cube_results[1]=errq

writefits,path+name+suff1_write+'_q.fits',cube_results;*mask_AGN;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_q.fits'

;tvscl,qim/errq,/nan

 if keyword_set(verify) then print,ohim[0,0]
 if keyword_set(verify) then print,minmax(ohim)

if keyword_set(verify) then print,'KK04 at center and error, Q error: ',minmax(ohim[15:25,15:25],/nan),minmax(e_ohim[15:25,15:25],/nan),minmax(qim[15:25,15:25],/nan),minmax(errq[15:25,15:25],/nan)

;OH_R23=12-2.939-0.2*R23-0.237*R23^2-0.305*R23^3-0.0283*R23^4-O32*(0.0047-0.0221*R23-0.102*R23^2-0.0817*R23^3-0.00717*R23^4)
;OH_R23(low_branch)=12-4.944+0.767*R23(low_branch)+0.602*R23(low_branch)^2-O32(low_branch)*(0.29+0.332*R23(low_branch)-0.331*R23(low_branch)^2)


if ohim eq 0 then ohim=!VALUES.F_NAN
if e_ohim eq 0 then e_ohim=!VALUES.F_NAN

cube_results[0]=ohim
cube_results[1]=e_ohim

writefits,path+name+suff1_write+'_OH_R23.fits',cube_results;*mask_AGN;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_OH_R23.fits',min=7.8,max=9.3

; Metallicity R23 Tremonti
cube_results[0]=(9.185-0.313*R23-0.264*R23^2-0.321*R23^3);*mask_AGN
cube_results[1]= abs((-0.313-2*0.264*R23-3*0.321*R23^2))*err_R23;*mask_AGN
writefits,path+name+suff1_write+'_OH_R23_T.fits',cube_results;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_OH_R23_T.fits'


;Pylugin
errados=f_OII/f_Hb
erratres=(f_OIIIa+f_OIIIb)/f_Hb
enados=(f_NIIa+f_NIIb)/f_Hb
essados=(f_SIIa+f_SIIb)/f_Hb
e_errados=sqrt((e_OII/f_Hb)^2+(f_OII*e_Hb/(f_Hb^2))^2)
e_erratres=sqrt((e_OIIIa/f_Hb)^2+(e_OIIIb/f_Hb)^2+((f_OIIIa+f_OIIIb)*e_Hb/(f_Hb^2))^2)
e_enados=sqrt((e_NIIa/f_Hb)^2+(e_NIIb/f_Hb)^2+((f_NIIa+f_NIIb)*e_Hb/(f_Hb^2))^2)
e_essados=sqrt((e_SIIa/f_Hb)^2+(e_SIIb/f_Hb)^2+((f_SIIa+f_SIIb)*e_Hb/(f_Hb^2))^2)

ll=size(enados,/dim)
if (ll[0] lt 31) then begin
   print,'I suppose you want the azimutal average mean of the PYLUGIN measurement'
   if (alog10(enados) lt -0.1) then begin
        cube_results[0]=8.642+0.077*alog10(erratres)+0.411*alog10(errados)+0.601*alog10(enados/errados)
        cube_results[1]=1/alog(10.)*sqrt((0.077*e_erratres/erratres)^2+(-0.190*e_errados/errados)^2+(0.601*e_enados/enados)^2)
   endif else begin
      cube_results[0]=8.606-0.105*alog10(erratres)-0.410*alog10(errados)-0.150*alog10(enados/errados) 
      cube_results[1]=1/alog(10.)*sqrt((-0.105*e_erratres/erratres)^2+(-0.260*e_errados/errados)^2+(-0.150*e_enados/enados)^2) 
   endelse
endif else begin
    print,'I suppose you want the common PYLUGIN measurement'
    if (alog10(enados) lt -0.1) then begin
          cube_results[i,j,0]=8.642+0.077*alog10(erratres)+0.411*alog10(errados)+0.601*alog10(enados/errados)
          cube_results[i,j,1]=1/alog(10.)*sqrt((0.077*e_erratres/erratres)^2+(-0.190*e_errados/errados)^2+(0.601*e_enados/enados)^2)
    endif else begin
        cube_results[i,j,0]=8.606-0.105*alog10(erratres)-0.410*alog10(errados)-0.150*alog10(enados/errados) 
        cube_results[i,j,1]=1/alog(10.)*sqrt((-0.105*e_erratres/erratres)^2+(-0.260*e_errados/errados)^2+(-0.150*e_enados/enados)^2) 
    endelse
endelse

writefits,path+name+suff1_write+'_OH_P10.fits',cube_results;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_OH_P10.fits'


ll=size(enados,/dim)
if (ll[0] lt 31) then begin
    print,'I suppose you want the azimutal average mean of the PYLUGIN-11 measurement'
    if (alog10(enados) lt -0.1) then begin
        cube_results[0]=8.456+0.082*alog10(erratres)+0.391*alog10(enados)+0.290*alog10(enados/essados)
        cube_results[1]=1/alog(10.)*sqrt((0.082*e_erratres/erratres)^2+(0.681*e_enados/enados)^2+(-0.290*e_essados/essados)^2)
    endif else begin
      cube_results[0]=8.454-0.216*alog10(erratres)-0.362*alog10(essados)-0.101*alog10(enados/essados) 
      cube_results[1]=1/alog(10.)*sqrt((-0.216*e_erratres/erratres)^2+(-0.101*e_enados/enados)^2+(-0.261*e_essados/essados)^2) 
    endelse
endif else begin
    print,'I suppose you want the common PYLUGIN-11 measurement'
    if (alog10(enados) lt -0.1) then begin
          cube_results[0]=8.456+0.082*alog10(erratres)+0.391*alog10(enados)+0.290*alog10(enados/essados)
          cube_results[1]=1/alog(10.)*sqrt((0.082*e_erratres/erratres)^2+(0.681*e_enados/enados)^2+(-0.290*e_essados/essados)^2)
    endif else begin
        cube_results[0]=8.454-0.216*alog10(erratres)-0.362*alog10(essados)-0.101*alog10(enados/essados) 
        cube_results[1]=1/alog(10.)*sqrt((-0.216*e_erratres/erratres)^2+(-0.101*e_enados/enados)^2+(-0.261*e_essados/essados)^2) 
    endelse
endelse

writefits,path+name+suff1_write+'_OH_P11.fits',cube_results;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_OH_P11.fits'

; Star-formation map
;f_Ha=f_Ha*10^(c1*(1.13-0.3))
EMpc=3.0857d24
c =3.e8
Lo=3.826d33
pi=3.1415926
DL=dl(redshift,/silent)* EMpc
L_Ha=double(4*!DPI*(DL^2)*scale*f_Ha)*(1+redshift)
cube_results[0]=(7.9e-42*L_Ha)
err_L_Ha=double(4*!PI*(DL^2)*scale*e_Ha)*(1+redshift)
cube_results[1]=(7.9e-42*err_L_Ha)
cube_results=cube_results                                                                 ;*1e-16
writefits,path+name+suff1_write+'_SFR.fits',cube_results;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_SFR.fits'

conver=arctokpc(1.,redshift)
cube_results=cube_results/(conver*conver)
writefits,path+name+suff1_write+'_SFRsqkpc.fits',cube_results;,header
if keyword_set(make_ps) then colorflux,path+name+suff1_write+'_SFRsqkpc.fits'





END
