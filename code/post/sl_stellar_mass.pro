function sl_stellar_mass,sl_fit_file,redshift,ini=ini

par=read_sl_fit(sl_fit_file)

     ; Star-formation map
     ;f_Ha=f_Ha*10^(c1*(1.13-0.3))
EMpc=3.0857d24
Lo=3.826d33
DL=dl(redshift,/silent)* EMpc


if keyword_set(ini) then mass=double(par.m_initot*4*!DPI*DL*DL/Lo) $  ;falta 10e-16 perque he ficat unitats de flux a l'espectre observat
                    else mass=double(par.m_cortot*4*!DPI*DL*DL/Lo)   ;falta 10e-16 perque he ficat unitats de flux a l'espectre observat

return,alog10(mass)    ;mass
END



; M⋆ = Mcor tot × 10−16 × 4πd2 × (3.826 × 1033)−1      to obtain the present mass in stars (in M⊙).

;Mini⋆ = Mini tot × 10−16 × 4πd2 × (3.826 × 1033)−1    how many M⊙’s have been processed into stars throughout the galaxies life.