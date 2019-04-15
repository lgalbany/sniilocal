!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Current version: PANcMEx_StarlightChains_v03.for

! (Called tmp_v022 & tmp_v02f.for during tests)

! BIG/IMPORTANT CHANGE! Optimize fn: I now allow for the scaling (normalization)
!   factor to be evaluated ANALITICALY, such that dchi2/dfn = 0. This
!   was something which should have been implemented long ago, but I
!   never did ... It seems that by letting fn free, as a parameter, most
!   pop vectors experimented during the chains had their chain energies
!   larger than they could be, just because of a bad/unlucky choice of
!   fn! I did not test this thoroughly, but the few tests I did suggest
!   this, and it makes perfect sense! With this new option, we may have
!   to evaluate many less steps than before!
!
!   The changes associated to this are marked by !OPTfn!. Apart from the
!   chane in arq_config and the common/Optimize_fn_OPT/, the main
!   changes are in routines SlowF_chi2, VeryFastF_chi2, and
!   MEX_ScubiduCadeVoce. The to-do-it-or-not-to-do-it flag
!   IsOptimize_fn_OPT (read from arq_config) is saved (in an odd
!   place...) in the output file

! Is1stLineHeader: I also removed the ugly Is1stLineHeader flag from arq_config. Now
!     arq_obs can have ANY number of header lines starting with a '#',
!     with no need for any previous information on that (as always did
!     for base spectra)

! Added i_SaveBestSingleCompFit in this output file to help reading
!     macros. This info determines whether a 5th spectral columns is
!     saved in the or not.

! Spectral flux format: Relucted to do this, but finally changed spectra
!     format in output file. (Sometimes, for crazy spectra, number would
!     be much larger than 1 (in units of fobs_norm...) and the format
!     would go correspondingly crazy. A harmless thing anyway.

! Flagged pix in output weight col: Instead of assigining -2 to every flagged
!	  pixel, I now keep the values of the input flag (possibly resampled) in
!     the output. If the input flag-spec has different values meaning different 
!	  situations (sky lines, hot pix, etc etc), this alows one to identify this
!	  in the output file itself. A purely cosmetical change!
!	  ElCid@Sanchhica - 10/Feb/2012


!2BEDONE! Some to-be-done things... (some marked in the code with !2BEDONE!)

!     * Implement save chi2_FIR , chi2_QHR... in the PDF-sampling option! (Some comments
!     were included in preparation for that.)
!
!     * Implement Jesus-Maiz law! (He does not have the equations yet!!)

!     * Introduce in config the meaning of flag-values? This may be
!     trickier than anticipated, since th eidea that flag > 1 signals a
!     bad pixel is more embedded in the code than I recalled.For
!     instance, in the PRE-PROCESSING I step I set Oflag_obs(i) = 99 for
!     pixels outside the Olsyn_ini to Olsyn_fin interval and those with
!     negative Oerr_obs(i). This only means that this change needs
!     care... since it is not really relavent, I'll forget it. The
!     flag-resampling issue is one of the reasons to worry abouyt
!     this...

!     * Think about wei-clip-filter since it will probably be used with
!     Bernd's errors...

!     * Minimum Nl_eff in clipping in 50 ... should be more smarly scaled!

!     * ATT: fobs_norm computation IGNORES bad (flaged/masked) pixels!!
!            This is dangerous and could/should be changed one day...

!     * Clean up the code!! This looks more like a diary!

! ElCid@Sanchica - 06--10/Feb/2012


! Minor & purely cosmetic changes wrt PANcMEx_StarlightChains_v02b.for:
! 	> idt & wdt_* times format changed.
!	> SSP_x format changed
!
! ElCid@Sanchica - 25/Nov/2001
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Was PANcMEx_StarlightChains_v01.for - started to tidy it up in GRX/Ago/2011!!
!	ElCid@Sanchica - 15/Ago/2011
!
!	The old compilation warning
!           do u=u_low,u_upp,du                                            
!           Warning: Deleted feature: Step expression in DO loop at (1) must be integer
!     was fixed with a new version of routine GaussSmooth which replaces the
!     (real) u-loop by an (integer) i_u-loop. 
!
!	Some testing was done, but not all things were implemented .... Anyway, decided to  
!	rename to PANcMEx_StarlightChains_v02.for to start CALIFA work!
!
!	ElCid@Sanchica - 08/Oct/2011

!FIX4CALIFA! ATT: In PANcMEx_StarlightChains_v01.for, used for the U/LIRGS work, 
!	the parameter n_pt_MaxBF = 10000, but this is TOO LARGE if one is dealing
!	with just one extinction, as I am now doing for the CALIFA data. In this case 
!	one has nBF_Dim = 1 and thus
!		n_pt_PerDim = (float(n_pt_MaxBF))**(1. / float(nBF_Dim)) = 10000
!		n_pt_BF     = n_pt_PerDim**nBF_Dim = 10000
!	For instance, in good-old StarlightChains_v04 the equivalent was NAV_try = 100, 
!	ie, 100 times less! 
!	I thus changed it to 100 here.... it's a 1 line change, marked with !FIX4CALIFA!
!	If I now try 2 extinctions, then there will be only 100**(1/2) = 10 pts per AV
!	Should rethink if this need to go to a config file or set it up on the code depending
!	on how many extinction I actually use in the fit...
!
!	ElCid@Sanchica - 10/Oct/2011

!	Found & fixed 2 BUGs in the same place: (1) In iEXOs_IsAve_x_of_lambda I used x_min, 
!	even if EXOs_PopVector_option = 'AVE'. (2) The AV_tot array needed in AverageXoverLambda
!	was NOT defined! This may have caused some damage in my MEx-fits of U/LIRGs...
!
!	I have relocated the READ MASK FILE block to BEFORE the read-observed-spectrum block, 
!	and used the masks in the estimation of a S/N in the absence of an error spectrum.
!	The motivation is that we are building CALIFA masks for Pre-Processing which mask away
!	main bugs, like sky-lines (besides EmLines, of course), and it is unfair to include any
!	such "faulty" region in an estimate of the rms which represents the error.
!
!	Added more blurb on the chains output headers (at Natalia's & Rafa's request)
!
!	Some experiments w/CALIFA suggest that EX0s produce x-solution which are more focused 
!	than those in the chains... need to understand this better! Is this related to the 
!	infamous bermuda triangle issue?? Is is a localization into a local minimum while fine-tuning?
!
!	The chains-output contains Mcor_j, NOT Mini_j. I can add that ~ easily, but
!	is is also easy to do the transformation using the information already provided in
!	the output, so I decided NOT to add more bytes to arq_out.

!	ElCid@Sanchica - 17/Oct/2011


!FIX-DAE2-BUG! ATT! Found a serious PROBLEM. In dae2.iaa.es, the gfortran version is 4.6.1,
!	while here in my humble MacBook it is (thank's God) older: 4.4.1. While doind some
!	fits of **horrible** ARP220 CALIFA spectra, the code stopped with a segmentation fault at dae2,
!	but not in my MacBook! Furungating in the code I realized that @ dae2 the FUNCTION fun_v0_chi2
!	was being called (from with those dammed NR routines .... mnbrack & golden)
!	with completrly absurd values... The input crappy-spectrum deserved such a bad treatment, but 
!	the code should not crash! I have thus changed both fun_v0_chi2 & fun_vd_chi2 FUNCTIONs,
!	such that they now return an "infinite" (= 3.0e38) chi2 everytime they are
!	called with vd or |v0| > 5000 km/s. For some reason gfortran 4.4.1 can cope with the infinities...
!	This should really not make any difference for minimally decent spectra (I've checked it)!
!	The changes are marked with !FIX-DAE2-BUG!
!
!	BTW, I have also fixed the gfortran (correct!) complaint/warning:
!   	call HYPERZ_red_spec(l_aux,f_int_aux,1.0,f_obs_aux,1,1,ilaw_rd)                             1
!		Warning: Actual argument contains too few elements for dummy argument 'wl' (1/1300) at (1)
!	This was a silly and harmless thing. Fix marked with !FIX-DAE2-WARNING-HYPERZ!
!
!	ElCid@Sanchica - 21/Oct/2011


!MORE-BUGS: Fixed an un-noticed bug in GaussSmooth (the new version which replaced do u=... 
!	by do i_u=... to avoid/skip compilation Warnings). See comments on SUBROUTINE GaussSmooth.
!       Many CALIFA galaxies were analised with the bugged version! No harm, though, since the bug
!	was spurred by the absurd vd's tried by NR's mnbrack). Nonetheless, some segmentation faults
!       did occur because of this (rgb@iaa.es, for instance, had several such cases in his 1st IAA-GRID runs).
!
!	While doing this I noticed another bug, now in FUNCTION fun_vd_chi2(vd). My dirty fix of the 
!	mnbrack-induced bug was wrongly implemented! I wrote  fun_v0_chi2 = 3.0e38 instead of 
!	fun_vd_chi2 = 3.0e38 as a result of a careless cut & paste... stupid! 
!	  
!   I therefore RENAMED this version from PANcMEx_StarlightChains_v02.for to PANcMEx_StarlightChains_v02b.for.
!
!	ElCid@Sanchica - 05/Nov/2011
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Minor fixes:

!      (1) Fixed Z_base format to f8.5 to satisfy Paula Coelho. She says
!      this is needed when Z is expressed in log units (eg, -1.5, as in
!      Vazdekis models)
  
!      (2) Natalia realized that when the SN-window has no info to compute
!      SN, then the trick to use all the spectrum to compute a fake SN
!      screws up all next galaxies in a grid! This happens because
!      llow_SN & llup_SN are read from the grid file, but changed
!      forever once this safety-trick is activated... To fix this bug I
!      simply renamed llow_SN & llupp_SN in the driver to the main code
!      to Ollow_SN & Ollupp_SN, and restore them everytime I call the
!      main code.

!     (3) One thing which alway anoyed me was that when repeating large
!     grids using the same grid-file, it to ~ 1 second per run (all in
!     the rest-stuff block) to reach the SkipExistingOutFiles-check. For
!     large grids, this implies a huge waste of time! I've changed the
!     order of things (read input info and test if file exist 1st, and
!     only then reset arrays) to fix this "bug".

!     (4) In Jan/2011 I had the impression the idum_orig = 0 trick was NOT
!     working, but now this bug seems to have vanished !

!     (5) One of the age-old compilation warnings was: 
!          if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'          
!          Warning: Deleted feature: PAUSE statement at (1)
!     This was "fixed" replacing pause by a full stop!



!	SHOULD IMPLEMENT AN AV(Y) = K * AV(ALL)  (K=2 WOULD BE STANDARD...)

!	SHOULD MAKE STRINGS LARGER TO COPE WIHT LARGER DIR & FILE NAMES!

! 	Shoud re-check PDF option, and maybe save only 1 out of "M" steps 
!	to avoid huge PDF files with little x-variations.... anyway, re-check it.

! 	DO I HAVE CONTROL OVER AV-RANGE REALLY? CAN I FIX AV_low = 1 and AV_upp =1 ou 2 for instance??

! 	IDEA: maybe i can impose limits in, say, EW(G-band) using PHO! .... think about it!

!	PHO stuff is still messy/incomplete/untested.

!   ElCid@Sanchica - 15/Ago/2011
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



c     ELR implemented:-)
c

! TO-DO's ... 11/Jan/2011

! - Abilio vai ver AB_mag's pra PHO-input...
! - MUST REVISE ALL NEW-FORMULAE & NUMBERS, particularly those for QHR...
! - The EX0s selection of relevant components does NOT take into account
!     non-optical data... we may throw out a population responsible
!     for, say, most Hbeta or FIR emission!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Read the comments on the output block called
!             Preparations for ETC [= FIR & QHR & PHO] output           !
!  Probably all my previous FIRc-runs [w/FIRcMEx_StarlightChains_v02.for, 
!  in Jan-Feb/2010] stored incorrect values of FIRModlogLFIR because I 
!  did not realize that until now!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

*************************************************************************
***                   PANcMEx_StarlightChains_v01.for                 ***
*************************************************************************

c     PANCRHOMATIC-contraints & MULTIPLE EXTINCTION STARLIGHT:
c     [FIR and/or QHR-ELR and/or PHO-constraints]

c     ABSTRACT: Expanding upon FIRcMEx_StarlightChains_v02.for [from
c     ~Jan-Feb/2010 work], this version now accounts for two other types
c     of input constraints: ionizing photon rate (QH) Related [hence
c     "QHR"] stuff and photometric data ("PHO"). Within QHR, one can use
c     emission line log-luminosities (logY) and/or use an
c     EmissionLineRatio (eg, Ha/Hb) to constrain the fit! (ELR-fits
c     should be particularly useful for MEx-fits, and help contraining
c     young & reddened populations.)  The notes below report some of the
c     things associated to this new version.


c     **> NEW ETC-file: The user now provides a **single** separate
c     input file (arq_ETCinfo) with FIR and/on QHR and/on PHO data. The
c     file is such that when reading it [in routine Read_ETC_Info] the
c     code realizes which FIR/QHR/PHO data are available and how to use
c     them. The scheme is pretty flexible. For instance, you do not need
c     to specify info for all FIR + QHR + PHO constraints, and even if
c     you do so, they can be accomodated in any order in arq_ETCinfo. See
c     comments in routine Read_ETC_Info for more on this.

c     The GRID FILE now contains IsFIRcOn, IsQHRcOn & IsPHOcOn flags.
c     These admit 3 options: 
c
c        1 = YES, use info to the constrain the fit. 
c        0 = NO, do not use the info in arq_ETCinfo at all.
c       -1 = do not use the ETC-info in the fit, but instead PREDICT IT
c
c     [eg, predict the FIR luminosity even without any FIR-data, or
c     predict PHOtometric luminosities given some input filters.]
c     Options 1 & -1 produce FIR/QHR/PHO-related output in arq_out,
c     while option 0 does not.


c     **> HANDLING 4 CHI2's! The way the FIR+QHR+PHO-info is actually
c     used in the fit is identical to that in earlier FIRc-work, ie, via
c     gaussian-flat-top likelyhood profiles which allow for range
c     fitting, normal fitting and scaling of the corresponding chi2 wrt
c     the optical-fit chi2. In practice, we now have
c
c     chi2 = chi2_Opt + chi2_ETC
c
c     where chi2_Opt is the good-old-chi2 of the optical fit (before all
c     new-stuff/ETC was implemented), and
c
c     chi2_ETC = chi2_FIR + chi2_QHR + chi2_PHO
c
c     The FIR, QHR & PHO-chi2 functions are all called from within the
c     old-chi2 functions/routines: Slow_Fchi2, FK_chi2, VeryFastF_chi2,
c     fun_vd_chi2 & fun_v0_chi2. [These are the 5 key chi2 routines, as
c     far as the code is concerned.] In fact, a single FUNCTION
c     ETC_CalcETCchi2(x,AV_tot,chi2_OPTICAL,WhereFrom) is called, which
c     joins the chi2's from the new FUNCTIONs FIR_CalcFIRchi2,
c     QHR_CalcQHRchi2 & PHO_CalcPHOchi2. Besides these functions, there
c     are initialization SUBROUTINES (FIRInitialization,
c     QHRInitialization, PHOInitialization, Read_ETCinfo) and a bunch of
c     commons (eg, QHR_BaseStuff, PHO_DataStuff, FIR_OutputStuff...)
c     involved in the implementation of these extra-info in the spectral
c     fits.
c
c     chi2 = chi2_Opt * [1 + chi2_FIR/chi2_Opt + chi2_QHR/chi2_Opt + chi2_PHO/chi2_Opt]

c     OBS: Curiously, with these changes it became somewhat tricky to
c     obtain a pure optical chi2... The reason is that all other
c     chi2-functions need Chi2_Opt for scaling purposes, while the
c     values returned by the old chi2-functions (eg, FK_chi2) are now
c     TOTAL chi2's... In the output section of the code there are a few
c     manouvers which illustrate how to go about this.
c

c     **> FIR/QHR/PHO CHI2-SCALE-FACTORS: FIRc-fits introduce only one
c     observable, FIR_logLFIR_obs, and thus only one FIR-to-Optical
c     chi2-scaling factor is needed: FIRChi2ScaleFactor. QHRc & PHOc
c     fits, however, may introduce more than one observable each
c     QHR_logYobs(iY), iY=1...NQHR_Ys, and PHO_logYobs(iY),
c     iY=1...NPHO_Ys. Each of this comes with its own scaling factor,
c     QHR_Chi2ScaleFactor(iY) or PHO_Chi2ScaleFactor(iY), but in order
c     to attribute a global meaning to these, there is also a
c     QHR_GlobalChi2ScaleFactor (and similarly for PHO) which says what
c     is the global weight that the QHR chi2 should have wrt to the
c     optical one. For instance, if you use two emission lines (Ha and
c     Hb), each with QHR_Chi2ScaleFactor(iY) = 1, and chose
c     QHR_GlobalChi2ScaleFactor = 1, in practice each of the lines will
c     be scaled by 1/2 of the optical chi2.


c     **> FRAC2MODEL, LOW <= OBS <= UPP - RANGE & FLAT-TOP-GAUSSIAN
c     LYKELYHOODS: For any of FIR/QHR/PHO, the user inform in
c     arq_ETCinfo a total luminosity and the fraction of it to be
c     modeled (a trick to account for, say, aperture effects), and this
c     defines the observed value used in the fit [eg, QHR_logY_obs(iY)
c     for the log of the iY emission line luminosity]. Then a RANGE is
c     specified, from which one defined *_low & *_upp limits for the
c     observable (eg, logLFIR). If RANGE = 0, then UPP = LOW and the
c     likelyhood wil be a gaussian with the given input error. If RANGE
c     NE 0, then we create a flat-top-cored gaussian, with gaussian
c     wings below *_LOW & above *_UPP, but flat behaviour (corresponding
c     to chi2 = 0!) in the middle, such that the fit will try to
c     constrain the value of the observable between *_LOW & *_UPP. For
c     RANGE > 0, we set *_UPP = *_OBS and *_LOW = *_OBS - RANGE
c     (suitable for upper-limit fitting), whereas (at Natalia's
c     suggestion) for RANGE < 0, we do the opposite (ie, *_LOW = *_OBS
c     & *_UPP = *_OBS + RANGE, suitable for lower limit fitting).
c     Setting a very small error, these ranges work like chi2 walls
c     ("likelyhhod abisses"). This scheme works for FIR, QHR & PHO info.


c     **> LUMINOSITIES & FLUXES: All galaxy-data in arq_ETCinfo is given
c     in **LUMINOSITIES** (Lsun for FIR-work and also for emission line
c     luminosities used in QHR-work, or Lsun/Angs), but the optical
c     spectrum is a flux... The user now inputs LumDistInMpc in the GRID
c     FILE, and this is used to define Lobs_norm and then scale things
c     to absolute quantities.

c     ATT: It is the user's responsability to make sure that the
c     luminosities in arq_ETCinfo were computed with the SAME
c     LumDistInMpc informed in the GRID file. Notice that this REVERTS
c     the double-checking policy in FIRcMEx_StarlightChains_v02.for...


c     **> NEW GRID-FILE ENTRIES:
c
c     etc_dir      = dir where to find arq_ETCinfo
c     IsFIRcOn     = 0/1/-1 = no FIRc-fit/FIRc-fit/only *predict* FIR things.
c     IsQHRcOn     = 0/1/-1 = no QHRc-fit/QHRc-fit/only *predict* QHR things.
c     IsPHOcOn     = 0/1/-1 = no PHOc-fit/PHOc-fit/only *predict* PHO things.
c     arq_ETCinfo  = file where all FIR/QHR/PHO info is provided. 
c     LumDistInMpc = same distance used to compute the luminosities in arq_ETCinfo


c     **> NEW CONFIG-FILE ENTRIES: If any FIR and/or QHR and/or PHO
c     stuff is to be used, then read ALL the following, regardless of
c     which of the ETC-features will be used in practice:
c
c     FIRAV_grid parameters (= PHOAV_grid parameters): NFIR_AV,
c     FIRAV_low, FIRAV_upp (= NPHO_AV, PHOAV_low, PHOAV_upp), 
c     QHRbeta_I, FIRbeta_I, FIRbeta_D


c   **> NEW ROUTINES introduced during FIR+QHR+PHO=ETC-implementation:
c
c     SUBROUTINE FIRInitialization(i_verbose,N_base,Nl,l,l_norm,red_law_option,FIRbeta_I,FIRbeta_D)
c     FUNCTION FIR_CalcFIRchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)
c     SUBROUTINE AverageXoverLambda(x,AV_tot,avelamb_x)
c
c     SUBROUTINE QHRInitialization(i_verbose,N_base,Nl,l,l_norm,red_law_option)
c     FUNCTION QHR_CalcQHRchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)
c
c     SUBROUTINE PHOInitialization(PHO_redshift,NPHO_Ys,i_verbose,N_base,NOl_base,Ol_base,l_norm,red_law_option,etc_dir)
c     FUNCTION PHO_CalcPHOchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)
c
c     SUBROUTINE Read_ETC_Info(i_verbose,arq)
c     FUNCTION ETC_CalcETCchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)
c
c     OBS: Apart from input & output changes, several changes were
c     needed in SUBROUTINE MEX_EX0s_Fit_Scubi to store & restore
c     ETC-related arrays.


c     **> FIR/QHR/PHO INPUTS in arq_ETCinfo: See also routine Read_ETCinfo.
c
c     FIR_logLFIR_TOT
c     FIR_LFIRFrac2Model
c     FIR_ErrlogLFIR
c     FIR_RangelogLFIR
c     FIRChi2ScaleFactor
c
c     NQHR_Ys     = Number of line-luminosities (Y) used for QHRc-fits
c     QHR_GlobalChi2ScaleFactor = to re-scale each QHR_Chi2ScaleFactor(iY)
c     ELR-info: IsELROn, ELR_lambdaA & B, ELR_ErrlogR, ELR_RangelogR, ELR_Chi2ScaleFactor.
c         These info are used inside QHR-routines to predict L_A / L_B and
c         compare it to the observed value.
c     QHR_lambda(iY)      = line wavelength, in Angs
c     QHR_frecomb(iY)     = # of recombination/line photons from Osterbrock-like theory
c     QHR_logY_TOT(iY)    = observed luminosity [in Lsun]
c     QHR_YFrac2Model(iY) = fration of Y_TOT which should be considered in fit
c     QHR_ErrlogY(iY)     = Error in logY
c     QHR_RangelogY(iY)   = Range in logY - defines flat-top core in likelyhood
c     QHR_Chi2ScaleFactor(iY) = weight wrt to optical fit
c     ...where iY=1,NQHR_Ys
c


c     PHO_redshift = redshift to be applied to a model spectrum such
c                    that the model-filter-luminosities can be compared
c                    with the observed ones (avoids need of K-correcting
c                    data).
c     NPHO_Ys      = # of photoetric data points to be used/considered.
c     PHO_name(iY)        = a char*10 string to identify the filter
c     PHO_Filter_file(iY) = the transmissin curve file
c     PHO_logY_TOT(iY)    = total luminosity in the filter [Lsun/Angs]
c     PHO_YFrac2Model(iY) = fraction of Y_TOT to be modeled
c     PHO_ErrlogY(iY)     = error in logY [dex]
c     PHO_RangelogY(iY)   = range in logY [dex]
c     PHO_Chi2ScaleFactor(iY) = weigth wrt to optical chi2
c     ...where iY=1,NPHO_Ys



c     **> GENERAL NOTES ON QHR-WORK! QHRc-fits denote fits which use
c     QH-Related constraints, QH being the rate of ionizing photons
c     predicted by a model. These fits use the observed luminosities of
c     NQHR_Ys recombination emission lines (e.g., Ha, Hb...) in
c     conjunction with the QH measured in the base spectra to predict
c     emission line luminosities and (hopefully) better constrain the
c     spectral fit. *QHR* marks changes in the code related to
c     QHRc-modeling.
c
c     See the abstracts to SUBROUTINE QHRInitialization and FUNCTION
c     QHR_CalcQHRchi2 for details on the modeling of QHR-things.

c     There is still quite a lot of testing to be done.... In any case,
c     I expect that by imposing consistency with LHa & LHb, and FIR-data
c     too, I can get more consistent/reliable resuts for MEx fits...
c
c     NOTE: All QHRelated things rely entirely in the input-base models
c     for lambda <= 912 Angs. If you do not trust these models, do not
c     use this code, or else, fiddle the <= Angs base spectra to make
c     them be more of your likening... This fiddling is enteriely left
c     for the user...

c     Within QHRc-fits, we can now fit an EmissionLineRatio (ELR)! This
c     is very cool, since we may not trust the, say, Hb & Ha
c     luminosities (due to claibration, apertures & other effects), but
c     we my still want the model to predict/fit the correct/observed
c     Ha/Hb! In practice, since the yous pops are the ones whichproduce
c     ionizing photos, fitting a given Ha/Hb is ~ equivalent of fixing
c     the AV_tot of these you populations! I always wanted to do this,
c     and it is now working :-)

c     **> GENERAL NOTES ON PHO-WORK! !AKI! I just wrote the code, but
c     had not tested it seriously yet! I'll wait for tests before
c     writing something here...



c     * MISCELANEOUS NOTES/TO-DO's: 

c     > While developing QHR, I have fixed a small number of relatively
c     minor bugs related to FIR-stuff in FIRcMEx_StarlightChains_v02.for.
c     The bugs were minor, but they affected all previous runs ...

c!!!> ATT: Probably all my previous FIRc-runs [with
c     FIRcMEx_StarlightChains_v02.for, in Jan-Feb/2010]
c     stored **incorrect** values of FIRModlogLFIR because I did not
c     realize that calling chi2-routines, which I do in the
c     output-phase, change FIR-related things...

c     Renamed FIR_x4kin to ETC_x4kin & same for AV &
c     common/FIR_xAndAV4Kinematics/... a harmless cosmetic change.

c     I could/should rename FIRlogLobs_norm to FIR_logLobs_norm [and
c     similarly for many other FIR-variales, like FIRq_norm ==>
c     FIR_q_norm, FIRChi2ScaleFactor ==> FIR_Chi2ScaleFactor] for
c     cosmetic elegance / coherent notation...

c     Removed all refs to FIQ_LumDist & QHR_LumDist, as well as double
c     check on LumDistInMpc!

c     Cleaned a few silly comments like "!Cid29Mar07-FIX! Defining base,
c     obs, mask & out directories strings!". Can clean much more after
c     real testing...

c     The iEXOs_IsAve_x_of_lambda does not care about info in other
c     bands in the definition of relevance ... not sure how to fix this,
c     neither if I want to!

c     Is1stLineHeader in arq_config = ugly & obsolete ... (can read
c     headers!)

c     Should completely remove fit-Power-Law option from the code!! It's
c     been turned off anyway!


c   Cid@Lagoa - 11/Jan/2011
*************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     The misterious "Illegal Instruction" Bug and how I fixed it ...   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     While implementing PHO-things, suddenly the code would compile
!     well but, when running it, this is what happend:
!
!     cid@MacBook-de-Roberto-Cid-Fernandes:~/IAA/LIRGS$ ./novo5.exe < qhtest.in5
!     Illegal instruction
!
!     It crashed immediately. This happened while I was reseting
!     PHO_3Mat(j,iY,i_AV) in SUBROUTINE MEX_EX0s_Fit_Scubi. It had to be
!     something about gfortran... Before this I had already written
!     notes to myslef saying that the FIRAV_grid array was stupidly
!     overdimensioned at Nl_max (= 14000) just for lazyness. The same
!     thing propagates to FIR (j,i_AV) matrices FIRR_opt, FIRR_Lya
!     FIRR_LCE & FIRRmat, to the new PHOAV_grid and
!     PHO_3Mat(j,iY,i_IAV). For the default values of Nb_max,
!     NPHO_Ys_max and Nl_max, the size of PHO_3Mat was 300 x 100 x 14000
!     = 4.2e8! No wonder it caused gfortran to crash...
!
!     The solution was to do what I should have done earlier, ie,
!     rationalize the size of AV-grid related arrays. I've introduced a
!     new parameter NAVgrid_max = 1001, and used it in the declarations
!     of:
!
!     FIRAV_grid(NAVgrid_max) , FIRR_opt(Nb_max,NAVgrid_max) ,
!     FIRR_Lya(Nb_max,NAVgrid_max) , FIRR_LCE(Nb_max,NAVgrid_max) ,
!     FIRRmat(Nb_max,NAVgrid_max)
!     PHOAV_grid(NAVgrid_max) ,PHO_3Mat(Nb_max,NPHO_Ys_max,NAVgrid_max)
!
!     and the bug was gone... These changes were marked with !IIB!.
!
!     Cid@Lagoa - 11/Jan/2011
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



*************************************************************************
***                  FIRcMEx_StarlightChains_v02.for                  ***
*************************************************************************

!     *ABSTRACT: STARLIGHT + Multiple Extinctions + FIR-constraints!

!     This version uses a simple (energy-balance/bolometric) dust model
!     to compute the luminosity expected in the FIR due to dust
!     reprocessing of extincted stellar radiation. This extinction- and
!     pop-vector-dependent prediction is (optionally) compared to the
!     observed FIR luminosity (eg, IRAS) to constrain the mixture of
!     populations & exctinctions, ie, the solution! Originally developed
!     for the U/LIRG work, this version may be useful for other projects
!     too. In particular, the changes introduced to account for the
!     FIR-constraints (alas "FIRc") are very similar to the ones I'll
!     need to introduce to include multi-band photometric constraints!

!     The dust reprocessing model is described in routine
!     FIRInitialization, and the computation and scaling of the FIR-chi2
!     is explained in the function FIR_CalcFIRchi2 (which allows for
!     both classical fitting and range/upper/lower-limit fittings).
!     FIR-related routines were basically added to this version of STARLIGHT!
!     Nonetheless, a number of adjustments had to be made in other parts
!     of the code. The changes are all marked with "*FIR*" and
!     commented. A new AverageXoverLambda(x,AV_tot,avelamb_x) routine
!     was also included (see below).

!     This version is supposed to run as the old one when no FIRc-fit is
!     intended, though this has not really been tested thoroughly as
!     yet...

!     There are IMPORTANT and MAJOR changes in the input grid file, as
!     well as minor aditions to arq_config! These were implemented
!     because they are needed in FIRc-fits, but they are useful **in
!     general**. For instance, the user must now provide the physical
!     UNITS of the input f_obs(l_obs) spectrum, and the DISTANCE! With
!     these changes the output masses are now in units M_sun! Things
!     like the bolometric lum of the model (in L_sun) and the predicted dust
!     reprocessed luminosity (also in L_sun) can also be computed. If no
!     distance is known (or if distance is irrelevant), the user needs
!     to cheat STARLIGHT by providing some entry for it;-)

!     An extra option to deal with how-to-EX0 the fit was added. The
!     user can now chose to decide which elements of the base are
!     excluded on the basis of an Lambda-Averaged pop-vector (LAx). I
!     found this useful when working with U/LIRGS, where a component may
!     be weak at lambda = l_norm but much stronger (and thus more
!     relevant) at other wavelengths. This is presumably a better option
!     in general! 

!     In addition, the user MUST provide an arq_FIRinfo file, which, as
!     name suggests, provides the FIR information needed for FIRc
!     fits. If no FIRc-fit is intended (and there is an entry in the
!     grid file to state this), just set arq_FIRinfo to nodata.stupid,
!     or something like that.

!     There are also MAJOR changes in the OUTPUT. I tried to change as
!     little as possible, but it was not possible to avoid changes which
!     imply adaptating the programs & macros which read STARLIGHT's
!     output...  The two major changes are a new block where
!     Lambda-Averaged pop-vectors (LAx*) are reported, and at the very
!     end of the file, where the FIR results are outputed.
!
!     Notice that correct modeling of the dust-reprocessing REQUIRES
!     bases which provide spectra from l << 912 Angs to the optical,
!     since ionizing photons and non-ionizing UV light are the major
!     contributors to dust emission in the Far IR... At the time of
!     writing, only BC03/CB07 provide this info!
!
!     Silly bugs in previous versions (like the bug in AV-estimation for
!     SSP fits and the saving-in-the-wrong-dir when fobs_norm = 0) were
!     also fixed (as in StarlightChains_v05.for). Also added the
!     save-best-ssp-spec option from v05.
!
!     OBS: There are important differences wrt
!     FIRcMEx_StarlightChains_v01.for (which was given a final-looking
!     name too early...).


!     * NEW GRID-FILE INPUTS:
!
!     IsFIRcOn     = 0/1/-1 = no FIRc-fit/FIRc-fit/only *predict* FIR lum.
!     flux_unit    = multiplicative factor to conver f_obs to ergs/s/cm2/A
!     LumDistInMpc = luminosity distance in Mpc.
!     arq_FIRinfo  = file where FIR info is provided. 
!
!
!     * arq_FIRinfo INPUTS:
!
!     FIR_logLFIR_TOT    = TOTal log FIR luminosity [Lsun]
!     FIR_LFIRFrac2Model = frac of LFIR_TOT to consider in fit
!     FIR_ErrlogLFIR     = error in log L_FIR [dex]
!     FIR_RangelogLFIR   = Range in log L_FIR [dex]
!     FIRChi2ScaleFactor = multiply chi2_FIR by this and by chi2_OPTical!
!
!
!     * NEW ROUTINES introduced during FIRc-implementation:
!
!     FIRInitialization(N_base,Nl,l,l_norm,..,FIRbeta_I,FIRbeta_D)
!     FIR_CalcFIRchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)
!     AverageXoverLambda(x,AV_tot,avelamb_x)
!
!     * NEW CONFIG-FILE ENTRIES:
!
!     > i_SaveBestSingleCompFit = inherited from v05 & MCclusters work:
!       = 1/0 = Y/N - save best single component fit on output
!     > iEXOs_IsAve_x_of_lambda = 0/1: (0) Make EX0 relevance-decision
!       with x(l_norm) or (1) with <x(lambda)>. The latter is the
!       Lambda-Averaged (weighting by weight_obs) pop-vector.
!     > FIRc technical parameters (needed to model dust reprocessing):
!       FIRbeta_I, FIRbeta_D, NFIR_AV, FIRAV_low & FIRAV_upp


!     * ARRAY-RESETINGS: During the implementation of FIR-constraints I
!     realized that some arrays were not reset to = 0. This causes a
!     harmless population of unused array elements with RAM-garbage. In
!     principle this does NOT affect the code, since the elements that
!     matter are correctly dealt with, but it is annoying and drove me
!     crazy for a few days in January/2010. Because of that I now reset
!     ALL output quantities, including ALL Nb_max elements of x() &
!     exAV(), among many other arrays. The reseting is done even in
!     routines which are called-up all the time, like
!     MEX_Unpack_par2xAVYAVfn, but I have verified that this does not
!     slow down the code significantly at all. These changes are marked
!     with "*FIR*", "Reset" or something similar throughout.


!     * MISCELANEOUS NOTES/TO-DO's: 

!     > Take care NOT to become FIR-dependent, or else protect against
!     problems when the base models do not cover the required range (as
!     BC03 do) for FIR-work...


!     Cid@Granada+LaLaguna - Jan+Feb/2010
*************************************************************************



!     The StarlightChains_v05.for abstract below was added to explain
!     what was introduced/fixed in that version.
!     ------------------------------------------------------------------
!     |                    StarlightChains_v05.for                     |
!     ------------------------------------------------------------------
!
!     !!Found a TRUE & +/- SERIOUS bug in StarlightChains_v04.for:( The
!     bug affects ONLY the single-component (SSP) fits carried out at
!     the end of each run. It MAY lead to wrong SSP-parameters,
!     particularly when they are already bad anyway. The bug was:
!
!         call mnbrak(ax,bx,cx,fa,fb,fc,F_SSP_chi2_AV)
!         chi2_AV_Gold = golden(ax,bx,cxf,F_SSP_chi2_AV,tol,AV_Gold)
!
!     where "cxf", instead of "cx", was passed to golden, leading to a
!     (not always) incorrect estimation of AV_Gold. Cid & Rosa noticed
!     this bug in the AV_SSP x SSP age curves while fiting star cluster
!     spectra. AV_SSP becomes constant at some values of age (but,
!     luckly, not in the range of best SSP fits). This version corrects
!     this bug.
!     
!     This version also reads a new i_SaveBestSingleCompFit flag from
!     the config file. If = 1, the best fit SSP spectrum is saved as a
!     5th column in arq_out, otherwise the output is as in v04.
!
!     Finally, I also take the chance to fix the harmless previously
!     knwon bug that when fobs_norm = 0 arq_out is saved in the running
!     dir, instead or in out_dir. Now such bugged output-files are saved
!     in out_dir.
!
!     All changes are marked with "!v05!".
!
!     Cid@Lagoa - 10/May/2008
!     ------------------------------------------------------------------



*************************************************************************
*                        MEX_StarlightChains_v01.for                    *
*************************************************************************
*
*     This is the 1st "official" version of the Multiple EXtinctions
*     (a.k.a. "MEX") version of Starlight. Based on
*     StarlightChains_v04.for, but with MEX-modifications. Except for
*     output details, the code is exactly = multext_draft9_com_v04.for,
*     which I used for initial tests.
*
*     There are still a "few" things to tidy up, but this will be left
*     for later ... Now I will concentrate on running this version for
*     the VLIRGS project with Rosa. Later on, I should consider making
*     this version public (in which case I should 1st check that it does
*     work exactly as the previous one in the simpler case of a single
*     AV).
*
*     Cid@Lagoa - 30/July/2007
*************************************************************************


*     Miscelaneous comments/notes...

*     retomando... versao 8 tinha um format bug imbecil...  percebi que
*     os PDF nao estavam sendo salvos no output dir e arrumei!
*     cid@lagoa - 21/July/2007
*     
*     Implemented CUMULATIVE AV_tot scheme. Most of what I had to do was
*     to adapt routine MEX_Update_AV_array(N_par,N_base,AV,exAV,AV_tot)
*     and a few loose ends.
*
*     eliminated useless input in arq_config: YAV_low & _upp +
*     PowerLaw-stuff.

*     Adicionei cousas no never-really-used chains-output...
*     devo decidr se vale a pena a ensuing alteracao do output size/format!!

!     tah indo bem, mas o adaptive-steop pros AV & exAV inda parece
!     estranho...  Passos em av parecem too large!! eg:
!     a*eps_AV=8.25E-01 Sao os "ultimos" !!UEPA!!'s antes de usar pra
!     valer com as VLIRGS...ou simuacoes ... ou SDSS galaxies...
!     aparentemente i_UpdateAVYAVStepSeparately = 0 dah um jeito...

*     Cid@Lagoa - 28/May/2007


*     Introduced "revert to slow chi2 in case of horror" scheme!
*     *BUT* It will only be used if NRC_AV_Default < 0 in arq_config!!
*     Ideally I should make it revert-to-slow-chi2 automathicaly when
*     the RC-scheme fails.... but I have NOT done it!
*     Cid@Lagoa - 28/May/2007

c     OBS: me dei conta de que ha uma maneira bem + simples de
c     implementar um EX0s scheme: eh soh definir um array com os indices
c     das componentes relevantes e soh perturbar elas! :(
*     ======================================================================


*     !!MEX!! Peguei o StarlightChains_v04.for e vou comecar a meter as
*     multiple extinctions nele...  Marcacoes com !!MEX!!  (Cid@Lagoa -
*     19/April/2007)

*     Repensar o NRC_AV ... usar >= default alto, tipo 10 ou 12?? Ou
*     aumentar o dchi2_threshold ???
* fitspec  2 ??MEX?? ATT: Increasing # of terms in AV-series...   1.67569E+02   1.40381E+02   1.93675E-01   1.00000E-03    0.970     0.970   -99.000   -99.000 11 ==> 12

*     Cid@Lagoa - 22/April/2007




!     ------------------------------------------------------------------
!     |                    StarlightChains_v04.for                     |
!     ------------------------------------------------------------------

!     v04 = 1st public version of STARLIGHT!

!     This new version is still similar to StarlightChains_v02.for, but
!     with a large suite of new features, including some which imply
!     changes in the grid file info, base-file (arq_base), config-file
!     (arq_config) and output file (arq_out) too!
!
!     The ~ minor changes are (some imply changes in the grid, config,
!     base & arq_out files):
!
!     (01) Fixed the "infinite loop" bug plus the silly i<->j mistake in
!     arrays-reseting (found by William). Also, changed Nl_max from 7000
!     to 14000. (Paula's base has ~ 10400 pixels!).
!     
!     (02) Added some time() & etime() commands to figure how long a run
!     takes (useful to test configs).
!
!     (03) Added protections against non-increasing l_obs() & l_base(),
!     Mstars() <= 0 & idum_orig > 0.
!
!     (04) If idum_orig = 0 in the grid.in file, then idum_orig will be
!     kept the same for all runs in the grid file. Also added protection
!     against idum_orig > 0.
!
!     (05) Added a nice & simple feature: If the output file arq_out
!     already exists, then we SKIP THE FIT altogether and go to the next
!     one on the list! Very useful to run grids in "modo-burro"! At
!     William's request, I added a Y/N flag in arq_config which decides
!     whether this trick is to be used or not.
!
!     (06) In preparation for alpha enhanced bases, I have increased to
!     character*15 the size of the component nickname (base_name), and
!     introduced a alphaFe_base() arrays which MUST be read from
!     arq_base. alphaFe_base() is echoed in the output, along with x &
!     etc.
!
!     (07) Mask base edges! The first & last dl_cushion Angs in the base
!     spectra will be masked to prevent edge effects in the kinematical
!     filter! 
!
!     (08) dl_cushion used to be fixed at 50 Angs in previous
!     versions. Now we read it from arq_config! (whose format thus
!     changed!)
!
!     (09) After STARLIGHT is done (ie, just before output) we now do
!     "Single-Component-Fits", where each of the base spectra is
!     compares to the data individually. We fit AV and the
!     norm. factor. The results are stored along with the main fit pop
!     vector (see SSP_* arrays). This may be used to define a "best
!     single component" representing a spectrum (potentially usefull for
!     star-clusters and elliptical galaxies.)
!
!     (10) Now reading fEX0_MinBaseSize from arq_config. This is used
!     to define N_base_min = nint(fEX0_MinBaseSize * N_base) in EX0-fit.
!
!     (11) Removed unnecessary common/YAV_stuff/ in ConvertLight2Mass.
!
!     (12) Now the EX0s routine will NOT be called when NEX0s_loop = 0.
!     Part of the output was computed & printed with EX0s* parameters
!     (like x, av, v0, vd...), so I had to fix that to allow for the
!     NEX0s_loop = 0 case. (changes marked with "!Cid24Mar07-FIX!")
!
!     (13) Added "censored weights" (Weight-control-filter-stuff) info
!     in the output file just for control-purposes.
!
!     (14) Added a scheme to count number of columns in arq_obs to check
!     consistency of IsErrSpecAvailable & IsFlagSpecAvailable...
!
!     (15) Added protection againts using i_FastBC03_FLAG = 1 and non
!     BC03/STELIB models. If this happens, reset i_FastBC03_FLAG to 0.
!
!     (16) The grid file now hase obs_dir, mask_dir & out_dir, which are
!     dire directorieds whre arq_obs & arq_masks will be read from and
!     arq_out will be written to.
!     
!     (17) Format of l_obs on arq_out changed from i4 to f8.2 to cope
!     with IR data and high-res too.
!
!
!     The ~ major changes (implying changes in the grid file) are:
!
!     (01) grid-file "header" info is now vertical! Ie, the previoulsy
!     horribly long 3rd line input is now broken into one thing per
!     line. Similarly, the config file has been restructured and
!     reordered for clarity! 
!
!     (02) Introduced Gordon et al (2003) SMC & LMC extinction curves
!     (red_law_option = GD1, GD2 & GD3). They are still interpolated
!     from a table, but better and different from HyperZ-curves.
!     Also, introduced a general "from-a-file" reddening law read from a
!     file named ExtinctionLaw.red_law_option! 
!
!     (03) Now the initial estimates for v0 & vd are read FOR EACH
!     GALAXY from the grid file! The previous version all used v0 = 0 &
!     vd = 150 km/s as a starting point. This implied a change to the
!     format of grid files...
!
!     (04) Introduced a new "FXK" option for the FitOrSamplePDF_option
!     string, which keeps vd & v0 FIXED throught the fits. In other
!     words, kinematics is NOT fitted at all, and v0 & vd are kept = to
!     their initial values (now read from the grid file)! The base is
!     shifted by v0 & vd before any call to "Scubidu" routines and never
!     shifted nor smoothed again. This may be useful in a number of
!     cases. This was the major change. Search for "FXK" to find where
!     the changes were made.
!
!     (05) While writing the manual and checking the code I found an
!     invisible bug! In the Refit-after-clipping and Burn-In calls to
!     Fit_Scubidu the variable i_RestartChains was UNDEFINED! That means
!     it was 0, and neither in the Refit nor Burn-In stages the chains
!     were reset to the best parameters... I now set i_RestartChains =
!     i_RestartChains_FF in both calls. This was most lilkey harmless,
!     though (in general) we should now get slightly better fits (if
!     i_RestartChains_FF = 1, of course).
!
!     Later (01April07) I found another related bug! The new bug is that
!     by setting the last values equal to the best ones (in Fit_Scubi
!     with i_RestartChains = 1), if one of these best chain values = the
!     global best so far then the trick of resetting the chain with
!     worst energy = to the best global model may result in this
!     chain-values appearing TWICE!! Ie, 2 of the chains would start
!     from IDENTICAL positions.... To see how I fixed this search for
!     "Cid01April07-FIX".  I also introduce a new option:
!     i_RestartChains = 2, which will reset ONLY the worst chain to the
!     best model, but leave the others wherever they are.
!
!
!     Cid@Lagoa - 07/February/2007 --> 11/April/2007

!

!     **Disorganized list of things to improve**

!     - Better/analytical reddening laws for SMC & LMC! & fix minor
!     "step" in CAL law. 

!     - The i_FastBC03_FLAG for Fast-rebinning should be changed to cope
!     with other bases...

!     - Write robust resampling routine...

!     -Option to deal with rectified spectra??? Dificult to see how this
!     can be done...

!     -COMO FIXAR AV E NAO FITA-LO?? ISTO SERAH OUTRO PROG, COM
!     RETIFICACAO E O ESCAMBAU...  UMA possibilidade eh inventar
!     red_law_option = 'NONE' e setar red_term(lambda) = 0 ou 1 ...
!     desse modo o av/yav variaria mas sem efeito pratico (?) e talvez
!     as chains convirjam anyway...  neste caso teriamos que mudar
! 100     AV_test  = 0.3
!         YAV_test = 0.3
!     onde se testa os chi2 slow x fast... que acontece se setar NRC_AV=0?
!     ------------------------------------------------------------------


!     ------------------------------------------------------------------
!     |                      StarlightChains_v02.for                   |
!     ------------------------------------------------------------------
!     This version fixes a few minor things in StarlightChains_v01.for.
!     Changes are essentially cosmetic or safety-measures.
!
!     (1) No 'Running model # ...' screen output anymore.
!
!     (2) Introduced Is1stLineHeader flag in config.
!
!     (3) If fobs_norm <= 0 the run is aborted.
!
!     (4) Added safety-check in the computation of S/N's in the
!     llow_SN--lupp_SN window...
!
!     (5) Changed a bit the way we constrain initial chain fn, AV & YAV
!     values to allow them to be constrained by the low->upp values in
!     the config file.
!
!     (6) Introduced IsScaleNstepsInEX0sFit flag in config. If turned
!     on, the numbers N_sim & Nmax_steps, which control the length of
!     Metropolis runs in Scubidu, will be scaled according to the actual
!     number of relevant parameters in the EX0s fits. This ONLY changes
!     things in the EX0s fits.

!     (7) Introduced a IsNoKin4LargeBaseInEX0sFits flag & a
!     frac_NoKin4LargeBaseInEX0sFits variable, which, multiplied by
!     N_base, defines what is the minimum size of a EX0s-condensed base
!     for which we turn-off kinematical fits. (In v01 this was fixed at
!     50 elements). This too ONLY changes things in the EX0s fits.
!
!     (8) Introduced the 5 reddening laws of HYPERZ ('HZ1' ... 'HZ5'),
!     of which HZ4 is the SMC (which I'd expect to be most relevant)
!
!     Only (6) & (7) may have some relevant effect wrt v01, particular
!     (6). Setting IsScaleNstepsInEX0sFit = 0 should yield nearly
!     identical results, but setting it to 1 will result in shorter
!     (faster) EX0s-fits, which may be slightly worse in chi2.
!
!     Cid@Lagoa - 09/October/2005
!     ------------------------------------------------------------------



!     This is my "1st-final" version of "ScubiStar" (built upon
!     ScubiStar7.for version produced at IAA in May/2005). The main
!     difference wrt previous versions of Starliht is that this version
!     uses N_chains parallel Markov chains, which run untill a Gelman &
!     Rubin-like Burn-In convergence criterion is satisfied. Full
!     documentation is not ready yet...
!
!     A few other cosmetic changes were made. Base files are now stored
!     in a user-defined directory (read from the input file), avoiding
!     multiple copies of base-spectra throughout the disk.  Another
!     change is that now the input file specifies two flags
!     (IsErrSpecAvailable & IsFlagSpecAvailable) which define whether
!     error & flag-spectra are read from the input file (as for SDSS
!     data) or are computed online with the old rms-in-a-window
!     method. Yet another change is that now the base, mask and config
!     files & extinction curve-option are specified for each galaxy!
!     These changes (and others) were made in the sense of having a
!     single version of the code useful for different purposes (eg, SDSS
!     and my LLAGN data). More details will be added later.
!
!     I have not run a full set of tests/simulations of this code...
!     However, quick & dirty tests of different configurations show
!     that results are very similar for all configurations I
!     tried. Hence, I'll distribute the code and see how it
!     performs. Changes will be made "on-the-fly", ie, as problems appear.
!
!     Sent to Abilio on 30/May/2005
!
!     Cid@IAA - 30/May/2005

c     voltei a mexer no codigo - ... abstracts-totalmente desatualizados....
c     Cid@IAA - 3/May/2005
c
c     Tentando botar Scubidu no Starlight05_Fast_v01.for...
c     melhor seguir testes sem YAV e depois resolver este bug novo...
c     Cid - 17-2-05


c     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c     This is my latest update of the Starlight spectral synthesis code.

c     BRIEF HISTORY of Starlight: So far there have been 2 main
c     versions.

c     (1) Starligh03_v02.for (September/2003) was the original code,
c     used in the Seyfert 2 paper (CF, Gu, Melnick ... et al 2004;
c     MNRAS). The same code was used in my LLAGN paper III (CF, Gonzalez
c     Delgado ... et al 2004; MNRAS), with a base of template
c     galaxies. Initial experiments with SDSS spectra were also
c     performed using a slightly modified version of the code
c     (Starlight03_v04SDSS.for; which dealt with errors and flags).

c     (2) In Starlight03_v04SDSS_Fast2.for I've introduce the trick of
c     expanding the extinction term in a truncated series (the
c     "RapidChi2" or RC-trick), which speeds things up by > 200 times!
c     That was used in SDSS-Paper-I (CF, Mateus ... et al MNRAS 2005). A
c     slightly modified version called Starlight03_v04_ind_Fast2.for was
c     produced to deal with individual spectra (which require galaxy by
c     galaxy input info); it also introduced the Calzetti reddening-law
c     option.

c     NEW FEATURES OF THIS VERSION: The main novelty in this 3rd version
c     is that we now allow for 2 different extinctions: AV & YAV. Each
c     population j (j=1...N_base) is reddened by AV_tot(j) = AV + YAV,
c     and only some of them are allowed to have the extra extinction
c     YAV. The array ind_YAV(j) says who can have this extra
c     extinction. The idea is to do something like a Charlot & Fall
c     (2000) extinction law, which has one (higher) AV for young
c     populations (<= 10 Myr) and a global AV affecting both young and
c     old stars. However, the SAME reddening law is assumed for both
c     populations, otherwise heavy changes would be required in the
c     RC-chi2 routines... The code was also re-organized, with
c     improvement structure and comments. Finally, I now read the base
c     DIRECTORY (string base_dir) from the input/grid file (2nd
c     line). This avoids duplication of base spectral files. Use "./"
c     for base_dir if the base files are in the current directory!
c     Note that symbolic links work:)

c     Other changes: Allowance for returned-mass in the light-to-mass
c     fractions conversion was introduced. The conjugate-gradient
c     routines were removed altogether (they never proved worthy). Both
c     config & base files now have different formats. The output file
c     too is now slightly different (so auxiliary table-making &
c     plotting codes will have to be updated accordingly).

c     FUTURE CHANGES: There is a major change in progress to make the
c     whole algorithm more efficient, using the estimates of the
c     covariance matrix to identify principal directions for movement in
c     the parameter space. The burn-in phase should also be mapped, and
c     full error estimates of the parameters will be done. Jean, William
c     & myself are working on it.

c     The code is heavily commented, but not yet fully comprehenssible
c     by people other than myself! Full documentation & tips will be
c     provided elsewhere, one day...
c
c     Cid@Lagoa - 23/January/2005
c     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



c ###########################################################################
c #                         DRIVER FOR MAIN CODE                            #
c ###########################################################################
      program Starlight_With_MEx_And_FIRc_And_QHRc

      implicit real (a-h,j-m,o-z)
      implicit integer (i,n)
      character*3 FitOrSamplePDF_option
      character*100 base_dir , obs_dir , mask_dir , out_dir , etc_dir

c     1st line in input/grid-file is the number of fits to be computed
      read (*,*) N_models

c     Input some global grid-parameters here & pass them on to main code
      read (*,'(a)') base_dir
      read (*,'(a)') obs_dir
      read (*,'(a)') mask_dir
      read (*,'(a)') etc_dir
      read (*,'(a)') out_dir
      read (*,*) idum_orig
!NAT-SN-BUG! Renamed ll*_SN to Oll_*SN - ElCid@Sanchica - 15/Ago/2011
      read (*,*) Ollow_SN
      read (*,*) Olupp_SN
      read (*,*) Olsyn_ini
      read (*,*) Olsyn_fin
      read (*,*) Odl_syn
      read (*,*) fscale_chi2
      read (*,*) FitOrSamplePDF_option
      read (*,*) IsErrSpecAvailable
      read (*,*) IsFlagSpecAvailable

*PHO* 1/0/-1 = Y/N/Predict-PHO Flag to turn PHOc On/Off/Predict-only. 
      read (*,*) IsPHOcOn
         write (*,'(70(a))') ('*',i=1,38)
         if (IsPHOcOn.EQ.1) then 
            write (*,'(a,a)')  '*  ATT:  You asked for PHOc fits!!!  *'
         else if (IsPHOcOn.EQ.0) then 
            write (*,'(a,a)')  '* ATT: You asked for NON-PHOc fits!! *'
         else if (IsPHOcOn.EQ.-1) then 
            write (*,'(a,a)')  '* ATT: PHO will just be predicted!   *'
         else
            write (*,'(a,i4,a)')
     $           '* ATT: IsPHOcOn should be -1/0/1, but you chose : ',
     $           IsPHOcOn , ' ... :-('
            stop
         end if
         write (*,'(70(a))') ('*',i=1,38)

*QHR* 1/0/-1 = Y/N/Predict-QHR Flag to turn QHRc On/Off/Predict-only. 
      read (*,*) IsQHRcOn
         write (*,'(70(a))') ('*',i=1,38)
         if (IsQHRcOn.EQ.1) then 
            write (*,'(a,a)')  '*  ATT:  You asked for QHRc fits!!!  *'
         else if (IsQHRcOn.EQ.0) then 
            write (*,'(a,a)')  '* ATT: You asked for NON-QHRc fits!! *'
         else if (IsQHRcOn.EQ.-1) then 
            write (*,'(a,a)')  '* ATT: QHR will just be predicted!   *'
         else
            write (*,'(a,i4,a)')
     $           '* ATT: IsQHRcOn should be -1/0/1, but you chose : ',
     $           IsQHRcOn , ' ... :-('
            stop
         end if
         write (*,'(70(a))') ('*',i=1,38)

*FIR* 1/0/-1 = Y/N/Predict-FIR Flag to turn FIRc On/Off/Predict-only. 
*     With -1 FIR-model-related-things are initialized but no actual
*     FIR-data is read. The purpose is just to predict the FIR
*     luminosity implied by the optical STARLIGHT-fit.
*     ATT: This flag defines what is read the ETC & config files!!
      read (*,*) IsFIRcOn
         write (*,'(70(a))') ('*',i=1,38)
         if (IsFIRcOn.EQ.1) then 
            write (*,'(a,a)')  '*  ATT:  You asked for FIRc fits!!  *'
         else if (IsFIRcOn.EQ.0) then 
            write (*,'(a,a)')  '* ATT: You asked for NON-FIRc fits!! *'
         else if (IsFIRcOn.EQ.-1) then 
            write (*,'(a,a)')  '* ATT: L_FIR will just be predicted! *'
         else
            write (*,'(a,i4,a)')
     $           '* ATT: IsFIRcOn should be -1/0/1, but you chose : ',
     $           IsFIRcOn , ' ... :-('
            stop
         end if
         write (*,'(70(a))') ('*',i=1,38)

*FIR* flux_unit is the factor by which one must multiply the input
*     observed spectrum to obtain it in ergs/s/cm2/Angs. Knowing the
*     physical units is absolutely ESSENTIAL for FIRc work (since fluxes
*     ate two bands will be analysed simultaneously), but ABSOLUTELY
*     IRRELEVANT for optical-fits alone! 
*     OBS: flux_unit is unnecessary for the U/LIRGS project with Rosa
*     (the spectra are already in the right units), but may be useful
*     later when I apply this code to other data (SDSS, for instance,
*     comes in units of flux_unit = 1e-17)
      read (*,*) flux_unit

!Cid14Jan07-FIX! idum-trick: Sometimes I may want to repeat the run with 
!     the same idum, but different config, or base, or red-law. To allow
!     for this possibility using a single grid-file, just set idum = 0
!     in the grid.in file and all runs will have the same idum:)
!     Also, added a protection against idum_orig > 0.
      i_ChangeIdum = 1
      if (idum_orig.EQ.0) then
         i_ChangeIdum = 0
         idum_orig = -123456
      end if
      if (idum_orig.GT.0) idum_orig = -idum_orig

c     Safety-check on some of the input info...
      if ((FitOrSamplePDF_option.NE.'PDF').AND.
     $     (FitOrSamplePDF_option.NE.'FIT').AND.
     $     (FitOrSamplePDF_option.NE.'FXK')) then
         write (*,'(a)') ' Problem with your grid file! '
         write (*,'(a,a,a)') ' What are you doing? ' ,
     $        FitOrSamplePDF_option ,' should be = FIT or FXK ...'
         stop
      end if 

      if ((IsErrSpecAvailable.NE.0).AND.(IsErrSpecAvailable.NE.1)) then
         write (*,'(a)') ' Problem with your grid file! '
         write (*,'(a,a,i3)') ' What are you doing? ' ,
     $        ' IsErrSpecAvailable = ' , IsErrSpecAvailable
         stop
      end if 

      if ((IsFlagSpecAvailable.NE.0).AND.(IsFlagSpecAvailable.NE.1))then
         write (*,'(a)') ' Problem with your grid file! '
         write (*,'(a,a,i3)') ' What are you doing? ' ,
     $        ' IsFlagSpecAvailable = ' , IsFlagSpecAvailable
         stop
      end if 

      if ((IsErrSpecAvailable.EQ.0).AND.(IsFlagSpecAvailable.EQ.1)) then
         write (*,'(a)') ' Problem with your grid file! '
         write (*,'(a,a,i3,a,i3)') ' What are you doing? ' ,
     $        ' IsFlagSpecAvailable = ' , IsFlagSpecAvailable ,
     $        ' but IsErrSpecAvailable = ' , IsErrSpecAvailable
         stop
      end if 


c     Run grid of N_models spectral fits
      do i_model=1,N_models

!NAT-SN-BUG! Reset ll*_SN = Oll_*SN - ElCid@Sanchica - 15/Ago/2011
         llow_SN = Ollow_SN
         lupp_SN = Olupp_SN

*ETC* Added IsFIRcOn, IsQHRcOn, IsPHOcOn flags, flux_unit & etc_dir.
         call FitSpectrum(idum_orig,llow_SN,lupp_SN,Olsyn_ini,Olsyn_fin
     $        ,Odl_syn,fscale_chi2
     $        ,base_dir,obs_dir,mask_dir,etc_dir,out_dir
     $        ,FitOrSamplePDF_option,IsErrSpecAvailable
     $        ,IsFlagSpecAvailable,IsFIRcOn,IsQHRcOn,IsPHOcOn,flux_unit)

c     changes idum-seed just to avoid doing always the same sequence of
c     random numbers. Also, idum-trick: only changes if i_ChangeIdum.NE.0!!
         if (i_ChangeIdum.NE.0) idum_orig = idum_orig - 1

      end do


      end
c ###########################################################################



c @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c @                                MAIN CODE                                @
c @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     ATT: This is an incomplete & rather disorganized ABSTRACT of the
!     code. Full documentation will be produced one day... (yeah, right...)

c     This code synthesizes a galaxy spectrum f_obs(l_obs) with a base
c     of N_base spectra f_base(l_base). The base can be made up of
c     observed template galaxies or theoretical spectra (eg, BC03,
c     Starburst99, etc). The ouput consists of the population vector
c     x(i); i=1...N_base, plus the V-band global extinction AV and the
c     extra-extinction YAV to selected components, plus kinematical
c     parameters v0 (residual velocity-shift) & vd (velocity
c     dispersion).

c     The code seeks to minimize the chi2:
c
c     chi2(x,AV,YAV,v0,vd) = SUM [(f_syn - f_obs) * weight_obs]**2
c
c     where the sum is over all lambdas, and f_syn =
c     f_syn(x,AV,YAV,v0,vd) is the model spectrum for the given
c     parameters:
c
c     f_syn(i) = Sum_j x(j) * f_base(j,i) * 10.**(AV_tot(j) * red_term(i))
c                CONVOLVED with G(v0,vd)
c
c     where f_base(j,i) is the base spectrum of component # j at lambda
c     l_obs(i), convolved with the Gaussian kinematical filter G(v0,vd).
c
c     Masked lambdas are given weight_obs = 0. For other lambdas,
c     weight_obs = 1 / error_spectrum, or mask_weight times larger if
c     specified on arq_masks.

c     * NORMALIZATION: The base spectra f_base are normalized by their
c     flux at l_norm, while the observed spectrum is normalized by
c     fobs_norm = the median flux between llow_norm & lupp_norm (which
c     must bracket l_norm).  We do NOT force the population vector to be
c     normalized to 100%! The returned value of x_j is the flux in the
c     model spectrum (without convolution by G!) in component j, in unit of
c     fobs_norm.

c     * ERROR & FLAG SPECTRA: The input spectrum MUST come with an
c     associated error spectrum (err_obs) and a flag-spectrum
c     (flag_obs). The coding for the flag spectrum is like this: flag =
c     0 or 1 is OK, while flag > 1 will ignore the points from the fit
c     (like the masks described below). Flags are useful to identify bad
c     pixels. If you do NOT have an error spectrum, then invent one, but
c     beware that the chi2's computed here depend on this error.

c     * MASKS: The synthesis is carried out in all Nl_obs observed
c     wavelengths l_obs, but some are "masked out". The user-supplied
c     file arq_masks indicates which intervals should be excluded from
c     the modeling (eg, those containing emission lines or bad
c     pixels). We further impose a limit f_obs >= f_cut (in units of the
c     normalization flux). Weaker fluxes are masked out assignig
c     weight_obs = 0 for f_obs <= f_cut. (Just set f_cut = 0 if you do
c     not want this safety feature.)  The file arq_masks can also
c     specify windows over which we want to give larger weight (for
c     instance, on key absorption features). The array mask_weight()
c     decides by how much the weight should be multiplied in a given
c     window. (Example: mask_weight = 0 ignores the region in the fit
c     and mask_weight = 2 gives is twice as much weight as it would have
c     otherwise.)

c     * SAMPLING: The input spectrum f_obs(l_obs) is resampled/rebinned in
c     wavelength to lsyn_ini ==> lsyn_fin, in steps of dl_syn
c     (subroutine RebinSpec). This is the range which will actually be
c     modeled (but note that we internaly extend the desired fitting
c     range by 50 Angs on the blue & red edges to circumvent problems in
c     the kinematical filter). The base spectra f_base are resampled to
c     the same wavelengths (ie, l_obs = lsyn_ini ==> lsyn_fin in dl_syn
c     steps). Note that vd depends on dl_syn, as well as on the actual
c     spectral resolution of the input spectra... To correct vd for the
c     different resolutions of observed & base spectra do a quadrature
c     correction. For instance, if your spectra have sigma_inst = 62
c     km/s and the base has 86 km/s resolution do:

c     vd_corrected = sqrt(vd_given_by_the_code**2 + 86.**2 - 62.**2)

c     ==> It is VERY STRONGLY recommended that you provide both observed
c     & base spectra already resampled to 1 Angs! (This makes things
c     easier.) One reason for this is that is that my Rebinning routine
c     does NOT handle the error spectrum properly. It just resamples
c     err_obs, whereas it should actually correct for the fact that
c     going from, say, 1 to 2 Angs sampling the error should be smaller...


c     * NUMERICAL-SCHEME: I've tried several schemes to estimate the
c     stellar population (x, AV & YAV) and kinematical parameters (v0 &
c     vd).
c
c     1 - Metropolis in the (x,AV,YAV,v0,vd) space - (slow method)
c     2 - Metropolis on x,AV followed by Metropolis on (x,AV,YAV,v0,vd).
c     3 - Annealing schemes
c     4 - etc (in previous versions I tried conjugate gradients, dicrete
c         grids, combinations of these,...)
c
c     I have settled for an adapted version of stimulated annealing,
c     which runs a series of Metropolis runs to find the best (x,AV,YAV)
c     for fixed (v0,vd) and, at each step, re-estimates (v0,vd) for
c     fixed (x,AV,YAV). [Also, after each Metropolis run (in the
c     annealing loop) I do a quick search for better AV & YAV for fixed
c     x, v0 & vd.]
c
c     At each step in the annealing routine (Fit_xAVYAVv0vd) I multiply
c     the weights (= 1 / error) by a factor wf which runs from wf_low to
c     wf_upp. The step size (in x & AV) is concomitantly reduced from
c     eps_upp to eps_low, and the number of steps is adjusted so that a
c     long Metropolis trip can be made. The last loop is the slowest
c     one, with more points. The maximum number of steps at any given
c     annealing interaction is limited by Nmet_max. If you give a very
c     large number for this parameter then all steps that "should" be
c     made are actualy made, but it may take a long time to run...


c     [CLIPPIG & EX0s & etc] The code allows other tricks: 

c     (A) Sigma-Clipping can be made, removing all points with more
c     than, say, sig_clip_threshold = 3 sigma residuals, and re-doing
c     the fit. This can be done Nclip_loops_max times.

c     (B) Quick-Clipping: If Nclip_loops_max < 0 on input, then a "fast"
c     clipping routine is performed (only once) before going for the
c     actual fit. This works prettry well and saves time, so we
c     recommend it.

c     (C) Condensed-Base-Fitting: Setting NEX0s_llops >= 1 on arq_config
c     will re-fit the spectrum, taking as initial point the values
c     obtained in the main call to Fit_xAVYAVv0vd, and reducing the base
c     just to the N_non_zeros components which have x_j > 0. This speeds
c     up computations (because of the smaller base), and serves to
c     refine the fit.

c     (D) The rebinning of BC03 spectra is now much faster (and less
c     general, since we assume dlambda = 1). Use i_FastBC03_FLAG = 1 in
c     arq_config.

c     (E) WEIGHT-CONTROL-FILTER: I have implemented a "safety-filter" which
c     prevents any given point from having a weight more than 2 sigma
c     above the mean weight of unmasked points. I did this after
c     realizing that some SDSS spectra come with very low errors, which
c     produce large weight_obs(i) and are thus given too much weight in
c     the synthesis.

c     (F) fscale_chi2: I've introduced a "cooking factor", fscale_chi2,
c     which should be set to fscale_chi2 ~ actual-spectral-resolution /
c     Odl_syn. When you synthesize a spectrum of say, resolution 3 Angs
c     but sampled at every 1 Angs then only ~ 1 out of every 3 points is
c     actually independent. The chi2 is thus ~ 3 times larger than what
c     it should be, so dividing it by 3 would fix the chi2 scale. This
c     is NOT really relevant insofar as the value of the reduced chi2 is
c     concerned. In fact, the whole quest for a minimum chi2 through the
c     annealing scheme used here makes the actual absolute scale of chi2
c     rather irrelevant. The only noticeable effect of this re-scaling
c     of chi2 is in the METROPOLIS algorithm, where one uses a
c     likelyhood ratio to decide whether upwards steps are taken. This
c     ratio is exp(-(Enew - Eold)), where E = 0.5 * chi2, which IS
c     sensitive to the absolute scale of chi2. I therefore implement the
c     chi2 scaling ONLY in this part of the metropolis routine, which
c     becomes: exp(-(Enew - Eold) / fscale_chi2). The expected effect
c     for fscale_chi2 > 1 is that more updwards metropolis-steps will be
c     accepted, increasing the overall metropolis-efficiency. In
c     practice, since we fiddle with chi2 with the weight-factors in the
c     annealing loop, this extra-parameter probably has the same effect
c     as changing the weighting factors! The bottom line is: If you
c     believe your error spectrum (particularly its amplitude), just set
c     fscale_chi2 = 1 and forget this feature. (I myself never used it
c     seriously).


!     * RAPID-CHI2 TRICK: This used to be Starlight03_v04SDSS.for, but
!     is now VERY different. There were several minor changes and one
!     HUGE change. The huge modification is the use of a truncated
!     series to represent the extinction term in the model
!     spectrum. This simple trick, which works fine as long as AV is not
!     too large, allows a MUCH-MUCH faster calculation of the chi2 by
!     circumventing the need to sum over all lambdas every time. See
!     notes for details! The bottom line is that we can now let the
!     Metropolis routine run for MUCH-MUCH longer than before in less
!     time! (Indeed, the advantage is only noticeable when we do force
!     huge Metropolis runs.) The routines involved are Calc_RC_terms &
!     VeryFastF_chi2. The screen output register info on the difference
!     between approximate & exact chi2's, but nothing is stored in the
!     ouput files! Note that we have not eliminated the full-exact-chi2
!     routines altogether; some are still used (for output, to minimize
!     over kinematic parameters v0 & vd, etc). The faster routines are
!     used essentially only in the (x,AV,YAV) Metropolis chain...

c     The RC-scheme involves the computation of 3 terms: RC1, RC2(n,j)
c     and RC3(n,j,k). RC2 & RC3 involve the base spectra f_base(j,i). My
c     implementation of this scheme uses f_base(j,i) which are already
c     convolved with the kinematical filter G(v0,vd). This means that I
c     first convolve with G and only then apply the extinction
c     factor. This is slightly different than appliyng the extinction
c     factor and then convolving correction (the difference is small
c     because the extinction curve is ~ constant over scales of
c     vd). This RC-chi2 is inly used in the estimation of x, AV &
c     YAV. The kinematical fitting of v0 & vd uses a Golden Search
c     method with a chi2 function which included the extinction term
c     inside the convolution. This note is just to recall myself that
c     there are two slightly different chi2's in the code...

!     Still to be implemented: 

!     * A bias-correcting trick... (for the AV-expansion series)
!
!     * There might be a reasonable method to define an initial
!     annealing weight-factor (<=> Temperature) based on the typical
!     value of chi2 for a few random models and what factor
!     (Temperature) would it take to make it go to a chi2 ~ 1/N_dof in
!     not-too-many steps... Think about it!
!
!     * The "MCMC-thing", with ~ 12 paralell chains, mean & std dev
!     estimates (with burn-in convergence test)

c     ATT: This abstract is permanently under construction...
c
c     Cid@Lagoa - 27/Jan/2005


      SUBROUTINE FitSpectrum(idum_orig,llow_SN,lupp_SN,Olsyn_ini
     $     ,Olsyn_fin,Odl_syn,aux_fscale_chi2
     $     ,base_dir,obs_dir,mask_dir,etc_dir,out_dir
     $     ,FitOrSamplePDF_option ,IsErrSpecAvailable
     $     ,IsFlagSpecAvailable,IsFIRcOn,IsQHRcOn,IsPHOcOn,flux_unit)


c ***************************************************************************
c *				DEFINITIONS				    *
c ***************************************************************************
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
      parameter (Nchains_max = 1000,Npar_max = 601)
      parameter (NQHR_Ys_max=10, NPHO_Ys_max=100)
!IIB! Changed Nl_max to NAVgrid_max!
      parameter (NAVgrid_max = 1001)

c     Spectral quantities: observed, synthetic, mask & auxiliar stuff
      real l_obs(Nl_max)       , f_obs(Nl_max)
      real err_obs(Nl_max)     , weight_obs(Nl_max)
      real l_syn(Nl_max)       , f_syn(Nl_max)
      real llow_mask(Nl_max)   , lupp_mask(Nl_max)
      real mask_weight(Nl_max) , f_mask(Nl_max)
      real l_aux(Nl_max)
      real f_aux1(Nl_max)      , f_aux2(Nl_max)   , f_aux3(Nl_max)

c     Base spectra, age, Z & "absolute"-luminosity @ l_norm (= fbase_norm)
      real f_base(Nb_max,Nl_max) , fbase_norm(Nb_max) 
      real age_base(Nb_max) , Z_base(Nb_max) , alphaFe_base(Nb_max)
      character*15 base_name(Nb_max)

c     Population vector arrays in flux (x) & mass (mu_ini) fractions
      real x(Nb_max) , x_min(Nb_max) , mu_ini(Nb_max)

c     Auxiliar arrays
      real aux1(Nb_max) , aux2(Nb_max) , aux3(Nb_max)

c     Mstars = array which says what fraction of the INITIAL stellar
c     mass in population # j is still in stars at age(j). The name
c     "Mstars" is inspired by the bc2003 nomenclature. Mstars(j) is only
c     used in the flux-to-mass fractions conversions (ie, for output
c     purposes). mu_cor() are the corrected mass-fractions.
      real Mstars(Nb_max) , mu_cor(Nb_max)

c     AV_array: AV_tot(j) = population dependent extinction
      real AV_tot(Nb_max)

c     Arrays to store Original (hence the "O") observed & base spectra
      real Ol_obs(Nl_max)   , Of_obs(Nl_max)
      real Ol_base(Nl_max)  , Of_base(Nb_max,Nl_max)
      real Oerr_obs(Nl_max) , Oflag_obs(Nl_max)

c     Chain-arrays!
      real ChainPar(Npar_max,Nchains_max)
      real Gpar_min(Npar_max) , Gpar_ave(Npar_max)
      real Cpar_min(Npar_max,Nchains_max) , CEne_min(Nchains_max)
      real EX0sChainPar(Npar_max,Nchains_max)
      real EX0sGpar_min(Npar_max) , EX0sGpar_ave(Npar_max)
      real EX0sCpar_min(Npar_max,Nchains_max) ,EX0sCEne_min(Nchains_max)

c     file names - including arq_ETCinfo
      character*100 arq_obs , arq_base , arq_masks , arq_config
      character*100 arq_out , arq_aux , arq_lix , arq_ETCinfo

c     Reddening-law option string 
      character*3 red_law_option

c     Option-String which defines which method will be used to define
c     who is relevant and who is not in the pop vector x_Test!
      character*5 EXOs_method_option
      character*3 EXOs_PopVector_option

c     Defining base, obs, mask, etc & out directories strings!
      character*100 base_dir , obs_dir , mask_dir , out_dir , etc_dir

c     Clipping-Method-opting-string
      character*6 clip_method_option

c     Option-string which defines whether we'll fit or sample the PDF
      character*3 FitOrSamplePDF_option

c     commons...
      common/Red_factor/ red_term(Nl_max)
      common/Obs/ Nl_obs , l_obs , f_obs , weight_obs
      common/Base/ N_base , f_base
      common/Obase/ Of_base
      common/N_Gauss/ N_int_Gauss
      common/Age_and_Met/ age_base , Z_base
      common/CalcMass/ q_norm , fbase_norm , i_FitPowerLaw , fobs_norm
      common/ScaleChi2_Factor/ fscale_chi2
      common/RC_AV_Stuff/ NRC_AV , i_RC_CRASH_FLAG , NRC_AV_Default

c     Lower & upper limits on parameters (x,AV,exAV,fn)
      common/ParLimits/ par_low(Npar_max) , par_upp(Npar_max)

c     ScubiduCadeVoce Technical Parameters
      common/Scubi_TechnicalParameters/ alpha , Falpha , IsGRTestHard ,
     $     eff_IDEAL , i_UpdateEps , Feps , i_UpdateAlpha ,
     $     i_MoveOneParOnly , i_HelpParWithMaxR , prob_jRmax ,
     $     Npar4GRtest ,i_UpdateAVYAVStepSeparately ,
     $     ind_HelpParMove2Average(Npar_max) ,
     $     prob_HelpParMove2Average(Npar_max)

c     Opened file units, for output
      common/OutputFileUnits/ i_UnitScubi

c     Willy time functions
      real willy_tiempo(2)

!Cid19Jan07-FIX! Added chi2, AV & norm. factor for best SSP fits!
!v05! Array to store best-fit SSP spectrum
      real SSP_AV(Nb_max) , SSP_x(Nb_max) , SSP_chi2(Nb_max) ,
     $     SSP_adev(Nb_max)
      real Best_f_SSP(Nl_max)
      common/SSP_chi2_stuff/ x_SSP , f_SSP(Nl_max)
      external F_SSP_chi2_AV

!Cid14Jan07-FXK! Added FXK = Fixed kinematics option!!
      common/FXK_fits/ IsThisAnFXKFit

*     !!MEX!!
      real exAV(Nb_max)
      real exAV_min(Nb_max) , exAV_low(Nb_max) , exAV_upp(Nb_max)
      common/exAV_stuff/ N_exAV , ind_exAV(Nb_max)

!!MEX!! IMPORTANT CHANGE: Introduced flag to use SLOW chi2 always!
      common/ToRCOrNotToRC_ThatIsTheQuestion/ i_UseSlowChi2Always

*ETC* Common to get FIR, QHR & PHO chi2's from function ETC_CalcETCchi2 for output 
      common/ETC_Chi2s/ chi2_FIR , chi2_QHR , chi2_PHO

*FIR* Arrays to store Lambda-Averaged pop-vectors: min, ave &
*     chain-solutions. OBS: They are ONLY used in the output-phase!!
      real LAx_min(Nb_max) , LAx_ave(Nb_max) , 
     &     LAx_chain(Nb_max,Nchains_max)

*FIR* Data-info-common for FIR-work
      common/FIR_DataStuff/ FIRlogLobs_norm , 
     $     FIR_logLFIR_TOT , FIR_LFIRFrac2Model ,
     $     FIR_logLFIR_obs , FIR_ErrlogLFIR ,
     $     FIR_logLFIR_low , FIR_logLFIR_upp , FIR_RangelogLFIR ,
     $     FIRChi2ScaleFactor , iFIR_IsFIRcOn

*FIR* Base-info-common. 
*     OBS: In the main code (ie here) these stuff (actually just part of
*     it!) is only used to produce FIR-related output
!IIB! Changed Nl_max to NAVgrid_max!
      common/FIR_BaseStuff/ NFIR_AV , FIRAV_grid(NAVgrid_max) ,
     $     FIRAV_low , FIRAV_upp , FIRdAV , 
     $     FIRLbol(Nb_max) , FIRFracLion(Nb_max) , FIRBolCor(Nb_max) ,
     $     FIRR_opt(Nb_max,NAVgrid_max) , FIRR_Lya(Nb_max,NAVgrid_max) , 
     $     FIRR_LCE(Nb_max,NAVgrid_max) ,
     $     FIRRmat(Nb_max,NAVgrid_max) , FIRq_norm , NFIR_base

*FIR* Common (used just for output) for predicted L_FIR & x_FIR & x_BOL
*     pop-vectors
      common/FIR_OutputStuff/ FIRModlogLFIR , FIRModlogLBOL , 
     &     x_FIR(Nb_max),x_BOL(Nb_max)

*QHR* Data-info-common for QHR-work
      common/QHR_DataStuff/ QHR_lambda(NQHR_Ys_max) ,
     $     QHR_frecomb(NQHR_Ys_max)  ,
     $     QHR_logY_TOT(NQHR_Ys_max) , QHR_YFrac2Model(NQHR_Ys_max) ,
     $     QHR_ErrlogY(NQHR_Ys_max)  , QHR_RangelogY(NQHR_Ys_max) ,
     $     QHR_Chi2ScaleFactor(NQHR_Ys_max) , 
     $     QHR_logY_obs(NQHR_Ys_max) ,
     $     QHR_logY_upp(NQHR_Ys_max) , QHR_logY_low(NQHR_Ys_max)  ,
     $     QHR_qlambda(NQHR_Ys_max) ,
     $     QHR_GlobalChi2ScaleFactor ,
     $     QHR_logLobs_norm , 
     $     NQHR_Ys , iQHR_IsQHRcOn

*QHR* Base-info-common. 
*     OBS: In the main code (ie here) these stuff (actually just part of
*     it!) is only used to produce FIR-related output
      common/QHR_BaseStuff/ QHR_qH40(Nb_max) , QHR_Q2Lnorm40(Nb_max) ,
     $     QHR_q_norm , QHRbeta_I , NQHR_base

*QHR* Common (used just for output) for predicted QHR-&-ELR-things
      common/QHR_OutputStuff/ QHR_x_QH0(Nb_max) , QHR_x_QHeff(Nb_max) ,
     $     QHR_x_Y(Nb_max,NQHR_Ys_max) , QHR_ModlogY(NQHR_Ys_max) ,
     $     QHR_chi2_Y(NQHR_Ys_max) , QHR_log_QH0 , QHR_log_QHeff ,
     $     ELR_chi2_logR , ELR_ModlogR

*ELR* Data-info for ELR-work - goes together w/QHR-things!
      common/ELR_DataStuff/ ELR_lambdaA , ELR_lambdaB , 
     $     ELR_logRint , ELR_logRobs , 
     $     ELR_ErrlogR , ELR_logR_low , ELR_logR_upp , 
     $     ELR_RangelogR , ELR_Chi2ScaleFactor , 
     $     iY_ELR_A , iY_ELR_B , IsELROn , 
     $     IsELRInputInconsistentWarningFlag

*PHO* filter name & file name common for PHO-work
      character*10  PHO_name(NPHO_Ys_max)
      character*100 PHO_Filter_file(NPHO_Ys_max)
      common/PHO_NamesStuff/ PHO_name , PHO_Filter_file

*PHO* Data-info-common for PHO-work
      common/PHO_DataStuff/ PHO_logY_TOT(NPHO_Ys_max) ,
     $     PHO_YFrac2Model(NPHO_Ys_max) , PHO_ErrlogY(NPHO_Ys_max)  ,
     $     PHO_RangelogY(NPHO_Ys_max) , 
     $     PHO_Chi2ScaleFactor(NPHO_Ys_max) ,
     $     PHO_logY_obs(NPHO_Ys_max) ,
     $     PHO_logY_low(NPHO_Ys_max) , PHO_logY_upp(NPHO_Ys_max) ,
     $     PHO_GlobalChi2ScaleFactor , 
     $     PHO_redshift , PHO_logLobs_norm , NPHO_Ys , iPHO_IsPHOcOn

*PHO* Base-info-common. 
!IIB! Changed Nl_max to NAVgrid_max!
      common/PHO_BaseStuff/ NPHO_AV , PHOAV_grid(NAVgrid_max) ,
     $     PHOAV_low , PHOAV_upp , PHOdAV , 
     $     PHO_MeanLambda(NPHO_Ys_max) , PHO_StdDevLambda(NPHO_Ys_max) ,
     $     PHO_q_MeanLambda(NPHO_Ys_max) ,
     $     PHO_3Mat(Nb_max,NPHO_Ys_max,NAVgrid_max) ,
     $     PHO_q_norm , NPHO_base

*PHO* Common (used just for output) for predicted PHO-things
      common/PHO_OutputStuff/ PHO_ModlogY(NPHO_Ys_max) , 
     $     PHO_chi2_Y(NPHO_Ys_max) , PHO_x_Y(Nb_max,NPHO_Ys_max)

!OPTfn! Optically optimal fn - common. Used to pass IsOptimize_fn_OPT
!     from here to wher it is needed.
      common/Optimize_fn_OPT/ IsOptimize_fn_OPT , Optimal_fn_OPT

c     Implemented more-than-one header-lines - ElCid@Sanchica - 29/Jan/2012
      character*1 cerquinha
c ***************************************************************************



c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c +                       INPUT & PRE-PROCESSING STEPS                      +
c +                                                                         +
c +   In this part we reset variables, read input data & perform all        +
c + necessary pre-processing, like re-sampling to a common lambda-scale,    +
c + normalization, definition of weights, etc.                              +
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


c ***************************************************************************
c *	                          INPUT INFO                                * 
c ***************************************************************************
*ATT* This block used to be **after** the reset-things block below, but
*     I moved it here to force the SkipExistingOutFiles-check to be done
*     before the slow (~ 1 second) reset-things block. This should speed
*     up the repetition of large grids ...
*
*     OBS: Note that I still read the do-it-or-not flag i_SkipExistingOutFiles 
*     from arq_config instead of the grid file. This is rather stupid,
*     but saves me the trouble of adapting files and scripts...
*
*     ElCid@Sanchica - 15/Ago/2011


*ETC* Input observed spectrum, ETCinfo file-name, luminosity-distance
*     (in Mpc), config-file, base, masks, initial v0 & vd, and output
*     filename on a galaxy by galaxy basis! 
*     ATT: Provide an arq_ETCinfo filename (eg 'nodata') even if no
*     ETCc-fit is being attempted. Likewise, if no distance is known,
*     provide something (e.g., 100 Mpc). 
      write (*,*)
      write (*,'(a)') '**> ... Reading input from grid file...'
      read (*,*) arq_obs , arq_ETCinfo , LumDistInMpc , 
     &     arq_config , arq_base , arq_masks , red_law_option , 
     &     v0_ini , vd_ini , arq_out 

*NEW* Skips runs already made: A 2007-implemented trick now in a new &
*     faster skin:-)
      open (unit=1,file=arq_config,status='old')
      do i=1,38 
        read (1,*)
      end do
      read (1,*) i_SkipExistingOutFiles
      close(1)

      if (i_SkipExistingOutFiles.EQ.1) then
         iTest = iCheckExistingOutFiles(out_dir,arq_out)
         if (iTest.NE.0) then
          write (*,'(a,a)')' OOPS! arq_out already exists!! Skipping ',
     $       arq_out
           return
         end if
      end if


*ETC* Store flags IsFIRcOn, IsQHRcOn & IsPHOcOn into variables
*     acessibles by the FIR/QHR/PHO routines via the respective
*     data-commons.
      iFIR_IsFIRcOn = IsFIRcOn
      iQHR_IsQHRcOn = IsQHRcOn
      iPHO_IsPHOcOn = IsPHOcOn

*ETC* Added this for a coherent output...
      if ((iPHO_IsPHOcOn.NE.1).AND.(iQHR_IsQHRcOn.NE.1).AND.
     $     (iFIR_IsFIRcOn.NE.1)) arq_ETCinfo = 'nodata'

*     Screen check...
      if (i_verbose.GT.0) then
         write (*,'(a,a)')    'arq_obs        = ', arq_obs
         write (*,'(a,a)')    'arq_ETCinfo    = ', arq_ETCinfo
         write (*,'(a,f8.3)') 'LumDistInMpc   = ', LumDistInMpc  
         write (*,'(a,a)')    'arq_config     = ', arq_config
         write (*,'(a,a)')    'arq_base       = ', arq_base
         write (*,'(a,a)')    'arq_masks      = ', arq_masks
         write (*,'(a,a)')    'red_law_option = ', red_law_option  
         write (*,'(a,f8.3)') 'v0_ini         = ', v0_ini
         write (*,'(a,f8.3)') 'vd_ini         = ', vd_ini
         write (*,'(a,a)')    'arq_out        = ', arq_out
      end if

*ETC* Set iETC_IsETCcOn = 1 if something needs to be read from
*     arq_ETCinfo (and/or arq_config), even if only for predictions.
      iETC_IsETCcOn = 0
      if  ((iFIR_IsFIRcOn.EQ.1).OR.
     $     (iQHR_IsQHRcOn.EQ.1).OR.
     $     (iQHR_IsQHRcOn.EQ.-1).OR.
     $     (iPHO_IsPHOcOn.EQ.1).OR.
     $     (iPHO_IsPHOcOn.EQ.-1)) then
         iETC_IsETCcOn = 1
      end if
c ***************************************************************************



c ***************************************************************************
c *			        RESET STUFF			            * 
c ***************************************************************************
!Cid01Apr07-FIX! Reseting all time counters (just in case)
      i_ini_time        = 0
      w_ini_TotTime     = 0
      w_ini_UsrTime     = 0
      w_ini_SysTime     = 0
      i_fin_time        = 0 
      w_fin_TotTime     = 0
      w_fin_UsrTime     = 0
      w_fin_SysTime     = 0

      i_ini_time_FF     = 0 
      w_ini_FF_TotTime  = 0 
      w_ini_FF_UsrTime  = 0 
      w_ini_FF_SysTime  = 0 
      i_fin_time_FF     = 0 
      w_fin_FF_TotTime  = 0 
      w_fin_FF_UsrTime  = 0 
      w_fin_FF_SysTime  = 0 

      i_ini_time_RF     = 0 
      w_ini_RF_TotTime  = 0 
      w_ini_RF_UsrTime  = 0 
      w_ini_RF_SysTime  = 0 
      i_fin_time_RF     = 0 
      w_fin_RF_TotTime  = 0 
      w_fin_RF_UsrTime  = 0 
      w_fin_RF_SysTime  = 0 

      i_ini_time_BI     = 0 
      w_ini_BI_TotTime  = 0 
      w_ini_BI_UsrTime  = 0 
      w_ini_BI_SysTime  = 0 
      i_fin_time_BI     = 0 
      w_fin_BI_TotTime  = 0 
      w_fin_BI_UsrTime  = 0 
      w_fin_BI_SysTime  = 0 

      i_ini_time_EX0    = 0 
      w_ini_EX0_TotTime = 0 
      w_ini_EX0_UsrTime = 0 
      w_ini_EX0_SysTime = 0 
      i_fin_time_EX0    = 0 
      w_fin_EX0_TotTime = 0 
      w_fin_EX0_UsrTime = 0 
      w_fin_EX0_SysTime = 0 


!Cid29Dec06-FIX! Introduced time-counters (to help in config-tests)
      i_ini_time    = time()
      w_ini_TotTime = etime(willy_tiempo)
      w_ini_UsrTime = willy_tiempo(1)
      w_ini_SysTime = willy_tiempo(2)

c     Reset All arrays
      do i=1,Nl_max
         f_aux1(i)      = 0.
         f_aux2(i)      = 0.
         f_aux3(i)      = 0.
         l_obs(i)       = 0.
         f_obs(i)       = 0.
         l_syn(i)       = 0.
         f_syn(i)       = 0.
         err_obs(i)     = 0.
         weight_obs(i)  = 0.
         llow_mask(i)   = 0.
         lupp_mask(i)   = 0.
         mask_weight(i) = 0.
         f_mask(i)      = 0.
         Ol_obs(i)      = 0.
         Of_obs(i)      = 0.
         Oerr_obs(i)    = 0.
         Oflag_obs(i)   = 0.
         Ol_base(i)     = 0.
         do j=1,Nb_max
            f_base(j,i)  = 0.
            Of_base(j,i) = 0.
         end do
         Best_f_SSP(i)  = 0.    !v05!
         f_SSP(i)       = 0.    !v05!
      end do
      
!Cid29Dec06-FIX! Fixed silly i->j bug found by William on Dec/23/2006!
      do j=1,Nb_max
         age_base(j)     = 0.
         Z_base(j)       = 0.
         fbase_norm(j)   = 0.
         x(j)            = 0.
         x_min(j)        = 0.
         AV_tot(j)       = 0.
         mu_ini(j)       = 0.
         mu_cor(j)       = 0.
         Mstars(j)       = 0.
         base_name(j)    = 'none'
         SSP_AV(j)       = 0.
         SSP_x(j)        = 0.
         SSP_chi2(j)     = 0.
         SSP_adev(j)     = 0.
      end do

*     !!MEX Reset x(), AV & exAV estimates to zero flag (Cid/Apr/19/2007)
      do j=1,Nb_max
         x(j)        = 0.
         x_min(j)    = 0.
         exAV_min(j) = 0.
         exAV(j)     = 0.
         exAV_low(j) = 0.
         exAV_upp(j) = 0.
         ind_exAV(j) = 0
         aux1(j)     = 0.
         aux2(j)     = 0.
         aux3(j)     = 0.
      end do
      AV_min = 0.
      N_exAV = 0

c     Reset clip-bug flag & CRASH flag for RC-scheme (RC = Rapid Chi2)
      iCLIPBUG_flag   = 0
      i_RC_CRASH_FLAG = 0

c     Reset sigma-clipping & global step counters
      i_clip_loop   = 0
      Ntot_cliped   = 0
      Nglobal_steps = 0

c     Set idum = seed for random-#-generator. idum_orig is preserved for output
      idum = idum_orig

c     Set chi2 scaling factor (used ONLY in Metropolis routine)
      fscale_chi2 = aux_fscale_chi2

c     Reset arrays ind_ & prob_HelpParMove2Average arrays
      do j=1,Npar_max
         ind_HelpParMove2Average(j)  = 0
         prob_HelpParMove2Average(j) = 0.0
      end do

c     Reset Chain-arrays!
      do j=1,Npar_max
         Gpar_min(j)     = 0
         Gpar_ave(j)     = 0
         EX0sGpar_min(j) = 0
         EX0sGpar_ave(j) = 0
         do i_chain=1,Nchains_max
            ChainPar(j,i_chain)     = 0
            Cpar_min(j,i_chain)     = 0
            EX0sCpar_min(j,i_chain) = 0
            CEne_min(i_chain)       = 0
            EX0sCEne_min(i_chain)   = 0
         end do 
      end do 

c     Reset SaveStep-flag = file unit!
      i_UnitScubi = 0

* !!MEX!! Reset desperate use-only-slow-chi2 option!
      i_UseSlowChi2Always = 0

*FIR* Reset FIR_BaseStuff arrays & other things for safety!
      do j=1,Nb_max
         x_FIR(j)       = 0.
         x_BOL(j)       = 0.
         FIRLbol(j)     = 0.
         FIRFracLion(j) = 0.
         FIRBolCor(j)   = 0.
!IIB! Changed Nl_max to NAVgrid_max!
         do i_AV=1,NAVgrid_max
            FIRR_opt(j,i_AV) = 0.
            FIRR_Lya(j,i_AV) = 0.
            FIRR_LCE(j,i_AV) = 0.
            FIRRmat(j,i_AV)  = 0.
         end do
         LAx_min(j) = 0.
         LAx_ave(j) = 0.
         do i_chain=1,Nchains_max
            LAx_chain(j,i_chain) = 0.
         end do
      end do
      FIRfobs_norm       = 0.
      FIRlogLobs_norm    = 0.
      FIRq_norm          = 0.
      NFIR_base          = 0
      FITRbeta_D         = 0.
      FITRbeta_I         = 0.
      FIR_logLFIR_TOT    = 0.
      FIR_logLFIR_obs    = 0.
      FIR_LFIRFrac2Model = 0.
      FIR_ErrlogLFIR     = 99.
      FIR_RangelogLFIR   = 0.
      FIRChi2ScaleFactor = 0.
      FIR_logLFIR_low    = -99.
      FIR_logLFIR_upp    = 99.
      NFIR_AV            = 0
      FIRAV_low          = 0.
      FIRAV_upp          = 0.
      FIRdAV             = 0.
!IIB! Changed Nl_max to NAVgrid_max!
      do i_AV=1,NAVgrid_max
         FIRAV_grid(i_AV) = 0.
      end do
*FIR* OBS: The 99's above are just to avoid inf's & nan's if I decide to
*     call FIR-routines in non-FIRc fits (ie, when IsFIRcOn = 0)!

*QHR* Reset all QHR stuff for safety!
      QHR_GlobalChi2ScaleFactor = 0.
      QHR_q_norm       = 0.
      QHRbeta_I        = 0.
      NQHR_base        = 0
      QHR_logLobs_norm = 0. 
      NQHR_Ys          = 0
      QHR_log_QH0      = 0.
      QHR_log_QHeff    = 0.
      do iY=1,NQHR_Ys_max
         QHR_lambda(iY)          = 0. 
         QHR_frecomb(iY)         = 0.
         QHR_logY_TOT(iY)        = 0.
         QHR_YFrac2Model(iY)     = 0.
         QHR_ErrlogY(iY)         = 0.
         QHR_RangelogY(iY)       = 0.
         QHR_Chi2ScaleFactor(iY) = 0.
         QHR_logY_obs(iY)        = 0.
         QHR_logY_upp(iY)        = 0.
         QHR_logY_low(iY)        = 0.
         QHR_qlambda(iY)         = 0.
         QHR_ModlogY(iY)         = 0.
         QHR_chi2_Y(iY)          = 0.
      end do
      do j=1,Nb_max
         QHR_qH40(j)      = 0.
         QHR_Q2Lnorm40(j) = 0.
         QHR_x_QH0(j)     = 0.
         QHR_x_QHeff(j)   = 0.
         do iY=1,NQHR_Ys_max
            QHR_x_Y(j,iY) = 0.
         end do
      end do
*ELR* reset ELR-things
      ELR_lambdaA          = 0.
      ELR_lambdaB          = 0.
      ELR_logRint          = 0.
      ELR_logRobs          = 0.
      ELR_ErrlogR          = 0.
      ELR_logR_low         = 0.
      ELR_logR_upp         = 0. 
      ELR_RangelogR        = 0.
      ELR_Chi2ScaleFactor  = 0.
      iY_ELR_A      = 0
      iY_ELR_B      = 0
      IsELROn       = 0
      ELR_chi2_logR = 0.
      ELR_ModlogR   = 0.
      IsELRInputInconsistent = 0

*PHO* Reset all PHO-stuff for safety!
      PHO_GlobalChi2ScaleFactor = 0. 
      NPHO_Ys          = 0
      PHO_redshift     = 0.
      PHO_logLobs_norm = 0.
      PHO_q_norm       = 0.
      NPHO_base        = 0
      do iY=1,NPHO_Ys_max
         PHO_name(iY)            = 'none'
         PHO_Filter_file(iY)     = 'none'
         PHO_logY_TOT(iY)        = 0.
         PHO_YFrac2Model(iY)     = 0.
         PHO_ErrlogY(iY)         = 0.
         PHO_RangelogY(iY)       = 0.
         PHO_Chi2ScaleFactor(iY) = 0.
         PHO_logY_obs(iY)        = 0.
         PHO_logY_upp(iY)        = 0.
         PHO_logY_low(iY)        = 0.
         PHO_MeanLambda(iY)      = 0.
         PHO_StdDevLambda(iY)    = 0.
         PHO_ModlogY(iY)         = 0.
         PHO_chi2_Y(iY)          = 0.
         do j=1,Nb_max
            PHO_x_Y(j,iY) = 0.
         end do
      end do
      NPHO_AV   = 0
      PHOAV_low = 0.
      PHOAV_upp = 0.
      PHOdAV    = 0.
!IIB! Changed Nl_max to NAVgrid_max!
      do i_AV=1,NAVgrid_max
         PHOAV_grid(i_AV) = 0.
      end do
      do j=1,Nb_max
         do iY=1,NPHO_Ys_max
            do i_AV=1,NAVgrid_max
               PHO_3Mat(j,iY,i_AV) = 0.
            end do
         end do
      end do

      Best_SSP_chi2 = 1.e36   !v05!
c ***************************************************************************



c ***************************************************************************
c *          READ "TECHNICAL-PARAMETERS" FROM CONFIGURATION FILE            *
c ***************************************************************************
      write (*,'(a,a)') '**> ... Reading arq_config = ',arq_config
      open (unit=1,file=arq_config,status='old')
      read (1,*)

c     Set normalization wavelength (for base & synthetic spectra) &
c     window limits for normalization of observed spectrum. Make sure
c     that llow_norm <= l_norm <= lupp_norm!
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*) l_norm
      read (1,*) llow_norm
      read (1,*) lupp_norm

      if ((l_norm.LT.llow_norm).OR.(l_norm.GT.lupp_norm)) then
         write (*,*) '[MAIN] l_norm outside llow_norm ==> lupp_norm!'
         write (*,*) llow_norm , l_norm , lupp_norm
         stop
      end if

c     Set a priori limits on extinction (AV), extra-extinction (YAV),
c     normalization factor (fn) and kinematical parameters (v0 & vd).
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*) AV_low
      read (1,*) AV_upp
!!MEX!! REMOVED YAV_low & _upp from arq_config!! Meaningless in MEX-fits!
      read (1,*) fn_low
      read (1,*) fn_upp
      read (1,*) v0_low
      read (1,*) v0_upp
      read (1,*) vd_low
      read (1,*) vd_upp

c     Set clip_method_option-string & sigma-clipping threshold
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*) clip_method_option
      read (1,*) sig_clip_threshold
 
c     Set weight control filter ... Points with weights larger than
c     wei_ave + wei_nsig_threshold * wei_sig will be "censored".
c     wei_nsig_threshold <= 0 means that no weight filter will be
c     applied.
      read (1,*) wei_nsig_threshold

!Cid20Jan07-FIX! Set safety margin cushion (in Angs). Previoulsy = 50 Angs.
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*) dl_cushion

c     Set minimum flux (in units of fobs_norm) to be modeled. Points
c     with f_obs < f_cut are masked out!
      read (1,*) f_cut

c     Set # of points in the brute-force convolution in GaussSmooth.
c     (Smarter integrators did not prove worthy!)
      read (1,*) N_int_Gauss

c     Set verbose... (i_verbose = 0 gives silent output)
c     iMassAgeZ_verbose lists estimates of Mass & Z just just for control...
      read (1,*) i_verbose
      read (1,*) i_verbose_anneal
      iMassAgeZ_verbose = i_verbose_anneal

c     Read fast-BC03-rebining option
      read (1,*) i_FastBC03_FLAG

c     Set if a POWER-LAW F_lambda ~ lambda**alpha_PowerLaw is to be
c     added to the base components from arq_base. Useful for AGN...  If
c     a power-law is added it will be the LAST component!
!!MEX!! Built-In Power-Law REMOVED from arq_config in MEX-version 
!     (but NOT removed completely from the code, so let's turn it off!!
      i_FitPowerLaw = 0
      alpha_PowerLaw = -9.99

!Cid08Fev07-FIX! Set if an existing arq_out file should be skiped (1) 
c     or if the fit should be redone and the file overwriten (0).
c     (Added this at William's request)
c     ATT: IRRELEVANT: It's been read above to seped up this trick!
      read (1,*) i_SkipExistingOutFiles

!v05! Decide whether or not to save the best fit SSP spectrum on arq_out
      read (1,*) i_SaveBestSingleCompFit

c     Read # of chains & max(x_j) for initial random chain pop-vectors.
      read (1,*) 
      read (1,*) 
      read (1,*) 
      read (1,*) 
      read (1,*) N_chains
      read (1,*) xinit_max

c     Read a bunch of technical parameters...
      read (1,*) i_UpdateEps
      read (1,*) i_UpdateAlpha
      read (1,*) Falpha
      read (1,*) i_MoveOneParOnly
      read (1,*) i_UpdateAVYAVStepSeparately
      read (1,*) i_HelpParWithMaxR
      read (1,*) prob_jRmax
      read (1,*) i_HelpPopVectorMove2Average
      read (1,*) prob_HelpPopVectorMove2Average
      read (1,*) i_HelpAVYAVMove2Average
      read (1,*) prob_HelpAVYAVMove2Average

c     Set NRC_AV_Default = initial # of terms in the RC AV-series
      read (1,*) NRC_AV_Default

* !!MEX!! Quick & Dirty way to turn-on use-only-slow-chi2 option!
      if (NRC_AV_Default.LT.0) then
         NRC_AV_Default = abs(NRC_AV_Default)
         i_UseSlowChi2Always = 1
         write (*,*)
         write (*,'(a)') '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
         write (*,'(a)') '@@> ATT: Will use only SLOW chi2!!'
         write (*,'(a)') '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
         write (*,*)
      end if

!OPTfn! Read choice of whether to (= 1) optimize fn or else (!= 1) fit
!     it as another par (ie, via Metropolis & Chains). A VERY important change!!
!     ElCid@Sanchica - 28/Jan/2012
      read (1,*) IsOptimize_fn_OPT

c     Read First-Fits (FF) technical parameters...
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*) Temp_ini_FF
      read (1,*) Temp_fin_FF
      read (1,*) fini_eps_FF
      read (1,*) ffin_eps_FF
      read (1,*) R_ini_FF
      read (1,*) R_fin_FF
      read (1,*) IsGRTestHard_FF
      read (1,*) N_loops_FF
      read (1,*) i_RestartChains_FF
      read (1,*) fN_sim_FF
      read (1,*) fNmax_steps_FF
      read (1,*) eff_IDEAL_FF

c     Read Burn-In stage technical parameters...
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*) R_BurnIn
      read (1,*) IsGRTestHard_BurnIn

c     Read EX0s-Fit (EX0s) technical parameters...
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*) EXOs_PopVector_option
      read (1,*) EXOs_method_option
*FIR* Added iEXOs_IsAve_x_of_lambda entry...
*     OBS: Motivated by, but independent of FIRc-stuff (Cid@IAA, 31/01/10)
      read (1,*) iEXOs_IsAve_x_of_lambda
      read (1,*) EX0s_Threshold
      read (1,*) Temp_ini_EX0s
      read (1,*) Temp_fin_EX0s
      read (1,*) fini_eps_EX0s
      read (1,*) ffin_eps_EX0s
      read (1,*) R_ini_EX0s
      read (1,*) R_fin_EX0s
      read (1,*) IsGRTestHard_EX0s
      read (1,*) N_loops_EX0s
      read (1,*) i_RestartChains_EX0s
      read (1,*) fN_sim_EX0s
      read (1,*) fNmax_steps_EX0s
      read (1,*) eff_IDEAL_EX0s
!     New stuff in v02!
      read (1,*) IsScaleNstepsInEX0sFits
      read (1,*) IsNoKin4LargeBaseInEX0sFits
      read (1,*) frac_NoKin4LargeBaseInEX0sFits

!Cid18Feb07-FIX! Define maximal base reduction factor for EX0-fits
!     such that at least NEX0s_base >= N_base * fEX0_MinBaseSize
!     (In previous versions this was fixed at 10%)
      read (1,*) fEX0_MinBaseSize

* ---------------------------------------------------------------------------
*ETC* If FIR and/or QHR and/or PHO stuff is to be used (either fits or
*     just predictions), then read ALL the following, regardless of
*     which of the ETC-features will be used in practice! This allows me
*     to keep arq_config with a rigid (yet general) format.
*
*     Cid@Lagoa - 09/Jan/2011
      if (iETC_IsETCcOn.EQ.1) then
         do i_skip_12_lines=1,12
            read (1,*)
         end do
*     FIRAV_grid parameters = PHOAV_grid parameters!
         read (1,*) NFIR_AV
         read (1,*) FIRAV_low
         read (1,*) FIRAV_upp
         NPHO_AV   = NFIR_AV
         PHOAV_low = FIRAV_low
         PHOAV_upp = FIRAV_upp
         read (1,*)
         read (1,*)
         read (1,*)
         read (1,*)
         read (1,*) QHRbeta_I
         read (1,*) FIRbeta_I
         read (1,*) FIRbeta_D

*     Basic sanity checks on ETC-things read from arq_config...
         if ((NFIR_AV.LT.10).OR.(FIRAV_low.LT.0.0).OR.(FIRAV_upp.LT.3.0)
     $        .OR.(FIRbeta_D.LT.0.0).OR.(FIRbeta_D.GT.1.0)
     $        .OR.(FIRbeta_I.LT.0.0).OR.(FIRbeta_I.GT.1.0)
     $        .OR.(QHRbeta_I.LT.0.).OR.(QHRbeta_I.GT.1.)) then
            write (*,'(a,a)')  '[ETC]>> Something wrong with the input', 
     &           ' read from arq_config = ',arq_config
            write (*,*) ' FIRbeta_I=',FIRbeta_I,'  FIRbeta_D=',FIRbeta_D
     $           ,' | QHRbeta_I = ', QHRbeta_I
            write (*,*) ' NFIR_AV=',NFIR_AV,'  FIRAV_low=',FIRAV_low, 
     &           '  FIRAV_upp=',FIRAV_upp
            stop
         end if

*     Reset things which won't be really used! (just for the sake of a
*     fully consistent output)
         if (iFIR_IsFIRcOn.EQ.0) then
            NFIR_AV   = 0
            FIRAV_low = 0.
            FIRAV_upp = 0.
            FIRbeta_I = 0.
            FIRbeta_D = 0.
         end if
         if (iQHR_IsQHRcOn.EQ.0) QHRbeta_I = 0.
         if (iPHO_IsPHOcOn.EQ.0) then
            NPHO_AV   = 0
            PHOAV_low = 0.
            PHOAV_upp = 0.
         end if
      end if
* ---------------------------------------------------------------------------

      close(1)
      write (*,'(a)') '**> ... finished reading arq_config!'

      if (i_verbose.GT.0) then 
         write (*,*)
         write (*,'(a,a)') '**> Welcome to PANcMEx_Starlight (v02c)!' ,
     $        ' Keep your fingers crossed...'
         write (*,*)
         write (*,'(a,a)') '**> Read config file: ' , arq_config
      end if

!     Check some of the input quantities...
      if ((EXOs_method_option.NE.'CUMUL').AND.(EXOs_method_option.NE.
     $     'SMALL')) then
         write (*,'(a,a)') ' [MAIN] OOPS! Bad EXOs_method_option ' ,
     $        EXOs_method_option
         stop
      end if

      if ((i_HelpPopVectorMove2Average.NE.0).AND.
     $     (i_HelpPopVectorMove2Average.NE.1)) then
         write (*,*) '[MAIN] Bad i_HelpPopVectorMove2Average! ' ,
     $        i_HelpPopVectorMove2Average
         stop
      end if
         
      if ((i_HelpAVYAVMove2Average.NE.0).AND.
     $     (i_HelpAVYAVMove2Average.NE.1)) then
         write (*,*) '[MAIN] Bad i_HelpAVYAVMove2Average! ' ,
     $        i_HelpAVYAVMove2Average
         stop
      end if

      if ((EXOs_PopVector_option.NE.'AVE').AND.
     $     (EXOs_PopVector_option.NE.'MIN')) then 
         write (*,*) '[MAIN] Bad EXOs_PopVector_option! ' ,
     $        EXOs_PopVector_option
         stop
      end if

      if (i_MoveOneParOnly.NE.1) then
         write (*,*) '[MAIN] Bad i_MoveOneParOnly: ' , i_MoveOneParOnly
         stop
      end if

      if ((clip_method_option.NE.'NOCLIP').AND.
     $     (clip_method_option.NE.'NSIGMA').AND.
     $     (clip_method_option.NE.'RELRES').AND.
     $     (clip_method_option.NE.'ABSRES')) then
         write (*,*) '[MAIN] Bad clip_method_option ',clip_method_option
         stop
      end if

c     Initialize # of terms in RC-AV series (passed by common)
      NRC_AV = NRC_AV_Default

c     Reset alpha to 1.
      alpha  = 1.0
c ***************************************************************************



c ***************************************************************************
c *	                          READ  MASK FILE   		        		    * 
c ***************************************************************************
!CALIFA! This part used to be locted some ~ 400 lines below, after the ETC block.
!	I am now reading the mask-file here, BEFORE the observed spectrum. The reason 
!	is that in this way I can disconsider masked-regions in the computation of a
!	pseudo-error when the input spectra come with NO errors and thus STARLIGHT needs
!	to guess one from a "S/N window". If this is not done then a S/N window informed 
!	in a grid file may inadvertedly include emission lines or sky-lines or other
!	things, leading to an overestimation of the actual rms in the S/N window.
      i_last_non_blank_char = 0
      do i=1,99
         if ((mask_dir(i:i+1).EQ.' ').AND.(i_last_non_blank_char.EQ.0))
     &        i_last_non_blank_char = i - 1
      end do
      arq_aux = mask_dir(1:i_last_non_blank_char) // arq_masks

c     Read MASK intervals: The array mask_weight() is used like this:
c     * mask_weight <= 0 on unwanted windows (eg emission lines)
c     * mask_weight > 0  is the factor by which we multiply the weight of the
c                        points within llow_mask ==> lupp_mask. This may be
c                        used to give more weight to key absorption features.
c     OBS: If any of the intervals overlap, then the one further down
c     the list will take precedence over the other, so be careful when
c     defining the intervals in arq_masks!
      write (*,'(a,a)') '**> ... Reading arq_masks = ',arq_aux
      open (unit=1,file=arq_aux,status='old')
      read (1,*) N_masks
      do i=1,N_masks
         read (1,*) llow_mask(i) , lupp_mask(i) , mask_weight(i)
      end do
      close(1)

!     Screen output
      if (i_verbose.GT.0) then
         write (*,'(a,a45,a,i3)') '**> Read masks file: ' , arq_masks , 
     &        ' N_masks = ' , N_masks
         write (*,'(a,a)') '    from dir = ' , 
     &        mask_dir(1:i_last_non_blank_char)
      end if
c ***************************************************************************



c ***************************************************************************
c *	                       READ OBSERVED SPECTRUM	          	            * 
c ***************************************************************************
c     OBS: The "O" arrays (where "O" stands for "Original") defined here
c     are renamed without the "O" later on, after pre-processing steps.
c
c     Implemented more-than-one header-lines - ElCid@Sanchica - 29/Jan/2012

!Cid29Mar07-FIX! Now read arq_obs from obs_dir, so figure what is the first
!     blank char in obs_dir & concatenate it: arq_obs --> arq_aux
      i_last_non_blank_char = 0
      do i=1,99
         if ((obs_dir(i:i+1).EQ.' ').AND.(i_last_non_blank_char.EQ.0))
     &        i_last_non_blank_char = i - 1
      end do
      arq_aux = obs_dir(1:i_last_non_blank_char) // arq_obs

      write (*,'(a,a)') '**> ... Reading arq_obs = ',arq_aux

c     First find-out how many header lines (those which start with '#') to skip
      open (unit=1,file=arq_aux,status='old')
      N_skip_arq_obs_header = 0
      do i=1,Nl_max
         read (1,'(a1)') cerquinha
         if (cerquinha.EQ.'#') then 
            N_skip_arq_obs_header = N_skip_arq_obs_header + 1
         else
            goto 2901
         end if
      end do
 2901 close(1)


!Cid24Mar07-FIX! Try to figure how many columns are in arq_obs to check
!     consistency of IsErrSpecAvailable & IsFlagSpecAvailable info given...
!     If a problem is detected issue a warning, but do NOT stop!!
!     OBS: Does NOT work if columns are separated by tabs or commas!
      write (*,'(a,a)') '**> ... Reading arq_obs = ',arq_aux
      open (unit=1,file=arq_aux,status='old')
      do i=1,N_skip_arq_obs_header
         read (1,*)
      end do
      read (1,*)
      read (1,'(a)') arq_lix
      close(1)
      icol = 0
      do i=1,99
         if ((arq_lix(i:i).EQ.' ').AND.(arq_lix(i+1:i+1).NE.' ')) 
     $        icol = icol + 1
      end do
      if (arq_lix(1:1).EQ.' ') icol = icol - 1

      if ((IsErrSpecAvailable.EQ.1).AND.(IsFlagSpecAvailable.EQ.1)
     $     .AND.(icol.NE.3)) then
 189  format(67('?'))
         write (*,189)
         write (*,'(a,a,a)') 'ATT: You said there are error & flag ' ,
     $        'columns in ' , arq_obs
         write (*,'(a,i3,a)') '     but I counted only ', icol + 1 ,
     $        ' columns!!'
         write (*,'(a)')      '     Hope you know what you are doing...'
         write (*,189)
      else if ((IsErrSpecAvailable.EQ.1).AND.(IsFlagSpecAvailable.EQ.0)
     $     .AND.(icol.LT.2)) then
         write (*,189)
         write (*,'(a,a,a)') 'ATT: You said there is an error ' ,
     $        'column in ' , arq_obs
         write (*,'(a,i3,a)') '     but I counted ', icol + 1 ,
     $        ' columns!!'
         write (*,'(a)')      '     Hope you know what you are doing...'
         write (*,189)
      end if


c     Read original/observed galaxy flux, error & flag spectra:
c     Of_obs(Ol_obs) , Oerr_obs(Ol_obs) & Oflag_obs(Ol_obs)
c     ATT: Implemented no error/flag-available possibilities...
c     Cid@IAA - 29/May/2005
      open (unit=1,file=arq_aux,status='old')

      do i=1,N_skip_arq_obs_header
         read (1,*)
      end do

      do i=1,Nl_max
      if ((IsErrSpecAvailable.EQ.1).AND.(IsFlagSpecAvailable.EQ.1)) then
         read (1,*,END=20) Ol_obs(i) , Of_obs(i) , Oerr_obs(i) ,
     $        Oflag_obs(i)
      else if ((IsErrSpecAvailable.EQ.1).AND.(IsFlagSpecAvailable.EQ.0))
     $        then
         read (1,*,END=20) Ol_obs(i) , Of_obs(i) , Oerr_obs(i)
         Oflag_obs(i) = 0
      else if ((IsErrSpecAvailable.EQ.0).AND.(IsFlagSpecAvailable.EQ.0))
     $        then
         read (1,*,END=20) Ol_obs(i) , Of_obs(i)
         Oflag_obs(i) = 0
      end if
      end do

      write (*,*) ' [MAIN] BUG reading input spectrum! N > Nl_max!! ' , 
     &     Nl_max
      stop
 20   NOl_obs = i - 1
      close(1)


!Cid07Feb07+Cid25Jan2010-FIX! Screen output 
      if (i_verbose.GT.0) then
         write (*,'(a,a30)') '**> Read data file: ',arq_obs
         write (*,'(a,a)') '    from dir = ' , 
     &        obs_dir(1:i_last_non_blank_char)
         write (*,'(4x,a,2(2x,f6.1,a),i5,a)') 'Data goes from: ' , 
     $        Ol_obs(1) , '  --> ' , Ol_obs(NOl_obs) , 
     $        ' Angs with N = ' , NOl_obs , ' pixels'
         if (IsErrSpecAvailable.EQ.1) write (*,'(4x,a,a)')
     $        'STARLIGHT has read an ERROR (column 3) from ' , arq_obs
         if (IsFlagSpecAvailable.EQ.1) write (*,'(4x,a,a)')
     $        'STARLIGHT has read an FLAG  (column 4) from ' , arq_obs
      end if

!Cid14Jan07-FIX! Check if Ol_obs is given in increasing order.
      do i=1,NOl_obs-1
         if (Ol_obs(i+1).LE.Ol_obs(i)) then
            write (*,*) ' [MAIN] PROBLEM w/spectrum: ',
     $           'lambda should be given in increasing order! : ' , i ,
     $           Ol_obs(i) , Ol_obs(i+1)
            stop
         end if
      end do

c     In the absence of a proper error spectrum, define Oerr_obs() using
c     the old method of rms in a window... (f_aux1 is just a
c     tmp-fake-weight-array to fool Stats4SN...). 
      if (IsErrSpecAvailable.EQ.0) then
         do i=1,NOl_obs
            f_aux1(i) = 100.
            if (Of_obs(i).LE.0.) f_aux1(i) = -100.
!CALIFA! f_aux1 now EXCLUDES windows masked away in the mask file
!	so not to overestimate the cooked-up error estimated here!
            do j=1,N_masks
         	   if ((Ol_obs(i).GE.llow_mask(j)).AND.
     &             (Ol_obs(i).LE.lupp_mask(j)).AND.
     &             (mask_weight(j).LE.0.)) f_aux1(i) = -100.
            end do
         end do

         call Stats4SN(NOl_obs,Ol_obs,Of_obs,f_aux1,llow_SN,lupp_SN
     $        ,aux_ave,aux_sig)

c     ATT: For safety, if the rms in the llow_SN -> lupp_SN is 0
c     (god-knows-why...) then recompute it using ALL the spectrum!
         if (aux_sig.LE.0.) then
            llow_SN = Ol_obs(1)
            lupp_SN = Ol_obs(NOl_obs)
            call Stats4SN(NOl_obs,Ol_obs,Of_obs,f_aux1,llow_SN,lupp_SN
     $           ,aux_ave,aux_sig)
         end if

         do i=1,NOl_obs
            Oerr_obs(i) = aux_sig
         end do

         if (i_verbose.GT.0) then
            write (*,'(a,a,f9.3)')
     $           '**> ATT: Defined error & flag spectra the old way!',
     $           ' S/N =' , aux_ave / aux_sig
         end if
      end if
c ***************************************************************************



c ***************************************************************************
c *	                READ BASE MASTER-FILE & BASE SPECTRA                    * 
c ***************************************************************************
c     Read base spectra: Of_base(i_base,Ol_base) (later renamed to
c     f_base). Notice that I'm now also reading Mstars & ind_exAV &
c     alphaFe_base from the base-file!!
c     ATT: All base spectra are read from base_dir, so we first find the
c     last non empty bit in string base_dir (to concatenate names below)
      i_last_non_blank_char = 0
      do i=1,99
         if ((base_dir(i:i+1).EQ.' ').AND.(i_last_non_blank_char.EQ.0))
     &        i_last_non_blank_char = i - 1
      end do

      write (*,'(a,a)') '**> ... Reading arq_base = ',arq_base
      open (unit=3,file=arq_base,status='old')
      read (3,*) N_base

      if (N_base.GT.Nb_max) then
         write (*,*) ' [MAIN] OOPS! N_base > Nb_max: ' , N_base , Nb_max
         write (*,*) ' [MAIN] Run aborted!'
         stop
      end if

      do j=1,N_base

c     Read base-info: filename, age, Z, code-name, Mstars, ind_exAV &
c     alphaFe_base
         read (3,*) arq_aux , age_base(j) , Z_base(j) , base_name(j) , 
     &        Mstars(j) , ind_exAV(j) , alphaFe_base(j)

c     Concatenate base spectrum filename with base directory & read it!
         arq_aux = base_dir(1:i_last_non_blank_char) // arq_aux
         call ReadSpecWithHeader(arq_aux,Nl_max,NOl_base,Ol_base,f_aux1)

!Cid29Mar07-FIX! Introduced tests to check if i_FastBC03_FLAG = 1 is 
!     being used with the original BC03/STELIB models. If not, I issue a
!     warning and revert i_FastBC03_FLAG to 0.
         if (i_FastBC03_FLAG.EQ.1) then
            i_shit = 0
            if (NOl_base.NE.6900) i_shit = 1
            if (int(Ol_base(1)).NE.91) i_shit = 1
            if (int(Ol_base(1000)).NE.3974) i_shit = 1
            if (int(Ol_base(2000)).NE.4974) i_shit = 1
            if (int(Ol_base(3000)).NE.5974) i_shit = 1
            if (int(Ol_base(4738)).NE.7712) i_shit = 1
            if (int(Ol_base(6361)).NE.9990) i_shit = 1
            if (int(Ol_base(6900)).NE.1600000) i_shit = 1
            if (i_shit.EQ.1) then
 188           format(80('?'))
               write (*,188)
               write (*,'(a,a)') 'OOPS! You have specified ' , 
     $              'i_FastBC03_FLAG = 1 (in arq_config) but the file: '
               write (*,'(a)') arq_aux
               write (*,'(a,a)') 'does NOT match the lambdas of ',
     $              'BC03/STELIB!! Reverting to i_FastBC03_FLAG = 0!'
               i_FastBC03_FLAG = 0
               write (*,188)
            end if
         end if

         do i=1,NOl_base
            Of_base(j,i) = f_aux1(i)
         end do

!Cid14Jan07-FIX! Check if Ol_base is given in increasing order.
         do i=1,NOl_base-1
            if (Ol_base(i+1).LE.Ol_base(i)) then
               write (*,*) ' [MAIN] PROBLEM w/base spectrum: ',
     $              'lambda should be given in increasing order! : ' , i 
     $              , Ol_base(i) , Ol_base(i+1)
               write (*,*) ' [MAIN] base-file = ', arq_aux
               stop
            end if
         end do
      end do
      close(3)

!Cid14Jan07-FIX!  Warn if Mstars(j) <= 0, but only fixes it if ALL are <= 0.
      IsAtLeastOneMstarsOk = 0
      do j=1,N_base
         if ((Mstars(j).GT.0.).AND.(IsAtLeastOneMstarsOk.EQ.0))
     $        IsAtLeastOneMstarsOk = 1
         if (Mstars(j).LE.0.) then
            write (*,*) ' [MAIN] ATT: Mstars <= 0!! ',
     $           'Fix it or disregard mass fractions on output...' , j
         end if
      end do
      if (IsAtLeastOneMstarsOk.EQ.0) then
         write (*,*) ' [MAIN] ATT: ALL Mstars <= 0!! WILL SET IT TO 1!'
         do j=1,N_base
            Mstars(j) = 1.
         end do
      end if

!     Screen output 
      if (i_verbose.GT.0) then
         write (*,'(a,a45,a,i3)') '**> Read base files: ' , arq_base , 
     &        ' N_base  = ' , N_base
         write (*,'(a,a)') '    from dir = ' , 
     &        base_dir(1:i_last_non_blank_char)
      end if
c ***************************************************************************



c ***************************************************************************
c *   !!MEX!!        SET UP EXTRA-EXTINCTIONS SCHEME (*exAV*)               *
c ***************************************************************************
*     Find Number of different extinctions to be fitted & Number of
*     components which will be allowed to have the extra-extinction exAV
*     (passed elsewhere by the common/exAV_stuff/).
      ind_exAV_max      = 0
      N_exAV_components = 0
      do j=1,N_base
         if (ind_exAV(j).GT.ind_exAV_max) ind_exAV_max = ind_exAV(j)
         if (ind_exAV(j).NE.0) N_exAV_components = N_exAV_components + 1
      end do
      N_exAV = ind_exAV_max

*     Check whether there are gaps in the ind_exAV array. It must be =
*     0, 1, 2, 3... with NO GAPS !
      do i_aux=1,N_exAV
         n_aux = 0
         do j=1,N_base
            if (ind_exAV(j).EQ.i_aux) n_aux = n_aux + 1
         end do
         if (n_aux.EQ.0) then
            write (*,*) '[MAIN] OOPS! Problem with ind_exAV array'
            write (*,*) '[MAIN] N_exAV = ', N_exAV ,
     $           ' but there is no pop with ind_exAV = ' , i_aux
            write (*,*) '[MAIN] The ind_exAV array must be = ', 
     $           '0, 1, ... Nex_AV, WITHOUT GAPS!'
            write (*,*) '[MAIN] Run aborted :('
            stop
         end if
      end do

*     Now read LIMITS to exAV from arq_base. These should be located at
*     the BOTTOM of the base file and AFTER 2 lines of comments...
      if (N_exAV.GT.0) then
         open (unit=3,file=arq_base,status='old')
         do i_skip=1,N_base+3
            read (3,*) 
         end do
         do i_aux=1,N_exAV
            read (3,*) i_exAV , z1 , z2
            exAV_low(i_exAV) = z1
            exAV_upp(i_exAV) = z2
         end do
      end if
      close(3)

c     Screen check exAV-stuff...
      if (i_verbose.GT.0) then
         write (*,'(a,i3)') '    Number of extra extinctions = ' , 
     $        N_exAV
         do i_exAV=1,N_exAV
            write (*,'(a,i3,2x,a,f5.2,2x,a,f5.2)') 
     $           '    * ind_exAV= ',i_exAV , 
     $           ' exAV_low= ', exAV_low(i_exAV) ,
     $           ' exAV_upp= ', exAV_upp(i_exAV)
         end do
         write (*,'(a,i3)')
     $        '    Number of components with different extinction = ' ,
     $        N_exAV_components
      end if
c ***************************************************************************



c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c ETC             READ ETC-info & initialize some ETC-stuff               ETC
c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c     Read ETC-data file, assumed to be in etc_dir
      if (iETC_IsETCcOn.EQ.1) then
         i_last_non_blank_char = 0
         do i=1,99
           if ((etc_dir(i:i+1).EQ.' ').AND.(i_last_non_blank_char.EQ.0))
     &           i_last_non_blank_char = i - 1
         end do
         arq_aux = etc_dir(1:i_last_non_blank_char) // arq_ETCinfo
         call Read_ETC_Info(i_verbose,arq_aux)
      end if

*FIR* Initialize model-related things: R-matrix, Lbol, & etc.
*     OBS1: This is done here because we need the full and original
*     Of_base spectra, which is chopped and changed below!
*     OBS2: This inialization is also necessary if we only want to
*     PREDICT the FIR-luminosity (iFIR_IsFIRcOn option = -1)
      if ((iFIR_IsFIRcOn.EQ.1).OR.(iFIR_IsFIRcOn.EQ.-1)) then
         call FIRInitialization(i_verbose,N_base,NOl_base,Ol_base,
     &        l_norm,red_law_option,FIRbeta_I,FIRbeta_D)
      end if

*QRH* Initialize QHR-related things
      if ((iQHR_IsQHRcOn.EQ.1).OR.(iQHR_IsQHRcOn.EQ.-1)) then

*QHR* Recover QHR_qlambda(iY) = A(lambda(iY)) / AV from a fake 1-element
*     red_term array... [just a trick]
         n_aux = 1
         do iY=1,NQHR_Ys
            l_aux(1) = QHR_lambda(iY)
            call Calc_red_term(red_law_option,l_norm,n_aux,l_aux,
     $           f_aux3,qaux_norm)
            QHR_qlambda(iY) = qaux_norm -2.5 * f_aux3(1)
         end do

*QHR* Initialize QHR-base-related things: QHR_qH40(), QHR_Q2Lnorm40(),
*     QHR_q_norm & NQHR_base, all in common/QHR_BaseStuff/.
         call QHRInitialization(i_verbose,N_base,NOl_base,Ol_base,
     $        l_norm,red_law_option)
      end if

*PHO* Initialize PHO-related things: PHO_MeanLambda(,), PHO_3Mat(,,),
*     etc., all in common/PHO_BaseStuff/
      if ((iPHO_IsPHOcOn.EQ.1).OR.(iPHO_IsPHOcOn.EQ.-1)) then
         call PHOInitialization(PHO_redshift,NPHO_Ys,i_verbose,N_base
     $        ,NOl_base,Ol_base,l_norm,red_law_option,etc_dir)
      end if
c ETC -ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC



c ***************************************************************************
c *	                          ..................		        		    * 
c ***************************************************************************
!CALIFA! The mask-file used to be read here, but not any more. See above.

!     Screen output (Added FXK info - Cid@Lagoa - 14/Jan/2007)
      if (i_verbose.GT.0) then
         write (*,'(a,a30,a,a)') '**> Modeling ',arq_obs,' ==> ',arq_out
         write (*,'(4x,3(a,f6.1),a,a,1p,e9.3)') 'Sampling: ',Olsyn_ini, 
     &        ' ==> ' , Olsyn_fin , '  with dl =' , Odl_syn , ' Angs' ,
     &        '  & fscale_chi2 = ' , fscale_chi2
         write (*,'(a,f6.1)') '    dl_cushion (A) = ' , dl_cushion
         write (*,'(a,a)')    '    Extinction Law = ' , red_law_option
         write (*,'(4x,a,f6.1,a,f6.1,a)')
     $        ' Initial kinematical parameters: v0 = ' , v0_ini ,
     $        ' & vd = ', vd_ini , ' km/s'

         if (FitOrSamplePDF_option.EQ.'PDF') then 
            write (*,'(a)')
     $           '    We are going to SAMPLE THE PDF for this galaxy'
         else if (FitOrSamplePDF_option.EQ.'FIT') then 
            write (*,'(a)')
     $           '    We are going to FIT the spectrum of this galaxy'
         else if (FitOrSamplePDF_option.EQ.'FXK') then 
            write (*,'(a,a)')
     $           '    We are going to FIT the spectrum of this galaxy' ,
     $           ' with FIXED KINEMATICAL PARAMETERS!'
         end if
      end if
c ***************************************************************************



c ***************************************************************************
c * PRE-PROCESSING I: RESAMPLING OBSERVED, BASE & WEIGHT/FLAG/MASK SPECTRA  *
c ***************************************************************************
c     Here I introduce the following trick: The wanted range is really
c     from Olsyn_ini ==> Olsyn_fin. I extend this range by ~ dl_cushion
c     Angstroms in both blue & red edges, but force the mask f_mask to
c     be = zero outside the Olsyn_ini ==> Olsyn_fin interval. Therefore
c     these extra points will NOT be synthesized. However, these
c     "safety" cushion will help minimizing edge effects in the
c     Gauss_smooth routine. In fact, these extra points are useful only
c     in the convolution of the synthetic spectrum with the kinematical
c     filter (in Gauss_smooth).

!Cid20Jan07-FIX! dl_cushion used to be fixed at 50 Angs. Now it is
!     read from arq_config!
c     Extend wavelength range by dl_cushion Angs to blue & red. 
      n_extra  = dl_cushion / Odl_syn
      lsyn_ini = Olsyn_ini - n_extra * Odl_syn
      lsyn_fin = Olsyn_fin + n_extra * Odl_syn
      dl_syn   = Odl_syn

c     Define lambda grid for spectral synthesis, l_syn (later renamed l_obs!)
      Nl_syn = (lsyn_fin - lsyn_ini) / dl_syn + 1

      if (Nl_syn.GT.Nl_max) then
         write (*,*) ' [MAIN] OOPS! Nl_syn > Nl_max: ' , Nl_syn , Nl_max
         stop
      end if

      do i=1,Nl_syn
         l_syn(i) = lsyn_ini + (i - 1) * dl_syn
      end do

c     Rebin observed spectra. I STRONGLY suggest you rebin your spectra
c     externally (eg, with IRAF's dispcor) and use Odl_syn = original
c     dl. The RebinSpec routine would then keep the spectrum untouched!
      call RebinSpec(0,NOl_obs,Ol_obs,Of_obs,Nl_syn,l_syn,f_aux1)

c     Rebin error & flag spectra too. 
c     **> ATT: THE REBINNING OF THE ERROR SPECTRUM SHOULD BE RE-THOUGHT!
c     (Binning adjacent pixels should reduce the noise, but my rebinning
c     routine simply "averages" it...). This is why I do NOT recommend
c     rebinning at all!
c     Rebinning the flag also changes the integer-flags to real
c     numbers... This, however, is not a problem: We simply flag-away
c     points with flag > 1.
!2BEDONE! ATT: Change here when/if flag-value-meaning is defined in arq_config!
      call RebinSpec(0,NOl_obs,Ol_obs,Oerr_obs,Nl_syn,l_syn,f_aux2)
      call RebinSpec(0,NOl_obs,Ol_obs,Oflag_obs,Nl_syn,l_syn,f_aux3)


c     OBS: This next step REDEFINES NOl_obs, Ol_obs(), Of_obs(),
c     Oerr_obs() & Oflag_obs()!!
      NOl_obs = Nl_syn
      do i=1,NOl_obs
!         write (*,'(i4,2x,4(1x,f8.3),10x,4(1x,f8.3))') i , 
!     &        Ol_obs(i) , Of_obs(i) , Oerr_obs(i) , Oflag_obs(i) ,
!     &        l_syn(i)  , f_aux1(i) , f_aux2(i)   , f_aux3(i)
         Ol_obs(i)    = l_syn(i)
         Of_obs(i)    = f_aux1(i)
         Oerr_obs(i)  = f_aux2(i)
         Oflag_obs(i) = f_aux3(i)

c     Mask/flag out points outside desired Olsyn_ini ==> Olsyn_fin range.
c     Also mask/flag out points with error <= 0! (Just a safety measure...)
         if ((Ol_obs(i).LT.Olsyn_ini).OR.
     &       (Ol_obs(i).GT.Olsyn_fin).OR.
     &       (Oerr_obs(i).LE.0.)) Oflag_obs(i) = 99.
!2BEDONE! ATT: Change here when/if flag-value-meaning is defined in arq_config!
      end do


c     Resample/Rebin base spectra to l_obs wavelength-scale.  I've
c     introduced a i_FastBC03_FLAG which greatly speeds up rebinning for
c     a BC03 base and dl = 1 and Ol_obs is within 3322-->9300 Angs.
c     i_FastBC03_FLAG is read from config file.
c     OBS: This step REDEFINES the Of_base() matrix!!
!
!     ATT: I've adapted the fast-rebin trick to work with Base-Rosa
!     spectra (br04*spec files), which cover the 3001 --> 6996 Angs
!     range. To use use fast-rebin with BR04 you must set
!     i_FastBC03_FLAG = 2!  This trick should be generalized in a more
!     elegant manner... (Cid@IAA - 30/May/2005)
!XXX! THIS SHOULD BE IMPROVED!
      do j=1,N_base
         do i=1,NOl_base
            f_aux1(i) = Of_base(j,i)
         end do

         call RebinSpec(i_FastBC03_FLAG,NOl_base,Ol_base,f_aux1,
     &        NOl_obs,Ol_obs,f_aux2)

         do i=1,NOl_obs
            Of_base(j,i) = f_aux2(i)
         end do
      end do


!Cid20Jan07-FIX! WARN if we will need to mask base edges to be safe
!     against edge effects in the kinematical filter! 
!     NOTE: I initially worried that masking my favorite norm-range
!     (4010-4060 A) in Pegase-HR models (which start at 4000 Ans) would
!     bring problems, but apparently not!
      if ((Ol_base(1).GT.lsyn_ini).OR.(Ol_base(NOl_base).LT.lsyn_fin))
     $     then
         write (*,'(a,a)') '@@@@@@@@@@@@@@@@@@@@@@@@@@@@> BASE ' , 
     $        'LIMITS WARNING <@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
         write (*,'(2(a,f9.1),a)') '@ Your base goes from lambda = ',
     $        Ol_base(1) ,' to ' , Ol_base(NOl_base) , ' Angs '
         write (*,'(2(a,f9.1),a)') '@ but the fit will run from    ' ,
     $        lsyn_ini , ' to ' , lsyn_fin ,
     $        ' Angs '
         write (*,'(a,f5.1,a)') '@ (including a +/-', dl_cushion ,
     $        ' Angs cushion for the kinematical filter)'
         write (*,'(2(a,f9.1),a)')
     $        '@ STARLIGHT will flag-out lambdas outside the ' ,
     $        Ol_base(1) + dl_cushion , ' to ' , 
     $        Ol_base(NOl_base) - dl_cushion , ' Angs range!'
         write (*,'(a,a)') '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
     $        ,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
      end if


c     Define a MASK-SPECTRUM f_mask() in the Ol_obs wavelength-scale.
c     The mask_weight array should be = 0 on the N_mask llow_mask ==>
c     lupp_mask windows which we do NOT want to model. If mask_weight >
c     0, this sets different/special weight for the corresponding
c     llow_mask ==> lupp_mask interval.  We also force the mask to be
c     zero outside the Olsyn_ini ==> Olsyn_fin interval.  Finally, we
c     mask-out points for which the (SDSS-inspired) quality flag is > 1.
c     These points are marked with f_mask = -2.
      do i=1,NOl_obs

         f_mask(i) = 1.

         do j=1,N_masks
            if ((Ol_obs(i).GE.llow_mask(j)).AND.
     &           (Ol_obs(i).LE.lupp_mask(j))) f_mask(i) = mask_weight(j)
         end do

         if ((Ol_obs(i).LT.Olsyn_ini).OR.
     &        (Ol_obs(i).GT.Olsyn_fin)) f_mask(i) = 0.

!Cid20Jan07-FIX! Introduced safety against people trying to fit data
!     within dl_cushion Angs of the base edges. The if below will
!     mask the first and last dl_cushon Angs from the base edges.
!     OBS: The possible/likely overlap with the if above is harmless!
         if ((Ol_obs(i).LT.(Ol_base(1)+dl_cushion)).OR.
     $        (Ol_obs(i).GT.(Ol_base(NOl_base)-dl_cushion))) then
            f_mask(i) = 0.
         end if

c     Accept points with flag <= 1! This means that we INCLUDE things
c     SDSS people call "emission lines" (flag = 1), but we've verified
c     these are not always real emission lines!
c     Cid&Jean@UFSC - 20/Nov/2003

!2BEDONE! ATT: Change here when/if flag-value-meaning is defined in arq_config!
!KEEP-FLAG-VALUES! This next line used to be: if (Oflag_obs(i).GT.1.) f_mask(i) = -2.
!	  I now keep the flag>1 values in f_mask. This is just to aid the identification
!	  of flagged points (which can mean sky-lines, hot pix, etc) in the output of
!	  STARLIGHT. Nothing else changes! But note that f_mask keeps the resampled
!	  flag, which will not be equal to the original when the data is rebinnedd..
!	  ElCid@Sanchhica - 10/Feb/2012
         if (Oflag_obs(i).GT.1.) f_mask(i) = -Oflag_obs(i)
      end do

c     OBS: After the steps above all arrays (Of_obs, Of_base & f_mask) 
c          are sampled at Ol_obs. (Ol_base is not used anymore).
c ***************************************************************************



c ***************************************************************************
c *                     PRE-PROCESSING II: NORMALIZATION                    *
c ***************************************************************************
c     Renormalize the observed spectrum by the MEDIAN flux between
c     llow_norm & lupp_norm, but the base spectra are normalized just by
c     the flux at l_norm! Normalizing the base at a single wavelength is
c     better (eg, it preserves normalization when applying the reddening
c     term to synthetic spectra) and should be Ok since the base has "no
c     noise".  The observed spectrum, however, has noise, so it is
c     better to normalize it by the mean or median flux in a window.  

c     Normalize Of_obs to MEDIAN flux between llow_norm & lupp_norm
c     Re-Normalize Oerr_obs too.  ATT: the median is computed regardless
c     of whether points are flagged or masked!!
!     ATT: NEW-Safety: If fobs_norm <= 0 then produce a silly output
!     file & abort the run!! The user should fix this before running.
!2BEDONE! ATT: fobs_norm computation IGNORES bad (flaged/masked) pixels!!
!         This is dangerous and could/should be changed one day...
      call Stats_F(NOl_obs,Ol_obs,Of_obs,llow_norm,lupp_norm,
     &     Of_ave,Of_sig,Of_med)
      fobs_norm = Of_med

      if (fobs_norm.LE.0.) then
         write (*,'(a,1p,e12.3)') '## ATT: fobs_norm = ' , fobs_norm
         write (*,'(a,a)') '## arq_obs = ' , arq_obs
         write (*,'(a)') '## Run aborted & silly file produced :('
         i_last_non_blank_char = 0
         do i=1,99
           if ((out_dir(i:i+1).EQ.' ').AND.(i_last_non_blank_char.EQ.0))
     &           i_last_non_blank_char = i - 1
        end do
        arq_aux = out_dir(1:i_last_non_blank_char) // arq_out
        open (unit=2,file=arq_aux,status='unknown')
        write (2,'(a,1p,e12.3)') '## ATT: fobs_norm = ' , fobs_norm
        write (2,'(a,a)') '## arq_obs = ' , arq_obs
        write (2,'(a)') '## Run aborted:('
        close(2)
        return
      end if

      do i=1,NOl_obs
         Of_obs(i)   = Of_obs(i) / fobs_norm
         Oerr_obs(i) = Oerr_obs(i) / fobs_norm
      end do

c     Normalize base spectra to their flux at l_norm.
c     OBS1: fbase_norm(j) is saved to compute M/L's & mass-fractions...
c     OBS2: fbase_norm(N_base) is meaningless if a power-law was added!
      do j=1,N_base
         do i=1,NOl_obs
            f_aux1(i) = Of_base(j,i)
         end do

         call calc_Fnorm(NOl_obs,Ol_obs,f_aux1,l_norm,fbase_norm(j))

         do i=1,NOl_obs
            Of_base(j,i) = f_aux1(i) / fbase_norm(j)
         end do
      end do
c ***************************************************************************



c ETC -ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c ETC               DEFINE Lobs_norm  (Cid@Lagoa - 07/Jan/2011)            ETC
c ETC -ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c     Compute the log(monochromatic-luminosity-@l_norm): Lobs_norm and
c     its log. The units are expected to be L_sun/Angs, as long as the
c     entries **flux_unit** and **LumDistInMpc** are Ok.  This was
c     introduced during the FIR-work, but is needed for any ETC-like
c     fit.  The physical flux-scale MUST be right for FIR/QHR/PHO-fits
c     (this is assumed in the related-routines), but is COMPLETELY
c     IRRELEVANT for the optical spectral fit itself!  In fact, up until
c     the very final output phase, ONLY the FIR/QHR/PHO routines know
c     and use the actual units of the observed optical spectrum.
      logLobs_norm = log10(fobs_norm * flux_unit) + 
     &     log10(4. * (4. * atan(1.))) +
     &     2. * log10(LumDistInMpc * 1.e6 * 3.086e18) -
     &     log10(3.826e33)
      Lobs_norm = 10.**logLobs_norm

c     Store logLobs_norm into variables * acessibles by the FIR/QHR/PHO
c     routines via the respective * data-commons.
      FIRlogLobs_norm = logLobs_norm
      QHR_logLobs_norm = logLobs_norm
      PHO_logLobs_norm = logLobs_norm
c ETC -ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC



c ***************************************************************************
c *      PRE-PROCESSING III: DEFINITION OF ERRORS/WEIGHTS & S/N ratios      *
c ***************************************************************************
c     Define errors & weights! Weight = 0 if f_mask <= 0 OR f_obs <=
c     f_cut. These points will thus be excluded from the synthesis.  For
c     other points, use INPUT error flux!  Points with mask_weight > 0
c     are given extra weight (= f_mask).
c     Nl_eff = the number of non-zero-weight points, ie, the effective
c              number of lambdas synthesized.
      Nl_eff = 0
      do i=1,NOl_obs
         if ((f_mask(i).LE.0.).OR.(Of_obs(i).LE.f_cut).OR.
     &        (Oerr_obs(i).LE.0.)) then
            err_obs(i)    = -99.
            weight_obs(i) = 0.
         else
            err_obs(i)    = Oerr_obs(i) / f_mask(i)
            weight_obs(i) = 1. / err_obs(i)
            Nl_eff        = Nl_eff + 1
         end if
      end do

c     Store Nl_eff onto NOl_eff ("Original Nl_eff"), since Nl_eff will
c     be modified if clipping is turned on.
      NOl_eff = Nl_eff


c     Apply a "WEIGHT-CONTROL-FILTER"! Compute mean and std dev of
c     weight_obs for the NOl_eff unmasked lambdas, and reset all weights
c     2-sigma above mean to weigh_obs = mean + 2 sigma! I've decided to
c     apply this "filter" as some of the SDSS spectra have very-low
c     errors on some lambdas, and these are then given a
c     disproportionately large weight in the synthesis.
c     Cid@lynx - 17/March/2004
c
c     I later realized that it is desirable to have more control over
c     this option. For instance, when err_obs() = constant (as in Joao's
c     data) but we want (through the mask file) to give more weight to a
c     particular window (say, Ca II K), we'd have wei_sig(i) = 0, which
c     would result in weight_obs(i) = wei_ave all accross the
c     spectrum... Hence, I introduced the parameter wei_nsig_threshold
c     in arq_config. Setting wei_nsig_threshold <= 0 will turn OFF the
c     weight filter. Otherwise, we censor points whose weights are >=
c     wei_nsig_threshold times larger than wei_sig.
c     Cid@Lagoa - 24/March/2005

!Cid24Mar07-FIX! Added counter of censored weights just for output purposes.
      n_censored_weights = 0
      wei_limit          = 0
      if (wei_nsig_threshold.GT.0.) then
         N_aux = 0
         do i=1,NOl_obs
            if (weight_obs(i).GT.0.) then 
               N_aux         = N_aux + 1
               f_aux1(N_aux) = weight_obs(i)
            end if
         end do

         call avevar(f_aux1,N_aux,wei_ave,wei_var)
         wei_sig   = sqrt(wei_var)
         wei_limit = wei_ave + wei_nsig_threshold * wei_sig

         do i=1,NOl_obs
            if (weight_obs(i).GT.wei_limit) then
               weight_obs(i)      = wei_limit
               n_censored_weights = n_censored_weights + 1
            end if
         end do
      end if


c     Compute the S/N ratio in the SN & normalization windows from the
c     the average-flux/rms-flux ratios (excluding masked points). Then
c     re-estimate these quantities using the average-error instead of
c     the rm-flux in the denominator ("SNerr_*"). The S/N in the
c     normalization window is just 1 over the noise flux in the
c     SN-window (a sensible, albeit silly, definition).  
c     OBS: ALL these S/N's are IRRELEVANT since we use an input error
c     spectrum! They are ONLY used for output purposes!
      call Stats4SN(NOl_obs,Ol_obs,Of_obs,weight_obs,
     &     llow_SN,lupp_SN,ave1,sig1)

c     ATT: For safety, if the rms in the llow_SN -> lupp_SN is 0
c     (god-knows-why...) then recompute it using ALL the spectrum! This
c     is just to avoid stupid NaN's in the output file! 
c     Cid@Lagoa - 08/Oct/2005
      if (sig1.LE.0.) then
         llow_SN = Ol_obs(1)
         lupp_SN = Ol_obs(NOl_obs)
         call Stats4SN(NOl_obs,Ol_obs,Of_obs,weight_obs,llow_SN,lupp_SN
     $        ,ave1,sig1)
      end if

      call Stats4SN(NOl_obs,Ol_obs,Oerr_obs,weight_obs,
     &     llow_SN,lupp_SN,ave2,sig2)
      SN_snwin      = ave1 / sig1
      SN_normwin    = 1.   / sig1
      SNerr_snwin   = ave1 / ave2
      SNerr_normwin = 1.   / ave2

      if (i_verbose.GT.0) then
         write (*,'(a,f8.3,2x,f8.3,a)')'    S/N     = ', 
     &     SN_snwin , SN_normwin , '     in S/N & norm. windows'
         write (*,'(a,f8.3,2x,f8.3,a)')'    S/N_err = ', 
     &     SNerr_snwin , SNerr_normwin , '     in S/N & norm. windows'
      end if
c ***************************************************************************



c ***************************************************************************
c *        PRE-PROCESSING IV: SET EXTINCTION LAW & FINAL DETAILS            *
c ***************************************************************************
c     Define reddening term = -0.4 * (A(lambda)/AV - A(l_norm)/AV).
c     The array red_term is passed by the common/Red_factor/ to routines
c     which need it!
c     OBS1: The user should apply his own correction for Galactic
c           extinction to the input observed spectrum before-hand...
!     (Removed the obsolete R_V parameter - Cid@Granada - 24/Jan/2010)
      call Calc_red_term(red_law_option,l_norm,NOl_obs,Ol_obs,red_term,
     $     q_norm)

c     Store O* variables onto proper names, passed by commons to other
c     routines.
      Nl_obs = NOl_obs
      do i=1,Nl_obs
         l_obs(i) = Ol_obs(i)
         f_obs(i) = Of_obs(i)
         do j=1,N_base
            f_base(j,i) = Of_base(j,i)
         end do
      end do
c ***************************************************************************


c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c +                DONE WITH INPUT & PRE-PROCESSING STEPS:)                 +
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+






c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c +                              SPECTRAL FIT                               +
c +                                                                         +
c +                                                                         +
c +   The pre-processing steps above produced the following main            +
c +   variables: l_obs() , f_obs() , weight_obs(), Of_base(,),              +
c +   f_base(,), age_base(), Z_base(), red_term(), N_base , Nl_obs,         +
c +   ind_YAV(), N_YAV_components. These are used in the spectral fit,      +
c +   which produces estimates of the following output variables: x(),      +
c +   AV_tot() [= AV + YAV], v0 & vd.                                       +
c +                                                                         +
c +   In this part we do the actual spectral fit. The steps are:            +
c +   (1) Initialize chains & set-up several details
c +   (2) First-Fits (FF), to make chains converge to a reasonable zone of  +
c +       parameter space. N_loops_FF annealing loops are done.
c +   (3) Identify & clip outliers.
c +   (4) Re-fit clipped spectrum. Only 1 loop.
c +   (5) Re-fit once more, this time with a more rigorous Burn-In 
c +       criterion. The chain & global-best parameters after this step 
c +       are stored for output...
c +   (6) Do a "Fine-tuning"-annealing-fit with a condense base.
c +   
c +   OBS: The code below also allows for sampling of the PDF, but this part
c +        part is not really ready ... (Jean should have used it...)
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


c ***************************************************************************
c * Initialize Kinematical Parameters v0 & vd (Cid/FXK@Lagoa - 14/Jan/2007) *
c ***************************************************************************
c     Set initial values for kinematical parameters v0 & vd (in km/s).
c     These are now read from the input for each galaxy!!
c     (Life would be much easier if we knew these beforehand...)
      v0_min = v0_ini
      vd_min = vd_ini

c     IsKinFixed = 0 allows v0 & vd to be fitted. IsKinFixed = 1 keep
c     them fixed. IsThisAnFXKFit is an auxiliar flag (associated to
c     common/FXK_fits/ which indicates that this is an FXK fit, so v0 &
c     vd will NEVER be changed wrt the input values, so the base
c     smoothing can be done only once! (In this case, I smooth both
c     f_base & Of_base).
      IsThisAnFXKFit = 0
      IsKinFixed     = 0

      if (FitOrSamplePDF_option.EQ.'FXK') then
         IsKinFixed     = 1
         IsThisAnFXKFit = 1

         if (i_verbose.GT.0) write (*,'(a,2(a,f9.3),a)')
     $        '[MAIN] FXK-fit: Smoothing base to' ,' v0 = ' , v0_min ,
     $        '  & vd = ' , vd_min , ' km/s... once and for all!'

         do i=1,Nl_obs
            l_aux(i)  = l_obs(i)
         end do
         Nl_aux = Nl_obs
         do j=1,N_base
            do i=1,Nl_obs
               f_aux1(i) = f_base(j,i)
            end do
            call GaussSmooth(Nl_aux,l_aux,f_aux1,v0_min,vd_min,
     &           N_int_Gauss,Nl_obs,l_obs,f_aux2)
            do i=1,Nl_obs
               f_base(j,i)  = f_aux2(i)
               Of_base(j,i) = f_aux2(i)
            end do
         end do
         if (i_verbose.GT.0) write (*,'(a)') '[MAIN] FXK-fit: ... done!'

      end if

c     Sanity check: Warn if initial v0 & vd are not within low->upp bounds!
      if ((v0_min.LT.v0_low).OR.(v0_min.GT.v0_upp).OR.
     $    (vd_min.LT.vd_low).OR.(vd_min.GT.vd_upp)) then
         write (*,*) ' [MAIN] Your initial v0/vd values exceed the' ,
     $        '_low -> _upp  values in the config file!! '
         write (*,*) '        I hope you know what you are doing!'
      end if
c ***************************************************************************



c ***************************************************************************
c * !!MEX!!                    Initialize Chains                            *
c ***************************************************************************
c     Define par_low & par_upp limits, where par =(x,AV,exAV,fn), in this
c     order!  Also define j_AV, j_exAV & j_fn indices used throughout!
c     OBS: If no extra AVs are being used (N_exAV = 0) then j_exAV = j_fn!
      N_par  = N_base + 1 + N_exAV + 1
      j_AV   = N_base + 1
      j_exAV = N_base + 2
      j_fn   = N_par

      do j=1,N_base
         par_low(j) = 0.
         par_upp(j) = fn_upp
      end do
      par_low(j_AV)  = AV_low
      par_upp(j_AV)  = AV_upp 
      do j=j_exAV,j_fn-1        !OBS: This loop only runs if Nex_AV > 0!
         j_aux = (j - j_exAV) + 1
         par_low(j) = exAV_low(j_aux)
         par_upp(j) = exAV_upp(j_aux)
      end do
      par_low(j_fn) = fn_low
      par_upp(j_fn) = fn_upp


c     Initialize N_chains with RANDOM initial parameters:
c     (x,AV,exAV,fn).  Store them onto ChainPar-matrix.  Note also that
c     Generate_PopVector returns a pop-vector normalized to unit sum, so
c     we renormalize it to Sum{x_j} = fn a posteriori.

!     ATT: To decrease initial chain energies, I'm using a more
!     constrained range for fn, AV & exAV so that chains do not start
!     from completely crazy positions allowed in the low ==> upp ranges!
!     Also, I force the 1st chain to start from x_j = 1/N_base, AV = exAV = 0.

!     ATT! Whereas in StarlightChains_v01.for the z*_low & *_upp values
!     for fn, AV & exAV were fixed, I now set them equal to the
!     corresponding input *_low & *_upp values and restrict to narrower
!     ranges them only if they exceed certain limits (0.9 -> 1.1 for fn;
!     0 -> 1 for AV). The exAV values are generated in a compressed range.
      zfn_low = fn_low
      zfn_upp = fn_upp
      if (fn_low.LT.0.9) zfn_low = 0.9
      if (fn_upp.GT.1.1) zfn_upp = 1.1

      zAV_low = AV_low
      zAV_upp = AV_upp
      if (AV_low.LT.0.0) zAV_low = 0.0
      if (AV_upp.GT.1.0) zAV_upp = 1.0

      ChainPar(j_fn,1) = 1.0
      do j=1,N_base
         ChainPar(j,1) = 1.0 / float(N_base)
      end do
      ChainPar(j_AV,1) = 0.0
      do j=j_exAV,j_fn-1        !OBS: This loop only runs if Nex_AV > 0!
         ChainPar(j,1) = 0.0
      end do

      do i_chain=2,N_chains
         ChainPar(j_fn,i_chain) = zfn_low + (zfn_upp - zfn_low)
     $        * ran2(idum)

         call Generate_PopVector(idum,N_base,xinit_max,x)
         do j=1,N_base
            ChainPar(j,i_chain) = x(j) * ChainPar(j_fn,i_chain)
         end do

         ChainPar(j_AV,i_chain)  = zAV_low + (zAV_upp - zAV_low)
     $        * ran2(idum)
*     ATT: Here I compress the initial exAV range to 1/3 of its allowed
*     range to avoid too crazy initial chain values...
         do j=j_exAV,j_fn-1     !OBS: This loop only runs if Nex_AV > 0!
            ChainPar(j,i_chain) = (par_low(j) + par_upp(j)) / 2. +
     $           (par_upp(j) - par_low(j)) * (ran2(idum) - 0.5) * 2./3.
         end do
      end do

* Screen check
      if (i_verbose.GT.0) then
         write (*,'(150a)') ('-',j=1,120)
         write (*,'(45x,a)') '- Initial Chain Parameters -'
         write (*,'(150a)') ('-',j=1,120)
         write (*,'(a,200(1x,i6))') 'Chain #    | par#',(j,j=1,N_par)
         do i_chain=1,N_chains
            write (*,'(a,i2,a,200(1x,f6.2))') 'Ch-Init ',i_chain ,
     $           ' | par=', (100*ChainPar(j,i_chain),j =1,N_base) ,
     $           (ChainPar(jj,i_chain),jj=j_AV,j_fn) 
         end do
         write (*,'(150a)') ('-',j=1,120)
      end if
c ***************************************************************************



c ***************************************************************************
c * !!MEX!!          Initialize help-arrays & other stuff                   *
c ***************************************************************************
c     Define arrays which say which component will be given an extra
c     help to converge and with what probability. Things are done
c     separately for the pop-vector & extinction parameters. These
c     arrays are passed by the /Scubi_TechnicalParameters/ common.
      do j=1,N_par
         ind_HelpParMove2Average(j)  = i_HelpPopVectorMove2Average
         prob_HelpParMove2Average(j) = prob_HelpPopVectorMove2Average
         if ((j.GE.j_AV).AND.(j.LT.j_fn)) then
            ind_HelpParMove2Average(j)  = i_HelpAVYAVMove2Average
            prob_HelpParMove2Average(j) = prob_HelpAVYAVMove2Average
         end if
      end do

c     Define Npar4GRtest = # of parameters used in the GR-convergence
c     (burn-in) test. The "- 1" is to exclude fn (the last entry in
c     the par vector)
      Npar4GRtest = N_par - 1

c     Set probability of picking the par_j with largest GR_R (used if
c     sampling is one-par-per-move), to help it move more & thus
c     hopefully converge faster... To turn-on this trick you must set
c     i_HelpParWithMaxR = 1. The value of prob_jRmax requires some
c     experimentation...  Note that prob_jRmax is actually set on the
c     input file. Here I just adjust it in the case of small N_par...
*      prob_jRmax = max((5.0 / float(Npar4GRtest)),prob_jRmax)
*     !!MEX!! Adjusted!! (Cid@Lagoa - 20/April/2007)
      z1 = max((5.0 / float(Npar4GRtest)),prob_jRmax)
      if ((z1.LE.0.5).AND.(prob_jRmax.LT.0.5)) prob_jRmax = z1
      if (prob_jRmax.GT.0.9) prob_jRmax = 0.5
c ***************************************************************************



c ***************************************************************************
c *                    FF: First Fits (w/full base)                         *
c ***************************************************************************
c     Here we perform the "First Fits" (FF). The goal is to take all
c     chains to a reasonable region of parameter space, using simulated
c     annealing plus a Gelman & Rubin burn-in verification criterion.

c     Set eff_IDEAL & IsGRTestHard (passed by common), N_sim &
c     Nmax_steps for First-Fits.  
      eff_IDEAL     = eff_IDEAL_FF
      IsGRTestHard  = IsGRTestHard_FF
      N_sim_FF      = fN_sim_FF * N_par
      Nmax_steps_FF = fNmax_steps_FF * N_par

!     Screen output
      if (i_verbose.GT.0) then
         write (*,'(80a)') ('*',j=1,80)
         write (*,'(a)') '* STARTING FIRST FITS ...'
         write (*,'(a,i3)') '* N_loops  = ' , N_loops_FF
         write (*,'(a,i3)') '* N_chains = ' , N_chains
         write (*,'(a,1p,e8.2,a,e8.2)') '* Temp     = ' , Temp_ini_FF ,
     $        ' ==> ' , Temp_fin_FF
         write (*,'(a,1p,e8.2,a,e8.2,0p,a,i2,a)') '* GR-R     = ' ,
     $        R_ini_FF ,' ==> ' , R_fin_FF , '  Method = ' ,
     $        IsGRTestHard , '  (0/1=Soft/Hard)'
         write (*,'(80a)') ('*',j=1,80)
      end if

!Cid29Dec06-FIX! Introduced time-counters (to help in config-tests)
      i_ini_time_FF    = time()
      w_ini_FF_TotTime = etime(willy_tiempo)
      w_ini_FF_UsrTime = willy_tiempo(1)
      w_ini_FF_SysTime = willy_tiempo(2)


c     Go for First Fits...
 123  call MEX_Fit_Scubidu(i_verbose_anneal,idum,N_loops_FF,v0_low
     $     ,v0_upp,vd_low,vd_upp,x_min,AV_min,exAV_min,v0_min,vd_min
     $     ,Ntot_FF,eff_FF,Nl_eff,iMassAgeZ_verbose,N_par,ChainPar
     $     ,N_chains,Temp_ini_FF,Temp_fin_FF,fini_eps_FF,ffin_eps_FF
     $     ,R_ini_FF,R_fin_FF,N_sim_FF,Nmax_steps_FF,i_RestartChains_FF
     $     ,Gpar_min,GE_min,Gpar_ave,Cpar_min,CEne_min,IsKinFixed
     $     ,IsBurInOver_FF)


!Cid29Dec06-FIX! Introduced time-counters (to help in config-tests)
      i_fin_time_FF    = time()
      w_fin_FF_TotTime = etime(willy_tiempo)
      w_fin_FF_UsrTime = willy_tiempo(1)
      w_fin_FF_SysTime = willy_tiempo(2)

c     Update global steps counter
      Nglobal_steps = Nglobal_steps + Ntot_FF
c ***************************************************************************



c ***************************************************************************
c *                        Clip-outliers & Re-Fit                           *
c ***************************************************************************
c     Clip deviant points in the spectrum! Use best model to compute
c     residuals...
      call MEX_Update_AV_array(N_par,N_base,AV_min,exAV_min,AV_tot)
      call Kcalc_Fsyn(Nl_obs,x_min,AV_tot,v0_min,vd_min,f_syn)
      call Clip_Spectrum(clip_method_option,f_syn,sig_clip_threshold
     $     ,f_mask,NOl_eff,Nl_eff,Ntot_cliped,iCLIPBUG_flag)


c     If clipping actually cliped something, then run a single re-fit
c     with the SAME conditions as in the last loop above. Otherwise
c     no-refitting is done!
      if (Ntot_cliped.LE.0) then

         if (i_verbose.GT.0) then
            write (*,'(a,i5,a,f3.1,a,a,a)') '* CLIP>  Clipped ' ,
     $           Ntot_cliped,' points with > ' , sig_clip_threshold
     $           ,' sigma residuals with method = ' , clip_method_option
     $           , '.  No re-fitting needed!'
         end if

      else

         Temp_ini     = Temp_fin_FF
         Temp_fin     = Temp_fin_FF
         fini_eps     = ffin_eps_FF
         ffin_eps     = ffin_eps_FF
         R_ini        = R_fin_FF
         R_fin        = R_fin_FF
         N_sim        = N_sim_FF
         Nmax_steps   = Nmax_steps_FF
         N_loops      = 1
         IsGRTestHard = IsGRTestHard_FF

!Cid22Jan07-FIX! i_RestartChains was UNDEFINED here! (The
!     implication is that chains where not restarted to their best
!     values before going to the next stage).
         i_RestartChains = i_RestartChains_FF

!     Screen output
         if (i_verbose.GT.0) then
            write (*,'(80a)') ('*',j=1,80)
            write (*,'(a)') '* STARTING RE-FIT AFTER CLIPPING ...'
            write (*,'(a,i5,a,f3.1,a,a)') '* Cliped ' , Ntot_cliped
     $           ,' points with > ' , sig_clip_threshold
     $           ,' sigma residuals with method = ' , clip_method_option
            write (*,'(a,i3)') '* N_loops = ' , N_loops
            write (*,'(a,1p,e8.2,a,e8.2)') '* Temp    = ' , Temp_ini ,
     $           ' ==> ' , Temp_fin
            write (*,'(a,1p,e8.2,a,e8.2,0p,a,i2,a)') '* GR-R    = ' ,
     $           R_ini ,' ==> ' , R_fin , '  Method = ' , IsGRTestHard ,
     $           '  (0/1=Soft/Hard)'
            write (*,'(80a)') ('*',j=1,80)
         end if

!Cid29Dec06-FIX! Introduced time-counters (to help in config-tests)
         i_ini_time_RF    = time()
         w_ini_RF_TotTime = etime(willy_tiempo)
         w_ini_RF_UsrTime = willy_tiempo(1)
         w_ini_RF_SysTime = willy_tiempo(2)


c     Go for re-fit
         call MEX_Fit_Scubidu(i_verbose_anneal,idum,N_loops,v0_low
     $        ,v0_upp,vd_low,vd_upp,x_min,AV_min,exAV_min,v0_min,vd_min
     $        ,Ntot_ClipReFit,eff_ClipReFit,Nl_eff,iMassAgeZ_verbose
     $        ,N_par,ChainPar,N_chains,Temp_ini,Temp_fin,fini_eps
     $        ,ffin_eps,R_ini,R_fin,N_sim,Nmax_steps,i_RestartChains
     $        ,Gpar_min,GE_min,Gpar_ave,Cpar_min,CEne_min,IsKinFixed
     $        ,IsBurInOver_ClipReFit)


!Cid29Dec06-FIX! Introduced time-counters (to help in config-tests)
         i_fin_time_RF    = time()
         w_fin_RF_TotTime = etime(willy_tiempo)
         w_fin_RF_UsrTime = willy_tiempo(1)
         w_fin_RF_SysTime = willy_tiempo(2)

c     Update global steps counter
         Nglobal_steps = Nglobal_steps + Ntot_ClipReFit

      end if
c ***************************************************************************



c ***************************************************************************
c *                            Burn-In Re-fit                               *
c ***************************************************************************
c     Re-fit (clipped) spectrum using the T = 1 & IsGRTestHard = 1
c     (meaning that every individual parameter must satisfy the GR-test)
c     & GR_R-limit <= 1.2. This is to ensure we're over the BURN-IN
c     phase before proceeding to further sampling or minimization. 
c     (OBS: Actually, the Hard/Soft GR-test and the R-Threshold used here
c     are define by the config values of IsGRTestHard_BurnIn &
c     R_BurnIn!)
c     If T = 1 & IsGRTestHard = 1 in the last loop above, then this
c     refit is unnecessary.
      if ((Temp_fin.LE.1.).AND.(IsGRTestHard.EQ.1).AND.(R_fin.LE.
     $     R_BurnIn)) then
         
         IsBurInOver_BurnIn = 1
         if (i_verbose.GT.0) then
            write (*,'(a,a)') '* No T = 1 & Hard-GR-test refit needed.'
     $           , ' Burn-in is already over!'
         end if

      else

!Cid11Feb07-FIX! Previsouly (StarlightChains_v02) I did the Burn-In 
!     at T = 1, seting Temp_ini = Temp_fin = 1.0. Now I do it at T =
!     Temp_fin_FF to be more general. All my fits with
!     StarlightChains_v02 used Temp_fin_FF = 1, so results should be
!     identical.
         Temp_ini     = Temp_fin_FF
         Temp_fin     = Temp_fin_FF
         fini_eps     = ffin_eps_FF
         ffin_eps     = ffin_eps_FF
         R_ini        = R_BurnIn
         R_fin        = R_BurnIn
         N_sim        = N_sim_FF
         Nmax_steps   = Nmax_steps_FF
         N_loops      = 1
         IsGRTestHard = IsGRTestHard_BurnIn

!Cid22Jan07-FIX! i_RestartChains was UNDEFINED here! (The
!     implication is that chains where not restarted to their best
!     values before going to the next stage).
         i_RestartChains = i_RestartChains_FF

!     Screen output
         if (i_verbose.GT.0) then
            write (*,'(80a)') ('*',j=1,80)
            write (*,'(a)') '* STARTING FINAL BURN-IN LOOP...'
            write (*,'(a,i3)') '* N_loops = ' , N_loops
            write (*,'(a,1p,e8.2,a,e8.2)') '* Temp    = ' , Temp_ini ,
     $           ' ==> ' , Temp_fin
            write (*,'(a,1p,e8.2,a,e8.2,0p,a,i2,a)') '* GR-R    = ' ,
     $           R_ini ,' ==> ' , R_fin , '  Method = ' , IsGRTestHard ,
     $           '  (0/1=Soft/Hard)'
            write (*,'(80a)') ('*',j=1,80)
         end if

!Cid29Dec06-FIX! Introduced time-counters (to help in config-tests)
         i_ini_time_BI    = time()
         w_ini_BI_TotTime = etime(willy_tiempo)
         w_ini_BI_UsrTime = willy_tiempo(1)
         w_ini_BI_SysTime = willy_tiempo(2)


c     Go for it 
         call MEX_Fit_Scubidu(i_verbose_anneal,idum,N_loops,v0_low
     $        ,v0_upp,vd_low,vd_upp,x_min,AV_min,exAV_min,v0_min,vd_min
     $        ,Ntot_BurnIn,eff_BurnIn,Nl_eff,iMassAgeZ_verbose,N_par
     $        ,ChainPar,N_chains,Temp_ini,Temp_fin,fini_eps,ffin_eps
     $        ,R_ini,R_fin,N_sim,Nmax_steps,i_RestartChains,Gpar_min
     $        ,GE_min,Gpar_ave,Cpar_min,CEne_min,IsKinFixed
     $        ,IsBurInOver_BurnIn)


!Cid29Dec06-FIX! Introduced time-counters (to help in config-tests)
         i_fin_time_BI    = time()
         w_fin_BI_TotTime = etime(willy_tiempo)
         w_fin_BI_UsrTime = willy_tiempo(1)
         w_fin_BI_SysTime = willy_tiempo(2)

c     Update global steps counter
         Nglobal_steps = Nglobal_steps + Ntot_BurnIn

      end if
c ***************************************************************************



c ***************************************************************************
c *                   Sample PDF option (4 Jean but never used...)          *
c ***************************************************************************
c     If FitOrSamplePDF_option = 'PDF', we'll sample the Probability
c     Distribution Function (PDF), otherwise proceed to optimization
c     with a shortened base. Instead of sampling the PDF of the actual
c     synthesis parameters we store in arq_out the values of some
c     "condensed-parameters" (such as mean age, mean Z, etc) for every
c     step of every chain.

      if (FitOrSamplePDF_option.EQ.'PDF') then 

c     To sample the PDF it is first necessary to fix all technical
c     sampling knobs/tricks, such that the proposal distribution does
c     NOT vary while sampling. In other words, stop helping parameters
c     move towards convergence, stop adjusting step-size to achieve
c     desired efficiency, etc. (All MCMC stuff I´ve read says that
c     changing/adapting the proposal distribution "on-the-fly" is no
c     good while sampling the PDF.)
         do j=1,N_par
            ind_HelpParMove2Average(j)  = 0.
            prob_HelpParMove2Average(j) = 0.
         end do

         i_HelpParWithMaxR           = 0
         prob_jRmax                  = 0.
         i_UpdateAVYAVStepSeparately = 0
         i_UpdateEps                 = 0
         i_UpdateAlpha               = 0
         Falpha                      = 1.0           
         i_MoveOneParOnly            = 1
         i_HelpAVYAVMove2Average     = 0
         prob_HelpAVYAVMove2Average  = 0
         i_RestartChains             = 0

c     Now fix, kinematics, T = 1 & step-sizes. 
         IsKinFixed   = 1
         Temp_ini     = 1.0
         Temp_fin     = 1.0
         fini_eps     = ffin_eps_FF
         ffin_eps     = ffin_eps_FF
         N_loops      = 1

c     To determine the number of samples, we first run Fit_Scubidu until
c     a Freeze-In condition is fullfilled, and figure out the number of
c     steps from the length of the Freeze-In period. Here we define the
c     GR-R criterion for the Freeze-in loop.
         R_ini        = R_BurnIn
         R_fin        = R_BurnIn
         IsGRTestHard = 1
         N_sim        = N_sim_FF
         Nmax_steps   = 2.e6

!     Screen output
         if (i_verbose.GT.0) then
            write (*,'(80a)') ('*',j=1,80)
            write (*,'(a)') '* STARTING FREEZE-IN LOOP...'
            write (*,'(a,1p,e9.3)') '* alpha      = ' , alpha
            write (*,'(a,i9)')      '* N_sim      = ' , N_sim
            write (*,'(a,i9)')      '* Nmax_steps = ' , Nmax_steps
            write (*,'(a,1p,e8.2,a,e8.2,0p,a,i2,a)') '* GR-R       = ' ,
     $           R_ini ,'==> ' , R_fin , '  Method = ' , IsGRTestHard ,
     $           '  (0/1=Soft/Hard)'
            write (*,'(80a)') ('*',j=1,80)
         end if


c     Initial sampling to define Freeze-in length: N_steps_FreezeIn
         call MEX_Fit_Scubidu(i_verbose_anneal,idum,N_loops,v0_low
     $        ,v0_upp,vd_low,vd_upp,x_min,AV_min,exAV_min,v0_min,vd_min
     $        ,Ntot_FreezeIn,eff_FreezeIn,Nl_eff,iMassAgeZ_verbose,N_par
     $        ,ChainPar,N_chains,Temp_ini,Temp_fin,fini_eps,ffin_eps
     $        ,R_ini,R_fin,N_sim,Nmax_steps,i_RestartChains,Gpar_min
     $        ,GE_min,Gpar_ave,Cpar_min,CEne_min,IsKinFixed
     $        ,IsBurInOver_FreezeIn)

         N_steps_FreezeIn = float(Ntot_FreezeIn) / float(N_chains)

c     Now set a ridiculously small value for GR-R, such that we control
c     the number of steps by Nmax_steps, which is defined as twice
c     the freeze in length (or 1e6, whichever is smaller).
         R_ini        = 0.1
         R_fin        = 0.1
         Nmax_steps   = min((2. * N_steps_FreezeIn),1.e5)
                                !JEAN: Limit/Define limit file-sizes here!
                                !JEAN: You may try 2, 3...10 times the 
                                !      Freeze in length!

!     Screen output
         if (i_verbose.GT.0) then
            write (*,'(80a)') ('*',j=1,80)
            write (*,'(a)') '* STARTING PDF SAMPLING...'
            write (*,'(a,i9)') '* N_FreezeIn = ' , N_steps_FreezeIn
            write (*,'(a,i9)') '* Nmax_steps = ' , Nmax_steps
            write (*,'(a,a)')  '* Steps saved in ' , arq_out
            write (*,'(80a)') ('*',j=1,80)
         end if


c     Open output steps-file (communicated though i_UnitScubi to actual
c     sampling routine) & Go for it: sample PDF!
         i_last_non_blank_char = 0
         do i=1,99
           if ((out_dir(i:i+1).EQ.' ').AND.(i_last_non_blank_char.EQ.0))
     &           i_last_non_blank_char = i - 1
         end do
         arq_aux = out_dir(1:i_last_non_blank_char) // arq_out

         i_UnitScubi = 2
         open (unit=i_UnitScubi,file=arq_aux,status='unknown')

         write (2,'(93(a1))') ('#',j=1,93)
         write (2,'(a,a)') '## OUTPUT of PANcMExStarlight_v03.for ',
     &        ' - SAMPLE-PDF option! - [ElCid@Sanchica 10/Feb/2012] ##'
         write (2,'(93(a1))') ('#',j=1,93)
         write (2,'(a)') '#'
         write (2,'(a,a)')  '# arq_obs     ' , arq_obs
         write (2,'(a,a)')  '# arq_base    ' , arq_base
         write (2,'(a,a)')  '# arq_config  ' , arq_config
         write (2,'(a,a)')  '# red_law     ' , red_law_option
         write (2,'(a,i8)') '# N_base      ' , N_base
         write (2,'(a,i8)') '# N_par       ' , N_par
         write (2,'(a,i8)') '# N_chains    ' , N_chains
         write (2,'(a,i8)') '# N_sim       ' , N_sim
         write (2,'(a,i8)') '# N_FreezeIn  ' , N_steps_FreezeIn
         write (2,'(a,i8)') '# Nmax_steps  ' , Nmax_steps
         write (2,'(a)') '#'
         write (2,'(a,a)') '# i_chain chi2 & 100 X par(1)...par(N_par)'

         call MEX_Fit_Scubidu(i_verbose_anneal,idum,N_loops,v0_low
     $        ,v0_upp,vd_low,vd_upp,x_min,AV_min,exAV_min,v0_min,vd_min
     $        ,Ntot_SamplePDF,eff_SamplePDF,Nl_eff,iMassAgeZ_verbose
     $        ,N_par,ChainPar,N_chains,Temp_ini,Temp_fin,fini_eps
     $        ,ffin_eps,R_ini,R_fin,N_sim,Nmax_steps,i_RestartChains
     $        ,Gpar_min,GE_min,Gpar_ave,Cpar_min,CEne_min,IsKinFixed
     $        ,IsBurInOver_SamplePDF)

         close (i_UnitScubi)
         i_UnitScubi = 0

c     We'll not optimize the solution, since the goal here was to Sample
c     the PDF). Hence we're done!
         return

      end if
c ***************************************************************************



c ***************************************************************************
c *                         Condensed-Base Fits                             *
c ***************************************************************************
c     If FitOrSamplePDF_option is not = 'PDF' we're thrown here, where
c     we'll try to find an "optimal" solution by (1) condensing the base
c     to "relevant" elements and (2) annealing down to low T to find a
c     best-chi2 solution.

!Cid24Mar07-FIX!  In previous versions, EX0s would run even with
!     N_loops_EX0s = 0! This is fixed now. But then I need to reset the 
!     things below, just in case N_loops_EX0s = 0.
      NEX0s_base = N_base
      do j=1,N_par
         EX0sGpar_min(j) = Gpar_min(j)
      end do
      EX0sv0_min        = v0_min
      EX0svd_min        = vd_min
      v0_min_BeforeEX0s = v0_min
      vd_min_BeforeEX0s = vd_min
      idt_EX0           = 0
      wdt_EX0_TotTime   = 0
      wdt_EX0_UsrTime   = 0
      wdt_EX0_SysTime   = 0


      if (N_loops_EX0s.GE.1) then

c     Store initial Chain-positions in EX0sChainPar (since I may want to
c     save the last ChainPar from the Burn-In loop above)
         do i_chain=1,N_chains
            do j=1,N_par
               EX0sChainPar(j,i_chain) = ChainPar(j,i_chain)
            end do
         end do

c ---------------------------------------------------------------------------
c     Define "x_Test". This is the pop-vector which will be used to
c     decide who stays in the condensed base and who is eliminated! If
c     EXOs_PopVector_option = 'AVE' we'll use x_Test = global average
c     pop-vector, which is a CONSERVATIVE strategy, since it has less
c     x_j = 0 components than any best-solution pop vector. For
c     EXOs_PopVector_option = 'MIN' we'll use x_Test = global best x so
c     far, which may be a better thing to do if we wanna refine this
c     minimum (This option should also produce shorter bases).

*FIR* (OBS: Done while implementing FIRc-stuff, but independent from it!)
*     Input flag (read from arq_config) iEXOs_IsAve_x_of_lambda decides
*     whether the "x_Test" above (which is either the AVE or the MIN
*     pop-vector) is maintained at l_norm (iEXOs_IsAve_x_of_lambda = 0)
*     or averaged over all lambdas (iEXOs_IsAve_x_of_lambda) before
*     calling the EXO-routine. Notice that x is just used to decide who
*     is relevant or not. Using iEXOs_IsAve_x_of_lambda = 0 you will
*     base this decision on the l_norm pop-fractions, whereas with
*     iEXOs_IsAve_x_of_lambda = 1 all lambdas in the fit (weigthing by
*     weight_obs(lambda)) will be considered. The latter option sounds
*     better, **specially for MEx-fits** (where one component can contribute
*	  little at l_norm but much more in other lambdas).
*
*     Cid@Granada - 31/Jan/2010
*
*FIX* There were 2 BUGs here: (1) The iEXOs_IsAve_x_of_lambda used x_min, even if
*     EXOs_PopVector_option = 'AVE'. (2) The AV_tot array needed in AverageXoverLambda
*     was NOT defined! I have now fixed both bugs and reorganized this part more compactly. 
*
*     ElCid@Sanchica - 17/Oct/2011
         if (EXOs_PopVector_option.EQ.'AVE') then 
          call MEX_Unpack_par2xAVYAVfn(N_par,N_base,Gpar_ave,
     &         x,AV,exAV,fn)
         else if (EXOs_PopVector_option.EQ.'MIN') then 
          call MEX_Unpack_par2xAVYAVfn(N_par,N_base,Gpar_min,
     &         x,AV,exAV,fn)
         end if
         call MEX_Update_AV_array(N_par,N_base,AV,exAV,AV_tot)

         if (iEXOs_IsAve_x_of_lambda.EQ.1) then
            call AverageXoverLambda(x,AV_tot,aux3)
            do j=1,N_base
               x(j) = aux3(j)
            end do
         end if
c ---------------------------------------------------------------------------


c     Set eff_IDEAL & IsGRTestHard (passed by common), N_sim &
c     Nmax_steps for EX0s-Fits
         eff_IDEAL       = eff_IDEAL_EX0s
         IsGRTestHard    = IsGRTestHard_EX0s
         N_sim_EX0s      = fN_sim_EX0s * N_par
         Nmax_steps_EX0s = fNmax_steps_EX0s * N_par

!     Screen output
         if (i_verbose.GT.0) then
            write (*,'(80a)') ('*',j=1,80)
            write (*,'(a)') '* STARTING EX0s FITS ...'
            write (*,'(a,a)') 
     $           '* EX0_method     =   ', EXOs_method_option

            if (iEXOs_IsAve_x_of_lambda.EQ.1) then
               write (*,'(a)') 
     &              '* x-vector for EX0s will be <x(lambda)>!'
            else if (iEXOs_IsAve_x_of_lambda.EQ.0) then
               write (*,'(a,a,f10.2)') 
     &              '* x-vector for EX0s will be <x(l_norm)>,' ,
     &              ' with l_norm = ',l_norm
            end if

            write (*,'(a,f7.3)')
     $           '* EX0s_Threshold = ' , EX0s_Threshold 
            write (*,'(a,a)')
     $           '* EXOs_PopVector =     ' , EXOs_PopVector_option
            write (*,'(a,i3)') '* N_loops = ' , N_loops_EX0s
            write (*,'(a,1p,e8.2,a,e8.2)') '* Temp    = ',Temp_ini_EX0s,
     $           ' ==> ' , Temp_fin_EX0s
            write (*,'(a,1p,e8.2,a,e8.2,0p,a,i2,a)')  '* GR-R    = ' ,
     $           R_ini_EX0s ,' ==> ' , R_fin_EX0s , '  Method = ' ,
     $           IsGRTestHard , '  (0/1=Soft/Hard)'
            write (*,'(80a)') ('*',j=1,80)
         end if

!Cid18Feb07-FIX! Define minimum acceptable base-size!
         N_base_min = nint(fEX0_MinBaseSize * N_base)
         if (N_base_min.LT.3) N_base_min = 3

!Cid29Dec06-FIX! Introduced time-counters (to help in config-tests)
         i_ini_time_EX0    = time()
         w_ini_EX0_TotTime = etime(willy_tiempo)
         w_ini_EX0_UsrTime = willy_tiempo(1)
         w_ini_EX0_SysTime = willy_tiempo(2)

c     Go for it!
         call MEX_EX0s_Fit_Scubi(i_verbose_anneal,idum,N_loops_EX0s
     $        ,v0_low,v0_upp,vd_low,vd_upp,x_min,AV_min,exAV_min
     $        ,EX0sv0_min,EX0svd_min,Ntot_EX0s,eff_NEX0s,Nl_eff
     $        ,iMassAgeZ_verbose,N_par,EX0sChainPar,N_chains
     $        ,Temp_ini_EX0s,Temp_fin_EX0s,fini_eps_EX0s,ffin_eps_EX0s
     $        ,R_ini_EX0s,R_fin_EX0s,N_sim_EX0s,Nmax_steps_EX0s
     $        ,i_RestartChains_EX0s,EX0sGpar_min,EX0sGE_min
     $        ,EX0sGpar_ave,EX0sCpar_min,EX0sCEne_min,x
     $        ,EXOs_method_option,EX0s_Threshold,N_base_min,IsKinFixed
     $        ,IsBurInOver_EX0s,NEX0s_base,IsScaleNstepsInEX0sFits
     $        ,IsNoKin4LargeBaseInEX0sFits
     $        ,frac_NoKin4LargeBaseInEX0sFits)


!Cid24Mar07-FIX!  Store EX0s_v0 & vd on v0_min & vd_min.
         v0_min = EX0sv0_min
         vd_min = EX0svd_min

!Cid29Dec06-FIX! Introduced time-counters (to help in config-tests)
         i_fin_time_EX0    = time()
         w_fin_EX0_TotTime = etime(willy_tiempo)
         w_fin_EX0_UsrTime = willy_tiempo(1)
         w_fin_EX0_SysTime = willy_tiempo(2)

c     Update global steps counter
         Nglobal_steps = Nglobal_steps + Ntot_EX0s

      end if
c ***************************************************************************


c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c +                       DONE WITH SPECTRAL FITTING:)                      +
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+





c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c +                     FITS WITH JUST ONE COMPONENT                        +
c +                                                                         +
c +    ~Quick fits with individual f_base spectra to see which individual   +
c +    base component would best fit the observed spectrum!                 +
c +    (These fits do NOT interfere with the composite fits done above!)    +
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

c ***************************************************************************
c *             Single-Component-fits - Cid@Lagoa - 20/Jan/2007             *
c ***************************************************************************
c     Just for the "fun" of it, before closing STARLIGHT will now fit
c     f_obs with each of the N_base individual SSPs in the base. With
c     just one component, there are only 2 parameters: AV and the
c     normalizatior factor. The best fit values are saved as SSP_AV()
c     and x_SSP() arrays, along with the corresponding SSP_chi2() and
c     SSP_adev() arrays. Minimization in SSP_AV is done with a Golden
c     Search, using FUNCTION F_SSP_chi2_AV(AV_SSP), and SSP_x is fitted
c     analytically inside that function.

c     The kinematical parameters used are the best ones at this stage,
c     and not changed from SSP to SSP fit. This exercise may prove
c     useful for, say, elliptical galaxies, to find the best
c     single-population fit. (Indeed, this option was inspired during
c     Ignacio de la Rosa's visit to Florianopolis on Jan/2007.)

c     Also useful to fit star-clusters, in which case it's a good idea
c     to chose an FXK fit and config parameters which minimize the time
c     spent on the fit-pop-mixtures part of the code (like N_loops_FF =
c     1 at high temperature, R-limits, no clipping , etc...), such that
c     the code gets to this point ASAP.

c     Cid@Lagoa - 20/Jan/2007

*ETC* ATT: These SSP-fits do NOT include FIR/QHR/PHO info, but that's Ok!

      do j=1,N_base

c     Apply kinematical filter to the j'th SSP, unless this is a FXK fit
c     (in which case f_base is already smoothed & shifted). The
c     resulting spectrum is stored in f_SSP().
         if (FitOrSamplePDF_option.EQ.'FXK') then
            do i=1,Nl_obs
               f_SSP(i) = f_base(j,i)
            end do
         else
            Nl_aux = Nl_obs
            do i=1,Nl_obs
               l_aux(i)  = l_obs(i)
               f_aux1(i) = f_base(j,i)
            end do
            call GaussSmooth(Nl_aux,l_aux,f_aux1,v0_min,vd_min,
     &           N_int_Gauss,Nl_obs,l_obs,f_SSP)
         end if

c     Golden Search for best SSP_AV (best norm.  factor x_SSP computed
c     analytically in F_SSP_chi2_AV and returned via common). Notice
c     tolerance = 1/1000, more than enough for this purpose!
         ax  = AV_low
         bx  = AV_upp
         tol = 1.e-3
         call mnbrak(ax,bx,cx,fa,fb,fc,F_SSP_chi2_AV)
!v05! ATT! "cfx" BUG-DETECTED BY CID&ROSA@Lagoa - 10/May/2008!!
!was this:  chi2_AV_Gold = golden(ax,bx,cxf,F_SSP_chi2_AV,tol,AV_Gold)
         chi2_AV_Gold = golden(ax,bx,cx,F_SSP_chi2_AV,tol,AV_Gold)
c     Do not allow SSP_AV outside the AV_low --> AV_upp limits.
         if (AV_Gold.LT.AV_low) AV_Gold = AV_low
         if (AV_Gold.GT.AV_upp) AV_Gold = AV_upp
         SSP_chi2(j) = F_SSP_chi2_AV(AV_Gold)
         SSP_x(j)    = x_SSP
         SSP_AV(j)   = AV_Gold

         do i=1,Nl_obs
            f_SSP(i) = 10.**(AV_Gold * red_term(i)) * x_SSP * f_SSP(i)
         end do
         SSP_adev(j) = F_adev(f_SSP)

!v05! Storing best-fit SSP spectrum
         if (j.EQ.1) Best_SSP_chi2 = SSP_chi2(j)
         if (SSP_chi2(j).LE.Best_SSP_chi2) then
            Best_SSP_chi2  = SSP_chi2(j)
            index_Best_SSP = j
            do i=1,Nl_obs
               Best_f_SSP(i) = f_SSP(i)
            end do
         end if

      end do
c ***************************************************************************
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+





c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c +                        POST-PROCESSING & OUTPUT                         +
c +                                                                         +
c +   Now that we have the final model, do a few final post-processing      +
c + (like converting flux to mass fractions) & produce output file arq_out. +
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

c ***************************************************************************
c *                              OUTPUT                                     *
c ***************************************************************************
!Cid29Dec06-FIX! Introduced time-counters (to help in config-tests)
      i_fin_time      = time()
      idt_all         = i_fin_time     - i_ini_time
      idt_FF          = i_fin_time_FF  - i_ini_time_FF
      idt_RF          = i_fin_time_RF  - i_ini_time_RF
      idt_BI          = i_fin_time_BI  - i_ini_time_BI
      idt_EX0         = i_fin_time_EX0 - i_ini_time_EX0
      w_fin_TotTime   = etime(willy_tiempo)
      w_fin_UsrTime   = willy_tiempo(1)
      w_fin_SysTime   = willy_tiempo(2)
      wdt_TotTime     = w_fin_TotTime     - w_ini_TotTime
      wdt_UsrTime     = w_fin_UsrTime     - w_ini_UsrTime
      wdt_SysTime     = w_fin_SysTime     - w_ini_SysTime
      wdt_FF_TotTime  = w_fin_FF_TotTime  - w_ini_FF_TotTime
      wdt_FF_UsrTime  = w_fin_FF_UsrTime  - w_ini_FF_UsrTime
      wdt_FF_SysTime  = w_fin_FF_SysTime  - w_ini_FF_SysTime
      wdt_RF_TotTime  = w_fin_RF_TotTime  - w_ini_RF_TotTime
      wdt_RF_UsrTime  = w_fin_RF_UsrTime  - w_ini_RF_UsrTime
      wdt_RF_SysTime  = w_fin_RF_SysTime  - w_ini_RF_SysTime
      wdt_BI_TotTime  = w_fin_BI_TotTime  - w_ini_BI_TotTime
      wdt_BI_UsrTime  = w_fin_BI_UsrTime  - w_ini_BI_UsrTime
      wdt_BI_SysTime  = w_fin_BI_SysTime  - w_ini_BI_SysTime
      wdt_EX0_TotTime = w_fin_EX0_TotTime - w_ini_EX0_TotTime
      wdt_EX0_UsrTime = w_fin_EX0_UsrTime - w_ini_EX0_UsrTime
      wdt_EX0_SysTime = w_fin_EX0_SysTime - w_ini_EX0_SysTime


*ETC* Now that the fit is over, I can inform the actual flux units used
*     and convert distance dependent quantities (like Mcor_tot masses)
*     to absolute ones using the input entries flux_unit and
*     LumDistInMpc. 
      fobs_norm = fobs_norm * flux_unit


c     Mark PL with an absurd L/M ratio... (for output-safety)
!!MEX!! No Power-Law
!      if (i_FitPowerLaw.EQ.1) fbase_norm(N_base) = 0.


!Cid29Mar07-FIX! Now write arq_out to out_dir, so figure what is the first
!     blank char in out_dir & concatenate it: arq_out --> arq_aux
      i_last_non_blank_char = 0
      do i=1,99
         if ((out_dir(i:i+1).EQ.' ').AND.(i_last_non_blank_char.EQ.0))
     &        i_last_non_blank_char = i - 1
      end do
      arq_aux = out_dir(1:i_last_non_blank_char) // arq_out

c     Open output file & write header
      open (unit=2,file=arq_aux,status='unknown')

      write (2,190)
 190  format(76('#'))
      write (2,'(a,3x,a)') 
     &     '## OUTPUT of PANcMExStarlight_v03.for ',
     &     '  [ElCid@Sanchica - 10/Feb/2012] ##'
      write (2,190)

      write (2,*)
      write (2,'(a)') '## Some input info'

*OBS: Added obs_dir because in CALIFA that's what identifies the galaxy!
      i_last_non_blank_char = 0
      do i=1,99
         if ((obs_dir(i:i+1).EQ.' ').AND.(i_last_non_blank_char.EQ.0))
     &        i_last_non_blank_char = i - 1
      end do
      write (2,'(a40,a,a)') arq_obs  , ' [arq_obs]  from ' , 
     &  obs_dir(1:i_last_non_blank_char)


      write (2,'(a40,a)')       arq_base    , ' [arq_base]'
      write (2,'(a40,a)')       arq_masks   , ' [arq_masks]'
      write (2,'(a40,a)')       arq_config  , ' [arq_config]'
      write (2,'(i5,35x,a)')    N_base      , ' [N_base]'
*ETC* Here added 3rd, 4th and 5th entries to inform whether this was a
*     FIR/QHR/PHO-fit!  If so, the FIR/QHR/PHO-results will be written
*     after all the rest of the standard MEX-Starlight output (ie, at
*     the end of the file).  This is quite ugly, but is the
*     change-the-least option!
      write (2,'(2(i4,4x),2x,3(i2,2x),10x,a,a)') N_exAV_components , 
     &     N_exAV , IsFIRcOn , IsQHRcOn , IsPHOcOn , 
     &    ' [N_exAV_components N_exAV    IsFIRcOn IsQHRcOn IsPHOcOn] ',
     &     '<=== FIRc-fit Y/N?  QHRc-fit Y/N?  PHOc-fit Y/N?'

      write (2,'(i5,35x,a)')  i_FitPowerLaw  , 
     &     ' [i_FitPowerLaw (1/0 = Yes/No) = NO in this version!!]'
      write (2,'(f5.2,35x,a)')  alpha_PowerLaw  , 
     &     ' [alpha_PowerLaw]'

      write (2,'(a3,37x,a)') red_law_option , ' [red_law_option]'
*ETC* Saving input-given flux_unit
      write (2,'(f9.6,5x,1p,e12.6,0p,14x,a)') q_norm , flux_unit ,
     &     ' [q_norm = A(l_norm)/A(V) & flux_unit]'

      write (2,*)
      write (2,'(a)') '## (Re)Sampling Parameters'
      write (2,'(f8.2,32x,a)') Olsyn_ini     , ' [l_ini (A)]'
      write (2,'(f8.2,32x,a)') Olsyn_fin     , ' [l_fin (A)]'
!Minor format change! ElCid@Sanchica - 12/Feb/2012
      write (2,'(f8.2,2x,f8.2,22x,a)') Odl_syn , dl_cushion ,
     $     ' [dl  &  dl_cushion (A)]'

      write (2,*)
      write (2,'(a)') '## Normalization info'
      write (2,'(f8.2,32x,a)')     l_norm   , ' [l_norm (A) - for base]'
      write (2,'(f8.2,32x,a)')  llow_norm   , 
     &     ' [llow_norm (A) - window for f_obs]'
      write (2,'(f8.2,32x,a)')  lupp_norm   , 
     &     ' [lupp_norm (A) - window for f_obs]'

*ETC* Inform fobs_norm AND Lobs_norm AND distance
      write (2,'(1p,2(e12.6,2x),0p,f8.3,4x,a)') 
     &     fobs_norm , Lobs_norm , LumDistInMpc ,
     &     ' [fobs_norm & Lobs_norm & LumDistInMpc]'
      write (2,*)
      write (2,'(a)') '## S/N'
      write (2,'(f8.2,32x,a)') llow_SN,' [llow_SN (A) - window for S/N]'
      write (2,'(f8.2,32x,a)') lupp_SN,' [lupp_SN (A) - window for S/N]'
      write (2,'(f8.3,32x,a)') SN_snwin   ,' [S/N in S/N window]'
      write (2,'(f8.3,32x,a)') SN_normwin ,' [S/N in norm. window]'
      write (2,'(f8.3,32x,a)')SNerr_snwin  ,' [S/N_err in S/N window]'
      write (2,'(f8.3,32x,a)')SNerr_normwin,' [S/N_err in norm. window]'
!OPTfn! Added IsOptimize_fn_OPT to output.
      write (2,'(f8.4,10x,i3,20x,a)') fscale_chi2 , IsOptimize_fn_OPT ,
     $    '[fscale_chi2 & IsOptimize_fn_OPT]'

      write (2,*)
      write (2,'(a)') '## etc...'
      write (2,'(i12,28x,a)')   idum_orig   , ' [idum_orig]'
      write (2,'(i6,34x,a)')   NOl_eff      , ' [NOl_eff]'
      write (2,'(i6,34x,a)')   Nl_eff       , ' [Nl_eff]'
      write (2,'(i6,4x,a6,25x,a)') Ntot_cliped , clip_method_option ,
     $     '[Ntot_cliped & clip_method]'
      write (2,'(i9,31x,a)')  Nglobal_steps , ' [Nglobal_steps]'
      write (2,'(i4,36x,a)')  N_chains , ' [N_chains]'
      write (2,'(i4,36x,a)')  NEX0s_base ,
     $     ' [NEX0s_base = N_base in EX0s-fits]'
!Cid24Mar07-FIX! Added censored weights info here
      write (2,'(3(i2,2x),2x,i5,2x,f4.1,2x,1p,e12.5,1x,a,a)') 
     $     iCLIPBUG_flag ,i_RC_CRASH_FLAG, IsBurInOver_BurnIn , 
     $     n_censored_weights , wei_nsig_threshold , wei_limit ,
     $     ' [Clip-Bug, RC-Crash & Burn-In warning-flags,' , 
     $     ' n_censored_weights, wei_nsig_threshold & wei_limit]'
!Cid29Dec06-FIX! Removed one empty line to write time-counters info!

! Changes format of idt & wdt_* (to avoid silly problems ocurred in my 1st CALIFA Runs,
! when time would sometime exceed the previous format!)
! ElCid@Sanchica - 25/Nov/2011
      write (2,'(i7,3(1x,f10.3),a)') 
     $     idt_all , wdt_TotTime , wdt_UsrTime , wdt_SysTime ,
     $     ' [idt_all, wdt_TotTime, wdt_UsrTime & wdt_SysTime (sec)]'

c     Synthesis results - Global Minimum!
      write (2,*)
      write (2,'(a)') '## Synthesis Results - Best model ##'
      write (2,*)

!Cid24Mar07-FIX! Before I mixed EX0sGpar_min, EX0sv0_min & etc with 
!     x_min, AV_min etc in this output stage. I've reorganized it so that
!     now everything is in terms of the *_min quantities.
      call MEX_Update_AV_array(N_par,N_base,AV_min,exAV_min,AV_tot)
      call Kcalc_Fsyn(Nl_obs,x_min,AV_tot,v0_min,vd_min,f_syn)
      chi2_min = FK_chi2(x_min,AV_tot,v0_min,vd_min) / float(Nl_eff)
      adev     = F_adev(f_syn)
      sum_x    = 100. * Sum_of_x(x_min,N_base)

      call ConvertLight2Mass(N_base,x_min,AV_tot,q_norm,Lobs_norm
     $     ,fbase_norm,i_FitPowerLaw,Mstars,Lum_tot,Mini_tot,Mcor_tot
     $     ,mu_ini,mu_cor)

*ETC* Added total chi2 to output & extra horizontal spaces here and below for elegance!
      write (2,'(1p,e12.5,5x,e12.5,38x,a)') chi2_min  , 
     &     chi2_min * float(Nl_eff) , ' [chi2/Nl_eff & chi2]'
      write (2,'(f12.5,55x,a)') adev         , ' [adev (%)]'
* ...........................................................................
*ETC* Adds a line with chi2_TOT, _Opt, _FIR, _QHR & PHO to the output!
*     OBS1: This line fills one which was previsouly empty in arq_out,
*           and hence does not screw up the output format.
*     OBS2: Notice a few bstruse manouvres to isolate chi2_Opt ...
*     OBS3: Notice the trick of calling ETC_CalcETCchi2 pretending
*           chi2_Optical = 1, such that the result is chi2_*/chi2_Opt
      chi2_TOT     = chi2_min * float(Nl_eff)

*ETC* Through common/ETC_Chi2s, recover chi2_FIR, chi2_QHR & chi2_PHO
*     per unit optical-chi2, and from this recover chi2_Opt.
      zz_nothing   = ETC_CalcETCchi2(x_min,AV_tot,1.0,'tmp')
      chi2_FIR2Opt = chi2_FIR
      chi2_QHR2Opt = chi2_QHR
      chi2_PHO2Opt = chi2_PHO
      chi2_Opt     = 
     $     chi2_TOT / (1. + chi2_FIR2Opt + chi2_QHR2Opt + chi2_PHO2Opt)

*ETC* Again through common/ETC_Chi2s, recovers the true absolute values
*     of chi2_FIR, chi2_QHR & chi2_PHO.
      zz_nothing  = ETC_CalcETCchi2(x_min,AV_tot,chi2_Opt,'tmp')
      
      write (2,'(1p,5(e12.5,1x),2x,a)') chi2_TOT , chi2_Opt , chi2_FIR ,
     $     chi2_QHR , chi2_PHO , 
     $   ' [chi2_TOT chi2_Opt chi2_FIR chi2_QHR chi2_PHO - ALL chi2s!]'
* ...........................................................................
      write (2,'(f12.5,55x,a)') sum_x        , ' [sum-of-x (%)]'
*ETC* OBS: If the base spectra are in Lsun/Angs/Msun and proper flux
*     units were given for the observed spectrum, this Lum_tot should be 
*     in Lsun/Angs and the Mcor_tot and Mini_tot masses should be in
*     Msun...
      write (2,'(1p,e12.5,55x,a)') Lum_tot , 
     &     ' [Lum_tot (Lsun/A if distance & flux_unit are Ok...)]'
      write (2,'(1p,e12.5,55x,a)') Mini_tot ,
     &     ' [Mini_tot (Msun  if distance & flux_unit are Ok...)]'
      write (2,'(1p,e12.5,55x,a)') Mcor_tot , 
     &     ' [Mcor_tot (Msun  if distance & flux_unit are Ok...)]'

      write (2,*)

*ETC* Added warnings to signal when vd/v0 are = initial values...
      if ((abs(v0_min-v0_ini).LT.0.01).AND.
     $    (FitOrSamplePDF_option.NE.'FXK')) then 
         write (2,'(f7.2,34x,a,f7.2,a)') v0_min ,
     $        '[v0_min  (km/s)] ... suspiciously close to v0_ini=' ,
     $        v0_ini , ' ... ?!'
      else
         write (2,'(f7.2,61x,a)') v0_min  ,'[v0_min  (km/s)]'
      end if

      if ((abs(vd_min-vd_ini).LT.0.01).AND.
     $    (FitOrSamplePDF_option.NE.'FXK')) then 
         write (2,'(f7.2,34x,a,f7.2,a)') vd_min ,
     $        '[vd_min  (km/s)] ... suspiciously close to vd_ini=' ,
     $        vd_ini , ' ... ?!'
      else
         write (2,'(f7.2,61x,a)') vd_min  ,'[vd_min  (km/s)]'
      end if

      write (2,'(f7.4,61x,a)') AV_min  ,'[AV_min  (mag)]'
      write (2,'(10(1x,f6.3),10x,a,i3)') (exAV_min(j),j=1,10)
     $     ,'[exAV_min (mag)] 1...10 of N_exAV = ', N_exAV
      write (2,*)

!Cid20Jan07-FIX! included an alpha/Fe array in pop-vector table, as well as
!     new single component fits SSP_* arrays!
*ETC*  Added yet another column (after AV_tot) to report LAx_min(j) = the
*     best pop-vector (x_min) averaged over all lambdas in the fit
*     (weighting by weight_obs), ie, <x_lambda>_j. This is just an ouput
*     thing!
*     Cid@Granada - 31/Jan/2010 ....
      call AverageXoverLambda(x_min,AV_tot,LAx_min)

      write (2,191)
  191  format ('# j     x_j(%)      Mini_j(%)     Mcor_j(%)' , 
     &     '     age_j(yr)     Z_j      (L/M)_j   exAV?  Mstars',
     &     '   component_j        a/Fe...',
     &     '    |  SSP_chi2r SSP_adev(%)   SSP_AV   SSP_x(%)',
     &     '    |  AV_tot   <LAx>_j(%)')

      do j=1,N_base
         write (2,192) j , (100.*x_min(j)) , 
     &        (100.*mu_ini(j)) , (100.*mu_cor(j)) , 
     &        age_base(j) , Z_base(j) , fbase_norm(j) , 
     &        ind_exAV(j) , Mstars(j) , base_name(j) , alphaFe_base(j) ,
     &        SSP_chi2(j) / float(Nl_eff) , SSP_adev(j) , 
     &        SSP_AV(j) , 100. * SSP_x(j) ,
     &        AV_tot(j) , 100. * LAx_min(j)
      end do
! Changed Z_base format to f8.5 at Paula's suggestion - 15/Ago/2011
! She says this is better form MILES-formatted log-Z [eg, -1.5]
!
! Changes format of SSP_x from f9.4 to f9.3 (to avoid silly problems ocurred in my 1st CALIFA Runs,
! when SSP_x would sometimes exceed 9999.9999!)
! ElCid@Sanchica - 25/Nov/2011
 192  format (i3,2x,f10.4,3x,
     &     1p,e12.4,2x,e12.4,4x,
     &     e12.6,1x,0p,f8.5,2x,1p,e9.3,3x,
     &     0p,i2,3x,f6.4,3x,a15,4x,f9.4,4x,
     &     1p,e11.4,0p,3x,f8.4,2x,
     &     f7.4,2x,f9.3,6x,
     &     f7.4,2x,f10.4)

!     Screen output Best model results
      if (i_verbose.GT.0) then
         write (*,*)
         write (*,'(a,a)') '**> RESULTS for ' , arq_out
         write (*,812) N_base
 812     format ('   N_tot chi2/Nl_eff    adev |     v0    vd | ',
     &        'sum       A_V |     x_j... [%], j = 1 ...',i3)
         write (*,813) Nglobal_steps , chi2_min , adev , 
     &        v0_min , vd_min , sum_x , AV_min ,
     &        (100.*x_min(j),j=1,min(10,N_base))
 813     format (i8,1x,1p,e10.4,0p,2x,f7.4,' | ' ,
     &        f6.1,1x,f5.1,' | ',f7.3,1x,f5.2,' | ',10(1x,f6.2))
            do j_out=2,int(N_base/10)+1
               j1 = (j_out - 1) * 10 + 1
               j2 = min(j1+9,N_base)
               write (*,'(62x,10(1x,f6.2))') (100.*x_min(j),j=j1,j2)
            end do

            if (N_exAV.GT.0) then
               write (*,'(62x,a)') '   AV_j... [mag]'
               do j_out=1,int(N_base/10)+1
                  j1 = (j_out - 1) * 10 + 1
                  j2 = min(j1+9,N_base)
                  write (*,'(62x,10(1x,f6.2))') (AV_tot(j),j=j1,j2)
               end do
            end if

            write (*,'(70x,a,1p,e12.4,a,0p,i6,a)') 'Total chi2 = ' ,
     $           chi2_min * Nl_eff , '  for Nl_eff =' , Nl_eff ,
     $           ' lambdas!'
         write (*,*)

      end if


c     Now output results for best (repeated), average & last-chain
c     parameters after Burn-In loop! (This should give us an idea of the
c     errors, and in any case it'd be a waste to throw all these numbers
c     away!)
      write (2,*)
      write (2,*)
      write (2,'(a)') '## Synthesis Results - Average & Chains ##'
      write (2,*)

c     Store chi2 & AV for each chain in f_aux variables, chain-mass
c     fractions in EX0sChainPar array & chain-masses in CEne_min (not
c     used anymore)...
      do i_chain=1,N_chains
         do j=1,N_par
            f_aux1(j) = ChainPar(j,i_chain)
         end do
!ATT: AV_tot is overwritten here! <====
         call MEX_Unpack_par2xAVYAVfn(N_par,N_base,f_aux1,x,AV,exAV,fn)
         call MEX_Update_AV_array(N_par,N_base,AV,exAV,AV_tot)
         chi2_red = FK_chi2(x,AV_tot,v0_min,vd_min) / float(Nl_eff)
         f_aux2(i_chain) = AV
         f_aux3(i_chain) = chi2_red
*ETC* Changed fobs_norm -> Lobs_norm & Flux_tot->Lum_tot
         call ConvertLight2Mass(N_base,x,AV_tot,q_norm,Lobs_norm
     $        ,fbase_norm,i_FitPowerLaw,Mstars,Lum_tot,chain_Mini_tot
     $        ,chain_Mcor_tot,mu_ini,aux1)
         do j=1,N_base
            EX0sChainPar(j,i_chain) = aux1(j) 
         end do
         aux2(i_chain) = chain_Mcor_tot
*ETC* Store in a LAx_chain matrix the Lambda-Averaged chain-pop-vectors, 
*     to be saved below. 
         call AverageXoverLambda(x,AV_tot,aux3)
         do j=1,N_base
            LAx_chain(j,i_chain) = aux3(j)
         end do
      end do

c     Now compute things for average solution
      call MEX_Unpack_par2xAVYAVfn(N_par,N_base,Gpar_ave,x,AV,exAV,fn)
      call MEX_Update_AV_array(N_par,N_base,AV,exAV,AV_tot)
      AV_ave = AV
      chi2_ave = FK_chi2(x,AV_tot,v0_min,vd_min) / float(Nl_eff)
*ETC* Changed fobs_norm -> Lobs_norm & Flux_tot->Lum_tot
      call ConvertLight2Mass(N_base,x,AV_tot,q_norm,Lobs_norm
     $     ,fbase_norm,i_FitPowerLaw,Mstars,Lum_tot,Mini_tot_ave
     $     ,Mcor_tot_ave,mu_ini,aux3)
*ETC* Store in LAx_ave the Lambda-Averaged average pop vector, to be
*     saved below
      call AverageXoverLambda(x,AV_tot,LAx_ave)

c     Print out x, mu, AV, chi2 & Masses for ave & chain-solutions
!!MEX!!  I now report ALL parameters, not just the x_j's, in this 1st
!     block. 
!     ATT: This CHANGES the size of the output file!!
!CALIFA! The "header" line was made more informative (I got complaints from Rafa & Natalia)
!	  ElCid@Sanchica - 18/Oct/2011
      ieAVl = N_base + 2
      ieAVu = N_par  - 1
      if (N_exAV.LE.0.) then
          ieAVl = -0
          ieAVu = -0
      endif
      write (2,'(a,i3,a,5(a,i3))')
     $     '# j    par_j: min, <> & last-chain-values for 1 ...'
     $     , N_chains , ' chains.'
     $     ,' x: 1 -> ',N_base, ' ; AV: ',N_base+1
     $     ,' ; exAV (if any): ',ieAVl,' -> ',ieAVu,' ; fn:',N_par
      do j=1,N_par
         aux = 100.
         if ((j.GT.N_base).AND.(j.NE.N_par)) aux = 1
         write (2,'(i3,40(2x,f10.4))') j , aux*EX0sGpar_min(j) ,
     $        aux*Gpar_ave(j) , ((aux*ChainPar(j,i_chain)),i_chain=1
     $        ,N_chains)
      end do

*ETC* Print out the LAx* Lambda-Averaged pop vector: min (repeated from above), 
*     average and chain-results.
!ATT> This block is NEW and it CHANGES the FILE-SIZE by N_base + 2 lines!!! 
!     Cid@Granada-31/Jan/2010
      write (2,*)
      write (2,'(a,a,i3,a)')
     $     '# j Lambda-Averaged pop-vectors <LAx_*>_j: ',
     $     'min, <> & last-chain-values for 1 ...',N_chains,' chains'
      do j=1,N_base
         write (2,'(i3,40(2x,f10.4))') j , 100.*LAx_min(j) ,
     $        100.*LAx_ave(j) , ((100.*LAx_chain(j,i_chain)),i_chain=1
     $        ,N_chains)
      end do

      write (2,*)
      write (2,'(a,i3,a)')
     $     '# j   Mcor_j: min, <> & last-chain-values for 1 ...'
     $     ,N_chains , ' chains'
      do j=1,N_base
         write (2,'(i3,1p,40(2x,e10.4))') j , 100.*mu_cor(j) , 100.
     $        *aux3(j), ((100.*EX0sChainPar(j,i_chain)),i_chain=1
     $        ,N_chains)
      end do

      write (2,*)
      write (2,'(a,a,i3,a)') '# chi2/Nl_eff & Mass for min, <> & ' ,
     $     'i_chain = 1 ...', N_chains , ' solutions'
      write (2,'(a5,1p,40(e10.4,2x))') 'chi2 ', chi2_min , chi2_ave ,
     $     (f_aux3(i_chain),i_chain=1,N_chains)
      write (2,'(a5,1p,40(e10.4,2x))') 'Mass ' , Mcor_tot ,
     $     Mcor_tot_ave ,(aux2(i_chain),i_chain=1,N_chains)
      write (2,'(f7.2,34x,a)') v0_min_BeforeEX0s ,
     $     '[v0_min (km/s) before EX0s...]'
      write (2,'(f7.2,34x,a)') vd_min_BeforeEX0s , 
     $     '[vd_min (km/s) before EX0s...]'
      write (2,*)


c     Observed, Synthetic & Weight Spectra: Save f_obs, f_syn &
c     weight_obs. Save only points inside originally input range:
c     Olsyn_ini ==> Olsyn_fin!  Points which have been clipped (those
c     with f_mask < 0) are saved with -1 in the weight_obs column just
c     for control.
      Nl_obs_save = 0
      do i=1,Nl_obs
         if ((l_obs(i).GE.Olsyn_ini).AND.(l_obs(i).LE.Olsyn_fin)) then
            Nl_obs_save = Nl_obs_save + 1
         end if
      end do

      write (2,*)
      write (2,*)

!v05! Saving Best_f_SSP() too (if wanted)
      if (i_SaveBestSingleCompFit.EQ.0) then
         write (2,'(a,a)') '## Synthetic spectrum (Best Model) ##',
     &        'l_obs f_obs f_syn wei'
      else 
         write (2,'(a,a)') '## Synthetic spectrum (Best Model) ##',
     &        'l_obs f_obs f_syn wei Best_f_SSP'
      end if
! Added i_SaveBestSingleCompFit in this line to help reading macros....
! ElCid@Sanchica - 29/Jan/2012
      write (2,'(i5,4x,i5,4x,i5,4x,a)') 
     $     Nl_obs_save , index_Best_SSP , i_SaveBestSingleCompFit ,
     $     '[Nl_obs & index_Best_SSP & i_SaveBestSingleCompFit]'

      do i=1,Nl_obs

         if ((l_obs(i).GE.Olsyn_ini).AND.(l_obs(i).LE.Olsyn_fin)) then
            aux = weight_obs(i)
            if  (f_mask(i).LT.0.) aux = f_mask(i)
            if (i_SaveBestSingleCompFit.EQ.0) then
! Relucted to do this, but finally changed spectra format... 
! ElCid@Sanchica - 29/Jan/2012
               write (2,'(f8.2,1p,3(1x,e12.5))') 
     &              l_obs(i) , f_obs(i) , f_syn(i) , aux
            else
               write (2,'(f8.2,1p,4(1x,e12.5))') 
     &              l_obs(i) , f_obs(i) , f_syn(i) , aux , Best_f_SSP(i)
            end if
         end if
      end do
c     Standard output ends here. ETC-related output follows below.
c ***************************************************************************



* ***************************************************************************
*             Preparations for ETC [= FIR & QHR & PHO] output               *
* ***************************************************************************
*ETC*  I have overwritten the best-fit AV_tot above (while saving
*     chains-results), so lets first restore it in case there is FIR
*     and/or QHR output to be saved below.
      call MEX_Update_AV_array(N_par,N_base,AV_min,exAV_min,AV_tot)

*     Every time any chi2 function [like FK_chi2, ...] is called, the
*     FIR_ & QHR_Calc*Chi2 functions are also called upp, and new
*     values are computed for things like FIRModlogLFIR and
*     QHR_chi2_Y(). In the savings above I have called up FK_chi2 with
*     Chain & Average parameters, and this screwed up these FIR & QHR
*     output quantities! This is why I recompute them below. Also, I use
*     the 'OUT' string to compute some final stuff for output.
      zz_nothing  = ETC_CalcETCchi2(x_min,AV_tot,chi2_Opt,'OUT')
* ***************************************************************************



c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
c FIR                        FIR-RELATED OUTPUT          .                FIR
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
*FIR* Produce/re-compute/compile FIR-related things for output
*     The if just below skips FIRc-output for non-FIRc fits (IsFIRcOn =
*     0). Note that for IsFIRcOn = -1 the predicted FIR luminosity is
*     given, even if no FIR-data was provided!
      if ((iFIR_IsFIRcOn.EQ.1).OR.(iFIR_IsFIRcOn.EQ.-1)) then
      write (2,*)
      write (2,*)
      write (2,*)
      write (2,'(35(a))') ('#',j=1,35)
      if (IsFIRcOn.EQ.-1) then
         write (2,'(a,a)') '##    FIRc-FIT RELATED OUTPUT    ##' ,
     &        ' (Predicted things only!) '
      else
         write (2,'(a)') '##    FIRc-FIT RELATED OUTPUT    ##'
      end if
      write (2,'(35(a))') ('#',j=1,35)

      write (2,*)
      write (2,'(a)') '## FIR: Input info'
      write (2,'(a40,a)')    arq_ETCinfo  , ' [arq_ETCinfo]'
*     OBS: Redundant LumDistInMpc included here just for backwards compatibility...
      write (2,'(f9.4,31x,a)') LumDistInMpc , 
     &     ' [LumDist (Mpc)]'
      write (2,'(f9.4,31x,a)') FIR_logLFIR_TOT , 
     &     ' [FIR_logLFIR_TOT (Lsun)]'
      write (2,'(f9.4,31x,a)') FIR_LFIRFrac2Model , 
     &     ' [FIR_LFIRFrac2Model]'
      write (2,'(f9.4,31x,a)') FIR_logLFIR_obs , 
     &     ' [FIR_logLFIR_obs (Lsun)]'
      write (2,'(f9.4,31x,a)') FIR_ErrlogLFIR  , 
     &     ' [FIR_ErrlogLFIR   (dex)]'
      write (2,'(f9.4,31x,a)') FIR_RangelogLFIR , 
     &     ' [FIR_RangelogLFIR (dex)]'
      write (2,'(f9.4,31x,a)') FIRChi2ScaleFactor , 
     &     ' [FIRChi2ScaleFactor]'
      write (2,'(2(f9.4,2x),18x,a)') FIR_logLFIR_low , FIR_logLFIR_upp , 
     &     ' [FIR_logLFIR_low -> FIR_logLFIR_upp (Lsun)]'
      write (2,'(2(f9.4,2x),18x,a)') FIRbeta_D , FIRbeta_I , 
     &     ' [FIRbeta_D & FIRbeta_I]'
*     Extra: Output a rough FIR/Optical "color"
      obs_FIR2Opt_ratio = FIR_logLFIR_obs - 
     &     (log10(l_norm) + FIRlogLobs_norm)
      write (2,'(f9.4,31x,a)') obs_FIR2Opt_ratio , 
     &     ' [log L_FIR/(l_norm * Lobs_norm) - a rough FIR/Opt ratio]'
      write (2,*)

*!!!* ATT: There was an ERROR here (discovered during QHR-development). 
*     This part was like this:
*      chi2_TOT = chi2_min * float(Nl_eff)
*      chi2_OPT = chi2_TOT - chi2_FIR
*     But chi2_FIR was NOT defined! This is now fixed, with the chi2's
*     computed above (just before this block). This ERROR/BUG only
*     affected chi2-related output on the FIR-output-block. This is
*     another reason to suspect that the Jan/Feb-2010 runs w/Rosa@IAA
*     where worng (even if only slightly so)...

*     Compute the predicted log LFIR, the FIR-population vector
*     (x_FIR) & chi2_FIR for the best model (x_min & AV_tot) by calling
*     the FIR_CalcFIRchi2 function with WhereFrom = 'OUT'. The predicted
*     bolometric luminosity and bolometric pop-vector are also computed
*     and outputed.
      write (2,'(a)') '## FIR: Results for best model'
      z1 = 0.
      if (IsFIRcOn.EQ.1) z1 = 10.**(FIRModlogLFIR - FIR_logLFIR_obs)
      write (2,'(f11.5,5x,1p,e9.3,0p,15x,a)') FIRModlogLFIR  , z1 ,
     $     ' [FIRModlogLFIR (Lsun) & Mod/Obs ratio]'
      z1 = 10.**(FIRModlogLFIR - FIRModlogLBOL)
      write (2,'(f11.5,5x,1p,e9.3,0p,15x,a)') FIRModlogLBOL , z1 ,
     $     ' [FIRModlogLBOL (Lsun) & FIR/BOL ratio]'
      write (2,'(1p,e11.5,5x,e11.5,0p,13x,a)') 
     $     chi2_FIR , chi2_Opt ,
     $     ' [chi2_FIR (already scaled) & chi2_Opt]'
      write (2,'(2(f11.5,2x),14x,a)') (100.*chi2_FIR/chi2_TOT) , 
     $     (100.*chi2_Opt/chi2_TOT) ,
     $     ' [%-chi2: FIR/TOT & optical/TOT]'
      write (2,*)

*     Screen output of FIR-results
      if (i_verbose.GT.0) then
         write (*,*)
         write (*,*)
         write (*,'(a)') '## FIR: Results for best model'
         z1 = 0.
         if (IsFIRcOn.EQ.1) z1 = 10.**(FIRModlogLFIR - FIR_logLFIR_obs)
         write (*,'(f11.5,5x,1p,e9.3,0p,15x,a)') FIRModlogLFIR  , z1 ,
     $        ' [FIRModlogLFIR (Lsun) & Mod/Obs ratio]'
         z1 = 10.**(FIRModlogLFIR - FIRModlogLBOL)
         write (*,'(f11.5,5x,1p,e9.3,0p,15x,a)') FIRModlogLBOL , z1 ,
     $        ' [FIRModlogLBOL (Lsun) & FIR/BOL ratio]'
         write (*,'(1p,e11.5,5x,e11.5,0p,13x,a)') 
     $        chi2_FIR , chi2_Opt , 
     $        ' [chi2_FIR (already scaled) & chi2_Opt]'
         write (*,'(2(f11.5,2x),14x,a)') (100.*chi2_FIR/chi2_TOT) , 
     $        (100.*chi2_Opt/chi2_TOT) ,
     $        ' [%-chi2: FIR/TOT & Opt/TOT]'
         write (*,*)
         write (*,1191)
      end if

*     OBS: Lbol(j) is the base-original per unit initial mass value!
      write (2,1191)
 1191 format ('# j  age_j(yr)     Z_j      exAV?  ',
     &     ' x_j(%)    AV_tot  | ',
     &     ' x_FIR(%)    x_BOL(%)  ',
     &     ' BolCor     FracLion   Lbol/M     ',
     &     ' Rmat       R_opt      R_Lya      R_LCE')
      do j=1,N_base
         ind_FIRAV = int(0.5 + (AV_tot(j) - FIRAV_low) / FIRdAV) + 1
         if (ind_FIRAV.LT.1) ind_FIRAV = 1
         if (ind_FIRAV.GT.NFIR_AV) ind_FIRAV = NFIR_AV
         write (2,1192) j , age_base(j) , Z_base(j) , ind_exAV(j) , 
     &        (100.*x_min(j)) , AV_tot(j) , 
     &        (100.*x_FIR(j)) , (100.*x_BOL(j)) ,
     &        FIRBolCor(j)    , FIRFracLion(j) , FIRLbol(j) ,
     &        FIRRmat(j,ind_FIRAV)  , FIRR_opt(j,ind_FIRAV) ,
     &        FIRR_Lya(j,ind_FIRAV) , FIRR_LCE(j,ind_FIRAV)
         if (i_verbose.GT.0) then
            write (*,1192) j , age_base(j) , Z_base(j) , ind_exAV(j) , 
     &           (100.*x_min(j)) , AV_tot(j) , 
     &           (100.*x_FIR(j)) , (100.*x_BOL(j)) ,
     &           FIRBolCor(j)    , FIRFracLion(j) , FIRLbol(j) ,
     &           FIRRmat(j,ind_FIRAV)  , FIRR_opt(j,ind_FIRAV) ,
     &           FIRR_Lya(j,ind_FIRAV) , FIRR_LCE(j,ind_FIRAV)
         end if
      end do
 1192 format (i3,2x,1p,e12.6,2x,0p,f7.5,2x,i3,2x,
     &     f10.4,2x,f7.4,2x,
     &     2(f10.4,2x),1p,
     &     3(2x,e9.3),1x,
     &     4(2x,e9.3))
      end if
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR



c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
c QHR                        QHR-RELATED OUTPUT          .                QHR
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
*QHR* Produce/re-compute/compile QHR-related things for output
*     The if just below skips QHRc-output for non-QHRc fits (IsQHRcOn =
*     0). Note that for IsQHRcOn = -1 predicted QHR things are saved,
*     even though no actual QHR-data was used.
      if ((iQHR_IsQHRcOn.EQ.1).OR.(iQHR_IsQHRcOn.EQ.-1)) then
         write (2,*)
         write (2,*)
         write (2,*)
         write (2,'(35(a))') ('#',j=1,35)
         if (IsQHRcOn.EQ.-1) then
            write (2,'(a,a)') '##    QHRc-FIT RELATED OUTPUT    ##' ,
     $           ' (Predicted things only!) '
         else
            write (2,'(a)') '##    QHRc-FIT RELATED OUTPUT    ##'
         end if
         write (2,'(35(a))') ('#',j=1,35)

         write (2,*)
         write (2,'(a)') '## QHR: Input info'
         write (2,'(f6.4,49x,a)') QHRbeta_I   , ' [QHRbeta_I]'
         write (2,'(a55,a)')     arq_ETCinfo  , ' [arq_ETCinfo]'
         write (2,'(f9.4,46x,a)') LumDistInMpc , ' [LumDist (Mpc)]'
         write (2,'(f9.4,46x,a)') QHR_GlobalChi2ScaleFactor , 
     $        ' [QHR_GlobalChi2ScaleFactor]'
         write (2,'(i3,52x,a)') NQHR_Ys , 
     $        ' [NQHR_Ys = number of em-lines included in fit]'

         write (2,*)
         write (2,'(a,a,a,a)') '# lambda[A]      frecomb     logY_TOT' , 
     $        '  YFrac2Model      ErrlogY     RangelogY' , 
     $        ' Chi2ScaleFactor' ,
     $        ' logY_obs     logY_low     logY_upp'
         do iY=1,NQHR_Ys
            write (2,'(i2,1x,f8.2,9(3x,f10.4))') iY , QHR_lambda(iY) ,
     $           QHR_frecomb(iY) ,QHR_logY_TOT(iY) , QHR_YFrac2Model(iY)
     $           ,QHR_ErrlogY(iY) , QHR_RangelogY(iY)
     $           ,QHR_Chi2ScaleFactor(iY) , QHR_logY_obs(iY)
     $           ,QHR_logY_low(iY) ,  QHR_logY_upp(iY)
         end do

*ELR* .......................................................................
         write (2,*) 
         if (IsELRInputInconsistentWarningFlag.EQ.0) then 
            write (2,'(a,a)') '## ELR input & results. ' , 
     $           ' OBS: chi2_ELR is **included** in chi2_QHR!'
         else
            write (2,'(4(a))') '## ELR input & results. ' , 
     $           ' OBS: chi2_ELR is **included** in chi2_QHR! ' ,
     $           'ATT: YOUR INPUT DATA IMPLY AV_neb(A,B) < 0 (!) ' , 
     $           'SO I **MODIFIED** IT!!!!! BEWARE!!! <==='
         end if

c     predicted AV_neb from ELR_lambdas A & B, & its error
         z1 = 0.
         z2 = 0.
         if ((iY_ELR_A.GT.0).AND.(iY_ELR_B.GT.0)) then
            z1 = -(2.5 / (QHR_qlambda(iY_ELR_A) -
     $           QHR_qlambda(iY_ELR_B))) * (ELR_logRobs - ELR_logRint)
            z2 = abs(2.5 / (QHR_qlambda(iY_ELR_A) -
     $           QHR_qlambda(iY_ELR_B))) * ELR_ErrlogR
         end if

         write (2,1357)
     $        IsELROn , ELR_lambdaA, ELR_lambdaB , iY_ELR_A , iY_ELR_B ,
     $        ELR_logRint , z1 , z2 ,
     $        ' [IsELROn lambda_A lambda_B ind_A ind_B' ,
     $        ' logRint AV_neb(A/B) +/- its error]'
 1357    format (i1,1x,f8.2,3x,f8.2,5x,i1,4x,i1,4x,f9.4,3x,2(1x,f9.4),
     $        1x,a,a)

         write (2,'(1x,2(f9.4,2x),3(f9.4,4x),7x,a,a)') 
     $        ELR_ErrlogR , ELR_RangelogR, 
     $        ELR_logR_low , ELR_logR_upp , ELR_Chi2ScaleFactor ,
     $        ' [Err_logR RangelogR logR_low logR_upp' ,
     $        ' Chi2ScaleFactor]' 

         z = 0.
         if (chi2_QHR.GT.0.) z = (ELR_chi2_logR/chi2_QHR)
         write (2,'(1x,2(f9.4,2x),1p,3(1x,e12.4),0p,7x,a,a)') 
     $        ELR_ModlogR , ELR_logRobs , ELR_chi2_logR ,
     $        z , (ELR_chi2_logR/chi2_TOT) , 
     $        ' [logRobs ModlogR' , 
     $        ' chi2_ELR chi2_ELR/chi2_QHR chi2_ELR/chi2_TOT]'
*ELR* .......................................................................

         write (2,*)
         write (2,'(a)') '## QHR: Results for best model'
         write (2,'(2(f10.5,2x),16x,a)') 
     $        QHR_log_QH0 , QHR_log_QHeff , 
     $        ' [log_QH0 & log_QHeff - in photons/sec]'
         write (2,'(2(f10.5,2x),16x,a)') 
     $        (100.*chi2_QHR/chi2_TOT) , 
     $        (100.*chi2_Opt/chi2_TOT) ,
     $        ' [%-chi2: QHR/TOT Opt/TOT]'
         write (2,'(1p,3(e11.5,2x),0p,1x,a,a)') 
     $        chi2_QHR , chi2_Opt , chi2_TOT , 
     $        ' [chi2_QHR chi2_Opt chi2_TOT]'
         write (2,*)
         write (2,'(a,a)') '# lambda[A]   q_lambda     logY_obs      ',
     $        'ModlogY       chi2_Y    chi2_Y/chi2_Opt  chi2_Y/chi2_TOT'
         do iY=1,NQHR_Ys
            write (2,2011) iY , QHR_lambda(iY) , QHR_qlambda(iY) ,
     $           QHR_logY_obs(iY) , QHR_ModlogY(iY) , QHR_chi2_Y(iY) , 
     $           (QHR_chi2_Y(iY)/chi2_Opt) , (QHR_chi2_Y(iY)/chi2_TOT)
         end do
 2011    format (i2,1x,f8.2,3x,f8.4,
     $        2(3x,f10.4),3x,1p,e12.4,
     $        2x,e12.4,5x,e12.4)
         write (2,*)

*     Repeat base info (age, Z, AV...) now with QHR-stuff too: 
*       qH        = 1e40 photons/sec/initial solar mass
*       QH/Lnorm  = 1e40 photons/sec per Lsun/Angs
*       x_QH0()   = %-contributions to total QH0
*       x_QHeff() = %-contributions to total QHeff
*       x_Y()     = %-contributions to each emission line
         write (2,2012)  (int(QHR_lambda(iY)),iY=1,NQHR_Ys)
 2012    format ('# j  age_j(yr)     Z_j      exAV?  ',
     &        ' x_j(%)    AV_tot  |  ',
     &        ' qH__40      QH2Lnorm__40  ',
     &        ' %-QH0 %-QHeff  |',
     &        ' %-Y:',10(i6))

         do j=1,N_base
            write (2,2013) j , age_base(j) , Z_base(j) , ind_exAV(j) , 
     &           (100.*x_min(j)) , AV_tot(j) , 
     &           QHR_qH40(j) ,  QHR_Q2Lnorm40(j),
     &           (100.*QHR_x_QH0(j)) , (100.*QHR_x_QHeff(j)) , 
     &           (100.*QHR_x_Y(j,iY),iY=1,NQHR_Ys)
         end do
 2013    format (i3,2x,1p,e12.6,2x,0p,f7.5,2x,i3,
     &        2x,f10.4,2x,f7.4,
     &        2x,1p,2(2x,e12.4),0p,
     &        2x,f6.2,2x,f6.2,
     &        9x,10(f5.1,1x))

*     Screen output of QHR-results
         if (i_verbose.GT.0) then
*ELR* .......................................................................
         write (*,*) 
         write (*,*) 
         if (IsELRInputInconsistentWarningFlag.EQ.0) then 
            write (*,'(a,a)') '## ELR input & results. ' , 
     $           ' OBS: chi2_ELR is **included** in chi2_QHR!'
         else
            write (*,'(4(a))') '## ELR input & results. ' , 
     $           ' OBS: chi2_ELR is **included** in chi2_QHR! ' ,
     $           'ATT: YOUR INPUT DATA IMPLY AV_neb(A,B) < 0 (!) ' , 
     $           'SO I **MODIFIED** IT!!!!! BEWARE!!! <==='
         end if
         z1 = 0.
         z2 = 0.
         if ((iY_ELR_A.GT.0).AND.(iY_ELR_B.GT.0)) then
            z1 = -(2.5 / (QHR_qlambda(iY_ELR_A) -
     $           QHR_qlambda(iY_ELR_B))) * (ELR_logRobs - ELR_logRint)
            z2 = abs(2.5 / (QHR_qlambda(iY_ELR_A) -
     $           QHR_qlambda(iY_ELR_B))) * ELR_ErrlogR
         end if
         write (*,1357)
     $        IsELROn , ELR_lambdaA, ELR_lambdaB , iY_ELR_A , iY_ELR_B ,
     $        ELR_logRint , z1 , z2 ,
     $        ' [IsELROn lambda_A lambda_B ind_A ind_B' ,
     $        ' logRint AV_neb(A/B) +/- its error]'
         write (*,'(1x,2(f9.4,2x),3(f9.4,4x),7x,a,a)') 
     $        ELR_ErrlogR , ELR_RangelogR, 
     $        ELR_logR_low , ELR_logR_upp , ELR_Chi2ScaleFactor ,
     $        ' [Err_logR RangelogR logR_low logR_upp' ,
     $        ' Chi2ScaleFactor]' 
         z = 0.
         if (chi2_QHR.GT.0.) z = (ELR_chi2_logR/chi2_QHR)
         write (*,'(1x,2(f9.4,2x),1p,3(1x,e12.4),0p,7x,a,a)') 
     $        ELR_ModlogR , ELR_logRobs , ELR_chi2_logR ,
     $        z , (ELR_chi2_logR/chi2_TOT) , 
     $        ' [logRobs ModlogR' , 
     $        ' chi2_ELR chi2_ELR/chi2_QHR chi2_ELR/chi2_TOT]'
*ELR* .......................................................................
            write (*,*)
            write (*,*)
            write (*,'(a)') '## QHR: Results for best model'
            write (*,'(f9.4,3x,f9.4,19x,a)') 
     $           QHR_log_QH0 , QHR_log_QHeff , 
     $           ' [log_QH0 & log_QHeff - in photons/sec]'
            write (*,'(f9.4,3x,f9.4,19x,a)') 
     $           (100.*chi2_QHR/chi2_TOT) , 
     $           (100.*chi2_Opt/chi2_TOT) ,
     $           ' [%-chi2: QHR/TOT & Opt/TOT]'
            write (*,'(1p,3(e11.5,2x),0p,1x,a)') 
     $           chi2_QHR , chi2_Opt , chi2_TOT , 
     $           ' [chi2_QHR chi2_Opt chi2_TOT]'
            write (*,*)
          write (*,'(a,a)') '# lambda[A]   q_lambda     logY_obs      ',
     $        'ModlogY       chi2_Y    chi2_Y/chi2_Opt  chi2_Y/chi2_TOT'
            do iY=1,NQHR_Ys
               write (*,2011) iY , QHR_lambda(iY) , QHR_qlambda(iY) ,
     $              QHR_logY_obs(iY) , QHR_ModlogY(iY) , QHR_chi2_Y(iY), 
     $              (QHR_chi2_Y(iY)/chi2_Opt),(QHR_chi2_Y(iY)/chi2_TOT)
            end do
            write (*,*)
            write (*,2012)  (int(QHR_lambda(iY)),iY=1,NQHR_Ys)
            do j=1,N_base
               write (*,2013) j , age_base(j) , Z_base(j) , ind_exAV(j), 
     &              (100.*x_min(j)) , AV_tot(j) , 
     &              QHR_qH40(j) ,  QHR_Q2Lnorm40(j),
     &              (100.*QHR_x_QH0(j)) , (100.*QHR_x_QHeff(j)) , 
     &              (100.*QHR_x_Y(j,iY),iY=1,NQHR_Ys)
            end do
         end if

      end if
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR



c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO
c PHO                        PHO-RELATED OUTPUT          .                PHO
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO

*PHO* Produce/re-compute/compile PHO-related things for output
*     The if just below skips PHOc-output for non-PHOc fits (IsPHOcOn =
*     0). Note that for IsPHOcOn = -1 predicted PHO things are saved,
*     even though no actual PHO-data was used.
      if ((iPHO_IsPHOcOn.EQ.1).OR.(iPHO_IsPHOcOn.EQ.-1)) then
         write (2,*)
         write (2,*)
         write (2,*)
         write (2,'(35(a))') ('#',j=1,35)
         if (IsPHOcOn.EQ.-1) then
            write (2,'(a,a)') '##    PHOc-FIT RELATED OUTPUT    ##' ,
     $           ' (Predicted things only!) '
         else
            write (2,'(a)') '##    PHOc-FIT RELATED OUTPUT    ##'
         end if
         write (2,'(35(a))') ('#',j=1,35)

         write (2,*)
         write (2,'(a)') '## PHO: Input info'
         write (2,'(a40,a)')     arq_ETCinfo  , ' [arq_ETCinfo]'
         write (2,'(f9.4,31x,a)') LumDistInMpc , ' [LumDist (Mpc)]'
         write (2,'(f9.6,31x,a)') PHO_redshift , ' [PHO_redshift]'
         write (2,'(f9.4,31x,a)') PHO_GlobalChi2ScaleFactor , 
     $        ' [PHO_GlobalChi2ScaleFactor]'
         write (2,'(i3,37x,a)') NPHO_Ys , 
     $        ' [NPHO_Ys = number of PHO-bands included in fit]'

         write (2,*)
         write (2,'(a,4x,a,a,a)') '#name/code  filter_file' ,
     $        ' logY_TOT  YFrac2Model  ErrlogY  RangelogY' ,
     $        ' Chi2ScaleFactor logY_obs    logY_low     logY_upp'
 2020    format(a10,2x,a15,5(f9.4,2x),3(2x,f9.4,2x))
         do iY=1,NPHO_Ys
            write (2,2020) PHO_name(iY) , 
     $           PHO_Filter_file(iY) , 
     $           PHO_logY_TOT(iY) , PHO_YFrac2Model(iY) ,
     $           PHO_ErrlogY(iY) , PHO_RangelogY(iY) ,
     $           PHO_Chi2ScaleFactor(iY) , PHO_logY_obs(iY) ,
     $           PHO_logY_low(iY) ,  PHO_logY_upp(iY)
         end do

         write (2,*)
         write (2,'(a)') '## PHO: Results for best model'
         write (2,'(f9.4,3x,f9.4,19x,a)') 
     $        (100.*chi2_PHO/chi2_TOT) , 
     $        (100.*chi2_Opt/chi2_TOT) ,
     $        ' [%-chi2: PHO/TOT & Opt/TOT]'
         write (2,'(1p,3(e11.5,2x),0p,1x,a)') 
     $        chi2_PHO , chi2_Opt , chi2_TOT , 
     $        ' [chi2_PHO chi2_Opt chi2_TOT]'

         write (2,*)
         write (2,'(a,a,a,a)') '#  name/code',
     $        '  MeanLamb StdDevLamb  q_MeanLamb', 
     $        '   logY_obs      ModlogY    ' ,
     $        '   chi2_Y   chi2_Y/chi2_Opt  chi2_Y/chi2_TOT'
         do iY=1,NPHO_Ys
            write (2,2021) iY , PHO_name(iY) , 
     $           PHO_MeanLambda(iY) , PHO_StdDevLambda(iY) ,
     $           PHO_q_MeanLambda(iY) ,
     $           PHO_logY_obs(iY) , PHO_ModlogY(iY) , PHO_chi2_Y(iY) , 
     $           (PHO_chi2_Y(iY)/chi2_Opt) , (PHO_chi2_Y(iY)/chi2_TOT)
         end do
 2021    format (i2,1x,a10,1x,f8.2,3x,f8.2,3x,f7.4,
     $        2(3x,f10.4),3x,1p,e12.4,2x,
     $        e12.4,5x,e12.4)

*     Repeat base info (age, Z, AV...) now with the PHO_x_Y pop vector:
*     PHO_x_Y(j,iY) = %-contribution of j-th component to filter iY
         write (2,*)
         write (2,2022)  (iY,iY=1,NPHO_Ys)
 2022    format ('# j  age_j(yr)     Z_j      exAV?  ',
     &        ' x_j(%)    AV_tot  |',
     &        ' %-Y:',100(i5,1x))

         do j=1,N_base
            write (2,2023) j , age_base(j) , Z_base(j) , ind_exAV(j) , 
     &           (100.*x_min(j)) , AV_tot(j) , 
     &           (100.*PHO_x_Y(j,iY),iY=1,NPHO_Ys)
         end do
 2023    format (i3,2x,1p,e12.6,2x,0p,f7.5,2x,i3,
     &        2x,f10.4,2x,f7.4,
     &        8x,100(f5.1,1x))


*     Screen output of PHO-results
         if (i_verbose.GT.0) then
            write (*,*)
            write (*,*)
            write (*,'(a)') '## PHO: Results for best model'
            write (*,'(f9.4,3x,f9.4,19x,a)') 
     $           (100.*chi2_PHO/chi2_TOT) , 
     $           (100.*chi2_Opt/chi2_TOT) ,
     $           ' [%-chi2: PHO/TOT & Opt/TOT]'
            write (*,'(1p,3(e11.5,2x),0p,1x,a)') 
     $           chi2_PHO , chi2_Opt , chi2_TOT , 
     $           ' [chi2_PHO chi2_Opt chi2_TOT]'
            write (*,*)
            write (*,'(a,a,a,a)') '#  name/code',
     $           '  MeanLamb StdDevLamb  q_MeanLamb', 
     $           '   logY_obs      ModlogY    ' ,
     $           '   chi2_Y   chi2_Y/chi2_Opt  chi2_Y/chi2_TOT'
            do iY=1,NPHO_Ys
               write (*,2021) iY ,  PHO_name(iY) , 
     $              PHO_MeanLambda(iY) , PHO_StdDevLambda(iY) ,
     $              PHO_q_MeanLambda(iY) ,
     $              PHO_logY_obs(iY) , PHO_ModlogY(iY) , PHO_chi2_Y(iY), 
     $              (PHO_chi2_Y(iY)/chi2_Opt) ,(PHO_chi2_Y(iY)/chi2_TOT)
            end do
            write (*,*)
            write (*,2022)  (iY,iY=1,NPHO_Ys)
            do j=1,N_base
               write (*,2023) j , age_base(j) , Z_base(j) , ind_exAV(j),
     &              (100.*x_min(j)) , AV_tot(j) , 
     &              (100.*PHO_x_Y(j,iY),iY=1,NPHO_Ys)
            end do
         end if

      end if
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO


c     Done with output :)
      close(2)


*     Moved this time-counters blablabla to the very end (= here!)
      if (i_verbose.GT.0) then
         write (*,*)
         write (*,'(5(a,i8,2x),a,a50)') 
     $        '*CidTIME>  FF= ' , idt_FF ,
     $        ' RF= ' , idt_RF  , ' BI= ' , idt_BI ,
     $        ' EX0= ', idt_EX0 , ' ALL= ', idt_all    ,
     $        ' sec - fit = ' , arq_out
         write (*,'(5(a,f8.3,2x),a,a50)') 
     $        '*TotTIME>  FF= ' , wdt_FF_TotTime ,
     $        ' RF= ' , wdt_RF_TotTime  , ' BI= ' , wdt_BI_TotTime ,
     $        ' EX0= ', wdt_EX0_TotTime , ' ALL= ', wdt_TotTime    ,
     $        ' sec - fit = ' , arq_out
         write (*,'(5(a,f8.3,2x),a,a50)') 
     $        '*UsrTIME>  FF= ' , wdt_FF_UsrTime ,
     $        ' RF= ' , wdt_RF_UsrTime  , ' BI= ' , wdt_BI_UsrTime ,
     $        ' EX0= ', wdt_EX0_UsrTime , ' ALL= ', wdt_UsrTime    ,
     $        ' sec - fit = ' , arq_out
         write (*,'(5(a,f8.3,2x),a,a50)') 
     $        '*SysTIME>  FF= ' , wdt_FF_SysTime ,
     $        ' RF= ' , wdt_RF_SysTime  , ' BI= ' , wdt_BI_SysTime ,
     $        ' EX0= ', wdt_EX0_SysTime , ' ALL= ', wdt_SysTime    ,
     $        ' sec - fit = ' , arq_out
         write (*,*)
      end if
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c +                  DONE WITH POST-PROCESSING & OUTPUT:)                   +
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+



      end
c @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c @                           END OF MAIN CODE                              @
c @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@





c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c +                       A BUNCH OF SIMPLE ROUTINES                        +
c +                                                                         +
c + iCheckExistingOutFiles - added by ElCid@Sanchica - 15/Ago/2011          +
c + ReadSpecWithHeader                                                      +
c + calc_Fnorm                                                              +
c + Stats_F                                                                 +
c + Sum_of_x                                                                +
c + F_adev                                                                  +
c + Interpol_Spec - not really used!!        <==> *PHO* IS NOW USING IT!    +
c + RebinSpec                                                               +
c + GaussSmooth                                                             +
c + Calc_Mean_Age_and_Z                                                     +
c + Update_AV_array                                                         +
c + Clip_Spectrum                                                           +
c + ConvertLight2Mass                                                       +
c + F_SSP_chi2_AV (added on 20/Jan/2007)                                    +
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+



c ###########################################################################
c     Checks if file out_dir/arq_out already exists (ih which case the
c     function returns 1), or not (function returns 0).
c     If the file does not exist them create a tmp one.
c
c     ElCid@Sanchica - 15/Ago/2011
c
c     Here is the old comment/abstract about this bit:
!
!Cid08Fev07-FIX! Skip runs already made!  New trick: Figure out if
!     arq_out already exists. If so, abort run & go to next one! Useful
!     to run grids in "modo-burro". If arq_out does NOT exit, then write
!     a tmp arq_out (to be overwritten at the end of the run) just to
!     avoid that another CPU takes this same job while the code is
!     running. 
!     OBS: Flag i_SkipExistingOutFile decides whether this test is done.

      FUNCTION iCheckExistingOutFiles(out_dir,arq_out)
      implicit real (a-h,j-m,o-z)
      implicit integer (i,n)
      character*100 out_dir , arq_out , arq_aux

      i_last_non_blank_char = 0
      do i=1,99
         if ((out_dir(i:i+1).EQ.' ').AND.(i_last_non_blank_char.EQ.0))
     &        i_last_non_blank_char = i - 1
      end do
      arq_aux = out_dir(1:i_last_non_blank_char) // arq_out

      open (unit=2,file=arq_aux,status='old',err=500)
      close(2)
      iCheckExistingOutFiles = 1
      return

 500  iCheckExistingOutFiles = 0
      open (unit=2,file=arq_aux,status='unknown')
      write (2,'(a)')
     $     'This file will be overwritten in a few minutes...'
      close(2)
      return
      end
c ###########################################################################


c ###########################################################################
c Reads a spectrum, but skips first N_skip lines if they start with a '#'.
c OBS: '#' must be in the 1st column!!
c Cid@lynx - 15/Aug/2003

      SUBROUTINE ReadSpecWithHeader(arq,N_max,N,l,f)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      real l(N_max) , f(N_max)
      character*100 arq
      character*1 cerquinha

c ** First find-out how many header lines (those which start with '#') to skip
      N_skip = 0
      open (unit=1,file=arq,status='old')
      do i=1,N_max
         read (1,'(a1)') cerquinha
         if (cerquinha.EQ.'#') then 
            N_skip = N_skip + 1
         else
            goto 10
         end if
      end do
 10   close(1)

c ** Now reads spectrum, skiping 1st N_skip lines
      open (unit=1,file=arq,status='old')
      do i=1,N_skip
         read (1,*)
      end do
      do i=1,N_max
         read (1,*,END=20) l(i) , f(i)
      end do
      write (*,*) ' [ReadSpecWithHeader] BUG! N > Nl_max!! ' , N_max
 20   close(1)
      N = i - 1

      end
c ###########################################################################


c ###########################################################################
c Interpolates flux at l_norm Angstrons in a f(l) flux vector. 
c Used for normalization purposes.
c Cid@Lynx - Aug/09/2003

      SUBROUTINE calc_Fnorm(N,l,f,l_norm,f_norm)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      real l(N) , f(N)

      do i=1,N-1
         if ((l(i).LE.l_norm).AND.(l(i+1).GE.l_norm)) then
            a = (f(i+1) - f(i)) / (l(i+1) - l(i))
            f_norm = a * (l_norm - l(i+1)) + f(i+1)
            return
         end if
      end do

      write (*,'(3(a,f7.1))')
     $     '[calc_Fnorm] Spectrum does NOT contain l_norm = ' , l_norm ,
     $     ' l goes from ' , l(1) , ' -> ' , l(N)
      stop

      end
c ###########################################################################


c ###########################################################################
c     Compute average, sigma & median flux between l1 & l2.
c     Cid@Lynx - Aug/24/2003
c     Added safety-check when i_tot < 2. Cid@Lagoa - 08/Oct/2005

      SUBROUTINE Stats_F(N,l,f,l1,l2,f_ave,f_sig,f_med)

      parameter (Nl_max = 14000)
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      real l(N) , f(N) , f_aux(Nl_max)

      i_tot = 0
      do i=1,N
         if ((l(i).GE.l1).AND.(l(i).LE.l2)) then
            i_tot = i_tot + 1
            f_aux(i_tot) = f(i)
         end if
      end do

      f_ave = 0.
      f_sig = 0.
      f_med = 0.
      if (i_tot.GT.1) then
         call avevar(f_aux,i_tot,f_ave,f_var)
         f_sig = sqrt(f_var)

         call MDIAN1(f_aux,i_tot,f_med)
      end if

      end
c ###########################################################################


c ###########################################################################
c     Compute average & sigma of f between l1 & l2 but considering only
c     points with w > 0. Silly routine used in the computation of S/N
c     ratios (for output purposes only!).
c     Cid@Lynx - Jan/27/2005
c     Added safety-check when i_tot < 2. Cid@Lagoa - 08/Oct/2005

      SUBROUTINE Stats4SN(N,l,f,w,l1,l2,f_ave,f_sig)

      parameter (Nl_max = 14000)
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      real l(N) , f(N) , w(N) , f_aux(Nl_max)

      i_tot = 0
      do i=1,N
         if ((l(i).GE.l1).AND.(l(i).LE.l2).AND.(w(i).GT.0.)) then
            i_tot = i_tot + 1
            f_aux(i_tot) = f(i)
         end if
      end do

      f_ave = 0.
      f_sig = 0.
      if (i_tot.GT.1) then
         call avevar(f_aux,i_tot,f_ave,f_var)
         f_sig = sqrt(f_var)
      end if

      end
c ###########################################################################


c ###########################################################################
c # Returns the sum of x
c # Cid@Lynx - Feb/21/2005

      FUNCTION Sum_of_x(x,N)

      parameter (Nb_max = 300)
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      real x(Nb_max)

      sum_x = 0.
      do j=1,N
         sum_x = sum_x + abs(x(j))
      end do
      Sum_of_x = sum_x

      end
c ###########################################################################


c ###########################################################################
c # Compute mean absolute relative deviation, in percent. 
c # OBS: only for the Nl_eff points actually used in the synthesis!
c # Note that most of the stuff is passed by a common!
c # Cid@Lynx - Mar/16/2003

c # Minor change: Nl_eff is now computed here.
c # Cid@Lynx - Aug/17/2003

      FUNCTION F_adev(f_syn)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000)
      real l_obs(Nl_max) , f_obs(Nl_max) , weight_obs(Nl_max)
      real f_syn(Nl_max)
      common/Obs/ Nl_obs , l_obs , f_obs , weight_obs

      adev   = 0.
      Nl_eff = 0
      do i=1,Nl_obs
         if (weight_obs(i).GT.0.) then
            adev = adev + abs((f_obs(i) - f_syn(i)) / f_obs(i))
            Nl_eff = Nl_eff + 1
         end if
      end do
      F_adev = 100. * adev / float(Nl_eff)

      return
      end
c ###########################################################################


c ###########################################################################
c # Simple linear interpolation of an input spectrum fo(lo) on to fs(ls),
c # where ls is given by the user.
c # Cid@Lynx - 03/Jan/2003

*PHO* Finally using this routine! I changed it to return a zero when ls
*     is outside the lo range!
*     Cid@Lagoa - 08/Jan/2011

      SUBROUTINE Interpol_Spec(No,lo,fo,Ns,ls,fs)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      real lo(No) , fo(No) , ls(Ns) , fs(Ns)

      do is=1,Ns
         fs(is) = 0.
         do io=1,No-1
            if ((ls(is).GE.lo(io)).AND.(ls(is).LE.lo(io+1))) then
               a = (fo(io+1) - fo(io)) / (lo(io+1) - lo(io))
               b = fo(io) - a * lo(io)
               fs(is) = a * ls(is) + b
               goto 10
            end if
         end do

!         write (*,*) ' [Interpol_Spec] BUG! ' , is , ls(is)
!         stop
 10      do_nothing_on_this_silly_line = 0.
      end do

      end
c ###########################################################################


c ###########################################################################
c # Rebin an input spectrum fo(lo) on to fs(ls) by direct integration...
c # Both lo & ls MUST BE sorted in increasing order.
c # This routine treats a spectrum as a "histogram": 
c #   fo(lo) is taken as constant within the lo bin
c #   lo(i) represents a bin from [lo(i-1)+lo(i)]/2 ==> [lo(i)+lo(i+1)]/2
c # Note that lo (nor ls) does not need to be uniformily spaced.
c # Flux density is preserved.
c # Cid@Lynx - 24/Aug/2003

c # I've introduced a i_FastBC03_FLAG which greatly speeds up rebinning
c # for a BC03 base and when dl = 1 and Ol_obs is within 3322-->9300 Angs. 
c # In this case no interpolation is needed! I just cut-out the piece of the
c # BC03 spectrum which matches the desired ls lambdas (assumed to be spaced
c # by 1 Angs, as the BC03 files).
c # Cid@Lagoa - 7/May/2004

c     I now test whether dls = ls(2) - ls(1) = 1 Angs, the BC03 sampling
c     in the high resolution optical range. If this is not 1 Angs, the
c     FastBC03 rebinning-trick should not be aplied!
c     Cid@Lagoa - 28/Jan/2005

c     I've adapted the fast-rebin trick to work with Base-Rosa spectra
c     (br04*spec files), which cover the 3001 --> 6996 Angs range. To
c     use fast-rebin with BR04 you must set i_FastBC03_FLAG = 2.
c
c     Cid@IAA - 27/May/2005

      SUBROUTINE RebinSpec(i_FastBC03_FLAG,No,lo,fo,Ns,ls,fs)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      real lo(No) , fo(No) , ls(Ns) , fs(Ns)

      if ((i_FastBC03_FLAG.GE.1).AND.((ls(2) - ls(1)).NE.1.0)) then
         write (*,'(a,a,i4,f8.3)') ' [RebinSpec] Tried FastBC03-' ,
     &        'rebinning but dl NE 1 Angs. Going for slow rebinning!' , 
     &        i_FastBC03_FLAG , ls(2) - ls(1)
      end if

      if ((i_FastBC03_FLAG.EQ.1).AND.((ls(2) - ls(1)).EQ.1.0).AND.
     $     (int(ls(1)).GE.3222).AND.(int(ls(Ns)).LE.9300)) then


c ** Just cut-out a BC03 spectrum to ls. See notes above.
         do is=1,Ns
            ind_in_BC03 = int(ls(is)) - 3322 + 348
            fs(is) = fo(ind_in_BC03)

c ** Paranoid check...
            if (int(ls(is)).NE.int(lo(ind_in_BC03))) then
               write (*,*) ' [RebinSpec] BUG in FastBC03 scheme! ' , 
     &              is , ls(is) , ind_in_BC03 , lo(ind_in_BC03)
               write (*,*) ' [RebinSpec] Try setting',
     $              ' i_FastBC03_FLAG = 0 in arq_config' 
               stop
            end if
         end do


! Same thing but for br04*spec, which cover 3001-->6996 Angs.
      else if ((i_FastBC03_FLAG.EQ.2).AND.((ls(2) - ls(1)).EQ.1.0).AND.
     $        (int(ls(1)).GE.3001).AND.(int(ls(2)).LE.6996)) then


c ** Just cut-out a BR04spectrum to ls. See notes above.
         do is=1,Ns
            ind_in_BR04 = int(ls(is)) - 3001 + 1
            fs(is) = fo(ind_in_BR04)

c ** Paranoid check...
            if (int(ls(is)).NE.int(lo(ind_in_BR04))) then
               write (*,*) ' [RebinSpec] BUG in FastBR04 scheme! ' , 
     &              is , ls(is) , ind_in_BR04 , lo(ind_in_BR04)
               write (*,*) ' [RebinSpec] Try setting',
     $              ' i_FastBC03_FLAG = 0 in arq_config' 
               stop
            end if
         end do


! Slow-rebinning.
      else


c ** Loop over ls-bins
      do is=1,Ns


c ** Define bin limits corresponding to ls(is)
c ** (1st & last bins are forced to be simmetric)
         if (is.EQ.1) then 
            ls_low = ls(is) - (ls(is+1) - ls(is)) / 2.
         else 
            ls_low = (ls(is) + ls(is-1)) / 2.
         end if

         if (is.EQ.Ns) then 
            ls_upp = ls(is) + (ls(is) - ls(is-1)) / 2.
         else 
            ls_upp = (ls(is) + ls(is+1)) / 2.
         end if


c ** Find lo points within current ls bin & integrate flux
         sum_bin = 0.
         do io=1,No


c ** Define bin limits corresponding to lo(io)
c ** (1st & last bins are forced to be simmetric)
            if (io.EQ.1) then 
               lo_low = lo(io) - (lo(io+1) - lo(io)) / 2.
            else 
               lo_low = (lo(io) + lo(io-1)) / 2.
            end if
            
            if (io.EQ.No) then 
               lo_upp = lo(io) + (lo(io) - lo(io-1)) / 2.
            else 
               lo_upp = (lo(io) + lo(io+1)) / 2.
            end if


c ** Check if ls & lo bins overlap. If not, skip this lo-bin. 
c ** If yes, integrate fo(lo) contribution to current ls bin.
            if ((ls_upp.LE.lo_low).OR.(ls_low.GE.lo_upp)) then

               ls_bin_does_not_overlap_with_current_lo_bin = 0.

            else

c ** set l_inf = MAX(lo_low,ls_low) and l_sup = MIN(lo_upp,ls_upp)
               l_inf = lo_low
               if (ls_low.GT.lo_low) l_inf = ls_low

               l_sup = lo_upp
               if (ls_upp.LT.lo_upp) l_sup = ls_upp

               sum_bin= sum_bin+ fo(io) * (l_sup - l_inf)
            end if 

         end do


c ** sum_bin= integrated flux in ls-bin & fs = sum_bin per delta-ls
         fs(is) = sum_bin / (ls_upp - ls_low)


c ** Screen check...
!         write (*,11) 'Rebin> ' , is , ls_low,ls(is),ls_upp,fs(is),sum
! 11      format (a,i5,2x,3(2x,f7.2),1p,2(2x,e9.3))
      end do


      end if


      end
c ###########################################################################




c ###########################################################################
c # Smooth an input spectrum fo(lo) with a Gaussian centered at velocity
c # v = v0 and with sig as velocity dispersion. The output is fs(ls).
c # The convolution is done with a BRUTE FORCE integration of n_u points
c # from u = -n_sig to u = +n_sig, where u = (v - v0) / sig.
c # OBS: Assumes original spectrum fo(lo) is EQUALLY SPACED in lambda!!
C #      Units of v0 & sig = km/s
c # Cid@Lynx - 05/Dec/2002

c # Added n_u as a parameter!
c # Cid@Lynx - 16/Mar/2003

*ETC* !AKI! gfortran does not like do u=u_low,u_upp,du! But works! Fix
*     after other tests are done! Cid@Lagoa - 11/Jan/2011

*     This version of GaussSmooth replaces the u-loop by one in an
*     integer index i_u. The differences wrt are minute (< 1e-5) and due
*     to precision...
*
*      ElCid@Sanchica - 16/Ago/2011

*ATT* Promoted to new oficial GaussSmooth without much testing!
*	  Previous version removed to skip compilation warnings...
*	ElCid@Sanchica - 08/Oct/2011

*ATT* Found an array bounds-bug when compiling with gfortran ... -fbounds-check
*	  Ex: Fortran runtime error: Index '3402' of dimension 1 of array 'fo' above upper bound of 3401
*	  Happens only when sig is very large. The old routine is kept below, and documents the bug briefly.
*	ElCid@Sanchica - 05/Nov/2011

      SUBROUTINE GaussSmooth(No,lo,fo,v0,sig,n_u,Ns,ls,fs)

** Definitions
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      real lo(No) , fo(No) , ls(Ns) , fs(Ns)
      pi = 4. * atan(1.)
      c  = 2.997925e5

c ** Check n_u
      if (n_u.LT.5) then
         write (*,*) ' [GaussSmooth] nu too small! ' , n_u
         stop
      end if

c ** Parameters for the brute-force integration in velocity
      n_sig = 6
      u_low = -n_sig
      u_upp = n_sig
      du    = (u_upp - u_low) / float(n_u - 1)
      d_lo  = lo(2) - lo(1)

c ** loop over ls, convolving fo(lo) with a gaussian
      do is=1,Ns

c ** reset integral of {fo[ll = ls/(1+v/c)] * Gaussian} & start integration
         sum_fg = 0.
         do i_u=1,n_u

c ** define velocity & lambda corresponding to u
            u  = u_low + float(i_u - 1) * du
            v  = v0 + abs(sig) * u
            ll = ls(is) / (1. + v/c)

c ** find fo flux for lambda = ll
            ind = (ll - lo(1)) / d_lo + 1
            if (ind.LE.1) then 
               ff = fo(1)
            else if (ind.GE.No) then 
               ff = fo(No)
            else
               a  = (fo(ind+1) - fo(ind)) / (lo(ind+1) - lo(ind))
               b  = fo(ind) - a * lo(ind)
               ff = a * ll + b
            end if

c ** smoothed spectrum
            sum_fg = sum_fg + ff * exp(-(u**2/2.)) 
            fs(is) = sum_fg * du / sqrt(2. * pi)
         end do

      end do

      end


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COPY OF THE OLD-BUGGED routine
! The ...else if (ll.GT.lo(No)) then... should be ...(ll.GE.lo(No))!
! This error caused array bounds bugs in some runs!
! ElCid@Sanchica - 05/Nov/2011

      SUBROUTINE OldAndBugged_GaussSmooth(No,lo,fo,v0,sig,n_u,Ns,ls,fs)

** Definitions
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      real lo(No) , fo(No) , ls(Ns) , fs(Ns)
      pi = 4. * atan(1.)
      c  = 2.997925e5

c ** Check n_u
      if (n_u.LT.5) then
         write (*,*) ' [GaussSmooth] nu too small! ' , n_u
         stop
      end if

c ** Parameters for the brute-force integration in velocity
      n_sig = 6
      u_low = -n_sig
      u_upp = n_sig
      du    = (u_upp - u_low) / float(n_u - 1)
      d_lo  = lo(2) - lo(1)

c ** loop over ls, convolving fo(lo) with a gaussian
      do is=1,Ns

c ** reset integral of {fo[ll = ls/(1+v/c)] * Gaussian} & start integration
         sum_fg = 0.
         do i_u=1,n_u

c ** define velocity & lambda corresponding to u
            u  = u_low + float(i_u - 1) * du
            v  = v0 + abs(sig) * u
            ll = ls(is) / (1. + v/c)

c ** find fo flux for lambda = ll
            if (ll.LT.lo(1)) then 
               ff = fo(1)
!!THIS IS THE BUGGED LINE! SHOULD BE ll.GE.lo(No)!!
            else if (ll.GT.lo(No)) then 
               ff = fo(No)
            else
               ind = (ll - lo(1)) / d_lo + 1
               a   = (fo(ind+1) - fo(ind)) / (lo(ind+1) - lo(ind))
               b   = fo(ind) - a * lo(ind)
               ff  = a * ll + b
            end if

c ** smoothed spectrum
            sum_fg = sum_fg + ff * exp(-(u**2/2.)) 
            fs(is) = sum_fg * du / sqrt(2. * pi)
         end do

      end do

      end
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c ###########################################################################


c ###########################################################################
c Silly mean age & Z - just for screen output/check purposes
c Cid@lynx - 04/Feb/04

      SUBROUTINE Calc_Mean_Age_and_Z(N_base,x,a_age,s_age,a_met,s_met)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nb_max = 300)
      real x(Nb_max) , age_base(Nb_max) , Z_base(Nb_max)
      common/Age_and_Met/ age_base , Z_base

      sum_t  = 0.
      sum_t2 = 0.
      sum_z  = 0.
      sum_z2 = 0.
      sum_x  = 0.
      do i=1,N_base
         sum_t  = sum_t  + x(i) * log10(age_base(i) + 1.)
         sum_t2 = sum_t2 + x(i) * (log10(age_base(i) + 1.))**2
         sum_z  = sum_z  + x(i) * Z_base(i)
         sum_z2 = sum_z2 + x(i) * (Z_base(i))**2
         sum_x  = sum_x  + x(i)
      end do

      a_age = sum_t / sum_x
      s_age = sqrt(abs((sum_t2/sum_x) - a_age**2))

      a_met = sum_z / sum_x
      s_met = sqrt(abs((sum_z2/sum_x) - a_met**2))

      end
c ###########################################################################


c ###########################################################################
c #    COMPUTE "TOTAL MASS" & MASS FRACTIONS IN COMPONENTS j=1...N_base     #
c ###########################################################################
c     Flux_tot is the total **dereddened** synthetic (WITHOUT
c     kinematics!) flux at l_norm, in the same units as in the input
c     spectrum.  Mass_tot is the total stellar mass. It's units must be
c     fixed by the user a posteriori!  For instance, if (a) the input
c     spectrum comes in units of 1e-17 erg/s/cm2/A, (2) the distance to
c     the galaxy (in cm) is d, and (3) the base spectra were given in
c     units of L_sun/A/M_sun (as is the case for the BC03 spectra), then
c     the actual mass is given by Mass_tot * 1e-17 * 4 * pi * d**2 /
c     3.826e33.
c
c     ATT: There are 2 versions of Mass_tot: Mini_tot & Mcor_tot! Why?
c     The original BC03 SSP spectra come normalized to a total INITIAL
c     mass of 1 M_sun. As time goes by, however, part of this stellar
c     mass goes back to the ISM. This implies that the mass (Mini_tot)
c     and mass-fractions (mu_ini(j)) are NOT the current ones! To
c     correct for the returned mass, I use the array Mstars (read fron
c     the base file), which says what fraction of the initial mass in
c     population # j is still in stars at the age of this same
c     component. The corrected total mass and mass-fractions are given
c     in Mcor_tot & mu_cor(j).
c     
c     OBS1: Note that a Power-Law component (which, if present, is
c     component # N_base) is excluded from the mass calculations (for
c     obvious reasons)!
c     OBS2: Masses & mass fractions may not make any sense at all if,
c     for instance, the base is made of template galaxies or normalized
c     models...
c     OBS3: Note that now AV is population dependent: AV_tot = AV + exAV!
c
c     Cid@Lagoa - 21/April/2007


      SUBROUTINE ConvertLight2Mass(N_base,x,AV_tot,q_norm,fobs_norm
     $     ,fbase_norm,i_FitPowerLaw,Mstars,Flux_tot,Mini_tot
     $     ,Mcor_tot,mu_ini,mu_cor)

c     Definitions
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nb_max = 300)

c     L2Mass (@ l_norm) in the base
      real fbase_norm(Nb_max) 

c     Population vector arrays in flux (x) & mass (mu_ini) fractions
      real x(Nb_max) , mu_ini(Nb_max)

c     Mstars = array which says what fraction of the INITIAL stellar
c     mass in population # j is still in stars at age(j). mu_cor() are
c     the corrected mass-fractions.
      real Mstars(Nb_max) , mu_cor(Nb_max)

c     AV_array: AV_tot(j) = population dependent extinction
      real AV_tot(Nb_max)


c     Computes Flux_tot & total initial & actual (corrected) masses.
      Flux_tot = 0.
      Mini_tot = 0.
      Mcor_tot = 0.
      do j=1,N_base
         flux_j = x(j) * fobs_norm * 10.**(0.4 * AV_tot(j) * q_norm)
         mini_j = flux_j / fbase_norm(j)
         mcor_j = flux_j / fbase_norm(j) * Mstars(j)
         if ((j.EQ.N_base).AND.(i_FitPowerLaw.EQ.1)) then 
            mini_j = 0.
            mcor_j = 0.
         end if
         Flux_tot = Flux_tot + flux_j
         Mini_tot = Mini_tot + mini_j
         Mcor_tot = Mcor_tot + mcor_j
      end do

c     Compute mass fractions (mu_ini & mu_cor) associated with each base
c     component.
      do j=1,N_base
         flux_j = x(j) * fobs_norm * 10.**(0.4 * AV_tot(j) * q_norm)
         mini_j = flux_j / fbase_norm(j)
         mcor_j = flux_j / fbase_norm(j) * Mstars(j)
         if ((j.EQ.N_base).AND.(i_FitPowerLaw.EQ.1)) then 
            mini_j = 0.
            mcor_j = 0.
         end if

         mu_ini(j)  = mini_j / Mini_tot
         mu_cor(j)  = mcor_j / Mcor_tot
      end do

      end
c ###########################################################################



c ###########################################################################
c #                      SIGMA-CLIPPING (excluding outliers)                #
c ###########################################################################
c     Given an observed (f_obs) and a model spectrum (f_syn)
c     mask-out/clip points with large residuals. 3 clipping methods are
c     offered: clip_method_option = "ABSRES", "RELRES" & "NSIGMA". In
c     the ABSRES method we clip points with |f_obs-f_syn| "ABSOLUTE"
c     residuals > sig_clip_threshold X the rms (f_obs - f_syn) over all
c     (non-masked) lambdas. [This was the only option in STARLIGHT up to
c     March/2005]. "RELRES" clipping consists of excluding points whose
c     RELATIVE residual |f_obs-f_syn|/error deviates by more than
c     sig_clip_threshold X the rms of this same quantity over all non
c     masked points. Finally, "NSIGMA" clipping consists of excluding
c     points whose |f_obs-f_syn| residual deviates by more than
c     sig_clip_threshold X the LOCAL error in f_obs. This scheme should
c     clip less points than ABSRES, since it may allow points with large
c     errors to be kept in f_obs for fitting. In all cases, exclusions
c     are added to weight_obs array, and clipped points are flaged as
c     f_mask(i) = -1 (just to mark them for output/control purposes).
c
c     Cid@Lagoa - 24/Mar/2005

      SUBROUTINE Clip_Spectrum(clip_method_option,f_syn
     $     ,sig_clip_threshold,f_mask,NOl_eff,Nl_eff,Ntot_cliped
     $     ,iCLIPBUG_flag)

               
c     Definitions
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
      real l_obs(Nl_max) , f_obs(Nl_max)  , weight_obs(Nl_max)
      real f_syn(Nl_max) , f_mask(Nl_max) , f_aux1(Nl_max)
      common/Obs/ Nl_obs , l_obs , f_obs , weight_obs
      character*6 clip_method_option


c     Clip-spectrum according to clip_method_option.
      if (clip_method_option.EQ.'NSIGMA') then


c     "NSIGMA" clipping! Compute # of points to be clipped, assign zero
c     weight & flag cliped points with f_mask = -1. Corrects
c     Nl_eff. Store original weight temporarely in f_aux1 (in case I
c     have to undo it below!)
         Ntot_cliped = 0
         do i=1,Nl_obs
            f_aux1(i) = weight_obs(i)
            if (weight_obs(i).GT.0.) then
               rel_resid = (f_syn(i) - f_obs(i)) * weight_obs(i)
               if (abs(rel_resid).GT.sig_clip_threshold) then
                  weight_obs(i) = 0.
                  f_mask(i)     = -1.
                  Ntot_cliped   = Ntot_cliped + 1
                  Nl_eff        = Nl_eff - 1
               end if
            end if
         end do


      else if (clip_method_option.EQ.'RELRES') then


c     "RELRES" clipping!  Store RELATIVE residuals (= residue/error for
c     each lambda) in f_aux1. Compute mean & std dev of relative
c     residuals (only for non masked points).
         i_aux = 0
         do i=1,Nl_obs
            if (weight_obs(i).GT.0.) then
               i_aux = i_aux + 1
               f_aux1(i_aux) = (f_syn(i) - f_obs(i)) * weight_obs(i)
            end if
         end do
         call avevar(f_aux1,i_aux,ave_rel_resid,var_rel_resid)
         sig_rel_resid = sqrt(var_rel_resid)

c     Compute # of points to be clipped, assign zero weight & flag
c     cliped points with f_mask = -1. Corrects Nl_eff. Store original
c     weight temporarely in f_aux1 (in case I have to undo it below!)
         Ntot_cliped = 0
         do i=1,Nl_obs
            f_aux1(i) = weight_obs(i)
            if (weight_obs(i).GT.0.) then
               rel_resid = (f_syn(i) - f_obs(i)) * weight_obs(i)
               if (abs(rel_resid).GT.(sig_clip_threshold*sig_rel_resid))
     $              then
                  weight_obs(i) = 0.
                  f_mask(i)     = -1.
                  Ntot_cliped   = Ntot_cliped + 1
                  Nl_eff        = Nl_eff - 1
               end if
            end if
         end do


      else if (clip_method_option.EQ.'ABSRES') then


c     "ABSRES" clipping!  Store residuals in f_aux1. Compute mean & std
c     dev of residuals (only for non masked points).
         i_aux = 0
         do i=1,Nl_obs
            if (weight_obs(i).GT.0.) then
               i_aux = i_aux + 1
               f_aux1(i_aux) = f_syn(i) - f_obs(i)
            end if
         end do
         call avevar(f_aux1,i_aux,ave_resid,var_resid)
         sig_resid = sqrt(var_resid)

c     Compute # of points to be clipped, assign zero weight & flag
c     cliped points with f_mask = -1. Corrects Nl_eff. Store original
c     weight temporarely in f_aux1 (in case I have to undo it below!)
         Ntot_cliped = 0
         do i=1,Nl_obs
            f_aux1(i) = weight_obs(i)
            if (weight_obs(i).GT.0.) then
               resid = f_syn(i) - f_obs(i)
               if (abs(resid).GT.(sig_clip_threshold*sig_resid)) then
                  weight_obs(i) = 0.
                  f_mask(i)     = -1.
                  Ntot_cliped   = Ntot_cliped + 1
                  Nl_eff        = Nl_eff - 1
               end if
            end if
         end do


      else


         if (clip_method_option.NE.'NOCLIP') then
            write (*,'(a,a)') '[Clip_Spectrum] Bad clip_method_option: '
     $           ,clip_method_option
            write (*,'(a)')   '[Clip_Spectrum] Assuming method = NOCLIP'
            clip_method_option = 'NOCLIP'
         end if


      end if


c     For safety, I first test whether Nl_eff has not become
c     ridiculously small due to excessive clipping. If that happens,
c     clipping is aborted (ie, I restore weight_obs & Nl_eff) and a
c     clip-bug-flag is turned on!
!2BEDONE! ATT: This minimum Nl_eff could be more smarly scaled ...
      if (Nl_eff.LT.50) then
         write (*,'(a,a,i5)') '[Clip_Spectrum] Too few points to model!'
     $        ,' clipping aborted ' , Nl_eff
         do i=1,Nl_obs
            weight_obs(i) = f_aux1(i)
         end do
         Nl_eff        = NOl_eff
         Ntot_cliped   = 0
         iCLIPBUG_flag = 1
      end if

c     Warn if > 20% of the data is being clipped...
      if (Ntot_cliped.GT.(0.2*NOl_eff)) then
         write (*,'(a,a,i5,3x,i5)')
     $        '[Clip_Spectrum] WARNING ! > 20% of original data',
     $        ' are being clipped-out... ' , Ntot_cliped , NOl_eff
      end if


      end 
c ###########################################################################


c ###########################################################################
c     Compute the chi2 between a spectrum in f_SSP (passes via common!)
c     and an observed spectrum f_obs, as a function of the
c     extinction. The normalization x_SSP of the fitted-spectrum is
c     computed analyticaly (it's the one satisfying dchi2/dx = 0 for the
c     given AV), and returned via common).  Useful to fit single
c     population spectra.
c
c     Cid@Lagoa - 20/Jan/2007

*ETC* ATT: These SSP-fits do NOT include ETC-constraints!
*     Cid@Lagoa - 06/Jan/2011

      FUNCTION F_SSP_chi2_AV(AV_SSP)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000)
      real l_obs(Nl_max) , f_obs(Nl_max) , weight_obs(Nl_max)
      common/Obs/ Nl_obs , l_obs , f_obs , weight_obs
      common/Red_factor/ red_term(Nl_max)
      common/SSP_chi2_stuff/ x_SSP , f_SSP(Nl_max)

      sum1 = 0.
      sum2 = 0.
      do i=1,Nl_obs
         if (weight_obs(i).GT.0.) then
            z1 = weight_obs(i)**2
            z2 = 10.**(AV_SSP * red_term(i)) * f_SSP(i)
            sum1 = sum1 + z1 * z2 * f_obs(i)
            sum2 = sum2 + z1 * z2 * z2
         end if
      end do
      x_SSP = sum1 / sum2

      chi2_SSP = 0.
      do i=1,Nl_obs
         if (weight_obs(i).GT.0.) then
            z1 = 10.**(AV_SSP * red_term(i)) * x_SSP * f_SSP(i)
            chi2_SSP = chi2_SSP + (weight_obs(i) * (f_obs(i) - z1))**2
         end if
      end do
 
      F_SSP_chi2_AV = chi2_SSP

      return
      end
c ###########################################################################
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+





c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c +                    ROUTINES FROM NUMERICAL-RECIPES                      +
c +                                                                         +
c + avevar, mdian1, sort, ran2, golden, mnbrak, gasdev, indexx              +
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar(data,n,ave,var)
      INTEGER n
      REAL ave,var,data(n)
      INTEGER j
      REAL s,ep
      ave=0.0
      do 11 j=1,n
        ave=ave+data(j)
11    continue
      ave=ave/n
      var=0.0
      ep=0.0
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        var=var+s*s
12    continue
      var=(var-ep**2/n)/(n-1)
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MDIAN1(X,N,XMED)
      DIMENSION X(N)
      CALL SORT(N,X)
      N2=N/2
      IF(2*N2.EQ.N)THEN
        XMED=0.5*(X(N2)+X(N2+1))
      ELSE
        XMED=X(N2+1)
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SORT(N,RA)
      DIMENSION RA(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION golden(ax,bx,cx,f,tol,xmin)
      REAL golden,ax,bx,cx,tol,xmin,f,R,C
      EXTERNAL f
      PARAMETER (R=.61803399,C=1.-R)
      REAL f1,f2,x0,x1,x2,x3
      x0=ax
      x3=cx
      if(abs(cx-bx).gt.abs(bx-ax))then
        x1=bx
        x2=bx+C*(cx-bx)
      else
        x2=bx
        x1=bx-C*(bx-ax)
      endif
      f1=f(x1)
      f2=f(x2)
1     if(abs(x3-x0).gt.tol*(abs(x1)+abs(x2)))then
        if(f2.lt.f1)then
          x0=x1
          x1=x2
          x2=R*x1+C*x3
          f1=f2
          f2=f(x2)
        else
          x3=x2
          x2=x1
          x1=R*x2+C*x0
          f2=f1
          f1=f(x1)
        endif
      goto 1
      endif
      if(f1.lt.f2)then
        golden=f1
        xmin=x1
      else
        golden=f2
        xmin=x2
      endif
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      REAL ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
      REAL dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
CU    USES ran2!
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran2
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran2(idum)-1.
        v2=2.*ran2(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
!MODIFIED-by-ElCid@Sanchica - 16/Ago/2011: "pause" causes silly compilation warnings...
!       if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(jstack.gt.NSTACK) then
           write (*,*) '[indexx] NSTACK too small in indexx... '
           write (*,*) '[indexx] run aborted :-( !!'
           stop
        end if
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+





c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c +                        REDDENING-RELATED ROUTINES                       +
c +                                                                         +
c +   Calc_red_term - computes -0.4 * (A(lambda)/A_V - A(l_norm)/A_V)       +
c +   Cardelli_RedLaw, Calzetti_RedLaw & HYPERZ_RedLaw                      +
c +                                                                         +
c +   To use HYPERZ_RedLaw use red_law_option = 'HZ1' ...'HZ5', where the   +
c +   number = 1 for Allen 76 MW, = 2 for Seaton 79 + Fitzpatrick 86 MW,    +
c +   = 3 for Fitzpatrick 86 LMC, = 4 for Prevot et al 84 and Bouchet et al +
c +  85 SMC & = 5 for Calzetti (astro-ph/9911459).                          +
c +                                                                         +
c +   For Gordon et al (2003, ApJ, 594, 279, Table 4), use red_law_option = +
c +   'GD1' => SMC Bar, 'GD2' => LMC2 SuperShell, 'GD3' => LMC  Average     +
c +                                                                         +
c +   The HZ1...HZ5 options are new to v02.       (Cid@Lagoa - 14/Oct/2005) + 
c +   And the GD1...GD3 options are new to v03.   (Cid@Lagoa - 20/Jan/2007) + 
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

c ###########################################################################
c     Routine Calc_red_term computes the "reddening term", defined by
c
c        red_term(lambda) = -0.4 * (A(lambda) / A_V  -  A(l_norm) / A_V)
c 
c     such that to redden a flux at a certain lambda and keep the
c     normalization at l_norm all we need to do is to multiply the
c     reddening-free flux by
c
c        10**(A_V * red_term(lambda))
c
c     The A(lambda) / A_V (== q) law is from of Cardelli et al (1989) or
c     Calzetti et al. (according to red_law_option)
c 
c     Cid@INAOE - July/06/2004
c
c     + Force red_term = 0 at l_norm (it may be NE 0 due to precision
c     issues...)
c     Cid@Lagoa - 06/Jan/2005
c
c     Introduced HYPERZ laws (red_law_option = 'HZ1' ... 'HZ5')

c     OBS: I didn't quite like their sloppy interpolation at optical
c          wavelengths of the SMC law (eg, a straight line between 5291
c          & 12500 Angs), bur could find nothing better!
c
c     Cid@Lagoa - 14/October/2005

c     Introduced Gordon SMC & LMC extinction curves, plus a general
c     "from-a-file" reddening law read from a file named
c     ExtinctionLaw.red_law_option! (Not fully tested yet!!) The
c     power-law option was included and then removed after implementing
c     the from-a-file option.
c
c     Cid@Lagoa - 07/February/2007

*FIR* Added CCC law (Calzetti + Claus + Cid). Better than CAL at UV and IR
*     Cid@Granada - 22/January/2010


      SUBROUTINE Calc_red_term(red_law_option,l_norm,nobs,lambda,
     &     red_term,q_norm)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)

      parameter (Nl_FromFileMax = 100000)
      real l_FromFile(Nl_FromFileMax) , q_FromFile(Nl_FromFileMax)

      real lambda(nobs) , red_term(nobs)
      real l_norm ,q_norm , R_V , q
      character*3 red_law_option

      if (red_law_option.EQ.'CCM') then 

!     (Removed the obsolete R_V parameter - Cid@Granada - 24/Jan/2010)
         R_V = 3.1
         q_norm = Cardelli_RedLaw(l_norm,R_V)
         do i=1,nobs
            q = Cardelli_RedLaw(lambda(i),R_V)
            red_term(i) = -0.4 * (q - q_norm)
            if (lambda(i).EQ.l_norm) red_term(i) = 0.
         end do

      else if (red_law_option.EQ.'CAL') then 

         q_norm = Calzetti_RedLaw(l_norm)
         do i=1,nobs
            q = Calzetti_RedLaw(lambda(i))
            red_term(i) = -0.4 * (q - q_norm)
            if (lambda(i).EQ.l_norm) red_term(i) = 0.
         end do

      else if (red_law_option.EQ.'CCC') then 

         q_norm = CCC_RedLaw(l_norm)
         do i=1,nobs
            q = CCC_RedLaw(lambda(i))
            red_term(i) = -0.4 * (q - q_norm)
            if (lambda(i).EQ.l_norm) red_term(i) = 0.
         end do

      else if (red_law_option.EQ.'HZ1') then 

         q_norm = HYPERZ_RedLaw(l_norm,1)
         do i=1,nobs
            q = HYPERZ_RedLaw(lambda(i),1)
            red_term(i) = -0.4 * (q - q_norm)
            if (lambda(i).EQ.l_norm) red_term(i) = 0.
         end do

      else if (red_law_option.EQ.'HZ2') then 

         q_norm = HYPERZ_RedLaw(l_norm,2)
         do i=1,nobs
            q = HYPERZ_RedLaw(lambda(i),2)
            red_term(i) = -0.4 * (q - q_norm)
            if (lambda(i).EQ.l_norm) red_term(i) = 0.
         end do

      else if (red_law_option.EQ.'HZ3') then 

         q_norm = HYPERZ_RedLaw(l_norm,3)
         do i=1,nobs
            q = HYPERZ_RedLaw(lambda(i),3)
            red_term(i) = -0.4 * (q - q_norm)
            if (lambda(i).EQ.l_norm) red_term(i) = 0.
         end do

      else if (red_law_option.EQ.'HZ4') then 

         q_norm = HYPERZ_RedLaw(l_norm,4)
         do i=1,nobs
            q = HYPERZ_RedLaw(lambda(i),4)
            red_term(i) = -0.4 * (q - q_norm)
            if (lambda(i).EQ.l_norm) red_term(i) = 0.
         end do

      else if (red_law_option.EQ.'HZ5') then 

         q_norm = HYPERZ_RedLaw(l_norm,5)
         do i=1,nobs
            q = HYPERZ_RedLaw(lambda(i),5)
            red_term(i) = -0.4 * (q - q_norm)
            if (lambda(i).EQ.l_norm) red_term(i) = 0.
         end do

!Cid29Dec06-FIX! Introduced Gordon SMC & LMC extinction curves
      else if ((red_law_option.EQ.'GD1').OR.
     $         (red_law_option.EQ.'GD2').OR.
     $         (red_law_option.EQ.'GD3')) then

         q_norm = GordonTab4_RedLaw(l_norm,red_law_option)
         do i=1,nobs
            q = GordonTab4_RedLaw(lambda(i),red_law_option)
            red_term(i) = -0.4 * (q - q_norm)
            if (lambda(i).EQ.l_norm) red_term(i) = 0.
         end do

      else

!Cid01Feb07-FIX! If none of the red_law_options above worked, then we'll
!     try reading a lambda q_lambda table from a file called:
!     'ExtinctionLaw'.red_law_option and interpolate from it! If that
!     fails, stop run!  The interpolation is in x = 10000 / lambda, and
!     done in a very non-elegant fashion!
         write (*,'(a,a,a,a)')
     $        '[Calc_red_term] --> Will try reading extinction law ' , 
     $        'from file ExtinctionLaw.' , red_law_option , ' ... '

         open (unit=1,file='ExtinctionLaw.'//red_law_option,
     $        status='old',err=501)
         do j=1,Nl_FromFileMax
            read (1,*,END=502) l_FromFile(j) , q_FromFile(j)
         end do
         close(1)
         write (*,'(a)') '[Calc_red_term] OOOPS! Should not be here!'
         write (*,'(a,a)') '[Calc_red_term] ' ,
     $      'Too many points in  your ExtinctionLaw file? Run aborted!'
         stop

 501     write (*,'(a,a)') '[Calc_red_term] OOPS! Bad red_law_option: ',
     $        red_law_option
         write (*,'(a,a,a)') '[Calc_red_term] File ExtinctionLaw.' ,
     $        red_law_option , ' does not exist either!'
         write (*,'(a,a,a)') '[Calc_red_term] Run aborted!'
         stop

 502     Nl_FromFile = j - 1
         close(1)

         q_norm = 0.
         do j=1,Nl_FromFile-1
            if ((l_FromFile(j).LE.l_norm).AND.
     $           (l_norm.LE.l_FromFile(j+1))) then
               x  = 10000. / l_norm
               x1 = 10000. / l_FromFile(j)
               x2 = 10000. / l_FromFile(j+1)
               y1 = q_FromFile(j)
               y2 = q_FromFile(j+1)
               q_norm  = y1 + (y2 - y1) / (x2 - x1) * (x - x1)
            end if
         end do

         do i=1,nobs
            q = 0.
            do j=1,Nl_FromFile-1
               if ((l_FromFile(j).LE.lambda(i)).AND.
     $              (lambda(i).LE.l_FromFile(j+1))) then
                  x  = 10000. / lambda(i)
                  x1 = 10000. / l_FromFile(j)
                  x2 = 10000. / l_FromFile(j+1)
                  y1 = q_FromFile(j)
                  y2 = q_FromFile(j+1)
                  q  = y1 + (y2 - y1) / (x2 - x1) * (x - x1)
               end if
            end do
            red_term(i) = -0.4 * (q - q_norm)
            if (lambda(i).EQ.l_norm) red_term(i) = 0.
         end do

         write (*,'(a,a)') '[Calc_red_term] --> ' ,
     $        'Ok! Read Extinction Law file!'
      end if

      end
c ###########################################################################


c ###########################################################################
c     q = A_lambda / A_V for Cardelli et al reddening law
c     l = lambda, in Angstrons
c     x = 1 / lambda in 1/microns
c     q = a + b / R_V; where a = a(x) & b = b(x)
c     Cid@INAOE - 6/July/2004

      FUNCTION Cardelli_RedLaw(l,R_V)

      real l , a , b , F_a , F_b , x , y , R_V , q

      a = 0.
      b = 0.
      x = 10000. / l


c     Far-Ultraviolet: 8 <= x <= 10 ; 1000 -> 1250 Angs 
      if ((x.GE.8.).AND.(x.LE.10.)) then

         a = -1.073 - 0.628 * (x - 8.) + 0.137 * (x - 8.)**2 - 
     &        0.070 * (x - 8.)**3
         b = 13.670 + 4.257 * (x - 8.) - 0.420 * (x - 8.)**2 +
     &        0.374 * (x - 8.)**3

c     Ultraviolet: 3.3 <= x <= 8 ; 1250 -> 3030 Angs 
      else if ((x.GE.3.3).AND.(x.LE.8.)) then

         F_a = 0
         F_b = 0
         if ((x.GE.5.9).AND.(x.LE.8)) then
            F_a = -0.04473 * (x - 5.9)**2 - 0.009779 * (x - 5.9)**3
            F_b =  0.2130 * (x - 5.9)**2 + 0.1207 * (x - 5.9)**3
         end if

         a =  1.752 - 0.316 * x - 0.104 / ((x - 4.67)**2 + 0.341) + F_a
         b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62)**2 + 0.263) + F_b

c     Optical/NIR: 1.1 <= x <= 3.3 ; 3030 -> 9091 Angs ; 
      else if ((x.GE.1.1).AND.(x.LE.3.3)) then

         y = 10000. / l - 1.82

         a = 1.+ 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3 + 
     &        0.72085 * y**4 + 0.01979 * y**5 - 0.77530 * y**6 + 
     &        0.32999 * y**7
         b = 1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3 - 
     &        5.38434 * y**4 - 0.62251 * y**5 + 5.30260 * y**6 -
     &        2.09002 * y**7  


c     Infrared: 0.3 <= x <= 1.1 ; 9091 -> 33333 Angs ; 
      else if ((x.GE.0.3).AND.(x.LE.1.1)) then

         a =  0.574 * x**1.61
         b = -0.527 * x**1.61

      end if


      q = a + b / R_V


c     Issue a warning if lambda falls outside 1000 -> 33333 Angs range of CCM89
      if ((x.LT.0.3).OR.(x.GT.10.)) then
         q = 0.
!         write (*,*) '[Cardelli_RedLaw] WARNING! ' , 
!     &        'lambda outside valid range: ' , l , x
      end if


      Cardelli_RedLaw = q

      return
      end
c ###########################################################################


c ###########################################################################
c     q = A_lambda / A_V for Calzetti et al reddening law (formula from
c     hyperz-manual).
c     l = lambda, in Angstrons
c     x = 1 / lambda in 1/microns
c     Cid@INAOE - 6/July/2004

      FUNCTION Calzetti_RedLaw(l)

      real l , x , R_V , q

      x   = 10000. / l
      R_V = 4.05


c     UV -> Optical: 1200 -> 6300 Angs
      if ((l.GE.1200.).AND.(l.LE.6300.)) then

         q = (2.659 / R_V) * (-2.156 + 1.509 * x - 0.198 * x**2 + 
     &        0.011 * x**3) + 1.

c     Red -> Infrared
      else if ((l.GE.6300.).AND.(l.LE.22000.)) then

         q = (2.659 / R_V) * (-1.857 + 1.040 * x) + 1.

c     Issue a warning if lambda falls outside 1200->22000 Angs range
      else 

         q = 0.
!         write (*,*) '[Calzetti_RedLaw] WARNING! ' , 
!     &        'lambda outside valid range: ' , l , x

      end if


      Calzetti_RedLaw = q

      return
      end
c ###########################################################################


c ###########################################################################
c     CCC = "Calzetti + Claus + Cid" - law: Returns $2 = q(l) = A(l) / AV; for
c     given l = (in Angs).  Strictly valid from 970 to 22000, but
c     extrapolated to 912 -> 31149 Angs.  UV: I switch from the
c     Leitherer to the Calzetti law at 1846 Angs. This is a bit longer
c     than the nominal limit of the Leitherer law, but is is where the
c     two equations meet and so the transition is smooth. In the IR,
c     extrapolating the Calzetti law to l > 22000 we find that it
c     becomes = 0 at 31149.77 Angs.  I thus set q = 0 for l > 3149.
c	
c     Cid@Granada - 22/Jan/2010

      FUNCTION CCC_RedLaw(l)

      real l , x , R_V , q

      x   = 10000. / l
      R_V = 4.05


c     FUV: 912 -> 1846.  Leitherer et al 2002 (ApJS, 140, 303), eq 14,
c     valid from 970 to 1800 Angs, but extended to 912 -> 1846. OBS:
c     They use E(B_V)_stars in their eq, so I divide by R_V, assumed to
c     be the same as for the Calzetti-law.
      if ((l.GE.912.).AND.(l.LE.1846.)) then

         q = (5.472 + 0.671 * x - 9.218e-3 * x**2 + 2.620e-3 * x**3)/R_V

c     UV -> optical: 1846 -> 6300. Calzetti law.
      else if ((l.GT.1846.).AND.(l.LE.6300.)) then

         q = (2.659 / R_V) * (-2.156 + 1.509 * x - 0.198 * x**2 + 
     &        0.011 * x**3) + 1.

c     Red -> Infrared: 6300 -> 31149. Calzetti law, originally valid up
c     to 22000 Angs, but stretched up to 31149 Angs, where it becomes
c     zero.
      else if ((l.GT.6300.).AND.(l.LE.31149.)) then

         q = (2.659 / R_V) * (-1.857 + 1.040 * x) + 1.

c     q = ZERO for l < 912 & l > 31149 Angs.
      else 

         q = 0.

      end if

      CCC_RedLaw = q

      return
      end
c ###########################################################################


c ###########################################################################
c     q = A_lambda / A_V for one of the 5 reddening laws in HYPERZ.
c     This is just a sloppy DRIVER for the HYPERZ fortran routine red_spec!
c     I call their routine with a wavelength array containing just one
c     lambda (l_aux(1) = l), and tell it to return me the reddened
c     spectrum f_obs_aux(1) for an intrinsic spectrum f_int_aux(1) = 1.0
c     and A_V = 1. From this we derive q = -2.5 * log10(f_obs_aux(1)).
c
c     l = lambda, in Angstrons
c     ilaw_rd=1 : the standard variation of absorption with wavelength 
c                 is used (Allen, 76) - MW
c     ilaw_rd=2 : the analytic fit of Seaton (79) curve by 
c                 Fitzpatrick (86) - MW 
c     ilaw_rd=3 : Fitzpatrick (86) - LMC
c     ilaw_rd=4 : Prevot et al. (1984) and Bouchet et al. (1985) - SMC
c     ilaw_rd=5 : Calzetti (astro-ph/9911459) - SB
c
c     Cid@Lagoa - 14/October/2005

!FIX-DAE2-WARNING-HYPERZ!
!	To avoid the compilation warning (and potential BUGS, like those reported
!	by Marianne Takamiya <takamiya@hawaii.edu> on her Wed, Jun 8, 2011 at 2:32 AM mail,
!	I've changed
!      real l_aux(1) , f_int_aux(1) , f_obs_aux(1)
!	to
!      real l_aux(1300) , f_int_aux(1300) , f_obs_aux(1300)
!	This compatibilizes this FUNCTION with the HYPER-one - at least dae2's gfortran 4.6.1 does
!	not issue any more warnings :-)
!
!	ElCid@Sanchica - 21/Oct/2011

      FUNCTION HYPERZ_RedLaw(l,ilaw_rd)

      real l_aux(1300) , f_int_aux(1300) , f_obs_aux(1300)
      real l , q

      l_aux(1)     = l
      f_int_aux(1) = 1.
      call HYPERZ_red_spec(l_aux,f_int_aux,1.0,f_obs_aux,1,1,ilaw_rd)
      q = -2.5 * log10(f_obs_aux(1))

      HYPERZ_RedLaw = q

      return
      end
c ###########################################################################


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This was "borrowed" from hyperz @ webast.ast.obs-mip.fr/hyperz/
C     I've just copied their red_spec.f (in zphot_src_1.1.tar.gz). I
C     have, though fixed their 'if lambda is in an interval' conditions,
C     as marked below.
C     Cid@Lagoa - 14/October/2005


      SUBROUTINE HYPERZ_red_spec(wl,fl,av,yob,ired,iwavl,law_rd)
c     Reads a spectrum and calculates the new one with internal
c     reddening. 
c     law_rd=1 : the standard variation of absorption with wavelength 
c                is used (Allen, 76) - MW
c     law_rd=2 : the analytic fit of Seaton (79) curve by 
c                Fitzpatrick (86) - MW 
c     law_rd=3 : Fitzpatrick (86) - LMC
c     law_rd=4 : Prevot et al. (1984) and Bouchet et al. (1985) - SMC
c     law_rd=5 : Calzetti (astro-ph/9911459) - SB
c    

!Cid! INCLUDE 'dimension.dec' I've included it here:)
      PARAMETER (mxz=141,mxage=51,mxtyp=15,mxfil=20,mxred=11,mxwl=1300)

      INTEGER ired,iwavl,law_rd
      REAL yob,alambda,wl,fl,wab,amab,aab,av,ff,ff11,ff12,slope,rv
      DIMENSION yob(mxwl),wl(mxwl),fl(mxwl),wab(30),amab(30),aab(30)

c ..... 1. Allen (1976) : wab = lambda in A, amab = k(lambda)/R, 
c          A(lambda) = k(lambda)*E(B-V) = amab*A_V ; R = 3.1

      IF (law_rd.eq.1) THEN
         rv=3.1
         wab(1)=1000.
         amab(1)=4.2
         wab(2)=1110.
         amab(2)=3.7
         wab(3)=1250.
         amab(3)=3.3
         wab(4)=1430.
         amab(4)=3.
         wab(5)=1670.
         amab(5)=2.7
         wab(6)=2000.
         amab(6)=2.8
         wab(7)=2220.
         amab(7)=2.9
         wab(8)=2500.
         amab(8)=2.30
         wab(9)=2850.
         amab(9)=1.97
         wab(10)=3330.
         amab(10)=1.69
         wab(11)=3650.
         amab(11)=1.58
         wab(12)=4000.
         amab(12)=1.45
         wab(13)=4400.
         amab(13)=1.32
         wab(14)=5000.
         amab(14)=1.13
         wab(15)=5530.
         amab(15)=1.00
         wab(16)=6700.
         amab(16)=0.74
         wab(17)=9000.
         amab(17)=0.46
         wab(18)=10000.
         amab(18)=0.38
         wab(19)=20000.
         amab(19)=0.11
         wab(20)=100000.
         amab(20)=0.0
         DO i=1,19
            aab(i)=amab(i)*av
         END DO
         DO ik=1,iwavl
            IF (wl(ik).le.wab(1)) then 
               ccc=(aab(1)-aab(2))/(wab(2)-wab(1))
               alambda=aab(1)+(wab(1)-wl(ik))*ccc
            ELSE IF(wl(ik).ge.wab(19))then 
               ccc=(aab(19)-aab(20))/(wab(20)-wab(19))
               alambda=aab(19)+(wab(19)-wl(ik))*ccc
            ELSE IF (wl(ik).gt.wab(1).and.wl(ik).lt.wab(19)) then
               DO i=1,19
!Cid-fix!         IF (wl(ik).ge.wab(i).and.wl(ik).lt.wab(i+1)) then
                  IF (wl(ik).ge.wab(i).and.wl(ik).le.wab(i+1)) then
                     ccc=(aab(i)-aab(i+1))/(wab(i+1)-wab(i))
                     alambda=aab(i)+(wab(i)-wl(ik))*ccc
                  END IF
               END DO
            END IF
            IF (alambda.lt.0.) alambda=0.
            IF (ired.eq.1) yob(ik)=fl(ik)*(10**(-0.4*alambda))
            IF (ired.eq.2) yob(ik)=fl(ik)*(10**(0.4*alambda))
         END DO
         RETURN

c ..... 2. Seaton (1979) - Fitzpatrick fit : 
c                       for lambda < 1200 A slope from 1100-1200
c                       for lambda > 3650 A points from Allen
c                       ff = E(lambda-V)/E(B-V) = k(lambda)-R
c                       A(lambda)=A_V*(ff/R +1) ; R = 3.1

      ELSE IF (law_rd.eq.2) THEN
         rv=3.1
         al0=4.595
         ga=1.051
         c1=-0.38
         c2=0.74
         c3=3.96
         c4=0.26
         call red(1100.,al0,ga,c1,c2,c3,c4,ff)
         ff11=ff
         call red(1200.,al0,ga,c1,c2,c3,c4,ff)
         ff12=ff
         slope=(ff12-ff11)/100.
         wab(11)=3650.
         amab(11)=1.58
         wab(12)=4000.
         amab(12)=1.45
         wab(13)=4400.
         amab(13)=1.32
         wab(14)=5000.
         amab(14)=1.13
         wab(15)=5530.
         amab(15)=1.00
         wab(16)=6700.
         amab(16)=0.74
         wab(17)=9000.
         amab(17)=0.46
         wab(18)=10000.
         amab(18)=0.38
         wab(19)=20000.
         amab(19)=0.11
         wab(20)=100000.
         amab(20)=0.0
         DO i=11,19
            aab(i)=amab(i)*av
         END DO
         DO ik=1,iwavl
            IF (wl(ik).ge.1200..and.wl(ik).le.3650.) then
               call red(wl(ik),al0,ga,c1,c2,c3,c4,ff)
               alambda=av*(ff/rv+1.) 
            ELSE IF (wl(ik).lt.1200.) then
               ff=ff11+(wl(ik)-1100.)*slope
               alambda=av*(ff/rv+1.)
            ELSE IF(wl(ik).ge.wab(19))then !extrapol at high wavelengths
               ccc=(aab(19)-aab(20))/(wab(20)-wab(19))
               alambda=aab(19)+(wab(19)-wl(ik))*ccc
            ELSE IF (wl(ik).gt.3650..and.wl(ik).lt.wab(19)) then
               DO i=11,19
!Cid-fix!         IF (wl(ik).gt.wab(i).and.wl(ik).lt.wab(i+1)) then
                  IF (wl(ik).ge.wab(i).and.wl(ik).le.wab(i+1)) then
                     ccc=(aab(i)-aab(i+1))/(wab(i+1)-wab(i))
                     alambda=aab(i)+(wab(i)-wl(ik))*ccc
                  END IF
               END DO
            END IF
            IF (alambda.lt.0.) alambda=0.
            IF (ired.eq.1) yob(ik)=fl(ik)*(10**(-0.4*alambda))
            IF (ired.eq.2) yob(ik)=fl(ik)*(10**(0.4*alambda))
         END DO
         RETURN
      
c ..... 3. Fitzpatrick (1986) : 
c                       for lambda < 1200 A slope derived from 1100-1200
c                       for lambda > 3330 A points from Allen
c                       ff = E(lambda-V)/E(B-V) = k(lambda)-R
c                       A(lambda) = A_V*(ff/R +1) ; R = 3.1

      ELSE IF (law_rd.eq.3) THEN
         rv=3.1
         al0=4.608
         ga=0.994
         c1=-0.69
         c2=0.89
         c3=2.55
         c4=0.50
         call red(1100.,al0,ga,c1,c2,c3,c4,ff)
         ff11=ff
         call red(1200.,al0,ga,c1,c2,c3,c4,ff)
         ff12=ff
         slope=(ff12-ff11)/100.
         wab(10)=3330.
         amab(10)=1.682
         wab(11)=3650.
         amab(11)=1.58
         wab(12)=4000.
         amab(12)=1.45
         wab(13)=4400.
         amab(13)=1.32
         wab(14)=5000.
         amab(14)=1.13
         wab(15)=5530.
         amab(15)=1.00
         wab(16)=6700.
         amab(16)=0.74
         wab(17)=9000.
         amab(17)=0.46
         wab(18)=10000.
         amab(18)=0.38
         wab(19)=20000.
         amab(19)=0.11
         wab(20)=100000.
         amab(20)=0.0
         DO i=10,19
            aab(i)=amab(i)*av
         END DO
         DO ik=1,iwavl
            IF (wl(ik).ge.1200..and.wl(ik).le.3330.) then
               call red(wl(ik),al0,ga,c1,c2,c3,c4,ff)
               alambda=av*(ff/rv+1.)
            ELSE IF (wl(ik).lt.1200.) then
               ff=ff11+(wl(ik)-1100.)*slope
               alambda=av*(ff/rv+1.)
            ELSE IF(wl(ik).ge.wab(19))then !extrapol at high wavelengths
               ccc=(aab(19)-aab(20))/(wab(20)-wab(19))
               alambda=aab(19)+(wab(19)-wl(ik))*ccc
            ELSE IF (wl(ik).gt.3330..and.wl(ik).lt.wab(19)) then
               DO i=10,19
!Cid-fix!         IF (wl(ik).gt.wab(i).and.wl(ik).lt.wab(i+1)) then
                  IF (wl(ik).ge.wab(i).and.wl(ik).le.wab(i+1)) then
                     ccc=(aab(i)-aab(i+1))/(wab(i+1)-wab(i))
                     alambda=aab(i)+(wab(i)-wl(ik))*ccc
                  END IF
               END DO
            END IF
            IF (alambda.lt.0.) alambda=0.
            IF (ired.eq.1) yob(ik)=fl(ik)*(10**(-0.4*alambda))
            IF (ired.eq.2) yob(ik)=fl(ik)*(10**(0.4*alambda))
         END DO
         RETURN

c ..... 4. Prevot (1984) and Bouchet (1985) :
c                        wab = lambda in A,  
c                        amab = E(lambda-V)/E(B-V) = k(lambda)-R
c                        A(lambda) = A_V*(amab/R +1)
c                        N.B. R = 2.72

         ELSE IF (law_rd.eq.4) THEN
            rv=2.72
            wab(1)=1275.
            amab(1)=13.54
            wab(2)=1330.
            amab(2)=12.52
            wab(3)=1385.
            amab(3)=11.51
            wab(4)=1435.
            amab(4)=10.80
            wab(5)=1490.
            amab(5)=9.84
            wab(6)=1545.
            amab(6)=9.28
            wab(7)=1595.
            amab(7)=9.06
            wab(8)=1647.
            amab(8)=8.49
            wab(9)=1700.
            amab(9)=8.01
            wab(10)=1755.
            amab(10)=7.71
            wab(11)=1810.
            amab(11)=7.17
            wab(12)=1860.
            amab(12)=6.90
            wab(13)=1910.
            amab(13)=6.76
            wab(14)=2000.
            amab(14)=6.38
            wab(15)=2115.
            amab(15)=5.85
            wab(16)=2220.
            amab(16)=5.30
            wab(17)=2335.
            amab(17)=4.53
            wab(18)=2445.
            amab(18)=4.24
            wab(19)=2550.
            amab(19)=3.91
            wab(20)=2665.
            amab(20)=3.49
            wab(21)=2778.
            amab(21)=3.15
            wab(22)=2890.
            amab(22)=3.00
            wab(23)=2995.
            amab(23)=2.65
            wab(24)=3105.
            amab(24)=2.29
            wab(25)=3704.
            amab(25)=1.81
            wab(26)=4255.
            amab(26)=1.00
            wab(27)=5291.
            amab(27)=0.00
            wab(28)=12500.
            amab(28)=-2.02
            wab(29)=16500.
            amab(29)=-2.36
            wab(30)=22000.
            amab(30)=-2.47

            DO i=1,30
               aab(i)=av*(amab(i)/rv+1.)
            END DO

            DO ik=1,iwavl
               IF (wl(ik).le.wab(1)) THEN !extrapol at short wavelengths
                  ccc=(aab(1)-aab(2))/(wab(2)-wab(1))
                  alambda=aab(1)+(wab(1)-wl(ik))*ccc
               ELSE IF(wl(ik).ge.wab(30))then !extrapol at high wavelengths
                  ccc=(aab(29)-aab(30))/(wab(30)-wab(29))
                  alambda=aab(29)+(wab(29)-wl(ik))*ccc  
               ELSE IF (wl(ik).gt.wab(1).and.wl(ik).lt.wab(30)) THEN
                  DO i=1,30
!Cid-fix!            IF (wl(ik).gt.wab(i).and.wl(ik).lt.wab(i+1)) THEN
                     IF (wl(ik).ge.wab(i).and.wl(ik).le.wab(i+1)) THEN
                        ccc=(aab(i)-aab(i+1))/(wab(i+1)-wab(i))
                        alambda=aab(i)+(wab(i)-wl(ik))*ccc
                     END IF
                  END DO
               END IF   
               IF (alambda.lt.0.) alambda=0.
               IF (ired.eq.1) yob(ik)=fl(ik)*(10**(-0.4*alambda))
               IF (ired.eq.2) yob(ik)=fl(ik)*(10**(0.4*alambda))
            END DO
            RETURN


c ..... Calzetti (1999) : 
c                   for lambda < 1200 A slope derived from 1100-1200 A
c                   for lambda > 22000A slope derived from 21900-22000 A
c                   ff = k(lambda); A(lambda) = A_V*k(lambda)/R ; R = 4.05

      ELSE IF (law_rd.eq.5) THEN
         rv=4.05
         p11=1/0.11
         ff11=2.659*(-2.156+1.509*p11-0.198*p11**2
     .        +0.011*p11**3)+rv
         p12=1/0.12
         ff12=2.659*(-2.156+1.509*p12-0.198*p12**2
     .        +0.011*p12**3)+rv
         slope1=(ff12-ff11)/100.
         ff99=2.659*(-1.857+1.040/2.19)+rv
         ff100=2.659*(-1.857+1.040/2.2)+rv
         slope2=(ff100-ff99)/100.
         DO ik=1,iwavl
            ala=wl(ik)*1.E-4    !wavelength in microns
            p=1./ala
            IF (ala.ge.0.63.and.ala.le.2.2) then
               ff=2.659*(-1.857+1.040*p)+rv
            ELSE IF (ala.ge.0.12.and.ala.lt.0.63) then
               ff=2.659*(-2.156+1.509*p-0.198*p**2+0.011*p**3)+rv
            ELSE IF (ala.lt.0.12) THEN
               ff=ff11+(wl(ik)-1100.)*slope1
            ELSE IF (ala.gt.2.2) THEN
               ff=ff99+(wl(ik)-21900.)*slope2
            END IF
            alambda=av*ff/rv
            IF (alambda.lt.0.) alambda=0.
            IF (ired.eq.1) yob(ik)=fl(ik)*(10**(-0.4*alambda))
            IF (ired.eq.2) yob(ik)=fl(ik)*(10**(0.4*alambda))
         END DO
         RETURN
      END IF

      WRITE(*,*) ' No reddening law is given.....'
      RETURN
      END

c ----------------------------------------------------------------------
      SUBROUTINE red(alan,al0,ga,c1,c2,c3,c4,ff)

c.......According to Fitzpatrick, 86

      ala=alan*1.E-4           !wavelength in microns
      p=1.0/ala
      IF (p.lt.5.9) c4=0.0
      t1=c3/(((p-(al0**2/p))**2)+ga*ga)
      t2=c4*(0.539*((p-5.9)**2)+0.0564*((p-5.9)**3))
      ff=c1+c2*p+t1+t2
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c ###########################################################################
c     Gordon SMC & LMC extinction curves
c     After failing to implement Gordon's LMC & SMC laws in analytical
c     form (their equation 5, with parameters from Table 3 and R_V from
c     Sec 3.4 & Table 2), I have just copied their Table 4 and
c     interpolate linearly (in x).
c
c     q = A_lambda / A_V from Gordon et al 2003, ApJ, 594, 279, Table 4.
c     l = lambda, in Angstrons
c     x = 1 / lambda in 1/microns
c     red_law_option = 'GD1' => SMC Bar
c                    = 'GD2' => LMC2 SuperShell
c                    = 'GD3' => LMC Average
c
c     Cid@Lagoa - 30/December/2006

      FUNCTION GordonTab4_RedLaw(l,red_law_option)

      implicit real (a-h,j-m,o-z)
      implicit integer (i,n)
      parameter (N_Gordon = 29)
      real x_Gordon(N_Gordon) , q_Gordon(N_Gordon)
      real q_SMCbar(N_Gordon) , q_LMC2(N_Gordon) , q_LMCave(N_Gordon)
      character*3 red_law_option

c     Table 4 of Gordon et al (Sample Average Curves)
c     OBS: The Table has 2 empty entries for LMC2 & LMCave, at x = 1.235
c     & 1.538. These entries were interpolated linearly in the
c     data-blocks below.
      data x_Gordon/ 0.455 , 0.606 , 0.800 , 1.235 , 1.538 , 1.818 , 
     $       2.273 , 2.703 , 3.375 , 3.625 , 3.875 , 4.125 , 4.375 ,
     $       4.625 , 4.875 , 5.125 , 5.375 , 5.625 , 5.875 , 6.125 ,
     $       6.375 , 6.625 , 6.875 , 7.125 , 7.375 , 7.625 , 7.875 , 
     $       8.125 , 8.375 /
      data q_SMCbar/ 0.016 , 0.169 , 0.131 , 0.567 , 0.801 , 1.000 , 
     $       1.374 , 1.672 , 2.000 , 2.220 , 2.428 , 2.661 , 2.947 , 
     $       3.161 , 3.293 , 3.489 , 3.637 , 3.866 , 4.013 , 4.243 , 
     $       4.472 , 4.776 , 5.000 , 5.272 , 5.575 , 5.795 , 6.074 , 
     $       6.297 , 6.436 /
      data q_LMC2/   0.101 , 0.097 , 0.299 , 0.5004, 0.7600, 1.000 , 
     $       1.349 , 1.665 , 1.899 , 2.067 , 2.249 , 2.447 , 2.777 , 
     $       2.922 , 2.921 , 2.812 , 2.805 , 2.863 , 2.932 , 3.060 , 
     $       3.110 , 3.299 , 3.408 , 3.515 , 3.670 , 3.862 , 3.937 , 
     $       4.055 , 3.969 /
      data q_LMCave/ 0.030 , 0.186 , 0.257 , 0.4705, 0.7457, 1.000 ,
     $       1.293 , 1.518 , 1.786 , 1.969 , 2.149 , 2.391 , 2.771 ,
     $       2.967 , 2.846 , 2.646 , 2.565 , 2.566 , 2.598 , 2.607 ,
     $       2.668 , 2.787 , 2.874 , 2.983 , 3.118 , 3.231 , 3.374 ,
     $       3.366 , 3.467 /

c     Define x & store appropriate q onto q_Gordon
      x = 10000. / l
      do i=1,N_Gordon
         if (red_law_option.EQ.'GD1') then
            q_Gordon(i) = q_SMCbar(i)
         else if (red_law_option.EQ.'GD2') then
            q_Gordon(i) = q_LMC2(i)
         else if (red_law_option.EQ.'GD3') then
            q_Gordon(i) = q_LMCave(i)
         else
            write (*,*) 'OOPS! Bad red_law_option = ',red_law_option
            stop
         end if
      end do

c     Linear interpolation of q = a x + b
      do i=1,N_Gordon-1
         if ((x.GE.x_Gordon(i)).AND.(x.LE.x_Gordon(i+1))) then
            q = q_Gordon(i) + (q_Gordon(i+1) - q_Gordon(i)) * 
     $           (x - x_Gordon(i)) / (x_Gordon(i+1) - x_Gordon(i))
            GordonTab4_RedLaw = q
            return
         end if
      end do

c     If lambda is outside Table 4 range, return zero,
      GordonTab4_RedLaw = 0.

      return
      end
c ###########################################################################
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+





c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c +        SYNTHETIC SPECTRUM, CHI2 & RAPID-CHI2 ("RC") ROUTINES            +
c +                                                                         +
c + calc_Fsyn                                                               +
c + Kcalc_Fsyn                                                              +
c + Slow_Fchi2                                                              +
c + FK_chi2                                                                 +
c + Calc_RC_terms                                                           +
c + VeryFastF_chi2                                                          +
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+



c ###########################################################################
c     Compute synthetic spectrum f_syn(l_syn) for given population &
c     extinction vectors, x(j) & AV_tot(j). Stellar kinematics (if any)
c     is inside the f_base matrix (as set elsewhere)!
c     OBS: Reddening term, base matrix & Nl_obs passed by common!  
c     Cid@Lagoa - Jan/11/2005

      SUBROUTINE calc_Fsyn(x,AV_tot,f_syn)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
      real l_obs(Nl_max) , f_obs(Nl_max) , weight_obs(Nl_max)
      real x(Nb_max) , AV_tot(Nb_max)
      real f_base(Nb_max,Nl_max) , f_syn(Nl_max)

      common/Red_factor/ red_term(Nl_max)
      common/Base/ N_base , f_base
      common/Obs/ Nl_obs , l_obs , f_obs , weight_obs

      do i=1,Nl_obs
         f_syn(i) = 0.
         do j=1,N_base
            f_syn(i) = f_syn(i) + f_base(j,i) * abs(x(j)) * 
     &           10.**(AV_tot(j) * red_term(i))
         end do
      end do

      end
c ###########################################################################



c ###########################################################################
c     Compute synthetic spectrum f_syn(l_syn) for given population &
c     extinction vectors, x(j) & AV_tot(j), AND kinematical parameters
c     v0 & vd. A synthetic spectrum f(x,AV_tot) is 1st computed & then
c     convolved with a Gaussian of velocity center = v0 & velocity
c     dispersion = vd (both in km/s).  
c     Cid@Lagoa - Jan/11/2005

      SUBROUTINE Kcalc_Fsyn(Nl_syn,x,AV_tot,v0,vd,f_syn)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
      real l_obs(Nl_max) , f_obs(Nl_max) , weight_obs(Nl_max)
      real l_aux(Nl_max) , f_aux(Nl_max) , f_syn(Nl_max)
      real x(Nb_max) , AV_tot(Nb_max)

      common/Obs/ Nl_obs , l_obs , f_obs , weight_obs
      common/N_Gauss/ N_int_Gauss

!Cid14Jan07-FXK! Added FXK = Fixed kinematics option!!
      common/FXK_fits/ IsThisAnFXKFit

!Cid14Jan07-FXK! If this is a FXK-fit then we should NOT apply the
!     kinematical filter, as it has already been applied to the base!!
!     Just call calc_Fsyn & return.
      if (IsThisAnFXKFit.EQ.1) then
         call calc_Fsyn(x,AV_tot,f_syn)
         return
      end if

c     Synthetic spectrum f_aux. Also define l_aux (used below). It is
c     presumed that the base spectra in f_base are in their original
c     form (ie, not convolved with a kinematical filter)!!
      call calc_Fsyn(x,AV_tot,f_aux)

c     Apply kinematical filter! (aux variables used to avoid repeating
c     variable names in the call to GaussSmooth). Returns the convolved
c     spectrum f_syn.
      Nl_aux = Nl_syn
      do i=1,Nl_syn
         l_aux(i) = l_obs(i)
      end do
      call GaussSmooth(Nl_aux,l_aux,f_aux,v0,vd,
     &     N_int_Gauss,Nl_syn,l_obs,f_syn)

      end
c ###########################################################################



c ###########################################################################
c     Compute chi2(x,AV_tot). Kinematical parameters are implicit in f_base!
c     I now call this function "Slow" because it is slow compared to the
c     Rapid-Chi2 tricks developed. Most of the code used rapid-chi2
c     functions.  This one is only used to check results.
c     Cid@Lagoa - Jan/11/2005

      FUNCTION SlowF_chi2(x,AV_tot)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
      real l_obs(Nl_max) , f_obs(Nl_max) , weight_obs(Nl_max)
      real x(Nb_max) , AV_tot(Nb_max)
      real f_base(Nb_max,Nl_max) , f_syn(Nl_max)

      common/Red_factor/ red_term(Nl_max)
      common/Base/ N_base , f_base
      common/Obs/ Nl_obs , l_obs , f_obs , weight_obs
!OPTfn! Optically optimal fn - computed here and passed on 
      common/Optimize_fn_OPT/ IsOptimize_fn_OPT , Optimal_fn_OPT

      call calc_Fsyn(x,AV_tot,f_syn)

      chi2 = 0.
      do i=1,Nl_obs
         chi2 = chi2 + (f_syn(i) - f_obs(i)) * weight_obs(i) *
     &        (f_syn(i) - f_obs(i)) * weight_obs(i)
      end do

*ETC* Adding all ETC-chi2's to optical-chi2!
      chi2_ETC = ETC_CalcETCchi2(x,AV_tot,chi2,'SFc')
      SlowF_chi2 = chi2 + chi2_ETC

!OPTfn! Added computation of Optimal_fn_OPT just in case... At the very least, the code 
!     will definitely need it in the unlikely case i_UseSlowChi2Always.Q.1.
!     ElCid@Sanchica - 28/Jan/2012
      chi2_2 = 0.
      chi2_3 = 0.
      do i=1,Nl_obs
         chi2_2 = chi2_2 - 2. *  f_syn(i) * weight_obs(i) * f_obs(i) *
     $        weight_obs(i)
         chi2_3 = chi2_3 + f_syn(i) * weight_obs(i) * f_syn(i) *
     $        weight_obs(i)
      end do
      Optimal_fn_OPT = -(chi2_2/2.) / chi2_3

      return
      end
c ###########################################################################



c ###########################################################################
c     Compute chi2(x,AV_tot,v0,vd), INCLUDING kinematical parameters!
c     This a "Slow" function, ie, no RC-tricks are used. In fact, it is
c     Very slow, since it convolves the spectrum with the G(v0,vd)
c     filter...
c     Cid@Lagoa - Jan/11/2005

      FUNCTION FK_chi2(x,AV_tot,v0,vd)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
      real l_obs(Nl_max) , f_obs(Nl_max) , weight_obs(Nl_max)
      real x(Nb_max) , AV_tot(Nb_max) , f_syn(Nl_max)

      common/Obs/ Nl_obs , l_obs , f_obs , weight_obs

      call Kcalc_Fsyn(Nl_obs,x,AV_tot,v0,vd,f_syn)

      chi2 = 0.
      do i=1,Nl_obs
         chi2 = chi2 + (f_syn(i) - f_obs(i)) * weight_obs(i) *
     &        (f_syn(i) - f_obs(i)) * weight_obs(i)
      end do

*ETC* Adding all ETC-chi2's to optical-chi2!
      chi2_ETC = ETC_CalcETCchi2(x,AV_tot,chi2,'FKc')
      FK_chi2 = chi2 + chi2_ETC


*FIR* adicionei este return (nada a ver com FIR, mas sempre é bom)
      return                    
      end
c ###########################################################################



c ###########################################################################
c     Compute RC-terms. These are the terms in chi2 which do NOT depend
c     on either x, AV nor YAV, and thus do not be recomputed every time
c     x, AV or YAV changes, allowing enormous speed gains! This trick
c     requires expanding the 10**(AV_tot(j) * red_term(lambda)) factor
c     in a series and truncating it to its first NRC_AV terms, which is
c     Ok as long as |AV| does not go much beyond ~ 1... The first term
c     (RC1) is always the same, whereas the 2nd and 3rd terms (RC2 &
c     RC3) need to be re-computed everytime something was done to v0
c     and/or vd.... The VeryFastF_chi2 routine below use RC1, RC2 & RC3,
c     as well as NRC_AV (the number of terms in the truncated
c     AV-series). See notes for more details!
c
c     This routine can be made faster by introducing some "memory" so
c     that only new terms in the AV-series need to be computed. However,
c     this should NOT be useful, since I expect it to be called mostly
c     when a full-recalculation is needed (eg, when f_base
c     changes). Re-calls to the routine due to precision bugs should be
c     only occasional. Hence I did not attempt to optimize it further.
c
c     Cid@Lagoa  - 9/May/2004

c     OBS: Some of the initial steps, like the definition of the
c     factorial and beta arrays, could be done only once in the code,
c     but they take so little time that I kept them here for
c     conciseness...

c     Now adapted to AV = AV_tot(j), ie, population dependent
c     extinctions!
c     Cid@Lagoa - 12/Jan/2005


      SUBROUTINE Calc_RC_terms


c ***************************************************************************
c *                                Definitions                              *
c ***************************************************************************
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
*     MEX-OBS: Changed NRCAV_max
      parameter (Nl_max = 14000,Nb_max = 300,NRCAV_max = 35)

      real l_obs(Nl_max) , f_obs(Nl_max) , weight_obs(Nl_max)
      real f_base(Nb_max,Nl_max) 

      double precision RC1 , RC2(NRCAV_max,Nb_max) , 
     &     RC3(NRCAV_max,Nb_max,Nb_max) , eta_AV
      double precision factorial(NRCAV_max) , rho(Nl_max) , beta(Nl_max)
      double precision soma , ln10

      common/Red_factor/ red_term(Nl_max)
      common/Base/ N_base , f_base
      common/Obs/ Nl_obs , l_obs , f_obs , weight_obs
      common/RapidChi2/ RC1 , RC2 , RC3 , eta_AV
      common/RC_AV_Stuff/ NRC_AV , i_RC_CRASH_FLAG , NRC_AV_Default
c ***************************************************************************


c ***************************************************************************
c *                              Initialization                             *
c ***************************************************************************
c     Test if NRC_AV exceed its limit! If so then set a CRASH_FLAG!
c     OBS: Once i_RC_CRASH_FLAG = 1, it'll stay so till the end of the
c          code! If a "crash" occurs we set NRC_AV = NRCAV_max and go on...
      if (NRC_AV.GT.NRCAV_max) then 
         write (*,9) 
 9       format (55('@'))
         write (*,*) '[calc_RC_terms] SERIOUS BUG! NRC_AV > NRCAV_max ',
     &        NRC_AV , NRCAV_max
         write (*,9) 
         NRC_AV =  NRCAV_max
         i_RC_CRASH_FLAG = 1
      end if


c     Define factorial of 0...NRC_AV-1. 
c     Note array off-set: factorial(i) = (i-1)!
      factorial(1) = 1.
      factorial(2) = 1.
      do i=3,NRC_AV
         factorial(i) = factorial(i-1) * (i - 1)
      end do


c     Define rho(lambda) = ln(10) * red_term(lambda) and compute its
c     mean.  No weighting is applied. This mean serves as the factor
c     "eta_AV".
      ln10 = log(10.)
      soma = 0.
      do i=1,Nl_obs
         rho(i)  = ln10 * red_term(i)
         soma = soma + rho(i)
      end do
      eta_AV = soma / Nl_obs


c     Define beta(lambda) = rho(lambda) - eta_AV
      do i=1,Nl_obs
         beta(i) = rho(i) - eta_AV
      end do
c ***************************************************************************


c ***************************************************************************
c *                    Compute RC terms RC1, RC2 & RC3                      *
c ***************************************************************************
c     RC1 = Sum-over-l  (O_l * w_l)**2. 
      RC1 = 0.
      do i=1,Nl_obs
         RC1 = RC1 + f_obs(i) * weight_obs(i) * f_obs(i) * weight_obs(i)
      end do


c     RC2 = -2 * Sum-over-l w_l**2 * O_l * gamma_l,j * beta_l**n
c     RC2 = RC2(iRC_AV,j)
      do j=1,N_base
         do iRC_AV=1,NRC_AV
            soma = 0.
            do i=1,Nl_obs
               soma = soma + f_obs(i) * weight_obs(i) * f_base(j,i) * 
     &              weight_obs(i) * beta(i)**(iRC_AV - 1)
            end do
            RC2(iRC_AV,j) = -(2 * soma) / factorial(iRC_AV)

! Screen check
!           write (*,2) '[calc_RC_terms] RC2: ' , j , iRC_AV , 
!    &           RC2(iRC_AV,j)
!2          format (a,2(2x,i4),3x,1p,3(4x,e12.4))
         end do
      end do


c     RC3 = Sum-over-l w_l**2 * gamma_l,j1 * gamma_l,j2 * beta_l**n
c     RC3 = RC3(iRC_AV,j1,j2)
c     OBS: Explore fact that RC3(iRC_AV,j1,j2) = RC3(iRC_AV,j2,j1)...
      do j1=1,N_base
         do j2=1,j1
            do iRC_AV=1,NRC_AV
               soma = 0.
               do i=1,Nl_obs
                  soma = soma + f_base(j1,i) * weight_obs(i) * 
     &                 f_base(j2,i) * weight_obs(i) * 
     &                 beta(i)**(iRC_AV - 1)
               end do
               RC3(iRC_AV,j1,j2) = soma / factorial(iRC_AV)
               RC3(iRC_AV,j2,j1) = RC3(iRC_AV,j1,j2)

! Screen check
!           write (*,3) '[calc_RC_terms] RC3: ' , j1 , j2 , iRC_AV ,
!    &           RC3(iRC_AV,j1,j2)
!3           format (a,3(2x,i4),3x,1p,e12.4)
            end do
         end do
      end do
c ***************************************************************************


      end
c ###########################################################################



c ###########################################################################
c     This version of chi2(x,AV_tot) uses the information in the flag
c     i_IsAVstillTheSame, which is = 1 ("YES") if AV is the same as in
c     the last call to this function. This allows us to skip the sum
c     over the NRC_AV terms in the truncated-series approximation to
c     chi2, thus saving quite a bit of time, hence allowing many more
c     Metropolis steps! See notes for more details!
c     Cid@Lagoa - 09/May/2004
c
c     Now adapted to AV = AV_tot(j), ie, population dependent
c     extinctions! Double precision variables were defined too!
c     Cid@Lagoa - 12/Jan/2005

!OPTfn! This function now returns (via common) the (single precision) optimal fn = -(chi2_2/2.) / chi2_3.
!       This considers ONLY the OPTICAL spectral fit!! (hence the _OPT subscript!)
!     ElCid@Sanchica - 28/Jan/2012

      FUNCTION VeryFastF_chi2(x,AV_tot,i_IsAVstillTheSame)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
*     MEX-OBS: Changed NRCAV_max
      parameter (Nl_max = 14000,Nb_max = 300,NRCAV_max = 35)

      real x(Nb_max) , abs_x(Nb_max) , AV_tot(Nb_max)
      real f_base(Nb_max,Nl_max)
      double precision RC1 , RC2(NRCAV_max,Nb_max) , 
     &     RC3(NRCAV_max,Nb_max,Nb_max) , eta_AV

      double precision SMAT1(Nb_max) , SMAT2(Nb_max,Nb_max)
      double precision chi2 , chi2_1 , chi2_2 , chi2_3 , soma , soma_j2

      SAVE SMAT1 , SMAT2

      common/Base/ N_base , f_base
      common/RapidChi2/ RC1 , RC2 , RC3 , eta_AV
      common/RC_AV_Stuff/ NRC_AV , i_RC_CRASH_FLAG , NRC_AV_Default

!!MEX!! IMPORTANT CHANGE: Introduced flag to use SLOW chi2 always!
      common/ToRCOrNotToRC_ThatIsTheQuestion/ i_UseSlowChi2Always

!OPTfn! Optically optimal fn - computed here and passed on
      common/Optimize_fn_OPT/ IsOptimize_fn_OPT , Optimal_fn_OPT

* !!MEX!! Use Slow Chi2 if told to do so!
      if (i_UseSlowChi2Always.EQ.1) then
         chi2           = SlowF_chi2(x,AV_tot)
         VeryFastF_chi2 = chi2
         return
      end if


c     Define auxiliar i_IsAVstillTheSame flag, which may be changed in
c     case of trouble (see bottom)
      i_aux_IsAVstillTheSame = i_IsAVstillTheSame


c     1st term = Sum_l (w_l O_l)**2
 10   chi2_1 = RC1 

c     If AV is NOT the same as in the last call then compute chi2 as in
c     the fast-chi2 routine, but saving the SMAT1 & SMAT2 arrays for
c     later use
      if (i_aux_IsAVstillTheSame.NE.1) then 


c     2nd term = -2 Sum_l w_l**2 O_l * M_l. 
c     Also, define abs_x() for convenience (= speed), as it's used below.
         chi2_2 = 0.
         do j=1,N_base
            abs_x(j) = abs(x(j))
            soma = 0.
            do iRC_AV=1,NRC_AV
               soma = soma + AV_tot(j)**(iRC_AV - 1) * RC2(iRC_AV,j)
            end do
            SMAT1(j) = soma * exp(eta_AV * AV_tot(j))
            chi2_2   = chi2_2 + abs_x(j) * SMAT1(j)
         end do


c     3rd term = Sum_l (w_l M_l)**2.
c     To speed up the triple loop, first sum over the off-diagnonal
c     elements j2 < j1, then use the (j1,j2) = (j2,j1) simmetry in RC3
c     (hence the chi2_3 = 2 * chi2_3 line below), and finally sum over
c     the j1 = j2 diagonal elements.
         chi2_3 = 0.
         do j1=1,N_base
            soma_j2 = 0.
            do j2=1,j1-1
               AV_aux = AV_tot(j1) + AV_tot(j2)
               soma   = 0.
               do iRC_AV=1,NRC_AV
                  soma = soma + AV_aux**(iRC_AV - 1) * RC3(iRC_AV,j1,j2)
               end do
               SMAT2(j1,j2) = soma * exp(eta_AV * AV_aux)
               SMAT2(j2,j1) = SMAT2(j1,j2)
               soma_j2      = soma_j2 + abs_x(j2) * SMAT2(j1,j2)
            end do
            chi2_3 = chi2_3 + soma_j2 * abs_x(j1)
         end do
         chi2_3 = 2 * chi2_3

         do j1=1,N_base
            AV_aux = AV_tot(j1) + AV_tot(j1)
            soma = 0.
            do iRC_AV=1,NRC_AV
               soma = soma + AV_aux**(iRC_AV - 1) * RC3(iRC_AV,j1,j1)
            end do
            SMAT2(j1,j1) = soma * exp(eta_AV * AV_aux)
            chi2_3       = chi2_3 + abs_x(j1) * abs_x(j1) * SMAT2(j1,j1)
         end do


      else


c     If AV IS the same as in the last call, then use previously defined
c     SMAT1 & SMAT2 arrays to speed up computation (skiping the
c     AV-series-loop)

c     2nd term ... Using SMAT1!!
         chi2_2 = 0.
         do j=1,N_base
            abs_x(j) = abs(x(j))
            chi2_2 = chi2_2 + abs_x(j) * SMAT1(j)
         end do

c     3rd term  ... Using SMAT2!!
         chi2_3 = 0.
         do j1=1,N_base
            soma_j2 = 0.
            do j2=1,j1-1
               soma_j2 = soma_j2 + abs_x(j2) * SMAT2(j1,j2)
            end do
            chi2_3 = chi2_3 + soma_j2 * abs_x(j1)
         end do
         chi2_3 = 2 * chi2_3

         do j1=1,N_base
            chi2_3 = chi2_3 + abs_x(j1) * abs_x(j1) * SMAT2(j1,j1)
         end do


      end if


c     Add the 3 terms to obtain the total chi2, test it & go back:)
      chi2 = chi2_1 + chi2_2 + chi2_3


!OPTfn! Computes the OPTically optimal fn
      Optimal_fn_OPT = -(chi2_2/2.) / chi2_3


*ETC* Adding all ETC-chi2's to optical-chi2!
*     OBS: aux_chi2 used to change to single precision in the call to
*     ETC_CalcETCchi2
      aux_chi2       = chi2
      chi2_ETC       = ETC_CalcETCchi2(x,AV_tot,aux_chi2,'VFc')
      VeryFastF_chi2 = chi2 + chi2_ETC

!!MEX!! IMPORTANT CHANGE: In case all else fails, revert to Slow chi2!
      if ((chi2.LT.0.).AND.((NRC_AV + 1).GT.NRCAV_max)) then 
         write (*,'(a,i3,a,f12.3,a,i3,a)')
     $        '[VeryFastF_chi2] ATT!! Tried to resize NRC_AV to' ,
     $        NRC_AV +1 , ' because chi2 = ' , chi2 , 
     $        ' but NRCAV_max = ' , NRCAV_max , ' !!'
         write (*,'(a)')
     $        '[VeryFastF_chi2] ..... Will revert to SLOW chi2!!'
         chi2           = SlowF_chi2(x,AV_tot)
         VeryFastF_chi2 = chi2

*FIR* Warning msg below was made more informative ...
         write (*,'(3(a,1p,e12.5,0p))')
     $        '[VeryFastF_chi2] ..... Done: chi2 = ', chi2 ,
     $        ' of which ' , chi2_FIR , 
     $        ' comes from FIR. Dif = ' ,chi2 - chi2_FIR      
      end if


c     Test whether series approximation produced chi2 < 0! If so, try
c     again with one further term [NRC_AV <== NRC_AV + 1]
      if (chi2.LT.0.) then
         NRC_AV = NRC_AV + 1
         write (*,'(a,i3,a,f12.3,a)')
     &        '[VeryFastF_chi2] ATT!! Resized NRC_AV to' , NRC_AV , 
     &        '  because chi2 = ' , chi2 , '!'

         call Calc_RC_terms

!CRASH: get out in case RC-crashes! If no-crash, re-do chi2 calculation.
         if (i_RC_CRASH_FLAG.NE.0) return

         i_aux_IsAVstillTheSame = 0
         go to 10
      end if

      return
      end
c ###########################################################################
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+





c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c +          CHI2 FUNCTIONS USED IN THE KINEMATICAL FITS (v0 & vd)          +
c +                                                                         +
c + fun_vd_chi2                                                             +
c + fun_v0_chi2                                                             +
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+



c ###########################################################################
c     This function is used in the Golden-Section-Search for a best vd.
c     It computes chi2 as a function of vd using the kinematics-free
c     (ie, v0 = vd = 0) synthetic spectrum "GSf_syn_NoKin", computed
c     elsewhere for the current best model. v0_min is passed by the
c     /GoldenSearch_aux/ common, as it is needed for GaussSmooth.
c     Nl_aux & l_aux are = Nl_obs & l_obs, but defined elesewhere
c     for convenience (= speed); their only role is to avoid calling
c     GaussSmooth with equally named variables...

c     Cid@Lagoa - Jan/13/2005

*ETC* Added common/ETC_xAndAV4KinFits/ to pass x & AV arrays to
*     fun_vd_chi2 & fun_v0_chi2. Also needed to add Nb_max in the
*     parameter list.
*     Cid@Lagoa - Jan/13/2010 + 07/Jan/2011

      FUNCTION fun_vd_chi2(vd)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
      real l_obs(Nl_max) , f_obs(Nl_max) , weight_obs(Nl_max)
      real l_aux(Nl_max) , f_syn(Nl_max)

      common/N_Gauss/ N_int_Gauss
      common/Obs/ Nl_obs , l_obs , f_obs , weight_obs
      common/GoldenSearch_aux/ v0_min , vd_min ,
     &     Nl_aux , l_aux , GSf_syn_NoKin(Nl_max)
      common/ETC_xAndAV4KinFits/ ETC_x4kin(Nb_max) , ETC_AV4kin(Nb_max)

!FIX-DAE2-BUG! Imposing a maximum vd! Set chi2 to "infty" is vd > 5000 km/s
!	  ElCid@Sanchica - 21/Oct/2011
!ATT! This FIX was wronly implemented in PANcMEx_StarlightChains_v02.for!
!	  The line which now says fun_vd_chi2 = 3.0e38 line was fun_v0_chi2 = 3.0e38
!	  as a result of a careless cut & paste... 
!	  ElCid@Sanchica - 05/Nov/2011
      if (abs(vd).GT.5000.0) then
      	write (*,*) '[fun_vd_chi2] ATT! Trying vd=',vd,'!!!'
      	write (*,*) '[fun_vd_chi2]      NOT allowed!'
        fun_vd_chi2 = 3.0e38
        return
      end if

      call GaussSmooth(Nl_aux,l_aux,GSf_syn_NoKin,v0_min,vd,
     &     N_int_Gauss,Nl_obs,l_obs,f_syn)

      chi2 = 0.
      do i=1,Nl_obs
         chi2 = chi2 + (f_syn(i) - f_obs(i)) * weight_obs(i) *
     &        (f_syn(i) - f_obs(i)) * weight_obs(i)
      end do

*ETC* Adding all ETC-chi2's to optical-chi2!
      chi2_ETC = ETC_CalcETCchi2(ETC_x4kin,ETC_AV4kin,chi2,'fvd')
      fun_vd_chi2 = chi2 + chi2_ETC

*FIR* adicionei este return (nada a ver com FIR, mas sempre é bom)
      return
      end
c ###########################################################################

c ###########################################################################
c     This function is used in the Golden-Section-Search for a best v0.
c     It computes chi2 as a function of v0 using the kinematics-free
c     (ie, v0 = vd = 0) synthetic spectrum "GSf_syn_NoKin", computed
c     elsewhere for the current best model. vd_min is passed by the
c     /GoldenSearch_aux/ common, as it is needed for GaussSmooth.
c     Nl_aux & l_aux are = Nl_obs & l_obs, but defined elesewhere
c     for convenience (= speed); their only role is to avoid calling
c     GaussSmooth with equally named variables...

c     Cid@Lagoa - Jan/13/2005

      FUNCTION fun_v0_chi2(v0)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
      real l_obs(Nl_max) , f_obs(Nl_max) , weight_obs(Nl_max)
      real l_aux(Nl_max) , f_syn(Nl_max)

      common/N_Gauss/ N_int_Gauss
      common/Obs/ Nl_obs , l_obs , f_obs , weight_obs
      common/GoldenSearch_aux/ v0_min , vd_min ,
     &     Nl_aux , l_aux , GSf_syn_NoKin(Nl_max)
      common/ETC_xAndAV4KinFits/ ETC_x4kin(Nb_max) , ETC_AV4kin(Nb_max)

!FIX-DAE2-BUG! Imposing a maximum |v0|! Set chi2 to "infty" is v0 > 5000 km/s
!	  ElCid@Sanchica - 21/Oct/2011
      if (abs(v0).GT.5000.0) then
      	write (*,*) '[fun_v0_chi2] ATT! Trying v0=',v0,'!!!'
      	write (*,*) '[fun_v0_chi2]      NOT allowed!'
        fun_v0_chi2 = 3.0e38
        return
      end if

      call GaussSmooth(Nl_aux,l_aux,GSf_syn_NoKin,v0,vd_min,
     &     N_int_Gauss,Nl_obs,l_obs,f_syn)

      chi2 = 0.
      do i=1,Nl_obs
         chi2 = chi2 + (f_syn(i) - f_obs(i)) * weight_obs(i) *
     &        (f_syn(i) - f_obs(i)) * weight_obs(i)
      end do

*ETC* Adding all ETC-chi2's to optical-chi2!
      chi2_ETC = ETC_CalcETCchi2(ETC_x4kin,ETC_AV4kin,chi2,'fv0')
      fun_v0_chi2 = chi2 + chi2_ETC

*FIR* adicionei este return (nada a ver com FIR, mas sempre é bom)
      return
      end
c ###########################################################################
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+





c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c +                      MARKOV-CHAINS RELATED ROUTINES                     +
c +             (This is where the "soul" of the code resides!!)            +
c +                                                                         +
c + Calc_GR_R                                                               +
c + Generate_PopVector                                                      +
c + Fpar_Slowchi2                                                           +
c + Fpar_VFchi2                                                             +
c + EX0s_Fit_Scubi                                                          +
c + Fit_Scubidu                                                             +
c + ScubiduCadeVoce                                                         +
c + Unpack_par2xAVYAVfn                                                     +
c + SaveStep                                                                +
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+



c ###########################################################################
c     Compute Gelman & Rubin "R" convergence-test-ratio GR_R(j) for each
c     parameter j = 1...N_par, given N_chains chains. The input matrices
c     Cpar_ave & Cpar_var must provide the average and variance of
c     parameter j in each chain i_chain. N_steps is the number of steps
c     taken by each chain. See Doran & Muller, astroph/0311331. Return
c     the R-tations for each parameter (GR_R) and the average, as well
c     as flags which verify if these R's are <= the given R_threshold.
c
c     Cid@Lagoa - 04/Feb/2005

      SUBROUTINE Calc_GR_R(N_steps,N_chains,N_par,Cpar_ave,Cpar_var
     $     ,R_threshold,GR_R,GR_R_ave,IsAveBurnInOver,IsIndBurnInOver)


c     Definitions
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Npar_max = 601 , Nchains_max = 1000)
      real Cpar_ave(Npar_max,Nchains_max)
      real Cpar_var(Npar_max,Nchains_max)
      real GR_R(Npar_max)


c     Start loop over all j = 1...N_par parameters
      sum_R_ave = 0.
      do j=1,N_par

c     Compute global average of parameter j over all chains! (= "mean of
c     means").
         sum_p = 0.
         do i_chain=1,N_chains
            sum_p = sum_p + Cpar_ave(j,i_chain)
         end do
         global_ave = sum_p / float(N_chains)

c     Define GR_B & GR_W: GR_B = "Between chains variance", while GR_W =
c     "Within chains variance" (eqs 5 & 6 of Doran & Muller).
         sum_B = 0.
         sum_W = 0.
         do i_chain=1,N_chains
            sum_B = sum_B + (Cpar_ave(j,i_chain) - global_ave)
     $           * (Cpar_ave(j,i_chain) - global_ave)
            sum_W = sum_W + Cpar_var(j,i_chain)
         end do
         GR_B = (float(N_steps) / float(N_chains - 1)) * sum_B
         GR_W = sum_W / float(N_chains)

c     Define GR_R(j) = convergence ratio for parameter # j (eq 7 of
c     Doran & Muller). Also, update summation for mean R.
!     Introduced "protection" against W = j (chains which do not walk!).
!     When this happens (eg, due to a small N_sim), I set GR_R(j) to a
!     very large value, such that neither <GR_R> not GR_R(j) will
!     satisfy the convergence criterion!
         z = N_steps
         if (GR_W.LE.0.) then
!            write (*,'(a,i3,1x,f9.6))')
!     $           '[Calc_GR_R] W_j <= 0! This should NOT happen!',j,GR_W
            GR_R(j) = 10. * N_par * R_threshold + 9999.
         else 
            GR_R(j) = (((z - 1.) / z) * GR_W + GR_B / z) / GR_W
         end if
         sum_R_ave = sum_R_ave + GR_R(j)

c     Check whether GR_R(j) <= 0, which may happen due to precision bugs!
         if (GR_R(j).LE.0.) then
            write (*,'(a,i3,1x,1p,e12.3)')
     $           '[Calc_GR_R] R_j <= 0! This should NOT happen!' , j ,
     $           GR_R(j)
         end if
      end do


c     Mean GR_R over all N_par parameters
      GR_R_ave = sum_R_ave / float(N_par)


c     Verify if Burn-In is over for the Average & Individual R's, given
c     the input threshold. Set 0/1 (No/Yes) flags accordingly.
      IsAveBurnInOver = 1
      if (GR_R_ave.GT.R_threshold) IsAveBurnInOver = 0

      IsIndBurnInOver = 1
      do j=1,N_par
         if (GR_R(j).GT.R_threshold) IsIndBurnInOver = 0
      end do


      end
c ###########################################################################


c ###########################################################################
c # Generate random pop vector x(j), j=1...N, with a scheme
c # that allows a few components to dominate. (Simple normalization
c # "a la Laerte" does not do this). x() is normalized to sum x = 1 
c # (to within roundoff errors...)
c # Cid@Lynx - Dec/15/2002

c # Introduced a maximum x_j amplitude (x_max)
c # Cid@Lynx - Nov/24/2004

! This may be re-written using sort!

      SUBROUTINE Generate_PopVector(idum,N,x_max,x)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nb_max = 300)
      real x(Nb_max) , ind(Nb_max)

      if (N.GT.Nb_max) then
         write (*,*) ' [Generate_PopVector] N > N_max! ' , N , Nb_Max
         stop
      end if

      do j=1,N
         ind(j) = -99
      end do

      do j=1,N
 1       ix = N * ran2(idum) + 1
         do i=1,N
            if (ind(i).EQ.ix) goto 1
         end do
         ind(j) = ix
      end do

      sum_x = 0.
      do j=1,N-1
         i = ind(j)
 2       x(i) = (1. - sum_x) * ran2(idum)
         if (x(i).GT.x_max) goto 2

         sum_x = sum_x + x(i)
      end do
      i = ind(N)
      x(i) = 1. - sum_x

      end
c ###########################################################################




c ###########################################################################
c #              SAVE STEP IN A FILE (used to sample PDF...)                #
c ###########################################################################
*     Heavily simplifyed for MEX-code...
*     Cid@Lagoa -22/April/2007

      SUBROUTINE SaveStep(i_sim,i_chain,N_par,par,Ene)

!2BEDONE! Should save chi2_FIR, chi2_QHR...

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nb_max = 300,Npar_max = 601)
      real par(Npar_max) 
      common/OutputFileUnits/ i_UnitScubi

      write (i_UnitScubi,1) i_chain , 2.*Ene ,(100.*par(j),j=1,N_par)
 1    format (i3,1x,1p,e11.4,2x,0p,300(1x,f6.1))

      end
c ###########################################################################


!     MISCELANEOUS NOTES...

!     Probelma da enegia muito baixa... como consertar??? Tem a ver com
!     a expansao... uma possib eh sempre usar SlowChi2 qdo se acha um
!     minimo...
c     o prob eh que as vezes da o bug FDP antes do VF_chi2 se dar conta
c     de q a coisa vai mal (ie, antes do chi2 ficar negativo)...
c
c     bug tem q ser resolvido calibrando qdo ele acontece em funcao de
c     AV e YAV (possibly only AV+AYV...) Simulations? William??
c     Talvez o melhor seja incluir elementos avermelhados na base! WILLIAM!

!??   SERAH que nao vale a pena fixar AV, estimar x e refazer o treco AV
!??   a AV? Tipo o q a gente faz pra v0,vd com a Golden-search .,.....

!     tenho a impressao que podemos ser mais relaxed com o GR_R pra AV,
!     pois os valores das cadeias estao BEM proximos (~ centesimos de
!     mag) embora GR_R(AV) > 1.2....

!     have to go ASAP to a condensed base! OR TO Fixed v0 & vd (such
!     that RC terms need not be recomputed often!)

!     se nao mudar kin entao nao precisa recalcular RC sempre! (soh qdo
!     AV crescer muito!) o q apressaria bastante a coisa!  ideia: 1 -
!     convuliar obs e base com vd = 300 km/s; 2 - fazer um fit rapido; 3
!     descnvoluit tudo; 4 - estimar vd com um fit usando o pop-vector
!     rapido; 5 - fixar kinematics thereafter!
c +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!     MEX-routines    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c ###########################################################################
c #     !!MEX!!    ANNEALING DRIVER FOR SCUBIDU .......................     #
c ###########################################################################


C????? ABSTRACT!!!      
c     ... if IsKinFixed = 1 then RCtems & base-convolution are done only
c     once, before the main loop!

c     Fit f_obs(l_obs) using a recipe like this:
c     1- Smooth & shift BASE spectra to current v0 & vd (Of_base ==> f_base)
c     2- Increase weights & decrease step-size, ajdusting number of steps
c     3- Run Metropolis for (x,AV,YAV,fn) only, using smothed base (f_base)
c     4- Un-smooth base & estimate (v0,vd) for current population parameters
c       (A Golden-Section-Search is used to find v0 & vd)
c     5- Goto 1, starting from new (x,AV,YAV,fn,v0,vd) initial point.
c     Tests so far proved that this works beautifully!!
c     Cid@Lynx - 21/Aug/2003
c     ....

      SUBROUTINE MEX_Fit_Scubidu(i_verbose,idum,N_loops,v0_low,v0_upp
     $     ,vd_low,vd_upp,x_min,AV_min,exAV_min,v0,vd,Ntot_anneal
     $     ,eff_anneal,Nl_eff,iMassAgeZ_verbose,N_par,Cpar_last
     $     ,N_chains,Temp_ini,Temp_fin,fini_eps,ffin_eps,R_ini,R_fin
     $     ,N_sim,Nmax_steps,i_RestartChains,Gpar_min,GE_min,Gpar_ave
     $     ,Cpar_min,CEne_min,IsKinFixed,IsBurInOver)


c ***************************************************************************
c *                                Definitions                              *
c ***************************************************************************
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
      parameter (Nchains_max = 1000,Npar_max = 601)

      real eps(Npar_max) , eps_ini(Npar_max) , eps_fin(Npar_max)
      real aux1(Npar_max) , aux2(Npar_max)

      real Cpar_last(Npar_max,Nchains_max) , CEne_last(Nchains_max)
      real Cpar_min(Npar_max,Nchains_max)  , CEne_min(Nchains_max)
      real Cpar_low(Npar_max) , Cpar_upp(Npar_max)
      real Gpar_min(Npar_max) , Gpar_ave(Npar_max)

      real l_obs(Nl_max) , f_obs(Nl_max) , weight_obs(Nl_max)
      real l_aux(Nl_max) , f_aux1(Nl_max) , f_aux2(Nl_max)
      real x_min(Nb_max) , AV_tot(Nb_max) 
      real f_base(Nb_max,Nl_max) , Of_base(Nb_max,Nl_max)
      real f_syn(Nl_max) , fbase_norm(Nb_max) 

      common/Obs/ Nl_obs , l_obs , f_obs , weight_obs
      common/Base/ N_base , f_base
      common/Obase/ Of_base
      common/N_Gauss/ N_int_Gauss
      common/GoldenSearch_aux/ v0_min , vd_min ,
     &     Nl_aux , l_aux , GSf_syn_NoKin(Nl_max)
      common/RC_AV_Stuff/ NRC_AV , i_RC_CRASH_FLAG , NRC_AV_Default
      common/CalcMass/ q_norm , fbase_norm , i_FitPowerLaw , fobs_norm

c     Lower & upper limits on parameters (x,AV,exAV,fn)
      common/ParLimits/ par_low(Npar_max) , par_upp(Npar_max)

      external fun_vd_chi2
      external fun_v0_chi2

*ETC* Added a common to pass x & AV arrays to fun_vd_chi2 & fun_v0_chi2
      common/ETC_xAndAV4KinFits/ ETC_x4kin(Nb_max) , ETC_AV4kin(Nb_max)

!Cid14Jan07-FXK! Added FXK = Fixed kinematics option!!
      common/FXK_fits/ IsThisAnFXKFit

*     MEX-arrays
      real exAV_min(Nb_max) , exAV_test(Nb_max)
      common/exAV_stuff/ N_exAV , ind_exAV(Nb_max)

*     MEX-BruteForce (AV,exAV) search stuff
!FIX4CALIFA! Was parameter (n_pt_MaxBF = 10000,n_Dim_MaxBF = 10,n_BF_Range = 5), but now is:
      parameter (n_pt_MaxBF = 100,n_Dim_MaxBF = 10,n_BF_Range = 5)
      integer ind_BFmat(n_pt_MaxBF,n_Dim_MaxBF)
      real d_BF(Nb_max) , par_BF(Npar_max) , BF_low(Npar_max) ,
     $     BF_upp(Npar_max)
c ***************************************************************************



c ***************************************************************************
c *                              Initialization                             *
c ***************************************************************************
*ETC* Overcaution/paranoid reset of several things ... 
      do j=1,Nb_max
         exAV_test(j)  = 0.
         ETC_x4kin(j)  = 0.
         ETC_AV4kin(j) = 0.
      end do
      do j=1,Npar_max
         aux1(j) = 0.
         aux2(j) = 0.
         eps_ini(j) = 0.
         eps_fin(j) = 0.
      end do
      do i=1,Nl_max
         l_aux(i)  = 0.
         f_aux1(i) = 0.
      end do

c     Reset best kinematic parameters to initial values
      v0_min = v0 
      vd_min = vd

c     Define AV, exAV & fn positions in par-vector
      N_par  = N_base + 1 + N_exAV + 1
      j_AV   = N_base + 1
      j_exAV = N_base + 2
      j_fn   = N_par

c     Set step-size schedule in terms of input parameters fini_eps &
c     ffin_eps. We restrict AV/exAV eps sizes to up to 1 mag / fini_eps
c     (initial) & 1 mag / ffin_eps (final).
      do j=1,N_par
         eps_ini(j) = (par_upp(j) - par_low(j)) / fini_eps
         eps_fin(j) = (par_upp(j) - par_low(j)) / ffin_eps
      end do
      eps_ini(j_AV)  = min(eps_ini(j_AV),(1.0/fini_eps))
      eps_fin(j_AV)  = min(eps_fin(j_AV),(1.0/ffin_eps))
      do j=j_exAV,j_fn-1        !OBS: This loop only runs if Nex_AV > 0!
         eps_ini(j) = min(eps_ini(j),(1.0/fini_eps))
         eps_fin(j) = min(eps_fin(j),(1.0/ffin_eps))
      end do

c     Reset total efficiency & number of steps counters. They will count
c     all steps for all chains and all annealing loops.
      eff_anneal  = 0.
      Ntot_anneal = 0


!Cid14Jan07-FXK! If IsThisAnFXKFit = 1 then base is **already smoothed**
!     and the steps below should be skipped! Start with a sanity check:
!     If IsThisAnFXKFit = 1 then IsKinFixed MUST be 1!
      if ((IsThisAnFXKFit.EQ.1).AND.(IsKinFixed.NE.1)) then
         write (*,'(a,2(a,f9.3),a)')
     $        '[Fit_Scubidu] BUG/PROBLEM! IsThisAnFXKFit = ' ,
     $        IsThisAnFXKFit , 'but IsKinFixed = ', IsKinFixed
         stop
      end if

      if (IsThisAnFXKFit.NE.1) then

c     Store l_obs & Nl_obs onto auxiliar variables (used in all calls to
c     GaussSmooth below).
         do i=1,Nl_obs
            l_aux(i)  = l_obs(i)
         end do
         Nl_aux = Nl_obs

c     Smooth & shift original base Of_base to f_base for input vd &
c     v0. This has to be done here (before the annealing loop) because
c     we need to compute chain energies. These energies MAY be used to
c     initialize the cooling schedule (but we're not really doing it!!)
         if (i_verbose.GE.1) write (*,'(a,2(a,f9.3),a)')
     $        '[Fit_Scubidu] Smoothing base to' ,' v0 = ' , v0 ,
     $        '  & vd = ' , vd , ' km/s...'
         do j=1,N_base
            do i=1,Nl_obs
               f_aux1(i) = Of_base(j,i)
            end do
            call GaussSmooth(Nl_aux,l_aux,f_aux1,v0,vd,
     &           N_int_Gauss,Nl_obs,l_obs,f_aux2)
            do i=1,Nl_obs
               f_base(j,i) = f_aux2(i)
            end do
         end do
         if (i_verbose.GT.1) write (*,'(a)') '[Fit_Scubidu] ... done!'

      end if


c     Initialize RC-terms (so that we can compute energies!) Note that
c     the 1st-guess kinematical filter has just been applied. Presumably
c     the default number of RC_AV terms will be used the first time the
c     code goes through this call to Calc_RC_terms...
      if (i_verbose.GE.1) write (*,'(a)')
     $     '[Fit_Scubidu] Computing RC-terms ...'
      call Calc_RC_terms
      if (i_verbose.GE.1) write (*,'(a)')
     $     '[Fit_Scubidu] ... RC-terms computed!'


c     Compute initial chain energies & reset initial chain-minima.
c     OBS: Note that if a previous best minimum was found and NOT stored
c     in one of the Cpar_last's then it will be lost!!
      do i_chain=1,N_chains
         do j=1,N_par
            aux1(j)             = Cpar_last(j,i_chain)
            Cpar_min(j,i_chain) = Cpar_last(j,i_chain)
         end do
         Ene_initial = 0.5 * MEX_Fpar_VFchi2(N_par,N_base,aux1,0)
         CEne_last(i_chain) = Ene_initial
         CEne_min(i_chain)  = Ene_initial
      end do


c     Reset initial GLOBAL best-model parameters & energy (Gpar_ & GE_min).
      GE_min = 2. * CEne_min(1)
      do i_chain=1,N_chains
         if (CEne_min(i_chain).LT.GE_min) then
            do j=1,N_par
               Gpar_min(j) = Cpar_min(j,i_chain)
            end do
            GE_min = CEne_min(i_chain)
         end if
      end do


c     Unpack Gpar_min to x_min, AV_min, exAV_min & fn_min, and then
c     update AV_tot with AV_min & exAV_min...
      call MEX_Unpack_par2xAVYAVfn(N_par,N_base,Gpar_min,x_min,AV_min,
     $     exAV_min,fn_min)
      call MEX_Update_AV_array(N_par,N_base,AV_min,exAV_min,AV_tot)
c ***************************************************************************


c ***************************************************************************
c *                             Annealing loop                              *
c ***************************************************************************
c     1- Smooth & shift f_base spectra to current v0 & vd
c     2- Recompute RC-terms & test NRC_AV
c     3- Recompute Chain Energies (min, last & global-min)
c     4- Adjust Temperature, R_Threshold & step-sizes according to
c        cooling schedule
c     5- Run Subidu-Metropolis for (x,AV,YAV,fn) only, using smothed base
c        until GR-convergence is reached.
c     6- Do Quick brute-force search for better AV & YAV for fixed x,fn,v0,vd 
c     7- Restart Chains: put each chain back to its min-position & replace
c        worst chain by global minimum
c     8- Un-smooth base & estimate (v0,vd) for current population parameters
c     9- Goto 1, starting from new (x,AV,YAV,fn,v0,vd) initial point.

      do i_loop=1,N_loops

c     Given that the base was smoothed & RC-terms were already computed
c     in the initialization above, we do not have to re-do these steps
c     it in the 1st loop....
         if ((i_loop.GT.1).AND.(IsKinFixed.NE.1)) then


c     -------------------> Smooth & shift BASE spectra <---------------------
c     Smooth & shift original base Of_base to f_base for current vd & v0
c     The smoothed f_base is passed by common to other routines!  ATT:
c     f_base is NOT renormalized after convolution with G(v0,vd). After
c     checking my formulae, I concluded that renormalization is not
c     desirable afterall.  
            if (i_verbose.GE.1) write (*,'(a,2(f9.3),a)')
     $           '[Fit_Scubidu] REsmoothing base to ',v0,vd ,' km/s...'
            do j=1,N_base
               do i=1,Nl_obs
                  f_aux1(i) = Of_base(j,i)
               end do
               call GaussSmooth(Nl_aux,l_aux,f_aux1,v0,vd,N_int_Gauss
     $              ,Nl_obs,l_obs,f_aux2)
               do i=1,Nl_obs
                  f_base(j,i) = f_aux2(i)
               end do
            end do
            if (i_verbose.GE.1) write (*,'(a)') '[Fit_Scubidu] ...done!'
c     -----------------------------------------------------------------------


c     -----------------------> Recompute RC-terms <--------------------------
c     Since the base-spectra changed (new v0,vd), we need to re-compute
c     the RC-terms.  Before doing so, I estimate the number NRC_AV of
c     terms in the AV-series as a function of AV_min, such that NRC_AV
c     increases linearly with |AV_min| up to a maximum value. This
c     relation needs to be calibrated through experiments...  In any
c     case, the test a few lines below fix any bad-guess for NRC_AV.
!!MEX!! Picking largest AV_tot for this ...
            if (i_verbose.GE.1) write (*,'(a)')
     $           '[Fit_Scubidu] Computing RC-terms ...'

            call MEX_Unpack_par2xAVYAVfn(N_par,N_base,Gpar_min,x_min,
     $           AV_min,exAV_min,fn_min)
            call MEX_Update_AV_array(N_par,N_base,AV_min,exAV_min,
     $           AV_tot)
            z1 = 0.
            do j=1,N_base
               if (abs(AV_tot(j)).GT.z1) z1 = abs(AV_tot(j))
            end do
            z1 = 5. + 6. * z1
            NRC_AV = min(z1,1.*NRC_AV_Default)
            call Calc_RC_terms
            if (i_verbose.GE.1) write (*,'(a,i3)')
     $           '[Fit_Scubidu] ... RC-terms re-computed: NRC_AV = '
     $           , NRC_AV
c     -----------------------------------------------------------------------


         end if


c     -------------> Test whether more RC terms  are needed <----------------
c     Compare the approximate (fast) and exact (slow) chi2's. If the
c     relative difference exceeds a threshold, then increase the number
c     of terms (NRC_AV) in the AV-truncated-series until convergence is
c     achieved. More terms are needed when AV increases... Since the two
c     terms always agree well for AV = 0, I carry out the test with AV =
c     0.3 (typical of SDSS galaxies...) unless |AV_min| is larger, and
c     likewise for YAV. Note also that the chi2 routines themselves may
c     change NRC_AV when they detect trouble (ie, negative chi2's)!
c     Hence, we've got plenty of protection. (ATT: we may fail to detect
c     a bad chi2 which, say, is going towards < 0 but hasn't yet gotten
c     there...)
!!MEX!! Checked this part and it seems Ok! It picks a large AV & exAV for
!     the test, which is what should be done.
 100     AV_test = max(0.3,(abs(AV_min)))
         aux_max_exAV_test = -1000.
         do jMEX=1,N_exAV
            exAV_test(jMEX) = max(0.3,(abs(exAV_min(jMEX))))
            if (exAV_test(jMEX).GT.aux_max_exAV_test) aux_max_exAV_test 
     $           = exAV_test(jMEX)
         end do
         call MEX_Update_AV_array(N_par,N_base,AV_test,exAV_test,AV_tot)

         z1 = VeryFastF_chi2(x_min,AV_tot,0)
         z2 = SlowF_chi2(x_min,AV_tot)
         rel_diff = (z1 - z2) / z2
         dchi2_threshold = 1.e-3

* !!MEX!! Slight change
         if (abs(z1 - z2).GT.dchi2_threshold * abs(z2)) then
            write (*,11) 'fitspec',i_loop,' ATT: Increasing'
     $           , ' # of terms in AV-series...' , z1 , z2 , rel_diff
     $           , dchi2_threshold , AV_test , AV_min ,
     $           aux_max_exAV_test , NRC_AV , ' ==> ' , NRC_AV +
     $           1
 11         format (a7,i3,a,a,1p,4(2x,e12.5),0p,3(f9.3,1x),2x,i2,a,i2)

            NRC_AV = NRC_AV + 1
            call Calc_RC_terms

!     CRASH: Get out in case RC-crashes! Otherwise repeat test
            if (i_RC_CRASH_FLAG.NE.1) goto 100
         end if

*     Screen output - controling difference between Slow & Fast - chi2!
         if (i_verbose.GT.0) write (*,12) i_loop , z1,z2, 100.* rel_diff
     $        , NRC_AV , AV_test , (exAV_test(jMEX),jMEX=1,N_exAV)
 12      format ('FitScub',i3,'  RC-test: Fast chi2=' ,1p,e12.4,2x,
     $        'Slow chi2=' ,e12.4,3x,' Dif=',e12.4,0p,
     $        '%  for NRC_AV=',i3,' w/AV_test=',f5.2,
     $        ' & exAV_test=',10(f5.2,1x))
c     -----------------------------------------------------------------------


c     -----------------------> Recompute Energies <--------------------------
c     Since v0 & vd may have changed, and maybe NRC_AV too, we must
c     recompute the minimum & last chain energies & the global minimum
c     before passing these values to ScubiduCadeVoce... Note that
c     energies are recomputed even if neither v0, vd nor NRC_AV changed,
c     which should do no harm.
         do i_chain=1,N_chains
            do j=1,N_par
               aux1(j) = Cpar_min(j,i_chain)
               aux2(j) = Cpar_last(j,i_chain)
            end do
            CEne_min(i_chain)  =0.5*MEX_Fpar_VFchi2(N_par,N_base,aux1,0)
            CEne_last(i_chain) =0.5*MEX_Fpar_VFchi2(N_par,N_base,aux2,0)
         end do
         GE_min = 0.5 * MEX_Fpar_VFchi2(N_par,N_base,Gpar_min,0)
c     -----------------------------------------------------------------------


c     ---------------------> Metropolis/Scubidu Run <------------------------
c     Update Temperature, R & step-sizes. Then run ScubiduCadeVoce.
         z1 = 0.
         if (N_loops.GT.1) z1 = float(i_loop - 1) / float(N_loops - 1)
         Temperature = Temp_ini * (Temp_fin / Temp_ini)**z1
         R_Threshold = R_ini * (R_fin / R_ini)**z1
         do j=1,N_par
            eps(j) = eps_ini(j) * (eps_fin(j) / eps_ini(j))**z1
         end do

         call MEX_ScubiduCadeVoce(idum,N_par,N_chains,N_sim,Nmax_steps
     $        ,Temperature,R_Threshold,i_verbose,IsBurInOver,Cpar_last
     $        ,CEne_last,eps,eff_Scubi,N_steps,Gpar_min,GE_min,Cpar_min
     $        ,CEne_min,Cpar_low,Cpar_upp,Gpar_ave)
c     -----------------------------------------------------------------------


*     -----------------------------------------------------------------------
*     -->  "Quick" brute force search for better AV/exAV's for each chain <--
*     BLABLABLA... BIG CHANGE IN SCHEME WITH ind_MatBF...
*     Cid@Lagoa - 22/April/2007
*     -----------------------------------------------------------------------

*     i_BruteForceFlag = just a flag to tell if the BF seach found a
*     better minimum for at least one chain! Reset is before we start.
         i_BruteForceFlag  = 0

*     nBF_Dim = number of dimensions to BF-search (1 for each of AV and
*     each of the exAV's)
         nBF_Dim = 1 + N_exAV
         if (nBF_Dim.GT.n_Dim_MaxBF) then
            write (*,*) '[Fit_Scubidu] OOPS!!!! too much!! ' , nBF_Dim
            stop
         end if

*     n_pt_PerDim = number of grid points in each dimension
*     n_pt_BF     = total number of points in BF-search grid
         n_pt_PerDim = (float(n_pt_MaxBF))**(1. / float(nBF_Dim))
         n_pt_BF     = n_pt_PerDim**nBF_Dim


*     echo BF numbers
         if (i_verbose.GT.0) then
            write (*,'(a,3(a,i5),a)')
     $           '[Fit_Scubidu] Starting BruteForce AV/exAV:' ,
     $           ' nBF_Dim = ' , nBF_Dim , '  n_pt_PerDim  = ' ,
     $           n_pt_PerDim , '  n_pt_BF = ' , n_pt_BF , ' ...'
         end if

*     Build ind_BFMat matrix
         do i_pt=1,n_pt_BF
            do i_Dim=1,nBF_Dim
               nS = n_pt_PerDim**(nBF_Dim - i_Dim)
               nL = n_pt_PerDim**(nBF_Dim - i_Dim + 1)
               ia = 1 + int((i_pt - 1) / nL) 
               ib = nL * (ia - 1) + 1
               ind_BFmat(i_pt,i_Dim) = (i_pt - ib) / nS + 1
            end do
!            write (*,'(a,i5,a,15i5)') 'HELP==BF ',i_pt,' : ',
!     $           ((ind_BFmat(i_pt,i_Dim)),i_Dim=1,nBF_Dim)
         end do


*     BF-search for each chain
         do i_chain=1,N_chains

*     Define BF step sizes for each of the dimensions
*     OBS: The parameter n_BF_Range defines the range of AV/exAV search
*     as +/- n_BF_range X the current corresponding step size. If this
*     range is < 0.1 mag then use 0.1 instead.
            do j=j_AV,j_fn-1
               BF_range  = max(n_BF_Range*eps(j),0.1)
               BF_low(j) = Cpar_min(j,i_chain) - BF_range
               BF_upp(j) = Cpar_min(j,i_chain) + BF_range
               if (BF_low(j).LT.par_low(j)) BF_low(j) = par_low(j)
               if (BF_upp(j).GT.par_upp(j)) BF_upp(j) = par_upp(j)
               d_BF(j)   = (BF_upp(j) - BF_low(j)) / (n_pt_PerDim - 1)
            end do

*     Set par_BF = best in this chain
            do j=1,N_par
               par_BF(j) = Cpar_min(j,i_chain)
            end do

*     BF-grid-search for better AV & exAV for current chain.
*     For each of the n_pt_BF points in the DCP-BF-grid, modify the AV &
*     exAV entries in par_BF. Then compute the energy for par_BF &
*     compate it to CE_min. If the new model is better, store it (and
*     update i_BruteForceFlag = 1 to inform this was all worth the
*     trouble!).
            do i_pt=1,n_pt_BF
               do i_Dim=1,nBF_Dim
                  j = j_AV + i_Dim - 1
                  par_BF(j) = BF_low(j) + (ind_BFmat(i_pt,i_Dim) - 1)
     $                 * d_BF(j)
               end do

               Ene_BF = 0.5 * MEX_Fpar_VFchi2(N_par,N_base,par_BF,0)

               if (Ene_BF.LT.CEne_min(i_chain)) then
!                  write (*,'(a,i5,i5,3(a,f9.3),a,30(1x,f7.3))') 'BF ',
!     $                 i_chain , i_pt , ' Eold=' , CEne_min(i_chain) ,
!     $                 ' Enew=' , Ene_BF ,
!     $                 ' Diff =' , Ene_BF - CEne_min(i_chain) ,
!     $                 ' | par_BF= ', (par_BF(j),j=j_AV,j_fn-1)
                  CEne_min(i_chain) = Ene_BF
                  do j=1,N_par
                     Cpar_min(j,i_chain) = par_BF(j)
                  end do
                  i_BruteForceFlag = 1
               end if
            end do

*     Goto BF-search for next chain
         end do


*     Now repeat the same brute-force search for the Global-best model
*     so far...
         do j=j_AV,j_fn-1
            BF_range  = max(n_BF_Range*eps(j),0.1)
            BF_low(j) = Gpar_min(j) - BF_range
            BF_upp(j) = Gpar_min(j) + BF_range
            if (BF_low(j).LT.par_low(j)) BF_low(j) = par_low(j)
            if (BF_upp(j).GT.par_upp(j)) BF_upp(j) = par_upp(j)
            d_BF(j)   = (BF_upp(j) - BF_low(j)) / (n_pt_PerDim - 1)
         end do

         do j=1,N_par
            par_BF(j) = Gpar_min(j)
         end do

         do i_pt=1,n_pt_BF
            do i_Dim=1,nBF_Dim
               j = j_AV + i_Dim - 1
               par_BF(j) = BF_low(j) + (ind_BFmat(i_pt,i_Dim) - 1)
     $              * d_BF(j)
            end do

            Ene_BF = 0.5 * MEX_Fpar_VFchi2(N_par,N_base,par_BF,0)

            if (Ene_BF.LT.GE_min) then
!               write (*,'(a,i5,3(a,f9.3),a,30(1x,f7.3))') 'GBF ',
!     $              i_pt , ' Eold=' , GE_min ,
!     $              ' Enew=' , Ene_BF , ' Diff =' , Ene_BF - GE_min ,
!     $              ' | par_BF= ', (par_BF(j),j=j_AV,j_fn-1)
               GE_min = Ene_BF
               do j=1,N_par
                  Gpar_min(j) = par_BF(j)
               end do
               i_BruteForceFlag = 2
            end if
         end do

*     Done with BF-searches
         if (i_verbose.GT.0) write (*,'(a,i2)')
     $        '[Fit_Scubidu] ... done! i_BruteForceFlag = ' ,
     $        i_BruteForceFlag
*     -----------------------------------------------------------------------


c     ---------------------> A few simple things ... <-----------------------
c     Re-find best Global minimum among the chains (which may have
c     changed after the brute-force AV/YAV minimization above).
         do i_chain=1,N_chains
            if (CEne_min(i_chain).LE.GE_min) then
               do j=1,N_par
                  Gpar_min(j) = Cpar_min(j,i_chain)
               end do
               GE_min = CEne_min(i_chain)
            end if
         end do

c     Update total annealing step-counter & efficiency summation. Note
c     that all chains are counted, but Brute-Force steps are not.
         Ntot_anneal = Ntot_anneal + N_chains * N_steps
         eff_anneal  = eff_anneal + eff_Scubi * N_chains * N_steps
c     -----------------------------------------------------------------------


c     ----------------> Replace Cpar_last by Cpar_min... <-------------------
c     Here we reset Cpar_last to Cpar_min, so that next time it will
c     start from its best location. Also, the chain with highest energy
c     is replaced by best Global model so far. This should improve their
c     convergence a tiny little bit. To use this feature you must set
c     i_RestartChains = "YES" (= 1).

!Cid01April07-FIX! Found a "BUG" here while preparing the final
!     v04! I did not realize this before because of another bug (namely,
!     that we did not use i_RestartChains = 1 in the BI-phase fits)! The
!     new bug is that by setting the last values equal to the best ones,
!     if one of these best chain values = the global best so far then
!     the trick of resetting the chain with worst energy = to the best
!     global model may result in this chain-values appearing TWICE!! Ie,
!     2 of the chains would start from IDENTICAL positions. This is not
!     a killer bug, but I never wanted to do this...  
!
!     To fix it I 1st figure out if any of the chain-minima differ from
!     the global-minimum by less than 1/1000 in energy (this is flagegd
!     with i_GEminIsAmongCEmins = 1 below). ONLY if NONE of them is this
!     close to the minimum I then apply the Worst2Best chain trick. (I
!     suspect that this is never gonna happen, but to be on the safe
!     side I decided to code it in this way).
!
!     I also introduce a new option: i_RestartChains = 2, which will
!     reset ONLY the worst chain to the best model, but leave the others
!     wherever they are.
         if (i_RestartChains.EQ.1) then
            i_GEminIsAmongCEmins = 0
            do i_chain=1,N_chains
               do j=1,N_par
                  Cpar_last(j,i_chain) = Cpar_min(j,i_chain)
               end do
               CEne_last(i_chain) = CEne_min(i_chain)

               test_diff = abs(GE_min - CEne_min(i_chain))
               if (test_diff.LE.(0.001*GE_min)) i_GEminIsAmongCEmins = 1
            end do

            if (i_GEminIsAmongCEmins.EQ.0) then
               WorstEne     = CEne_min(1)
               i_WorstChain = 1
               do i_chain=2,N_chains
                  if (CEne_min(i_chain).GT.WorstEne) then
                     WorstEne     = CEne_min(i_chain)
                     i_WorstChain = i_chain
                  end if
               end do
               do j=1,N_par
                  Cpar_last(j,i_WorstChain) = Gpar_min(j)
                  Cpar_min(j,i_WorstChain)  = Gpar_min(j)
               end do
               CEne_last(i_WorstChain) = GE_min
               CEne_min(i_WorstChain)  = GE_min
               if (i_verbose.GT.1) write (*,'(a,a)')
     $              '[Fit_Scubidu] GEmin NOT among the CEne_mins!' ,
     $              ' Worst2Best trick applied.'
            end if
         end if

!     New option: With i_RestartChains = 2 only the worst chain is
!     reset, and to the best global minimum.
         if (i_RestartChains.EQ.2) then
            WorstEne     = CEne_min(1)
            i_WorstChain = 1
            do i_chain=1,N_chains
               if (CEne_min(i_chain).GT.WorstEne) then
                  WorstEne     = CEne_min(i_chain)
                  i_WorstChain = i_chain
               end if
            end do

            do j=1,N_par
               Cpar_last(j,i_WorstChain) = Gpar_min(j)
               Cpar_min(j,i_WorstChain)  = Gpar_min(j)
            end do
            CEne_last(i_WorstChain) = GE_min
            CEne_min(i_WorstChain)  = GE_min
         end if
c     -----------------------------------------------------------------------


c     -----------------------------------------------------------------------
c     Unpack Gpar_min onto x_min, AV_min, YAV_min & fn_min
         call MEX_Unpack_par2xAVYAVfn(N_par,N_base,Gpar_min,x_min,
     $        AV_min,exAV_min,fn_min)
         call MEX_Update_AV_array(N_par,N_base,AV_min,exAV_min,AV_tot)
c     -----------------------------------------------------------------------


c     -----> Golden Search for better kinematical parameters: v0 & vd <------
         if (IsKinFixed.NE.1) then


c     1st RESTORE f_base to original un-smoothed & un-shifted base spectra!!
            if (i_verbose.GT.1) write (*,'(a)')
     $           '[Fit_Scubidu] Starting Re-fit of v0 & vd...'
            do i=1,Nl_obs
               do j=1,N_base
                  f_base(j,i) = Of_base(j,i)
               end do
            end do


c     Define "Kinematics-free" spectrum GSf_syn_NoKin, passed by the
c     common/GoldenSearch_aux/. 
            call calc_Fsyn(x_min,AV_tot,GSf_syn_NoKin)


*ETC* The external functions fun_vd_chi2 & fun_v0_chi2 called bellow do
*     NOT admit x and AV_tot as arguments (they must have just 1
*     argument to be used with NR's mnbrak & golden routines). I thus
*     store x_min & AV_tot in ETC_x4kin() & ETC_AV4kin(), now acessible
*     to fun_vd_chi2 & fun_v0_chi2 via common/ETC_xAndAV4Kinematics/.
            do j=1,N_base
               ETC_x4kin(j)  = x_min(j)
               ETC_AV4kin(j) = AV_tot(j)
            end do


c     Compute chi2 of current model - USES fun_v0_chi2, which convolves
c     GSf_syn_NoKin through the kinematical filter! 
c     OBS1: This is the same as chi2 = FK_chi2(x_min,AV_tot,v0_min,vd_min)
c     OBS2: This should differ slightly from the chi2_min returned from
c           the Metropolis routine because in fun_v0_chi2 & fun_vd_chi2 
c           the reddening term is convolved with the Gaussian filter too.
c     OBS3: Note that here NO AV-truncated-series trick is used in neither 
c           fun_v0_chi2 nor fun_vd_chi2!
            chi2 = fun_v0_chi2(v0_min)

c     Find (with a Golden-Section-Search method) the vd & v0 which
c     produce the best chi2-match using the current best (x,AV,YAV)
c     estimated above. The searches are done separately, and in this
c     order. Values of v0 & vd outside their _low ==> _upp limits are
c     rejected. Seems to work fine, but it's rather slow. (Can be
c     improved with FFT-convolution...). OBS: In Starlight05 I
c     introduced the chi2 = chi2_vd & chi2 = chi2_v0 resets which were
c     not present before.
            ax  = 0.
            bx  = 500.
            tol = 1.e-6
            call mnbrak(ax,bx,cx,fa,fb,fc,fun_vd_chi2)
            chi2_vd = golden(ax,bx,cx,fun_vd_chi2,tol,vd_Gold)
            vd_Gold = abs(vd_Gold)
            if ((chi2_vd.LT.chi2).AND.(vd_Gold.GE.vd_low).AND.
     &           (vd_Gold.LE.vd_upp)) then
               vd_min = vd_Gold
               chi2   = chi2_vd
            end if

            ax  = -300.
            bx  = 300.
            tol = 1.e-6
            call mnbrak(ax,bx,cx,fa,fb,fc,fun_v0_chi2)
            chi2_v0 = golden(ax,bx,cx,fun_v0_chi2,tol,v0_Gold)
            if ((chi2_v0.LT.chi2).AND.(v0_Gold.GE.v0_low).AND.
     &           (v0_Gold.LE.v0_upp)) then
               v0_min = v0_Gold
               chi2   = chi2_v0
            end if              

            if (i_verbose.GT.1) write (*,'(a)') '[Fit_Scubidu] ...done!'
         end if


c     Store v0_min & vd_min onto v0 (we could use only the *_min variables...)
         v0 = v0_min
         vd = vd_min
c     -----------------------------------------------------------------------


c     ------------------------> Screen output <------------------------------
!     Screen check! 
         if (i_verbose.GT.0) then
            write (*,912)
 912        format ('#   i_loop    Temp  eps           N ',
     &           'chi2         adev    eff  |   v0    vd   | ',
     &           'sum      A_V   |    x_j... [%]')

            if (IsKinFixed.EQ.1) then
               chi2 = 2. * GE_min / float(Nl_eff)
               call calc_Fsyn(x_min,AV_tot,f_syn)
            else 
               chi2 = chi2 / float(Nl_eff)
               call Kcalc_Fsyn(Nl_obs,x_min,AV_tot,v0_min,vd_min,f_syn)
            end if
            adev  = F_adev(f_syn)
            sum_x = 100. * Sum_of_x(x_min,N_base)

            write (*,913) 'FitScub' , i_loop , Temperature, eps(1) ,
     &           Ntot_anneal , chi2 , adev , 100.*eff_Scubi ,
     &           v0_min , vd_min ,
     &           sum_x , AV_min , 
     &           (100.*x_min(j),j=1,min(10,N_base))
 913        format (a7,i3,1x,1p,e7.1,0p,2x,f5.3,1x,i9,1x,1p,e10.4,0p,2x,
     $           f7.4,2x,f4.1,' | ',f6.1,1x,f5.1,' | ', f7.3,1x,f5.2,1x,
     $           ' | ', 10(1x,f5.2))
            do j_out=2,int(N_base/10)+1
               j1 = (j_out - 1) * 10 + 1
               j2 = min(j1+9,N_base)
               write (*,'(96x,10(1x,f5.2))') (100.*x_min(j),j=j1,j2)
            end do

            write (*,'(96x,a)') '   AV_j... [mag]'
            do j_out=1,int(N_base/10)+1
               j1 = (j_out - 1) * 10 + 1
               j2 = min(j1+9,N_base)
               write (*,'(96x,10(1x,f5.2))') (AV_tot(j),j=j1,j2)
            end do

            write (*,'(97x,a,1p,e12.4)') 'Total chi2 = ', chi2 * Nl_eff

!     Screen output mean age & z to check convergence. 
!     OBS: Gives stupid results when called from EXOs-routine!
            if (iMassAgeZ_verbose.GT.0) then
               call Calc_Mean_Age_and_Z(N_base,x_min,a_age,s_age,
     &              a_met,s_met)
               write (*,914) a_age , s_age , a_met , s_met
 914           format(97x,'<log t> = ',f7.4,' +/- ',f7.4,
     &              '   <Z> = ',f7.4,' +/- ',f7.4)


!     Screen output of "total mass". Actually, Mass_tot is just
!     proportional to the true mass. I just output it here to have a
!     coarse idea of how robust the total mass is...  
!     OBS: Gives stupid results when called from EXOs-routines!
!     Cid@lynx 19/Mar/2004
               N_aux = N_base
               if (i_FitPowerLaw.EQ.1) N_aux = N_base - 1
               Mass_tot = 0.
               do j=1,N_aux
                  Mass_tot = Mass_tot + x_min(j) * 
     &                 (10.**(0.4 * AV_tot(j) * q_norm)) / fbase_norm(j)
               end do
               Mass_tot = Mass_tot * 383171.17
               write (*,915) Mass_tot , log10(Mass_tot)
 915           format(97x,'Mass = ',1p,e12.3,6x,'log Mass = ',0p,f8.3,
     &           '  [ref units]')
            end if

!     Inform whether Brute-Force AV/exAV search helped...'
               if (i_BruteForceFlag.GT.0) write (*,'(97x,a)') 
     &           'OBS: Brute-Force AV/exAV search helped in this loop!'
         end if
c     -----------------------------------------------------------------------


c     Goto next annealing-loop iteration
      end do
c ***************************************************************************


c     Compute final annealing efficiency (just for output) & go home:)
      eff_anneal = eff_anneal / Ntot_anneal


c     RESTORE f_base to original un-smoothed & un-shifted base spectra
c     if no kinematical fitting was atempted above.
      if (IsKinFixed.EQ.1) then
         do i=1,Nl_obs
            do j=1,N_base
               f_base(j,i) = Of_base(j,i)
            end do
         end do
      end if


      end
c ###########################################################################




c ###########################################################################
c # !!MEX!!                     SCUBIDUCADEVOCE                             #
c ###########################################################################
c     Evolve N_chains with Metropolis N_sim steps at a time, iterating
c     with more steps until GR-convenrged is achieved or more than
c     Nmax_steps are taken.
c
c     LOTS of details to comment!!!!!!!
c
c     Cid@Lagoa - 25/Feb/2005

!OPTfn! Implemented optimal fn! See more on this below.
!     OBS: When IsOptimize_fn_OPT.EQ.1 the normalization factor fn will be
!     computed analyticaly. It is therefore pointless to (1) allow it to
!     be perturbed as any other par, and (2) to verify its low <--> upp
!     limits. We do both these pointless but HARMELESS things just for
!     fear of changing the code... 
!
!     ElCid@Sanchica - 28/Jan/2012


      SUBROUTINE MEX_ScubiduCadeVoce(idum,N_par,N_chains,N_sim
     $     ,Nmax_steps,Temperature,R_Threshold,i_verbose,IsBurInOver
     $     ,Cpar_last,CEne_last,eps_in,eff_Scubi,N_steps,Gpar_min,GE_min
     $     ,Cpar_min,CEne_min,Cpar_low,Cpar_upp,Gpar_ave)


c ***************************************************************************
c *				DEFINITIONS				    *
c ***************************************************************************
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nb_max = 300)
      parameter (Nchains_max = 1000,Npar_max = 601)

c     Current & tentative new parameter, and step-size-scales for each
c     parameter (for built-in Metropolis)
      real par(Npar_max) , par_try(Npar_max)
      real eps_in(Npar_max) , eps(Npar_max)

c     Global minimum & average & std. dev.
      real Gpar_min(Npar_max) , Gpar_ave(Npar_max) , Gpar_sig(Npar_max)

c     "Cpar" denotes Chain-parameter (ie, for each chain...) matrices...
      real Cpar_ave(Npar_max,Nchains_max)
      real Cpar_var(Npar_max,Nchains_max)
      real Cpar_last(Npar_max,Nchains_max)
      real Cpar_min(Npar_max,Nchains_max)
      real Cpar_low(Npar_max) , Cpar_upp(Npar_max)

c     Chain energies (current/last & minimum)
      real CEne_last(Nchains_max) , CEne_min(Nchains_max)

c     Summations of par & par**2 "Memory" arrays. Using double precision!
      double precision S_par(Npar_max,Nchains_max)
      double precision S_par2(Npar_max,Nchains_max)

c     Scaling factor & efficiency arrays for each chain
      real Chain_alpha(Nchains_max) , Chain_eff(Nchains_max)

c     GR_R(j) = Gelman & Rubin covergence ratio for parameter j.
      real GR_R(Npar_max)

c     Auxiliar array(s)
      real aux1(Npar_max) , aux2(Npar_max)

c     We don't change NRC_AV directly here, but it's good to know if it
c     changes on the way...
      common/RC_AV_Stuff/ NRC_AV , i_RC_CRASH_FLAG , NRC_AV_Default

c     Lower & upper limits on parameters (x,AV,YAV,fn)
      common/ParLimits/ par_low(Npar_max) , par_upp(Npar_max)

c     ScubiduCadeVoce Technical Parameters
      common/Scubi_TechnicalParameters/ alpha , Falpha , IsGRTestHard ,
     $     eff_IDEAL , i_UpdateEps , Feps , i_UpdateAlpha ,
     $     i_MoveOneParOnly , i_HelpParWithMaxR , prob_jRmax ,
     $     Npar4GRtest ,i_UpdateAVYAVStepSeparately ,
     $     ind_HelpParMove2Average(Npar_max) ,
     $     prob_HelpParMove2Average(Npar_max)

      common/OutputFileUnits/ i_UnitScubi

*     !!MEX!!
      common/exAV_stuff/ N_exAV , ind_exAV(Nb_max)

!OPTfn! Optically optimal fn - computed elsewhere and used here to adjust fn!
      common/Optimize_fn_OPT/ IsOptimize_fn_OPT , Optimal_fn_OPT
c ***************************************************************************


c ***************************************************************************
c *                             INITIALIZATION                              *
c ***************************************************************************
c     Reset what needs to be reset'd! 
      N_steps = 0
      do j=1,N_par
         GR_R(j) = 0.
         do i_chain=1,N_chains
            S_par(j,i_chain)    = 0.
            S_par2(j,i_chain)   = 0.
         end do
      end do

      do i_chain=1,N_chains
         Chain_alpha(i_chain) = alpha
      end do


c     Counter of number of chain-minima found throughout this Scubi run
c     (just for control).
      N_min  = 0


c     Store input eps (eps_in) in eps (which may change during the run).
      do j=1,N_par
         eps(j) = eps_in(j)
      end do


c     Define who's who in par() vector.
      N_base = N_par - 1 - N_exAV - 1
      j_AV   = N_base + 1
      j_exAV = N_base + 2
      j_fn   = N_par


c     Reset "Burn-In" flag & global-efficiency, as well as the index of
c     parameter with largest GR_R, and reset NRC_AV_last to its initial
c     value.
      IsBurInOver = 0
      eff_Scubi   = 0.
      jR_max      = 0
      NRC_AV_last = NRC_AV


c     Defines how many pars will enter the GR-test (exclude fn).  It's
c     necessary to define it here in case we call this routine using a
c     condensed base...
      Npar4GRtest = N_par - 1 
c ***************************************************************************


c ***************************************************************************
c *               BURN-IN CONVERGENCE LOOP (at fixed Temperature)           *
c ***************************************************************************
c     Define N_steps = number of steps taken in EACH chain AFTER the
c     next call to Metropolis is over. (Will come back here & add more
c     points if GR-convergence test fails...)
 100  N_steps = N_steps + N_sim 


c     Reset counters of total & rejected AV/YAV steps (valid only for
c     the next N_sim * N_chains moves).
      N_AVYAV_steps = 0
      N_AVYAV_skip  = 0


c     -------------------- Run N_chains Metropolis chains -------------------
c     Run N_chains Metropolis chains, each N_sim steps long.
      do i_chain=1,N_chains


c     ++++++++++++++++++ Metropolis run for chain # i_chain +++++++++++++++++
c     (0) Reset initial par-vector, Energy & N_skip. Start Metropolis run.
c     (1) Modify par => par_try (within low => upp boundaries).
c         Move can be in all par_j's or in just one-at-a-time, 
c         according to i_MoveOneParOnly = 0/1.
c     (2) Compute new energy (Ene_try); then accept or not the change.
c     (3) Verify if new state is an energy minimum (for current chain).
c     (4) Updates summations of p_j & p_j**2 & goto next step (i_sim).
c     (5) When done, compute chain-efficiency.
         do j=1,N_par
            par(j) = Cpar_last(j,i_chain)
         end do

!     OBS: Ene MUST be computed with a direct call to Fpar_VFchi2 with
!     i_IsAVstillTheSame = 0 = "NO" here, since AV/YAV varies from chain
!     to chain.
         Ene    = 0.5 * MEX_Fpar_VFchi2(N_par,N_base,par,0)
         N_skip = 0


!Cid06Dec05-FIX! Introduce counter of repeated steps (ie those which
!     don't vary)...
         N_Retry = 0

         do i_sim=1,N_sim


c     For safety, recompute Ene & CEne_min here if NRC_AV changed since
c     last time an energy was computed for this chain. NRC_AV_last
c     itself is only updated when all chains have run N_sim steps.
            if (NRC_AV.NE.NRC_AV_last) then
               Ene = 0.5 * MEX_Fpar_VFchi2(N_par,N_base,par,0)
               do j=1,N_par
                  aux1(j) = Cpar_min(j,i_chain)
               end do
               CEne_min(i_chain) = 0.5 * MEX_Fpar_VFchi2(N_par,N_base,
     $              aux1,0)
            end if


            if (i_MoveOneParOnly.EQ.1) then

c     ........... Modifies par --> par_try: One par_j per move!..............
c     First set par_try = par. Then chose ONE of the N_par variables &
c     perturb it. If the chosen variable doesn't really vary (eg, it was
c     and still is = to its lower limit), then try again, otherwise this
c     step would be wasted. This ensures that all N_sim steps actually
c     attempt to change the state. If i_HelpParWithMaxR = 1 then the
c     parameter with largest R_j will have an extra chance of beeing
c     picked.

c     par_try <= par
               do j=1,N_par
                  par_try(j) = par(j)
               end do

c     Pick jj 
 10            if ((i_HelpParWithMaxR.EQ.1).AND.(jR_max.GT.0).AND.
     $              (ran2(idum).LE.prob_jRmax)) then
                  jj = jR_max
               else
                  jj = N_par * ran2(idum) + 1
               end if 
               if (jj.GT.N_par) goto 10
!OPTfn! Pointless (but harmless) to pick jj = j_fn when IsOptimize_fn_OPT = 1.


c     Perturb jj-th component of par(jj) => par_try(jj).
               par_try(jj) = par(jj) + Chain_alpha(i_chain) * eps(jj)
     $              * gasdev(idum)


c     Here we give a little extra help to par_j, inverting the sign of
c     the perturbation (with probability = prob_HelpParMove2Average) if
c     the move took par_j away from the current global-average
c     value. This makes par_j converge faster, and is analogous to using
c     a smarter parameter-change proposal-scheme. Note that Gpar_ave is
c     only defined after N_steps > N_sim (ie, after the 1st loop).
               if ((ind_HelpParMove2Average(jj).EQ.1).AND.(N_steps.GT.
     $              N_sim)) then
                  z1    = par(jj) - (par_try(jj) - par(jj))
                  dist0 = abs(par_try(jj) - Gpar_ave(jj))
                  dist1 = abs(z1 - Gpar_ave(jj))
                  if ((dist1.LT.dist0).AND.(ran2(idum).LE.
     $                 prob_HelpParMove2Average(jj))) par_try(jj) = z1
               end if


c     Limit par_try(jj) to its boundaries & try again if nothing changed.
!OPTfn! Pointless (but harmless) if jj = j_fn and IsOptimize_fn_OPT = 1.
               if (par_try(jj).LT.par_low(jj)) par_try(jj) = par_low(jj)
               if (par_try(jj).GT.par_upp(jj)) par_try(jj) = par_upp(jj)
!INFINITE-LOOP-BUG!          if (par_try(jj).EQ.par(jj)) goto 10
!Cid06Dec05-FIX! Update N_Retry & only retry the step if N_Retry <= N_sim!
!     (Otherwise we may get into an infinite loop...)
               if (par_try(jj).EQ.par(jj)) then
                  N_Retry = N_Retry + 1
                  if (N_Retry.LE.N_sim) goto 10
               end if


c     Perturbed x(jj) or fn ==> so renormalize x() to sum = fn.
c     Also, AV & YAV are the same ==> so set i_IsAVstillTheSame to "YES"!
               if ((jj.LE.N_base).OR.(jj.EQ.j_fn)) then

                  i_IsAVstillTheSame = 1
                  fn_try = par_try(j_fn)
                  sum_x = 0.
                  do j=1,N_base
                     sum_x = sum_x + par_try(j)
                  end do
                  do j=1,N_base
                     par_try(j) = (par_try(j) / sum_x) * fn_try
                  end do

c                 fffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!OPTfn! If chosen by the user (in arq_config), we now rescale the pop vector (and thus fn) such 
!       that it produces the best possible chi2. Optimal_fn_OPT comes via common from the 
!       VeryFastF_chi2 function (called from within MEX_Fpar_VFchi2)
                  if (IsOptimize_fn_OPT.EQ.1) then
                     zilt = MEX_Fpar_VFchi2(N_par,N_base, par_try,1)
                     do j=1,N_base
                        par_try(j) = par_try(j) * Optimal_fn_OPT
                     end do
                     par_try(j_fn) = fn_try * Optimal_fn_OPT

!OPTfn! Screen test
!                     zolt = MEX_Fpar_VFchi2(N_par,N_base, par_try,1)
!                     write (*,'(a,5(a,1x,f10.4,1x))') 
!     &                    '[MEX_ScubiduCadeVoce] OPTfn: ' , 
!     &                    'fn changed from',fn_try,'to',par_try(j_fn),
!     &                    'and energy decreased(?) from', 0.5*zilt , 
!     &                    'to',0.5*zolt,'=> OLD/NEW energy=',zilt/zolt
                  end if
c                 fffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

c     Perturbed AV or exAV ==> so set i_IsAVstillTheSame to "NO". Also,
c     update counter of AV/exAV steps.
               else if ((jj.GE.j_AV).AND.(jj.LT.j_fn)) then

                  i_IsAVstillTheSame = 0
                  N_AVYAV_steps = N_AVYAV_steps + 1


c     Safety precaution (shouldn't really get here anyway)
               else

                  print *, ' ??? What am I doing here??? ' , 
     $                 i_IsAVstillTheSame , jj , par(jj) , par_try(jj) 
                  par_try(jj) = par(jj)
                  goto 10

               end if
c     .......................................................................

            else

c     ........ Modifies par --> par_try: Shakes ALL par_j's per move!........
c     Perturb/move par => par_try (limiting to low <=> upp boundaries).
c     Then renormalize population vector x_try = par_try(1...N_base) to
c     fn_try = par_try(j_fn). Also turn-off i_IsAVstillTheSame if AV/YAV
c     changed (which's very likely in this scheme).
               do j=1,N_par
                  par_try(j) = par(j) + Chain_alpha(i_chain) * eps(j)
     $                 * gasdev(idum)
                  if (par_try(j).LT.par_low(j)) par_try(j) = par_low(j)
                  if (par_try(j).GT.par_upp(j)) par_try(j) = par_upp(j)
               end do

               fn_try = par_try(j_fn)
               sum_x = 0.
               do j=1,N_base
                  sum_x = sum_x + par_try(j)
               end do
               do j=1,N_base
                  par_try(j) = (par_try(j) / sum_x) * fn_try
               end do

               i_IsAVstillTheSame = 1
               do j=j_AV,j_fn-1
                  if (par_try(j).NE.par(j)) i_IsAVstillTheSame = 0
               end do

c              ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!OPTfn! As above, but not really tested since I never used this shake-it-all-pars option anyway!
               if (IsOptimize_fn_OPT.EQ.1) then
                  zilt = MEX_Fpar_VFchi2(N_par,N_base,par_try,
     &                 i_IsAVstillTheSame)
                  do j=1,N_base
                     par_try(j) = par_try(j) * Optimal_fn_OPT
                  end do
                  par_try(j_fn) = fn_try * Optimal_fn_OPT
               end if
c              ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

c     .......................................................................

            
            end if


c     Compute energy at par_try. Then accept new position or not
c     according to Metropolis rule. Also, verify if par_try is a new
c     minimum for this chain.
            Ene_try = 0.5 * MEX_Fpar_VFchi2(N_par,N_base,par_try,
     $           i_IsAVstillTheSame)
            prob    = 1. - exp(-(Ene_try - Ene) / Temperature)

            if ((Ene_try.LE.Ene).OR.(ran2(idum).GT.prob)) then

               do j=1,N_par
                  par(j) = par_try(j)
               end do
               Ene = Ene_try
               if (Ene.LT.CEne_min(i_chain)) then
                  do j=1,N_par
                     Cpar_min(j,i_chain) = par(j)
                  end do
                  CEne_min(i_chain) = Ene
                  N_min = N_min + 1
               end if

            else

c     Step was rejected. Since we are using the fast skip-over-AV-loop
c     trick, the saved-summations in the fast-chi2-routine must be
c     restored to correct values if a change in AV (or YAV) was
c     rejected, as is the case in the "if" below. A simple call to the
c     routine using the previous values will fix it. In this case, we
c     also update the counter of rejected AV/YAV steps.
               N_skip = N_skip + 1
               if ((i_MoveOneParOnly.EQ.1).AND.
     $              (i_IsAVstillTheSame.EQ.0))
     $              rubish = MEX_Fpar_VFchi2(N_par,N_base,par,0)
               if ((jj.GE.j_AV).AND.(jj.LT.j_fn)) N_AVYAV_skip =
     $              N_AVYAV_skip + 1

            end if


c     Update par & par**2 summations for this chain (used in GR-test).
            do j=1,N_par
               S_par(j,i_chain)  = S_par(j,i_chain)  + par(j)
               S_par2(j,i_chain) = S_par2(j,i_chain) + par(j) * par(j)
            end do

!2BEDONE! AQUI teria que passar todos os chi2s pra salvar e ver a influencia deles nos fits!!
!2BEDONE!!tvz a linha que diz rubish = MEX_Fpar_VFchi2(N_par,N_base,par,0) atrapalhe... se houverem commmons...

c     Save Step if i_UnitScubi > 0
            if (i_UnitScubi.GT.0) call SaveStep(i_sim,i_chain,N_par,par
     $           ,Ene)


!!MEX!! Just testing AV & exAV steps...
!            if ((jj.GE.j_AV).AND.(jj.LT.j_fn)) write (*
!     $           ,'(4(i5,2x),20(f6.3,2x))') i_chain , i_sim , jj ,N_skip
!     $           ,(par(j),j=1,N_par)

c     Goto next step i_sim for chain # i_chain
         end do
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!     Screen check of running-AV/YAV-efficiency...
!         AVYAV_eff = 1. - float(N_AVYAV_skip) / float(N_AVYAV_steps)
!         write (*,'(a,3(1x,i5),3(a,f10.5))') ' AVYAVeff ' , i_chain ,
!     $        N_AVYAV_steps , N_AVYAV_skip ,  '  eff%=' , 100.
!     $        * AVYAV_eff , ' eps_AV=', eps(j_AV) , '  last_AV='
!     $        ,par(j_AV)


c     Compute chain-efficiency (latest N_sim steps only) & update global
c     efficiency summation (eff_Scubi).
         Chain_eff(i_chain) = 1. - float(N_skip) / float(N_sim)
         eff_Scubi          = eff_Scubi + N_skip

c     Store final chain parameters in Cpar_last (used to restart it the
c     next time), and compute Cpar_ave = mean & Cpar_var = variance of
c     par(j), used in the GR test. Note the abs in the variance to avoid
c     precision bugs...
         do j=1,N_par
            Cpar_last(j,i_chain) = par(j)
            Cpar_ave(j,i_chain)  = S_par(j,i_chain) / float(N_steps)
            Cpar_var(j,i_chain)  = abs(S_par2(j,i_chain) -
     $           float(N_steps) * Cpar_ave(j,i_chain) * Cpar_ave(j
     $           ,i_chain)) / float(N_steps - 1)
         end do


c     Goto next chain
      end do
c     -----------------------------------------------------------------------


c     -------------- Re-evaluate CEne_min's if NRC_AV changed ---------------
c     If NRC_AV changed since last loops through the chains it's a good
c     idea to recompute all minimum chain-energies. Also, reset
c     NRC_AV_last to current value.
      if (NRC_AV.NE.NRC_AV_last) then 
         do i_chain=1,N_chains
            do j=1,N_par
               aux1(j) = Cpar_min(j,i_chain)
            end do
            CEne_min(i_chain) = 0.5*MEX_Fpar_VFchi2(N_par,N_base,aux1,0)
         end do
         NRC_AV_last = NRC_AV
      end if
c     -----------------------------------------------------------------------


c     -- Brief break to search 4 a better Global minimum & define Gpar_ave --
c     Look up which chain has the smallest minimum & store it in the
c     Global minimum arrays if it's better than GE_min.
      do i_chain=1,N_chains
         if (CEne_min(i_chain).LT.GE_min) then
            do j=1,N_par
               Gpar_min(j) = Cpar_min(j,i_chain)
            end do
            GE_min = CEne_min(i_chain)
         end if
      end do


c     Verify if one of the N_chains-Cpar_ave's is a better Global-minimum...
      do i_chain=1,N_chains
         do j=1,N_par
            aux1(j) = Cpar_ave(j,i_chain)
         end do
         Ene_ave = 0.5 * MEX_Fpar_VFchi2(N_par,N_base,aux1,0)
         if (Ene_ave.LT.GE_min) then 
            do j=1,N_par
               Gpar_min(j) = aux1(j)
            end do
            GE_min = Ene_ave
         end if
      end do


c     Define global-par-average & verify if it is a better Global-minimum...
      do j=1,N_par
         Gpar_ave(j) = 0.
         do i_chain=1,N_chains
            Gpar_ave(j) = Gpar_ave(j) + Cpar_ave(j,i_chain)
         end do
         Gpar_ave(j) = Gpar_ave(j) / float(N_chains)
      end do
      GEne_ave = 0.5 * MEX_Fpar_VFchi2(N_par,N_base,Gpar_ave,0)
      if (GEne_ave.LT.GE_min) then 
         do j=1,N_par
            Gpar_min(j) = Gpar_ave(j)
         end do
         GE_min = GEne_ave
      end if
c     -----------------------------------------------------------------------



c     --------------- Update alpha & eps(AV/YAV) (if wanted) ----------------
c     Update alpha based on efficiency using an exponential recipe.
c     Updates can be made for all chains at once (i_UpdateAlpha = 1), in
c     which case the mean efficiency is used.  Using i_UpdateAlpha = 2
c     will update alpha for each chain using its own efficiency. ATT:
c     Alpha-updates are only done if the resulting step size is not too
c     big (> 20% of low==>upp range)!!
      if (i_UpdateAlpha.EQ.1) then

         call avevar(Chain_eff,N_chains,ave_eff,var_eff)
         alp_factor = Falpha**((ave_eff - eff_IDEAL) / eff_IDEAL)
         if (ave_eff.GE.eff_IDEAL) alp_factor = Falpha**((ave_eff -
     $        eff_IDEAL) / (1. - eff_IDEAL))
         test_alpha = Chain_alpha(1) * alp_factor

         do j=1,N_par
            z1 = test_alpha * eps(j)
            z2 = (par_upp(j) - par_low(j)) / 5.
            if (z1.GT.z2) goto 50
         end do

         do i_chain=1,N_chains
            Chain_alpha(i_chain) = test_alpha
         end do

      else if (i_UpdateAlpha.EQ.2) then

         do i_chain=1,N_chains
            alp_factor = Falpha**((Chain_eff(i_chain) - eff_IDEAL) /
     $           eff_IDEAL)
            if (Chain_eff(i_chain).GE.eff_IDEAL) alp_factor = 
     $           Falpha**((Chain_eff(i_chain) - eff_IDEAL) / (1. -
     $           eff_IDEAL))
            test_alpha = Chain_alpha(i_chain) * alp_factor

            do j=1,N_par
               z1 = test_alpha * eps(j)
               z2 = (par_upp(j) - par_low(j)) / 5.
               if (z1.GT.z2) goto 40
            end do
            Chain_alpha(i_chain) = test_alpha

 40         DoNothingOnThisSillyLine = 0.
         end do

      end if
 50   DoNothingOnThisSillyLine = 0.


c     Adjust the step-size for AV & YAV according to how efficient were
c     the moves in these 2 variables in the past N_sim steps of all
c     chains. This is equivalent to the alpha-update above, but here we
c     operate on eps itself instead of on alpha. Again, precautions are
c     taken against too large step-sizes. Also, we consider that the
c     effective step size depends on alpha (updated above).
      if (i_UpdateAVYAVStepSeparately.EQ.1) then
         AVYAV_eff = 1. - float(N_AVYAV_skip) / float(N_AVYAV_steps)

         eps_factor = Falpha**((AVYAV_eff - eff_IDEAL) / eff_IDEAL)
         if (AVYAV_eff.GE.eff_IDEAL) eps_factor = Falpha**((AVYAV_eff -
     $        eff_IDEAL) / (1. - eff_IDEAL))

         alp_ave = Chain_alpha(1)
         if (i_UpdateAlpha.EQ.2) call avevar(Chain_alpha,N_chains
     $        ,alp_ave,alp_var)

         z1 = alp_ave * eps_factor * eps(j_AV)
         z2 = (par_upp(j_AV) - par_low(j_AV)) / 5.
         if (z1.LT.z2) eps(j_AV) = eps_factor * eps(j_AV)

         do j=j_exAV,j_fn-1
            z1 = alp_ave * eps_factor * eps(j)
            z2 = (par_upp(j) - par_low(j)) / 5.
            if (z1.LT.z2) eps(j) = eps_factor * eps(j)
         end do

!!UEPA!!
!     Screen check
!        write (*,'(a,2(1x,i5),2(a,f10.6),a,9(2x,f10.6))')
!    $        '  *> AVYAV-step-update: ',N_AVYAV_steps , N_AVYAV_skip , 
!    $        '  eff%=' , 100.* AVYAV_eff, ' eps_AV=',eps(j_AV),
!    $        '  eps_exAV=',(eps(j),j=j_exAV,j_fn-1)
      end if
c     -----------------------------------------------------------------------


c     ----------------------- GR Convergence Test ---------------------------
c     Perform Gelman & Rubin R-test. Note that we EXCLUDE fn from the
c     test! This is justifiable because fn is already inside the
c     pop-vector... Since fn is the last entry in par(), we just pretend
c     there are N_par - 1 parameters. 
      call Calc_GR_R(N_steps,N_chains,Npar4GRtest,Cpar_ave,Cpar_var,
     $     R_Threshold,GR_R,GR_R_ave,IsAveBurnInOver,IsIndBurnInOver)


c     Verify if GR-R converged. If so, we are close enough to the end of
c     the burn-in to get out of the convergence loop and start a proper
c     mapping of the parameter space with other routines. Otherwise,
c     iterate again adding more points to the chains, until either R <
c     R_Threshold or we exceed our limit on N_steps... (<= Nmax_steps).
      i_aux = IsAveBurnInOver
      if (IsGRTestHard.EQ.1) i_aux = IsIndBurnInOver


c     Burn-In over:) print something & go home.
      if (i_aux.EQ.1) then 

         IsBurInOver = 1
         if (i_verbose.GE.1) write (*,'(a,i8,a,1p,e9.3)')
     $        '[Scubidu] Burn In over after ' ,N_steps ,
     $        ' steps :)    @ T = ' , Temperature
         goto 200

c     GR-convergence was NOT reached after Nmax_steps:( Issue a Warning
c     and give up!
      else if ((i_aux.EQ.0).AND.(N_steps.GT.Nmax_steps)) then

         if (i_verbose.EQ.1) write (*,'(a,i8,1p,2(a,e9.3),a)')
     $        '[Scubidu] ATT: NO CONVERGENCE after ' ,N_steps ,
     $        ' steps:(    @ Temperature = ' , Temperature, '  <R> ='
     $        ,GR_R_ave,' :(   ...but will keep going...'
         goto 200

      end if
c     -----------------------------------------------------------------------


c     ------------------- Burn In NOT over: Iterate More --------------------
c     Burn-In NOT over:(  Do a few things, print stuff & iterate more...
 300  alp_ave = Chain_alpha(1)
      alp_var = 0.
      if (i_UpdateAlpha.EQ.2) call avevar(Chain_alpha,N_chains,alp_ave
     $     ,alp_var)
      call avevar(Chain_eff,N_chains,ave_eff,var_eff)


c     Find out which parameter has the largest R. par(jR_max) will be
c     given an extra help in the next loop through the chains.
      R_max  = GR_R(1)
      jR_max = 1
      do j=2,Npar4GRtest
         if (GR_R(j).GT.R_max) then
            R_max  = GR_R(j)
            jR_max = j
         end if
      end do


c     Increase step size for parameters with large R_j (compared to <R>)
c     NEED TESTS!!! Maybe ABANDON this feature????
      if (i_UpdateEps.EQ.1) then
         do j=1,Npar4GRtest
            if ((GR_R(j)/GR_R_ave).GT.3.0) then
               test_eps = eps(j) * Feps
               z1 = (par_upp(j) - par_low(j)) / 5.
               if ((alp_ave * test_eps).LT.z1) then 
                  write (*,*) ' ????Re-scaled eps_j: ', j , eps(j)
     $                 ,test_eps , GR_R(j) , GR_R_ave
                  eps(j) = test_eps
               end if
            end if
         end do
      end if 


c     Print stuff...
      if (i_verbose.GT.1) then
         min_pRmax = Cpar_last(jR_max,1)
         max_pRmax = Cpar_last(jR_max,1)
         do i_chain=1,N_chains
            z1 = Cpar_last(jR_max,i_chain)
            if (z1.LT.min_pRmax) min_pRmax = z1
            if (z1.GT.max_pRmax) max_pRmax = z1
            aux1(i_chain) = z1
         end do
         call avevar(aux1,N_chains,ave_aux,var_aux)

         write (*,1) '[Scubidu] N_steps=',N_steps ,' 2Emin=',2*GE_min ,
     $        '  T=',Temperature , '  a=',alp_ave , '  a*eps_x='
     $        ,alp_ave*eps(1), '  a*eps_AV=',alp_ave*eps(j_AV),
     $        '  eff%=',100*ave_eff ,'  <R>=' ,GR_R_ave , ' <' ,R_max
     $        ,jR_max, ' d=', max_pRmax - min_pRmax ,  ' sig=',
     $        sqrt(var_aux)
 1       format (a,i9,1p,a,e10.4,4(a,e8.2),0p,3(a,f6.2),1x,i3,2(a,f7.4))
      end if


c     Iterate all chains another N_sim steps
      goto 100
c     -----------------------------------------------------------------------
c ***************************************************************************


c ***************************************************************************
c *                             FINALIZATION...                             *
c ***************************************************************************
c     Just check if <par> yields a better model than the last point for
c     each chain. If so, reset Cpar_last to Cpar_ave (& maybe Cpar_min
c     too) before returning. Take the oportunity to verify if a new
c     global minimum is found on the way.
 200  do i_chain=1,N_chains
         do j=1,N_par
            par(j)  = Cpar_last(j,i_chain)
            aux1(j) = Cpar_ave(j,i_chain)
         end do
         Ene_last = 0.5 * MEX_Fpar_VFchi2(N_par,N_base,par,0)
         Ene_ave  = 0.5 * MEX_Fpar_VFchi2(N_par,N_base,aux1,0)

         CEne_last(i_chain) = Ene_last

         if (Ene_ave.LT.Ene_last) then 
            do j=1,N_par
               Cpar_last(j,i_chain) = Cpar_ave(j,i_chain)
            end do
            CEne_last(i_chain) = Ene_ave
         end if

         if (Ene_ave.LT.CEne_min(i_chain)) then
            do j=1,N_par
               Cpar_min(j,i_chain) = Cpar_ave(j,i_chain)
            end do
            CEne_min(i_chain) = Ene_ave
         end if

         if (Ene_ave.LT.GE_min) then 
            do j=1,N_par
               Gpar_min(j) = Cpar_ave(j,i_chain)
            end do
            GE_min = Ene_ave
         end if
      end do


c     Store low-->upp range among last parameters of all chains in
c     Cpar_low & _upp, just for control...
      do j=1,N_par
         Cpar_low(j) = Cpar_last(j,1)
         Cpar_upp(j) = Cpar_last(j,1)
         do i_chain=2,N_chains
            z1 = Cpar_last(j,i_chain)
            if (z1.LT.Cpar_low(j)) Cpar_low(j) = z1
            if (z1.GT.Cpar_upp(j)) Cpar_upp(j) = z1
         end do         
      end do


c     Define mean alpha (to be returned as "alpha")
      alp_ave = Chain_alpha(1)
      alp_var = 0.
      if (i_UpdateAlpha.EQ.2) call avevar(Chain_alpha,N_chains
     $     ,alp_ave,alp_var)
      alpha = alp_ave


c     Define overall efficiency (to be returned)
      eff_Scubi = 1. - eff_Scubi / float(N_steps * N_chains)


!     Screen output <chi2>, <eff>, <alpha>, ... chi2(<p>)...Gpar_ave &
!     p_last-stats...
      if (i_verbose.GT.1) then
         call avevar(Chain_eff,N_chains,ave_eff,var_eff)
         call avevar(CEne_last,N_chains,ave_Ene,var_Ene)
         call avevar(eps,N_par,ave_eps,var_eps)

         do j=1,N_par
            aux1(j) = 0.
            do i_chain=1,N_chains
               aux1(j) = aux1(j) + Cpar_last(j,i_chain)
            end do
            aux1(j) = aux1(j) / float(N_chains)
         end do
         chi2_ave_last = MEX_Fpar_VFchi2(N_par,N_base,aux1,0)

         write (*,'(1p,3(a,e12.5))')
     $        '[Scubidu] Leaving routine with <chi2> = ',2*ave_Ene,
     $        ' +/- ' , 2*sqrt(var_Ene) , '   |  chi2(<p_last>) = '
     $        , chi2_ave_last
         zlix = MEX_Fpar_Slowchi2(N_par,N_base,Gpar_min)
         write (*,'(28x,1p,2(a,e12.5),a)') 'Gchi2_min = ' , 2. * GE_min
     $        , ' (FAST) &  ' , zlix , ' (SLOW)'
         write (*,'(30x,2(a,f12.5))') '<eff %> = ' , 100* ave_eff ,
     $        ' +/- ', 100 * sqrt(var_eff)
         write (*,'(30x,1p,2(a,e12.5))') '<alpha> = ' , alp_ave,
     $        ' +/- ' ,sqrt(alp_var)
         write (*,'(30x,1p,4(a,e12.5))') '<eps_j> = ' , ave_eps, ' +/- '
     $        , sqrt(var_eps) , '   |  <a> * <eps_j>  = ',alp_ave
     $        * ave_eps , '   |  <a> * eps_AV  = ',alp_ave * eps(j_AV)
         write (*,'(34x,2(a,f12.5),a,i1,a)') '<R> = ' , GR_R_ave ,
     $        '  | R_Threshold = ',R_Threshold , '  (',IsGRTestHard,')'
         write (*,'(30x,a,i8,a,i8)') 'N_steps = ' , N_steps ,
     $        '   <= Nmax_steps = ' , Nmax_steps

         write (*,'(28x,a,1p,e12.5)') 'eff_min % = ' , 100.
     $        * float(N_min) / float(N_steps * N_chains)

         write (*,*)
         write (*,'(a,a)') '[Scubidu] Gpar_ave & par_last(j) ' ,
     $        'chain-stats after Scubi-run'
         do j=1,N_par

c     Define std-dev among the Cpar_ave's for each parameter. Also
c     define std-dev & etc among the Cpar_last's for each parameter
            do i_chain=1,N_chains
               aux1(i_chain) = Cpar_ave(j,i_chain) - Gpar_ave(j)
               aux2(i_chain) = Cpar_last(j,i_chain)
            end do

            call avevar(aux1,N_chains,ave_diff,var_diff)
            Gpar_sig(j) = sqrt(abs(var_diff))
            if (var_diff.LT.0.) Gpar_sig(j) = - Gpar_sig(j)

            call avevar(aux2,N_chains,ave_p,var_p)
            sig_p  = sqrt(var_p)
            min_p  = Cpar_low(j)
            max_p  = Cpar_upp(j)
            p_best = Gpar_min(j)

            z1 = 1.
            if (j.LE.N_base) z1 = 100.
            write (*,'(a,i3,7(a,f6.2),a,f7.3,a,i3,a,f7.3)')
     $           '[Scubidu] j =' ,j , 
     $           '  <Gp>= ',z1*Gpar_ave(j),' +/- ',z1 *Gpar_sig(j) ,
     $           ' | <p_last> = ',z1*ave_p,' +/- ',z1*sig_p , 
     $           '  |  Range=',z1*min_p,' --> ',z1*max_p ,
     $           ' |  Delta= ',z1*max_p - z1*min_p ,
     $           '  |  Gp_min= ',z1*p_best , ' ~ ',nint(z1*p_best) , 
     $           ' | R_j = ' , GR_R(j)
         end do

         write (*,*)
         write (*,'(a)') '[Scubidu] Chains Minima after Scubi-run'
         write (*,'(a,1p,40(1x,e8.2))') '[Scubidu] Chain-chi2_min = '
     $        , (2.* CEne_min(i_chain),i_chain=1,N_chains)
         do j=1,N_par
            z1 = 1.
            if (j.LE.N_base) z1 = 100.
            write (*,'(a,i3,a,40(2x,f7.3))') '[Scubidu] j =' ,j ,
     $           '  |  p_j = ', (z1*Cpar_min(j,i_chain),i_chain=1
     $           ,N_chains)
         end do


c     Prints best chain energy & compare it to global minimum...
         z_min = CEne_min(1)
         i_min = 1
         do i_chain=2,N_chains
            if (CEne_min(i_chain).LT.z_min) then
               z_min = CEne_min(i_chain)
               i_min = i_chain
            end if
         end do
         write (*,'(a,i3,1p,3(a,e8.2))') '[Scubidu] Best Chain = ' ,
     $        i_min,'  2*CEne_min = ', 2.*z_min , ' | 2*GE_min = ',2.
     $        *GE_min,'  | Diff = ' , 2.*(z_min - GE_min)

         write (*,*)
      end if
c ***************************************************************************



      end
c ###########################################################################






c ###########################################################################
c     Unpack par-vector to (x,fn,AV,exAV) & compute chi2 using
c     SlowF_chi2. Assumes x is already normalized to sum x_j = fn!

      FUNCTION MEX_Fpar_Slowchi2(N_par,N_base,par)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nb_max = 300,Npar_max = 601)
      real par(Npar_max) , x(Nb_max) , AV_tot(Nb_max) ,exAV(Nb_max)

      call MEX_Unpack_par2xAVYAVfn(N_par,N_base,par,x,AV,exAV,fn)
      call MEX_Update_AV_array(N_par,N_base,AV,exAV,AV_tot)
      MEX_Fpar_Slowchi2 = SlowF_chi2(x,AV_tot)

      return
      end
c ###########################################################################



c ###########################################################################
c     Unpack par-vector to (x,fn,AV,exAV) & compute chi2 using
c     VeryFastF_chi2 Assumes x is already normalized to sum x_j = fn!
c     This version use the i_IsAVstillTheSame-trick!

      FUNCTION MEX_Fpar_VFchi2(N_par,N_base,par,i_IsAVstillTheSame)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nb_max = 300,Npar_max = 601)
      real par(Npar_max) , x(Nb_max) , AV_tot(Nb_max) , exAV(Nb_max)

      call MEX_Unpack_par2xAVYAVfn(N_par,N_base,par,x,AV,exAV,fn)
      call MEX_Update_AV_array(N_par,N_base,AV,exAV,AV_tot)
      MEX_Fpar_VFchi2 = VeryFastF_chi2(x,AV_tot,i_IsAVstillTheSame)

      return
      end
c ###########################################################################


c ###########################################################################
c     Unpack par-vector to x(), AV, exAV() & fn
c     Cid@Lagoa - 19/April/2007

      SUBROUTINE MEX_Unpack_par2xAVYAVfn(N_par,N_base,par,x,AV,exAV,fn)

      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nb_max = 300,Npar_max = 601)
      real par(Npar_max) , x(Nb_max) , exAV(Nb_max)

*FIR* Reset outputs for safety/caution/paranoia! - Cid@Lagoa 16/Jan/2010
      AV = 0.
      fn = 0.
      do j=1,Nb_max
         x(j)    = 0.
         exAV(j) = 0.
      end do

      if ((N_par.LE.0).OR.(N_base.LE.0)) then
         write (*,*)
     $        '[MEX_Unpack_par2xAVYAVfn] N_base <= 0!! N_base = '
     $        ,N_base , '  N_par = ' , N_par
         stop
      end if

      j_AV   = N_base + 1
      j_exAV = N_base + 2
      j_fn   = N_par

      do j=1,N_base
         x(j) = par(j)
      end do
      AV = par(j_AV)
      do j=j_exAV,j_fn-1        !OBS: This loop only runs if Nex_AV > 0!
         jMEX       = j - j_exAV + 1
         exAV(jMEX) = par(j)
      end do
      fn = par(j_fn)

      end
c ###########################################################################


c ###########################################################################
c     Simple routine just to store AV & YAV onto an AV_tot array. The
c     array ind_exAV(), defined elsewhere, says which components have
c     which of the exAV(1...N_exAV) N_exAV extra extinctions allowed.
c     Routines which compute a synthetic spectrum and the chi2-routines
c     combine AV & exAV to build the population dependent total
c     extinction AV_tot(j)

c     Cid@Lagoa - 21/April/2007

!!MEX!! VERY IMPORTANT! I have introduced a CUMULATIVE AV-scheme!
!     For instance, if a component has ind_exAV = 2, its AV will be
!     the global one PLUS exAV(1) AND exAV(2)..
!     Cid@Lagoa - 27/May/2007

      SUBROUTINE MEX_Update_AV_array(N_par,N_base,AV,exAV,AV_tot)
 
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nb_max = 300)
      real AV_tot(Nb_max) , exAV(Nb_max)
      common/exAV_stuff/ N_exAV , ind_exAV(Nb_max)

*FIR* Reset outputs for safety/caution/paranoia! - Cid@Lagoa 16/Jan/2010
      do j=1,Nb_max
         AV_tot(j) = 0.
      end do

!!MEX!! Cumulative AV_tot...
      do j=1,N_base
         AV_tot(j) = AV
         do jMEX = 1,ind_exAV(j)
            AV_tot(j) = AV_tot(j) + exAV(jMEX)
         end do
      end do

      end
c ###########################################################################




c ###########################################################################
c # !!MEX!!     CONDENSE BASE (excluding irrelevant x_j) & RUN Fit_Scubi    #
c ###########################################################################
c     This routine reduces/condenses the spectral base (Of_base &
c     f_base, used in all routines), excluding all components in
c     x_Test() which fail and EX0s-criterion which depends on the
c     EX0s_method_option-string. In the 'CUMUL' method we exclude the
c     weaker components which add up to < EX0s_Threshold. In the 'SMALL'
c     method we throw away any component which contributes with less
c     flux than EX0s_Threshold / N_base. This reduction is limited to a
c     smaller base dimension of N_base_min. x_Test may be, for instance,
c     the best, or the global-average population vector.
c
c     All "Irrelevant" components are zero'ed. We then call Fit_Scubidu
c     to obtain new fits.  Before exiting, the base matrices & other
c     arrays are restored to their original format, otherwise things
c     would become a total mess!  
c
c     Cid@Lagoa - 27/Mar/2005

      SUBROUTINE MEX_EX0s_Fit_Scubi(i_verbose,idum,N_loops,v0_low,v0_upp
     $     ,vd_low,vd_upp,x_min,AV_min,exAV_min,v0_min,vd_min,Ntot_EX0s
     $     ,eff_EX0s,Nl_eff,iMassAgeZ_verbose,N_par,ChainPar,N_chains
     $     ,Temp_ini,Temp_fin,fini_eps,ffin_eps,R_ini,R_fin,N_sim
     $     ,Nmax_steps,i_RestartChains,Gpar_min,GE_min,Gpar_ave,Cpar_min
     $     ,CEne_min,x_Test,EXOs_method_option,EX0s_Threshold,N_base_min
     $     ,IsKinFixed_EX0s,IsBurInOver_EX0s,N_base_Relevant
     $     ,IsScaleNstepsInEX0sFits,IsNoKin4LargeBaseInEX0sFits
     $     ,frac_NoKin4LargeBaseInEX0sFits)

c ***************************************************************************
c *                               DEFINITIONS                               *
c ***************************************************************************
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
      parameter (Nchains_max = 1000,Npar_max = 601)
!IIB! Changed Nl_max to NAVgrid_max!
      parameter (NAVgrid_max = 1001)
*PHO*
      parameter (NPHO_Ys_max=100)

c     Observed spectra & etc. (Not really needed here, but I'm using the
c     Obs-common to read Nl_obs)
      real l_obs(Nl_max) , f_obs(Nl_max) , weight_obs(Nl_max)

c     Population vector & base arrays
      real x_Test(Nb_max) , x_min(Nb_max)
      real f_base(Nb_max,Nl_max) , Of_base(Nb_max,Nl_max)
      real FULL_Of_base(Nb_max,Nl_max)

c     Auxiliar pop vector & its sorting-index-array
      real x_aux(Nb_max) 
      integer ind_x_aux(Npar_max) 

c     Auxiliar array which map components with x(j) > x_RelevantThreshold!
      integer ind_Relevant(Npar_max) 

c     Original Full-Base parameter limits
      real Opar_low(Npar_max) , Opar_upp(Npar_max)

c     Chain-stuff-arrays!
      real ChainPar(Npar_max,Nchains_max)
      real Gpar_min(Npar_max) , Gpar_ave(Npar_max)
      real Cpar_min(Npar_max,Nchains_max)  , CEne_min(Nchains_max)

c     Original Full-Base HelpParMove2Average arrays
      integer iO_ind_HelpParMove2Average(Npar_max)
      real Oprob_HelpParMove2Average(Npar_max)

c     Short-arrays!
      real Shortx_min(Nb_max) , ShortChainPar(Npar_max,Nchains_max)
      real ShortGpar_min(Npar_max) , ShortGpar_ave(Npar_max)
      real ShortCpar_min(Npar_max,Nchains_max)
      real ShortCEne_min(Nchains_max)

c     Option-String which defines which method will be used to define
c     who is relevant and who is not in the pop vector x_Test!
      character*5 EXOs_method_option

c     commons...
      common/Obs/ Nl_obs , l_obs , f_obs , weight_obs
      common/Base/ N_base , f_base
      common/Obase/ Of_base

c     Lower & upper limits on parameters (x,AV,YAV,fn)
      common/ParLimits/ par_low(Npar_max) , par_upp(Npar_max)

c     ScubiduCadeVoce Technical Parameters - All we actually need here
c     are the ind_HelpParMove2Average & prob_HelpParMove2Average arrays!
      common/Scubi_TechnicalParameters/ alpha , Falpha , IsGRTestHard ,
     $     eff_IDEAL , i_UpdateEps , Feps , i_UpdateAlpha ,
     $     i_MoveOneParOnly , i_HelpParWithMaxR , prob_jRmax ,
     $     Npar4GRtest ,i_UpdateAVYAVStepSeparately ,
     $     ind_HelpParMove2Average(Npar_max) ,
     $     prob_HelpParMove2Average(Npar_max)

*     MEX-arrays
      real exAV_min(Nb_max)
      common/exAV_stuff/ N_exAV , ind_exAV(Nb_max)

*     Original/FULL ind_exAV_array
      integer IFULL_ind_exAV(Nb_max)

*FIR* FIR-Base-common + Original FIRRmat matrix & FIRBolCor
*     ATT1: ONLY FIRRmat & FIRBoldCor are condensed & then restored!
*     ATT2: We condense-&-restore these FIR-things even for
*     non-FIR-fits, in which case these arrays may be empty/undefined!!
!IIB! Changed Nl_max to NAVgrid_max!
      common/FIR_BaseStuff/ NFIR_AV , FIRAV_grid(NAVgrid_max) ,
     $     FIRAV_low , FIRAV_upp , FIRdAV , 
     $     FIRLbol(Nb_max) , FIRFracLion(Nb_max) , FIRBolCor(Nb_max) ,
     $     FIRR_opt(Nb_max,NAVgrid_max) , FIRR_Lya(Nb_max,NAVgrid_max) , 
     $     FIRR_LCE(Nb_max,NAVgrid_max) ,
     $     FIRRmat(Nb_max,NAVgrid_max) , FIRq_norm , NFIR_base

      real OFIRRmat(Nb_max,NAVgrid_max) , OFIRBolCor(Nb_max)

*QHR* QHR-Base-common + Original OQHR_Q2Lnorm40 array.
*     ATT1: ONLY OQHR_Q2Lnorm40 is condensed & then restored, as it is
*          the only one actually used in function QHR_CalcQHRchi2.
*     ATT2: As above, this is done even for non-QHR-fits!
      common/QHR_BaseStuff/ QHR_qH40(Nb_max) , QHR_Q2Lnorm40(Nb_max) ,
     $     QHR_q_norm , QHRbeta_I , NQHR_base

      real OQHR_Q2Lnorm40(Nb_max)

*PHO* Base-info-common. 
!IIB! Changed Nl_max to NAVgrid_max!
      common/PHO_BaseStuff/ NPHO_AV , PHOAV_grid(NAVgrid_max) ,
     $     PHOAV_low , PHOAV_upp , PHOdAV , 
     $     PHO_MeanLambda(NPHO_Ys_max) , PHO_StdDevLambda(NPHO_Ys_max) ,
     $     PHO_q_MeanLambda(NPHO_Ys_max) ,
     $     PHO_3Mat(Nb_max,NPHO_Ys_max,NAVgrid_max) ,
     $     PHO_q_norm , NPHO_base

      real OPHO_3Mat(Nb_max,NPHO_Ys_max,NAVgrid_max)

*PHO* OBS: This routine does not have access to NPHO_Ys, and yet need
*     to set/reset PHO_3Mat & OPHO_3Mat in the following lines:     
*              PHO_3mat(j,iY,i_AV) = OPHO_3mat(j_in_full_base,iY,i_AV)
*              PHO_3mat(j,iY,i_AV) = OPHO_3mat(j,iY,i_AV)
*              OPHO_3Mat(j,iY,i_AV) = PHO_3Mat(j,iY,i_AV)

*     For simplicity, I just use NPHO_Ys_max in the iY set/reset do-loops,
*     which is fine.
c ***************************************************************************


c ***************************************************************************
c *                  INITIALIZATION: STORE ORIGINAL STUFF                   *
c ***************************************************************************
c    Store original N_base, N_par & AV/exAV/fn indices.
      NO_base = N_base
      NO_par  = N_base + 1 + N_exAV + 1
      jO_AV   = N_base + 1
      jO_exAV = N_base + 2
      jO_fn   = NO_par
      
c     Store original full base Of_base ==> FULL_Of_base, since Of_base
c     will be modified below!  Do the same thing for the ind_exAV array.
      do j=1,NO_base
         do i=1,Nl_obs
            FULL_Of_base(j,i) = Of_base(j,i)
         end do
         IFULL_ind_exAV(j) = ind_exAV(j)
      end do

c     Store original par-limits & arrays which say whether to help
c     parameters to converge and with what probability.
      do j=1,NO_par
         Opar_low(j) = par_low(j)
         Opar_upp(j) = par_upp(j)

         iO_ind_HelpParMove2Average(j) = ind_HelpParMove2Average(j)
         Oprob_HelpParMove2Average(j)  = prob_HelpParMove2Average(j)
      end do

*FIR* Store original FIRRmat matrix & FIRBolCor array
*QHR* Store original QHR_Q2Lnorm40 array
*PHO* Store original PHO_3Mat 3D array
*     ATT: using NPHO_Ys_max since NPHO_Ys is nor available to this
*     routine!
      do j=1,NO_base
         do i_AV=1,NFIR_AV
            OFIRRmat(j,i_AV) = FIRRmat(j,i_AV)
         end do
!IIB! This is where the Illegal Instruction bug appeared...
c Rafa Faltaba un bucle para NPHO_AV. Solo estaba para NFIR_AV
c pero si no ajustamos FIR entonces NFIR_AV = 0
c Aparece dos veces mas en MEX_EX0s_Fit_Scubi, ya esta cambiado :)

         do i_AV=1,NPHO_AV
            do iY=1,NPHO_Ys_max
               OPHO_3Mat(j,iY,i_AV) = PHO_3Mat(j,iY,i_AV)
            end do
         end do
         OFIRBolCor(j) = FIRBolCor(j)
         OQHR_Q2Lnorm40(j) = QHR_Q2Lnorm40(j)
      end do


c ***************************************************************************


c ***************************************************************************
c *            CONSTRUCT INDEX-ARRAY OF "RELEVANT" COMPONENTS               *
c ***************************************************************************
c     Zero-in ind_Relevant() array for precaution
      do j=1,Npar_max
         ind_Relevant(j) = 0
      end do
      N_base_Relevant = 0
      sum_x_Relevant  = 0

c     Figure out the normalization of x_Test
      sum_x_Test = 0
      do j=1,NO_base
         sum_x_Test = sum_x_Test + x_Test(j)
      end do

c     Now create x_aux, the normalized version of x_Test, and indexx it
c     in ind_x_aux, such that x_aux(ind_x_aux) is sorted in a increasing
c     order (stronger components last).
      do j=1,N_base
         x_aux(j) = x_Test(j) / sum_x_Test
      end do
      call indexx(N_base,x_aux,ind_x_aux)

!     Screen output
      if (i_verbose.GT.0) then 
         sum_x = 0
         do j=1,NO_base
            j_aux = ind_x_aux(j)
            sum_x = sum_x + x_aux(j_aux)
            write (*,'(2(a,i3),2(a,f7.3),a)') '[EX0s_Fit_Scubi] ', j,
     $           ' : x_Test(',j_aux,' ) = ' ,100.*x_aux(j_aux),
     $           '%   ==>  Cumulative Sum = ' ,100*sum_x ,'%'
         end do
      end if


c     CUMULATIVE-EX0s-METHOD: In this case EX0s_Threshold means the
c     fraction of light in x_Test which will be thrown away!
      if (EXOs_method_option.EQ.'CUMUL') then

c     Construct array ind_Relevant() which contains the indices of the
c     N_base_Relevant components in x_Test(j) which add-upp to >= 1 -
c     EXOs_Threshold. (Actually, the loop below works the other way
c     around: it discards from ind_Relevant() the smallest-x elements,
c     which add up to a fraction <= EX0s_Threshold.)
         sum_x  = 0
         i_flag = 0
         do j=1,NO_base
            j_aux = ind_x_aux(j)
            sum_x = sum_x + x_aux(j_aux)
            if ((i_flag.EQ.0).AND.(sum_x.GT.EX0s_Threshold)) i_flag = 1
            if (i_flag.EQ.1) then
               sum_x_Relevant  = sum_x_Relevant + x_aux(j_aux)
               N_base_Relevant = N_base_Relevant + 1
               ind_Relevant(N_base_Relevant) = j_aux
               if (i_verbose.GT.1) then 
                  write (*,'(a,4(i4),3(f10.3))')
     $                 '[EX0s_Fit_Scubi] Relevant comps.: ', j ,j_aux,
     $                 N_base_Relevant ,ind_Relevant(N_base_Relevant),
     $                 100 *x_aux(j_aux),100*sum_x,100*sum_x_Relevant
               end if
            end if
         end do

c     SMALL-X-EX0s-METHOD: In this case EX0s_Threshold / N_base defines
c     the minimum value of x for a component to be kept in the
c     base. Weaker components will be thrown away!
      else if (EXOs_method_option.EQ.'SMALL') then

c     Set Threshold for elimination!
         x_RelevantThreshold = EX0s_Threshold / float(N_base)

c     Construct array ind_Relevant() which contains the indices of the
c     N_Relevant components in x_aux (= normalized x_Test) which are
c     stronger (ie, >) than x_RelevantThreshold.
         do j=1,NO_base
            if (x_aux(j).GT.x_RelevantThreshold) then
               N_base_Relevant = N_base_Relevant + 1
               ind_Relevant(N_base_Relevant) = j
               sum_x_Relevant = sum_x_Relevant + x_aux(j)
               if (i_verbose.GT.1) then 
                  write (*,'(a,3(i4),2(f10.3))')
     $                 '[EX0s_Fit_Scubi] Relevant comps.: ', j,
     $                 N_base_Relevant ,ind_Relevant(N_base_Relevant),
     $                 100 *x_aux(j),100*sum_x_Relevant
               end if
            end if
         end do

      end if


c     PROTECTION against too large reductions in the base! If less than
c     N_base_min of the original base elements survived the reduction
c     above then just pick the N_base_min stonger elements. (Use the
c     sorting/index-array for this.)
      if (N_base_Relevant.LT.N_base_min) then

         write (*,'(a,i3,a,i3,a)')
         write (*,'(a,i5,a,i5,a)')
     $        '[EX0s_Fit_Scubi] ATT: N_base_Relevant = ' ,
     $        N_base_Relevant , ' < N_base_min = ' ,N_base_min
     $        ,' elements in the reduced base! '
         write (*,'(a,f7.3,a)') '[EX0s_Fit_Scubi] They they add up to '
     $        ,100.*sum_x_Relevant , '% of the light.'
         write (*,'(a,i3,a)')
     $        '[EX0s_Fit_Scubi] Will fix this picking the ',N_base_min ,
     $        ' strongest ones!!'
         write (*,*)

         N_base_Relevant = 0
         sum_x_Relevant  = 0
         i_flag = 0
         do j=NO_base,(NO_base-N_base_min+1),-1
            j_aux = ind_x_aux(j)
            sum_x_Relevant  = sum_x_Relevant + x_aux(j_aux)
            N_base_Relevant = N_base_Relevant + 1
            ind_Relevant(N_base_Relevant) = j_aux
!            write (*,'(a,3(i4),2(f10.3))')
!     $           '[EX0s_Fit_Scubi] Relevant comps.: ', j ,j_aux
!     $           ,N_base_Relevant ,100 *x_aux(j_aux),100 *sum_x_Relevant
         end do

      end if


!     Screen output
      if (i_verbose.GT.0) then 
         write (*,'(a,i3,a,f8.3,a)') '[EX0s_Fit_Scubi] There are '
     $        , N_base_Relevant ,
     $        ' elements in the reduced base, and they add up to ' ,
     $        100.*sum_x_Relevant , '% of the light.'
         write (*,'(a,a)') '[EX0s_Fit_Scubi] EX0-method     = '
     $        ,EXOs_method_option
         write (*,'(a,f7.3)')
     $        '[EX0s_Fit_Scubi] EX0s_Threshold = ' , EX0s_Threshold 
      end if


c     Freeze kinematics for EX0-fits if condensed base is still "too
c     large": > frac_NoKin4LargeBaseInEX0sFits * N_base. OBS: In
c     StarlightChains_v01 we froze the kinematics for > 50
c     components. Now we set this with the *NoKin4LargeBaseInEX0sFits
c     config variables.
!Cid14Jan07-FXK! IsKinFixed_EX0s may already be 1 if doing FXK fits!
      N_aux = nint(frac_NoKin4LargeBaseInEX0sFits * N_base)
      if ((IsKinFixed_EX0s.NE.1).AND.(IsNoKin4LargeBaseInEX0sFits.EQ.1)
     $     .AND.(N_base_Relevant.GT.N_aux)) IsKinFixed_EX0s = 1


!     Screen output - to freeze or not to freeze kinematical fits...
      if (i_verbose.GT.0) then 
         if (IsKinFixed_EX0s.EQ.1) then
!     Note: If IsThisAnFXKFit = 1, the message below will be shown but
!     means nothing which we did not already know!
            write (*,'(a,a)') '[EX0s_Fit_Scubi] ATT: ' ,
     $           ' We will freeze in the kinematical parameters!'
         else
            write (*,'(a,a)') '[EX0s_Fit_Scubi] ATT: ',
     $           ' We will NOT freeze in the kinematical parameters!'
         end if
      end if
c ***************************************************************************

      

c ***************************************************************************
c *              SCALE N_sim & Nmax_steps for reduced base!                 *
c *                      Cid@Lagoa - 09/October/2005                        *
c ***************************************************************************
!     I realized that in StarlightChains_v01.for I fixed the numbers of
!     steps N_sim & Nmax_steps for the EX0s runs as multiples of N_par
!     (as in other parts of the code). As a consequence runs where only
!     a few components are used in the EX0s fits do many more
!     steps/per-parameter than runs whereas galaxies with many EX0s
!     components do much less. A perhaps more reasonable thing to do is
!     to scale the number of steps in terms of the number of Relevant
!     parameters, which is only known once the code gets here. To use
!     this new feature I introduced the IsScaleNstepsInEX0sFits = Y/N
!     (1/0) flag.
      if (IsScaleNstepsInEX0sFits.EQ.1) then
         N_par_Relevant = N_base_Relevant + 1 + N_exAV + 1
         scale_down_fac = float(N_par_Relevant) / float(NO_par)

         N_sim      = nint(N_sim * scale_down_fac)
         Nmax_steps = nint(Nmax_steps * scale_down_fac)

!     Safety: Lower limits to N_sim & Nmax_steps
         if (N_sim.LT.(10*N_par_Relevant)) N_sim = 10 * N_par_Relevant
         if (Nmax_steps.LT.(10*N_sim)) Nmax_steps = 10 * N_sim

!     Screen output
         if (i_verbose.GT.0) then 
            write (*,'(2(a,i8))')
     $           '[EX0s_Fit_Scubi] N_sim      scaled from ' ,
     $           nint(N_sim/scale_down_fac), ' to ' , N_sim
            write (*,'(2(a,i8))')
     $           '[EX0s_Fit_Scubi] Nmax_steps scaled from ' ,
     $           nint(Nmax_steps/scale_down_fac) , ' to ',Nmax_steps
         end if
      end if
c ***************************************************************************



c ***************************************************************************
c *              CONDENSE BASE TO ITS "RELEVANT" COMPONENTS                 *
c ***************************************************************************
c     Reset all Short* arrays
      do j=1,Nb_max
         Shortx_min(j) = 0.
      end do
      do j=1,Npar_max
         ShortGpar_min(j) = 0.
         ShortGpar_ave(j) = 0.
         do i_chain=1,Nchains_max
            ShortChainPar(j,i_chain) = 0.
         end do
      end do


*     Redefine N_base, N_par & par-indices of AV, exAV's & fn
      N_base = N_base_Relevant
      N_par  = N_base + 1 + N_exAV + 1
      j_AV   = N_base + 1
      j_exAV = N_base + 2
      j_fn   = N_par


*     Fill in last N_exAV + 2 elements of ind_Relevant() with AV, exAV's
*     & fn indices.
      ind_Relevant(j_AV) = jO_AV
      do j=j_exAV,j_fn-1        !OBS: This loop only runs if Nex_AV > 0!
         ind_Relevant(j) = jO_exAV + (j - j_exAV)
      end do
      ind_Relevant(j_fn)  = jO_fn

*     Condense/shorten base N_base_Relevant components. Also condense
*     ind_exAV() array & recompute N_exAV_components for the smaller
*     base (both passed by the commom /exAV_stuff/ to other routines).
*     OBS: Everything must be restored before exiting!!
      do j=1,N_base
         j_in_full_base = ind_Relevant(j)

         do i=1,Nl_obs            
            Of_base(j,i) = FULL_Of_base(j_in_full_base,i)
            f_base(j,i)  = FULL_Of_base(j_in_full_base,i)
         end do

         ind_exAV(j) = IFULL_ind_exAV(j_in_full_base)
      end do


c     Condense par-limits & ind_HelpParMove2Average &
c     prob_HelpParMove2Average arrays (passed by common...)
      do j=1,N_par
         j_in_full_base = ind_Relevant(j)
         par_low(j) = Opar_low(j_in_full_base)
         par_upp(j) = Opar_upp(j_in_full_base)
         ind_HelpParMove2Average(j) =
     $        iO_ind_HelpParMove2Average(j_in_full_base)
         prob_HelpParMove2Average(j) =
     $        Oprob_HelpParMove2Average(j_in_full_base)
      end do


*     Reset exAV's NOT used in the condensed base to = 0!  This must be
*     done after the Scubi run too, since irrelevant exAV's are NOT
*     excluded form the parameter list for simplicity!

!!MEX!! NO! Now that we are using exAVs to compute CUMULATIVE total
!     extinctions, we should NOT reset any exAV to zero, even if the
!     components associated to it have zero flux! I have thus commented 
!     out the fragment below.
!
!     Cid@Lagoa - 27/May/2007

!      do jMEX=1,N_exAV
!         i_flag = 0
!         do j=1,N_base
!            if (ind_exAV(j).EQ.jMEX) i_flag = 1
!         end do
!         if (i_flag.EQ.0) then
!            write (*,'(a,a,i2,a)')
!     $           '[MEX_EX0s_Fit_Scubi] ATT: ',
!     $           'None of the Relevant components uses exAV (',jMEX,
!     $           ') so I will reset it to ZERO!'
!            j_aux = jMEX + j_AV
!            do i_chain=1,N_chains
!               ChainPar(j_aux,i_chain) = 0.
!               Cpar_min(j_aux,i_chain) = 0.
!            end do
!            Gpar_min(j_aux) = 0.
!            Gpar_ave(j_aux) = 0.
!         end if
!      end do


c     Condense pop vectors ChainPar() ==> ShortChainPar() for each
c     chain, plus the global average & minimum.
      do j=1,N_par
         j_in_full_base = ind_Relevant(j)
         do i_chain=1,N_chains
            ShortChainPar(j,i_chain) = ChainPar(j_in_full_base,i_chain)
         end do
         ShortGpar_min(j) = Gpar_min(j_in_full_base)
         ShortGpar_ave(j) = Gpar_ave(j_in_full_base)
      end do


c     Since we are throwing away all "irrelevant" components, it is a
c     good idea renormalize all chain-pop-vectors to their respective fn
c     (the irrelevant x_j's may still carry some flux if
c     x_RelevantThreshold > 0). 1st renomalize the chain-pop-vectors.
      do i_chain=1,N_chains
         sum_x = 0.
         do j=1,N_base
            sum_x = sum_x + ShortChainPar(j,i_chain)
         end do
         fn = ShortChainPar(j_fn,i_chain)
         do j=1,N_base
            ShortChainPar(j,i_chain) = ShortChainPar(j,i_chain)*sum_x/fn
         end do
      end do


c     Now renomalize the Global Average & Minimum pop-vectors.
c     (This is probably not necessary but it should do no harm!)
      sum1_x = 0.
      sum2_x = 0.
      do j=1,N_base
         sum1_x = sum1_x + ShortGpar_min(j)
         sum2_x = sum2_x + ShortGpar_ave(j)
      end do
      fn1 = ShortGpar_min(j_fn)
      fn2 = ShortGpar_ave(j_fn)
      do j=1,N_base
         ShortGpar_min(j) = ShortGpar_min(j) * sum1_x / fn1
         ShortGpar_ave(j) = ShortGpar_ave(j) * sum2_x / fn2
      end do


*     Define x_min, AV_min & exAV_min (probably un-necessary!)
!!MEX!! This IS indeed un-necessary, and so I comment it out.
!      do j=1,N_base
!         Shortx_min(j) = ShortGpar_min(j)
!      end do
!      AV_min  = ShortGpar_min(j_AV)
!      do j=j_exAV,j_fn-1        !OBS: This loop only runs if Nex_AV > 0!
!         jMEX           = j - j_exAV + 1
!         exAV_min(jMEX) = ShortGpar_min(j)
!      end do

ccccccccc

c     Empty x_min, ChainPar, Gpar_min & Gpar_ave arrays. These arrays
c     will all be filled back in after the fit with a condensed base is
c     over. Here we empty them just to make sure we are zero-ing the
c     "irrelevant" components!
      do j=1,Nb_max
         x_min(j) = 0.
      end do
      do j=1,Npar_max
         Gpar_min(j)      = 0.
         Gpar_ave(j)      = 0.
         do i_chain=1,Nchains_max
            ChainPar(j,i_chain) = 0.
         end do
      end do


*FIR* Condense FIR-related stuff: NFIR_base, FIRRmat(j,i_AV) &
*     FIRBolCor(j) OBS: Other FIR-arrays need not be condended (as
*     they're not used in the FIR_CalcFIRchi2 function!)
*QHR* Condense QHRelated stuff: NQHR_base & QHR_Q2Lnorm40(j)
*PHO* Condense PHO-related stuff: NPHO_base & PHO_3Mat(j,iY,i_AV)
*     ATT: using NPHO_Ys_max since NPHO_Ys is nor available to this
*     routine!
      NFIR_base = N_base
      NQHR_base = N_base
      NPHO_base = N_base
      do j=1,N_base
         j_in_full_base = ind_Relevant(j)
         do i_AV=1,NFIR_AV
            FIRRmat(j,i_AV) = OFIRRmat(j_in_full_base,i_AV)
         end do
c Rafa
         do i_AV=1,NPHO_AV
            do iY=1,NPHO_Ys_max
               PHO_3Mat(j,iY,i_AV) = OPHO_3Mat(j_in_full_base,iY,i_AV)
            end do
         end do
         FIRBolCor(j) = OFIRBolCor(j_in_full_base)
         QHR_Q2Lnorm40(j) = OQHR_Q2Lnorm40(j_in_full_base)
      end do

c     Screen output (Added on Jan/2010, but not a FIR-thing)
      if (i_verbose.GT.1) then 
         do j=1,N_base
            write (*,'(3(a,i5,3x))') 
     $           '[EX0s_Fit_Scubi] j=', j , 
     $           'ind_Relevant(j)=' , ind_Relevant(j) ,
     $           'ind_exAV(j) =' ,  ind_exAV(j) 
         end do
      end if
c ***************************************************************************


c ***************************************************************************
c *              FIT SPECTRUM USING SHORTENED/CONDENSED BASE                *
c ***************************************************************************
c     Call MEX_Fit_Scubidu with condensed base
      call MEX_Fit_Scubidu(i_verbose,idum,N_loops,v0_low,v0_upp,vd_low
     $     ,vd_upp,Shortx_min,AV_min,exAV_min,v0_min,vd_min,Ntot_EX0s
     $     ,eff_EX0s,Nl_eff,0,N_par,ShortChainPar,N_chains,Temp_ini
     $     ,Temp_fin,fini_eps,ffin_eps,R_ini,R_fin,N_sim,Nmax_steps
     $     ,i_RestartChains,ShortGpar_min,GE_min,ShortGpar_ave
     $     ,ShortCpar_min,ShortCEne_min,IsKinFixed_EX0s,
     $     IsBurInOver_EX0s)
c ***************************************************************************


c ***************************************************************************
c *             RESTORE ORIGINAL BASE & REARRANGE POP-VECTORS               *
c ***************************************************************************
c     Rearrange all parameters onto correct full-base order. (Irrelevant
c     components have been zero'ed above!)
      do j=1,N_par
         j_in_full_base = ind_Relevant(j)
         do i_chain=1,N_chains
            ChainPar(j_in_full_base,i_chain) = ShortChainPar(j,i_chain)
            Cpar_min(j_in_full_base,i_chain) = ShortCpar_min(j,i_chain) 
            CEne_min(i_chain)                = ShortCEne_min(i_chain)
         end do
         Gpar_min(j_in_full_base) = ShortGpar_min(j)
         Gpar_ave(j_in_full_base) = ShortGpar_ave(j)
      end do

*     Reset exAV's NOT used in the condensed base to = 0! 
!!MEX!! NO! Now that we are using exAVs to compute CUMULATIVE total
!     extinctions, we should NOT reset any exAV to zero, even if the
!     components associated to it have zero flux! I have thus commented 
!     out the fragment below.
!
!     Cid@Lagoa - 27/May/2007

!     do jMEX=1,N_exAV
!        i_flag = 0
!        do j=1,N_base
!           if (ind_exAV(j).EQ.jMEX) i_flag = 1
!        end do
!        if (i_flag.EQ.0) then
!           write (*,'(a,a,i2,a)')
!    $           '[MEX_EX0s_Fit_Scubi] ATT: ',
!    $           'None of the Relevant components uses exAV (',jMEX,
!    $           ') so I will reset it to ZERO!'
!            j_aux = jMEX + jO_AV
!            do i_chain=1,N_chains
!               ChainPar(j_aux,i_chain) = 0.
!               Cpar_min(j_aux,i_chain) = 0.
!            end do
!            Gpar_min(j_aux) = 0.
!            Gpar_ave(j_aux) = 0.
!        end if
!     end do

*     Restore Original base: FULL_Of_base ==> Of_base & f_base
      do i=1,Nl_obs
         do j=1,NO_base
            Of_base(j,i) = FULL_Of_base(j,i)
            f_base(j,i)  = Of_base(j,i)
         end do
      end do

*     Restore Original ind_exAV array
      do j=1,NO_base
         ind_exAV(j) = IFULL_ind_exAV(j)
      end do

*     Restore original par-limits & arrays which say whether to help
*     parameters to converge and with what probability.
      do j=1,NO_par
         par_low(j) = Opar_low(j)
         par_upp(j) = Opar_upp(j)
         ind_HelpParMove2Average(j)  = iO_ind_HelpParMove2Average(j)
         prob_HelpParMove2Average(j) = Oprob_HelpParMove2Average(j)
      end do

*     Restore original base dimension, number of exAV components, number
*     of parameters & indices of AV, exAV's & fn
      N_base = NO_base
      N_par  = NO_par
      j_AV   = jO_AV
      j_exAV = jO_exAV
      j_fn   = jO_fn


*     Unpack Gpar_min onto x_min,AV_min & YAV_min explicitly.
      call MEX_Unpack_par2xAVYAVfn(N_par,N_base,Gpar_min,x_min,AV_min
     $     ,exAV_min,fn_min)


*FIR* Restore Original values of FIRRmat(j,ind_FIRAV), FIRBolCor(j) &
*     NFIR_base
*QHR* Restore original values of QHR_Q2Lnorm40(j) & NQHR_base
*PHO* Restore originalvalues of NPHO_base & PHO_3Mat(j,iY,i_AV)
*     ATT: using NPHO_Ys_max since NPHO_Ys is nor available to this
*     routine!
      NFIR_base = NO_base
      NQHR_base = NO_base
      NPHO_base = NO_base
      do j=1,NO_base
         do i_AV=1,NFIR_AV
            FIRRmat(j,i_AV) = OFIRRmat(j,i_AV)
         end do
c Rafa
         do i_AV=1,NPHO_AV
            do iY=1,NPHO_Ys_max
               PHO_3Mat(j,iY,i_AV) = OPHO_3Mat(j,iY,i_AV)
            end do
         end do
         FIRBolCor(j) = OFIRBolCor(j)
         QHR_Q2Lnorm40(j) = OQHR_Q2Lnorm40(j)
      end do
c ***************************************************************************



      end
c ###########################################################################



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!                             !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!  END of MEX-routines block  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!                             !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
c FIR                                                                     FIR
c FIR                     Start of FIR-Routines Block                     FIR
c FIR                                                                     FIR
c FIR                                                                     FIR
c FIR FIRInitialization(i_verbose,N_base,Nl,l,..,FIRbeta_I,FIRbeta_D)     FIR
c FIR FIR_CalcFIRchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)                    FIR
c FIR AverageXoverLambda(x,AV_tot,avelamb_x)                              FIR
c FIR                                                                     FIR
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR


c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
!     Initialize FIR-related stuff.
!
!     For each of the j=1...N_base base spectra in Of_base:
!
!     (1) integrate the spectrum to obtain the bolometric luminosity:
!     FIRLbol(j) and the bolometric correction FIRBolCor(j), defined as
!     the monocromathic base luminosity at l_norm divided by Lbol;

!     (2) Compute the fraction of Lbol which comes from ionizing
!     photons: FIRFracLion(j)
!
!     (3) For each AV in a predefined grid FIRAV_grid(), compute the
!     fraction of Lbol which is absorbed via "optical"
!     extinction/attenuation [R_opt(j,i_AV)], through reprocessing of
!     Ly-alpha [R_Lya(j,i_AV)] and through Lyman Continuum Extinction
!     [R_LCE(j,i_AV)]. These things add up to the total reprocessing
!     matrix Rmat(j,i_AV), which is used to compute predicted FIR
!     luminosities.
!
!     The computation of R_opt is straight-forward: It comes from the
!     integral in lambda of L_SSP(j,l) * (1 - exp(-tau(l)) * dl, divided
!     by Lbol(j). Note that tau(l) = tau_V * q(l), with q(l) being the
!     reddening-law given by the input choice of red_law_option (at this
!     point I strongly recommend CCC!)

!     The conversion of ionizing luminosities to dust-reprocessed
!     radiation is done with a model based on Inoue (2002) ideas, but
!     with adaptations/corrections.  Inoue (2002) proposes (actually,
!     postulates!) that tau_LC^HII (see his paper and my notes) can be
!     scaled from tau_V^ISM, and says that his data shows that qsi =
!     tau_LC^HII/tau_V^ISM is ~ 1. Since he also assumes that
!     tau_V^total = tau_V^SF = tau_V^BC + tau_V^ISM and that tau_V^BC =
!     2 * tau_V^ISM (following Charlot & Fall 2000), we can estimate
!     tau_LC^HII/tau_V^tot is ~ 1/3. This 1/3 factor is what I call
!     FIRbeta_I, read from the config file, passed as an input argument
!     to this function, and used below to estimate tau_LC^HII from
!     AV_tot(j).
!
!     A key ingreedient in Inoue's formulation is f_INOUE = the fraction
!     of ionizing luminosity which is absorbed by gas (photoionization),
!     as a function of the dust optical depth to ionizing photons
!     (tau_LC^HII). The formula adopted by Inoue (2002) is eq 8 of
!     Petrosian et al (1972), but that formula gives WRONG results for
!     low tau_LC^HII! (Besides having the wrong assymptotic behaviour,
!     it also suffers from numerical instability). Luckly, I have
!     verified that exp(-0.785 * tau_LC^HII) gives very-very similar
!     results, so I adopt this expression here.
!
!     Besides LCE, Inoue postulates (quite reasonably) that ALL Ly-alpha
!     flux is thermalized by dust (after scattering around in the HII
!     region). We adopt Dopita's suggestion that a fraction beta_D (of
!     the order of 30%) of the ionizing luminosity is converted into
!     Ly-alpha, and ultimately transformed (by dust) into FIR
!     emission. This is Ok, but there is a caveat: We turn all Ly-alpha
!     flux/luminosity into FIR EVEN WHEN AV= 0!! Beware of this caveat.
!
!
!     All these was done when implementing the FIR-constraints in our
!     MEX-STARLIGHT analysis of LIRGs.
!
!
!     Misc-notes 1: Notice that the array FIRAV_grid() is actually
!     superfluous, as all that is actually needed to locate an entry in
!     this grid is FIRAV_low, FIRAV_upp and FIRdAV. Still, as it does no
!     harm and helps debugging, FIRAV_grid() is defined.
!
!     Misc-notes 2: FIRAV_grid(Nl_max) and FIRRmat(Nb_max,Nl_max) are
!     defined using Nl_max simply because this is a suitably large
!     number available thoughout the code.
!IIB! I now changed Nl_max to NAVgrid_max, in an effort to reduce memory
!     cost and fix the IIB! Cid@Lagoa - 11/Jan/2011
!
!     Misc-notes 3: Since the whole base spectrum (from "lambda = 0") is
!     needed, this routine MUST be called BEFORE chopping the original
!     base spectra to the optical range being synthesized! 
!
!     Cid@Granada - 29/Jan/2010

      SUBROUTINE FIRInitialization(i_verbose,N_base,Nl,l,l_norm,
     $     red_law_option,FIRbeta_I,FIRbeta_D)

c ***************************************************************************
c     Defs
      
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
!IIB! Changed Nl_max to NAVgrid_max!
      parameter (NAVgrid_max = 1001)
      real l(Nl_max) , dl(Nl_max) , q(Nl_max) , f_aux(Nl_max)
      real Of_base(Nb_max,Nl_max)
      character*3 red_law_option
      common/Obase/ Of_base
!IIB! Changed Nl_max to NAVgrid_max!
      common/FIR_BaseStuff/ NFIR_AV , FIRAV_grid(NAVgrid_max) ,
     $     FIRAV_low , FIRAV_upp , FIRdAV , 
     $     FIRLbol(Nb_max) , FIRFracLion(Nb_max) , FIRBolCor(Nb_max) ,
     $     FIRR_opt(Nb_max,NAVgrid_max) , FIRR_Lya(Nb_max,NAVgrid_max) , 
     $     FIRR_LCE(Nb_max,NAVgrid_max) ,
     $     FIRRmat(Nb_max,NAVgrid_max) , FIRq_norm , NFIR_base

c     Reset things not zero'ed in the MAIN routine.
      do i=1,Nl_max
         dl(i)         = 0.
         q(i)          = 0.
         f_aux(i)      = 0.
      end do

c     Warn & STOP (!) if the wavelength array seems too restricted for
c     FIR-modeling...
      if ((l(1).GT.800.0).OR.(l(Nl).LT.7000.0)) then
         write (*,'(a,2(a,1x,f9.1))')  '[FIRInitialization] ATT! ' ,
     &        '  lambda covers only from ',l(1),' to = ', l(Nl)
         stop
      end if

c     Warn if reddening-law being used is not = CCC. The CCC law is NOT
c     really mandatory, but it is the only one I tested!
      if (red_law_option.NE.'CCC') then 
         write (*,'(70(a))') ('*',i=1,70)
         write (*,'(4(a))')  '[FIRInitialization] ATT!' ,
     &        ' red_law_option is = ' , red_law_option , 
     &        ' but I recommend CCC!!'
         write (*,'(70(a))') ('*',i=1,70)
      end if

c     Setup q-law: Use Calc_red_term routine & invert red_term(l) = -0.4
c     * (q(l) - q(l_norm)) to obtain q(l), zeroing small things by hand
c     to circumvent numerical residuals. Just in case, also SET q = 0
c     FOR l < 912 Angs.
      call Calc_red_term(red_law_option,l_norm,Nl,l,f_aux,FIRq_norm)
      do i=1,Nl
         q(i) = FIRq_norm - 2.5 * f_aux(i)
         if ((l(i).LT.912.0).OR.(q(i).LE.1.e-5)) q(i) = 0.
      end do

c     Setup dl-array
      dl(1)  = (l(2)  - l(1)   ) / 2.
      dl(Nl) = (l(Nl) - l(Nl-1)) / 2.
      do i=2,Nl-1
         dl(i) = (l(i+1) - l(i-1)) / 2.
      end do

c     Find the index of lambda >= 912 Angs (for convenience/speed)
c     ind_l912 is the index of the **last lambda BEFORE 912 Angs**.
      ind_l912 = 0
      do i=1,Nl-1
         if ((l(i).LT.912.0).AND.(l(i+1).GE.912.0)) ind_l912 = i
      end do
      if (i_verbose.GT.0) write (*,'(a,i4,2(a,f7.2))')
     $     '[FIRInit..] index of last l < 912 Angs = ', ind_l912 ,
     $     ' ==> l(i)='   , l(ind_l912), ' & l(i+1)=' , l(ind_l912+1)
c ***************************************************************************


c ***************************************************************************
c     Store N_base in NFIR_base, such that it is available through the
c     FIR_BaseStuff common above.
      NFIR_base = N_base

c     Sanity check on NFIR_AV
      if (NFIR_AV.GT.NAVgrid_max) then
         write (*,'(a,2(a,1x,i5))')  '[FIRInitialization] ATT! ' ,
     &        '  NFIR_AV = ', NFIR_AV , ' too large! > ', NAVgrid_max
         stop
      end if

c     Setup FIRAV_grid for R_matrix(i_AV=1...NFIR_AV;j=1,N_base).
c     NFIR_AV , FIRAV_low & FIRAV_upp were read from arq_config.
      FIRdAV = (FIRAV_upp - FIRAV_low) / (NFIR_AV - 1)
      do i_AV=1,NFIR_AV
         FIRAV_grid(i_AV) = FIRAV_low + (i_AV - 1) * FIRdAV
      end do

c     Screen output
      if (i_verbose.GT.0) write (*,'(a,i4,3(2x,a,f6.3))') 
     &     '[FIRInit..] AVgrid: NFIR_AV=',NFIR_AV,
     &     'FIRAV_low=',FIRAV_low , 'FIRAV_upp=',FIRAV_upp ,
     &     'FIRdAV=',FIRdAV
c ***************************************************************************


c ***************************************************************************
c     For each base spectrum, compute Lbol(j), BolCor(j) (=
c     monocromathic base luminosity at l_norm divided by Lbol),
c     FracLion(j), & Reprocessed-Lbol-fraction matrices R_opt(j,i_AV),
c     R_Lya(j,i_AV), R_LCE(j,i_AV) & the total one Rmat(j,i_AV).
      do j=1,N_base

c     -----------------------------------------------------------------------
c     Lbol(j) & BolCor(j)
         Lbol = 0.
         do i=1,Nl
            Lbol = Lbol + Of_base(j,i) * dl(i)
            f_aux(i) = Of_base(j,i)
         end do
         call calc_Fnorm(Nl,l,f_aux,l_norm,f_norm)
         FIRLbol(j)   = Lbol
         FIRBolCor(j) = f_norm / Lbol

c     FracLion(j) = ionizing luminosity / bolometric luminosity
         Lion = 0.
         do i=1,ind_l912
            Lion = Lion + Of_base(j,i) * dl(i)
         end do
         FIRFracLion(j) = Lion / Lbol
c     -----------------------------------------------------------------------

c     -----------------------------------------------------------------------
c     R_opt + R_Lya + R_LCE & Rmat reprocessing (j,i_AV) matrices.
         do i_AV=1,NFIR_AV
            tauV = FIRAV_grid(i_AV) * 0.4 / log10(exp(1.))

c     R_opt(j,i_AV) = "optical" (l >= 912 Angs, non-ionizing)
c     reprocessing matrix, given as a fraction of Lbol.
            s = 0.
            do i=ind_l912+1,Nl
               s = s + Of_base(j,i) * dl(i) * (1. - exp(-tauV * q(i)))
            end do
            FIRR_opt(j,i_AV) = s / FIRLbol(j)

c     R_Lya(j,i_AV) & R_LCE(j,i_AV) reprocessing matrices (for ionizing
c     photons). Both give fractions of Lbol. Notice the use of the
c     beta_D & beta_I parameters (this is the only place they are
c     used!), as well as my approximate (and more correct in the small
c     tau limit!) version of eq 8 (for the factor f_INOUE) of Petrosian
c     et al (1972). This factor measures/estimates the fraction of Lion
c     which is used up by the gas (photoionization), so that 1 - f_INOUE
c     remains for the dust-LCE.
            tau_dust_LC_HII = FIRbeta_I * tauV
            f_Inoue = exp(-0.785 * tau_dust_LC_HII)
            FIRR_Lya(j,i_AV) = FIRFracLion(j) * FIRbeta_D * f_Inoue
            FIRR_LCE(j,i_AV) = FIRFracLion(j) * (1. - f_Inoue)

c     TOTAL reprocessing matrix: Rmat(j,i_AV) = R_opt + R_Lya + R_LCE
            FIRRmat(j,i_AV) = FIRR_opt(j,i_AV) + FIRR_Lya(j,i_AV) +
     $           FIRR_LCE(j,i_AV) 
         end do
c     -----------------------------------------------------------------------

      end do
c ***************************************************************************


c ***************************************************************************
c     Screen output
      if (i_verbose.GT.0) then
         istep = float(NFIR_AV) / float(21)
         ilast = 21 * istep
 11      format ('[FIRInit..]',i4,2(1x,f7.4),2x,f5.2,2x,a,201(1x,i4))
 12      format (153('-'))
 13      format ('[FIRInit..]   j    BC_j  Lbol_j  %-ion  ',
     &        'AV_grid=',21(1x,f4.2))
         write (*,12)
         write (*,13) (FIRAV_grid(i_AV),i_AV=1,ilast,istep)
         write (*,12)
         do j=1,N_base
            write (*,11) j , log10(FIRBolCor(j)) , log10(FIRLbol(j)) ,
     &           100.*FIRFracLion(j) , '%-R_opt=' ,
     &           (int(0.5+100*FIRR_opt(j,i_AV)),i_AV=1,ilast,istep)
            write (*,11) j , log10(FIRBolCor(j)) , log10(FIRLbol(j)) ,
     &           100.*FIRFracLion(j) , '%-R_Lya=' ,
     &           (int(0.5+100*FIRR_Lya(j,i_AV)),i_AV=1,ilast,istep)
            write (*,11) j , log10(FIRBolCor(j)) , log10(FIRLbol(j)) ,
     &           100.*FIRFracLion(j) , '%-R_LCE=' , 
     &           (int(0.5+100*FIRR_LCE(j,i_AV)),i_AV=1,ilast,istep)
            write (*,11) j , log10(FIRBolCor(j)) , log10(FIRLbol(j)) ,
     &           100.*FIRFracLion(j) , '%-Rmat =' , 
     &           (int(0.5+100*FIRRmat(j,i_AV)),i_AV=1,ilast,istep)
            if (j.NE.N_base) write (*,*)
         end do
         write (*,12)
         write (*,*)
      end if
c ***************************************************************************

      end
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR


c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
!     Computes the FIR-part of chi2...
!
!     Given x(j) and AV_tot(j) population and extinction vectors from my
!     MEx-STARLIGHT fits, and the stuff initialized in routine
!     FIRInitialization (passed though the FIR_BaseStuff common), this
!     function computes the predicted FIR luminosity (actually the total
!     bolometric dust-attenuated luminosity) and, comparing it to the
!     observed FIR luminosity, returns a chi2.
!
!     The value of chi2 can be a simple [(data-model)/error]**2, but it
!     can be more complicated/interesting/powerful than this, to allow
!     for upper-limits and range-fitting. If the predicted L_FIR
!     [FIRModlogLFIR, in L_sun if BC03-like bases are being used AND if
!     the input optical spectrum is given in units of ergs/s/cm2/Angs -
!     THIS IS A MUST for FIR-fits!!] lies within the FIR_logLFIR_low -->
!     FIR_logLFIR_upp range, then we assign chi2_FIR = 0. Outside this
!     flat core, chi2 grows quadraticaly, emulating a gaussian
!     likelyhood. Setting low = upp limits (they both come from the
!     input file arq_FIRinfo), we obtain a simple gaussian, but I
!     anticipate that this is not what we will be more interested in
!     the LIRGS-work with Rosa@IAA/2010.
!
!     In order that chi2_FIR, which is based on a single datum, competes
!     in absolute value with, chi2_OPTICAL, the folllowing scaling is
!     adopted: chi2_FIR --> chi2_FIR * FIRChi2ScaleFactor * chi2_OPTICAL
!     (This recipe is more general than scaling by Nl_eff, as done in
!     initial FIRc-versions.)
!
!     WhereFrom is just a string which tells from where this function
!     was called. When it is = 'OUT', then the FIR-pop vector is
!     computed (for output purposes only). 
!
!     Misc notes 1: In this version of the code EVERY time chi2 is
!     evaluated this function is called up! When the fit is NOT a
!     FIR-fit (i.e., no FIR data provided), then this function returns 0
!     immediately after being called.
!
!     Cid@La-Laguna - 07/Feb/2010
      
      FUNCTION FIR_CalcFIRchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)

c     Definitions
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
!IIB! Changed Nl_max to NAVgrid_max!
      parameter (NAVgrid_max = 1001)
      real x(Nb_max) , AV_tot(Nb_max)
      character*3 WhereFrom

!IIB! Changed Nl_max to NAVgrid_max!
      common/FIR_BaseStuff/ NFIR_AV , FIRAV_grid(NAVgrid_max) ,
     $     FIRAV_low , FIRAV_upp , FIRdAV , 
     $     FIRLbol(Nb_max) , FIRFracLion(Nb_max) , FIRBolCor(Nb_max) ,
     $     FIRR_opt(Nb_max,NAVgrid_max) , FIRR_Lya(Nb_max,NAVgrid_max) , 
     $     FIRR_LCE(Nb_max,NAVgrid_max) ,
     $     FIRRmat(Nb_max,NAVgrid_max) , FIRq_norm , NFIR_base

      common/FIR_DataStuff/ FIRlogLobs_norm , 
     $     FIR_logLFIR_TOT , FIR_LFIRFrac2Model ,
     $     FIR_logLFIR_obs , FIR_ErrlogLFIR ,
     $     FIR_logLFIR_low , FIR_logLFIR_upp , FIR_RangelogLFIR ,
     $     FIRChi2ScaleFactor , iFIR_IsFIRcOn

      common/FIR_OutputStuff/ FIRModlogLFIR , FIRModlogLBOL , 
     &     x_FIR(Nb_max),x_BOL(Nb_max)

c     Reset output chi2
      FIR_CalcFIRchi2 = 0.

c     If FIR is not turned on, just set the FIR-chi2 to zero & go back!
c     OBS: In fits without FIR, I may temporarely "turn FIR On" just to
c     call this routine to obtain the model FIR luminosity...
      if (iFIR_IsFIRcOn.EQ.0) return

c     Compute predicted L_FIR. See notes for formulae!
c     OBS1: NO interpolation of AV_tot() in the Rmat grid is done :( 
c           This **assumes/works-well-if** that the Rmat grid is 
c           sufficiently finely spaced (~0.1 mag should be enough).
c     OBS2: Skip negative extinctions, in case you're daring to allow
c           so. In principle, AV=0 could also be skipped (no dust =>
c           no FIR-reprocessing) but Ly-alpha is assumed to be eaten
c           by dust even when AV=0!
      s = 0.
      do j=1,NFIR_base
         if (AV_tot(j).GE.0.0) then

c     Find index of entry in FIRAV_grid which best matches AV_tot
            ind_FIRAV = int(0.5 + (AV_tot(j) - FIRAV_low) / FIRdAV) + 1
            if (ind_FIRAV.LT.1) ind_FIRAV = 1
            if (ind_FIRAV.GT.NFIR_AV) ind_FIRAV = NFIR_AV

c     Computes this element's addition to the integral s
            z1 = 10.**(0.4 * FIRq_norm * AV_tot(j))
            z2 = FIRRmat(j,ind_FIRAV) / FIRBolCor(j)
            s  = s + x(j) * z1 * z2

!OOPA! Testing... 
!           if (NFIR_base.LT.12) then
!            write (*,34) j,100.*x(j),AV_tot(j),ind_FIRAV,
!     $           FIRAV_grid(ind_FIRAV),FIRRmat(j,ind_FIRAV),
!     $           FIRBolCor(j),s
! 34         format ('[f] ',i2,' x=',f5.1,' AV=',f7.4,' iAV=',i4
!     $           ,' AVgrid=',f7.4,' Rmat=',f7.5,'  BC=',1p,e12.5,0p
!     $           ,' s=',1p,e9.3)
!           end if
         end if
      end do

c     Before taking the log of the sum (s), add 1e-30 for safety!
      FIRModlogLFIR = FIRlogLobs_norm + log10(s + 1.e-30)

c     chi2-value for a Gaussian-Smooth-Top-Hat Likelyhood function,
c     which gives parabolic chi2-profiles with a core @ chi2 = zero.
*ETC* Added if iFIR_IsFIRcOn.EQ.1 to cope with -1 (predict-only) option!
*     (I did this change during QHR+ETC-development, @ Jan/2010!)
      chi2_FIR = 0.
      if (iFIR_IsFIRcOn.EQ.1) then
        if (FIRModlogLFIR.GE.FIR_logLFIR_upp) then
          chi2_FIR = ((FIRModlogLFIR-FIR_logLFIR_upp)/FIR_ErrlogLFIR)**2
        else if (FIRModlogLFIR.LE.FIR_logLFIR_low) then
          chi2_FIR = ((FIRModlogLFIR-FIR_logLFIR_low)/FIR_ErrlogLFIR)**2
        end if
      end if

c     Scale chi2_FIR by FIRChi2ScaleFactor * chi2_OPTICAL
      FIR_CalcFIRchi2 = chi2_FIR * FIRChi2ScaleFactor * chi2_OPTICAL

c     FINAL OUTPUT: Calling this function with WhereFrom = 'OUT' means
c     we are done, and all that is left to do is (for output purposes
c     only), to compute the x_FIR population vector (the fractional
c     contribution of each base element to the dust-reprocessed
c     FIR-luminosity, normalize by the s-integral from above). This
c     vector, along with the predicted FIR-luminosity, is passed back to
c     the calling routine via the common FIR_OutputStuff.
      if (WhereFrom.EQ.'OUT') then
         do j=1,NFIR_base
            if (AV_tot(j).GE.0.0) then
               ind_FIRAV = int(0.5 + (AV_tot(j) - FIRAV_low)/FIRdAV) + 1
               if (ind_FIRAV.LT.1) ind_FIRAV = 1
               if (ind_FIRAV.GT.NFIR_AV) ind_FIRAV = NFIR_AV
               z1 = 10.**(0.4 * FIRq_norm * AV_tot(j))
               z2 = FIRRmat(j,ind_FIRAV) / FIRBolCor(j)
               x_FIR(j) = x(j) * z1 * z2 / s
            end if
         end do

c     Bolometric pop vector x_BOL(j) & luminosity
         sBOL = 0.
         do j=1,NFIR_base
            z1 = (10.**(0.4 * FIRq_norm * AV_tot(j))) / FIRBolCor(j)
            sBOL = sBOL + x(j) * z1
            x_BOL(j) = x(j) * z1
         end do
         do j=1,NFIR_base
            x_BOL(j) = x_BOL(j) / sBOL
         end do
         FIRModlogLBOL = FIRlogLobs_norm + log10(sBOL + 1.e-30)

      end if

!     Screen output/test
!     if (WhereFrom.EQ.'FKc') then
!      write (*,'(a,a,a,2(a,f7.3,1x),4(a,1p,e10.3,0p,2x))') 
!     $     '[FIR_CalcFIRchi2 from ' , WhereFrom , '] : ' ,
!     $     'LFIR=' , FIRModlogLFIR ,' X ' , FIR_logLFIR_obs ,
!     $     'chi2_FIR=' , FIR_CalcFIRchi2 , 
!     $     'chi2_Opt=' , chi2_OPTICAL ,
!     $     'FIR/Opt=', FIR_CalcFIRchi2 / (1.e-30+chi2_OPTICAL) ,
!     $     's=' , s
!      end if

      end
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR


c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
c     This routine uses a given input pop & extinction vector (x &
c     AV_tot) to compute the average of the pop-vector across all
c     lambdas. The motivation for this came during the FIR-work
c     (@Granada/Jan/2010), and is the fact that sometimes a base
c     component contributes little at lambda = l_norm bu may be much
c     more sigificant at another lambda. Hence, to decide which
c     components are relevant (as we don in the EX0s-routines), it may
c     be important to considered all lambdas in the fit, instead of just
c     l_norm. The output is avelamb_x(J=1...N_base), which give you the
c     weigthed average of x(lambda) across all lambdas; the weigthing is
c     done with weight_obs.
c
c     Cid@Granada - 31/Jan/2010

      SUBROUTINE AverageXoverLambda(x,AV_tot,avelamb_x)

c     Defs
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)

      real x_lambda(Nb_max,Nl_max) , aux(Nb_max) , avelamb_x(Nb_max)
      real x(Nb_max) , AV_tot(Nb_max)
      real f_base(Nb_max,Nl_max)
      real l_obs(Nl_max) , f_obs(Nl_max), weight_obs(Nl_max)
      common/Red_factor/ red_term(Nl_max)
      common/Base/ N_base , f_base
      common/Obs/ Nl_obs , l_obs , f_obs , weight_obs

c     1st compute x_lambda(j,i) = matrix of lambda-dependent pop vectors
c     OBS: Keep normalization equal to that of the input x.
      s_x = Sum_of_x(x,N_base)
      do i=1,Nl_obs
         do j=1,N_base
            aux(j) = (abs(x(j)) * f_base(j,i) * 
     &           10.**(AV_tot(j) * red_term(i)))
         end do
         s_aux = Sum_of_x(aux,N_base)
         do j=1,N_base
            x_lambda(j,i) = aux(j) * s_x / s_aux
         end do
      end do

c     Now, for each j, average x_lambda over lambda, weigthing by
c     weight_obs
      s_w = 0.
      do i=1,Nl_obs
         s_w = s_w + weight_obs(i)
      end do

      do j=1,N_base
         s = 0.
         do i=1,Nl_obs
            s = s + x_lambda(j,i) * weight_obs(i)
         end do
         avelamb_x(j) = s / s_w
      end do

      end 
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR


c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
c FIR                                                                     FIR
c FIR                      End of FIR-Routines Block                      FIR
c FIR                                                                     FIR
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR






c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
c QHR                                                                     QHR
c QHR                     Start of QHR-Routines Block                     QHR
c QHR                                                                     QHR
c QHR                                                                     QHR
c QHR QHRInitialization(i_verbose,N_base,Nl,l,l_norm,red_law_option)      QHR
c QHR QHR_CalcQHRchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)                    QHR
c QHR                                                                     QHR
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR


c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
!     Initialize QHR-related stuff.
!
!     For each base spectrum, compute qH in units of 1e40
!     photons/sec/initial-solar-mass[QHR_qH40(j)], and the ratio of QH
!     to the monochromatic luminosity at l_norm in units of 1e40
!     (photons/s) / (Lsun/Angs) [QHR_Q2Lnorm40(j)]
!     OBS: The units of 1e40 came about to allow single precision work.

!     This routine produces the following output, all of which are
!     passed through common/QHR_BaseStuff/:
!        QHR_qH40()
!        QHR_Q2Lnorm40() ,
!        QHR_q_norm
!        NQHR_base
!
!     Cid@Lagoa - 31/Dec/2010

      SUBROUTINE QHRInitialization(i_verbose,N_base,Nl,l,l_norm,
     $     red_law_option)

c ***************************************************************************
c     Defs
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
      real l(Nl_max) , dl(Nl_max) , f_aux(Nl_max)
      real Of_base(Nb_max,Nl_max)
      character*3 red_law_option
      common/Obase/ Of_base
      common/QHR_BaseStuff/ QHR_qH40(Nb_max) , QHR_Q2Lnorm40(Nb_max) ,
     $     QHR_q_norm , QHRbeta_I , NQHR_base

c     Reset things not zero'ed in the MAIN routine.
      do i=1,Nl_max
         dl(i)    = 0.
         f_aux(i) = 0.
      end do

c     Warn & STOP (!) if the wavelength array seems too restricted for
c     QHRc-modeling...
      if (l(1).GT.800.0) then
         write (*,'(a,2(a,1x,f9.1))')  '[QHRInitialization] ATT! ' ,
     &        '  lambda covers only from ',l(1),' to = ', l(Nl)
         stop
      end if

c     Defines QHR_q_norm, passed on to other routines via the
c     QHR_BaseStuff common.
      call Calc_red_term(red_law_option,l_norm,Nl,l,f_aux,QHR_q_norm)

c     Store N_base in NQHR_base, such that it is available through the
c     QHR_BaseStuff common above.
      NQHR_base = N_base

c     Setup dl-array for trapezoidal integrations.
      dl(1)  = (l(2)  - l(1)   ) / 2.
      dl(Nl) = (l(Nl) - l(Nl-1)) / 2.
      do i=2,Nl-1
         dl(i) = (l(i+1) - l(i-1)) / 2.
      end do

c     For each base spectrum, compute QHR_qH(j) in units of 1e40 photons/sec/initial-solar-mass
c     OBS: The units of 1e40 allow single precision work: 1e2 = 10.**(33 + 27 - 18 - 40)
c     Also compute the ratio of QH to the monochromatic luminosity at
c     l_norm, QHR_Q2Lnorm40(j), in units of (1e40 photons/s) /
c     (Lsun/Angs). In fact, this is what I really need for QHRc-fits.
      if (i_verbose.GT.0) write (*,12)
      K_Lsun_hc = (3.826 / (6.626 * 2.997925)) * 1.e2
      do j=1,N_base
         s = 0.
         do i=1,Nl
            if (l(i).LT.912.0) s = s + l(i) * Of_base(j,i) * dl(i)
            f_aux(i) = Of_base(j,i)
         end do
         call calc_Fnorm(Nl,l,f_aux,l_norm,f_norm)
         QHR_qh40(j)      = K_Lsun_hc * s
         QHR_Q2Lnorm40(j) = K_Lsun_hc * s / f_norm
         if (i_verbose.GT.0)  write (*,11) j , QHR_qh40(j) ,
     $        QHR_Q2Lnorm40(j)
 11      format ('[QHRInit..]',i4,' qH__40 = ',1p,e9.3,0p,3x,
     $        ' QH/Lnorm__40 = ',1p,e9.3)

c     RADICAL SAFETY MEASURE: Stop the code if any population has no
c     ionizing radiation! 
         if (s.LE.0.) then
            write (*,'(a,i4,a,e12.3)') '[QHRInit...] ATT population ' , 
     $           j , ' has no ionizing flux: s = ', s
            write (*,'(a)') '[QHRInit...] Run aborted:('
            stop
         end if
         
      end do
      if (i_verbose.GT.0) then
         write (*,12)
         write (*,*)
      end if
 12   format (62('-'))
c ***************************************************************************

      end
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR


c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
!     Computes the QHR-part of chi2...

!     Given x(j) and AV_tot(j) population and extinction vectors from my
!     MEx-STARLIGHT fits, and the stuff initialized in routine
!     QHRInitialization (passed though the QHR_BaseStuff common), this
!     function computes the predicted emission line luminosities (Y) for
!     the NQHR_Ys line informed in the arq_ETCinfo and passed on via
!     common/QHR_DataStuff/. The main output is the chi2 of the
!     comparison of these predicted Ys with the observed ones.

!     The predicted Ys take into account the contribution of each base
!     element to the ionizing photon rate, accounting for it's own
!     extinction [AV_tot(j)] both in the predicted Y (which must be
!     reddedned to be compared to the observed one) and in the reduction
!     of the population's QH due to dust [a la Inoue, as in the FIR-work].
!     An AV-dependent fraction of QH is presumed to be eaten up by dust.
!     The fraction which actually goes into photoionization is f_Inoue =
!     exp(-0.785*beta_I*tau_V) [see MY FIRc-related notes]. In routine
!     FIRInitialization this is coded in as follows:
!            tauV = FIRAV_grid(i_AV) * 0.4 / log10(exp(1.))
!            tau_dust_LC_HII = FIRbeta_I * tauV
!            f_Inoue = exp(-0.785 * tau_dust_LC_HII)
!     Here I rewrite this as 10**(-0.4*AV_tot(j)*0.785*QHRbeta_I). This
!     factor is included in the emission-line predictions below (as a
!     reduction of the actually available number of ionizing photons).
!
!     NOTE: When doing this, I realized there is a subtle detail:
!     f_INOUE, as defined in my FIRc-work, is the fraction of the
!     ionizing luminosity [hence energy] which goes into
!     photoionization, whereas what we'd need here is a fraction of
!     ionizing photons! I will IGNORE this difference, which tantamounts
!     to assuming no big differences in mean ionizing photon-energies
!     exists...

!     As in the FIRc-routines, the value of chi2 is = 0 between lower &
!     upper limits, rising parabolically beyond these limits, such that
!     the likelyhood is a flat-top with gaussian wings. The chi2 of each
!     emission line is scaled such that it can compete with the optical
!     fit (same scheme as for the FIR work).

!     WhereFrom is just a string which tells from where this function
!     was called. When it is = 'OUT', then some quantities are computed
!     for output purposes only!
      
!     Cid@Lagoa - 03/Jan/2011

!ELR! I've now included the option to fit an EmissionLineRatio (ELR)!
!     The ELR info is read in the QHR-block of arq_ETCinfo. If you only
!     want to fit the ELR, then you must still provide the QHR info, but
!     may set the chi2-scale factors of the QHR-lines to 0!  Notice that
!     chi2_ELR is included in the calculation of the QHR-chi2.
!
!     Cid@Granada - 29/Kan/2011


      FUNCTION QHR_CalcQHRchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)

c ***************************************************************************
c     Definitions
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nb_max = 300, NQHR_Ys_max=10)
      real x(Nb_max) , AV_tot(Nb_max)
      character*3 WhereFrom

      common/QHR_BaseStuff/ QHR_qH40(Nb_max) , QHR_Q2Lnorm40(Nb_max) ,
     $     QHR_q_norm , QHRbeta_I , NQHR_base

      common/QHR_DataStuff/ QHR_lambda(NQHR_Ys_max) ,
     $     QHR_frecomb(NQHR_Ys_max)  ,
     $     QHR_logY_TOT(NQHR_Ys_max) , QHR_YFrac2Model(NQHR_Ys_max) ,
     $     QHR_ErrlogY(NQHR_Ys_max)  , QHR_RangelogY(NQHR_Ys_max) ,
     $     QHR_Chi2ScaleFactor(NQHR_Ys_max) , 
     $     QHR_logY_obs(NQHR_Ys_max) ,
     $     QHR_logY_upp(NQHR_Ys_max) , QHR_logY_low(NQHR_Ys_max)  ,
     $     QHR_qlambda(NQHR_Ys_max) ,
     $     QHR_GlobalChi2ScaleFactor ,
     $     QHR_logLobs_norm , 
     $     NQHR_Ys , iQHR_IsQHRcOn

      common/QHR_OutputStuff/ QHR_x_QH0(Nb_max) , QHR_x_QHeff(Nb_max) ,
     $     QHR_x_Y(Nb_max,NQHR_Ys_max) , QHR_ModlogY(NQHR_Ys_max) ,
     $     QHR_chi2_Y(NQHR_Ys_max) , QHR_log_QH0 , QHR_log_QHeff ,
     $     ELR_chi2_logR , ELR_ModlogR

*ELR* Data-info for ELR-work - goes together w/QHR-things!
      common/ELR_DataStuff/ ELR_lambdaA , ELR_lambdaB , 
     $     ELR_logRint , ELR_logRobs , 
     $     ELR_ErrlogR , ELR_logR_low , ELR_logR_upp , 
     $     ELR_RangelogR , ELR_Chi2ScaleFactor , 
     $     iY_ELR_A , iY_ELR_B , IsELROn , 
     $     IsELRInputInconsistentWarningFlag

c     Reset output chi2
      QHR_CalcQHRchi2 = 0.

c     If QHR is not turned on, just set the QHR-chi2 to zero & go back!
c     OBS: In non-QHRc-fits, I may temporarely "turn QHR On" just to
c     call this routine to obtain the model QHRelated quantities...
      if (iQHR_IsQHRcOn.EQ.0) return
c ***************************************************************************


c ***************************************************************************
c     For each line, compute the predicted line luminosity, which is the
c     sum of luminosities coming from each population, each one with its
c     own total extinction, and accounting for the fraction of ionizing
c     radiation which is eaten up by dust (and thus not used in
c     photoionization). The log of each of these NQHR_Ys
c     line-luminosities is used to compute the chi2, in the fashion
c     described above.
c     OBS: AV < 0 is reset to zero here!
      s_chi2_Ys = 0.
      do iY=1,NQHR_Ys

c     Add upp contributions of all base elements to line-Y.
         s = 0.
         do j=1,NQHR_base
            aux_AV_tot = AV_tot(j)
            if (AV_tot(j).LT.0.0) aux_AV_tot = 0.
            z = 0.4 * aux_AV_tot * (QHR_q_norm - QHR_qlambda(iY) - 
     $           0.785 * QHRbeta_I)
            s = s + QHR_Q2Lnorm40(j) * x(j) * 10.**z
         end do

         if (s.LE.0.) then
            write (*,*) '[QHR_CalcQHRchi2] HELP!  s = ',s, ' iY=',iY
            write (*,*) '                  Run aborted :('
            do j=1,NQHR_base
               write (*,*) j , x(j) , AV_tot(j) ,  QHR_Q2Lnorm40(j) ,
     $              QHR_q_norm , QHR_qlambda(iY) , 0.785 * QHRbeta_I ,
     $              ' from : ', WhereFrom
            end do
            stop
         end if
         QHR_ModlogY(iY) = log10(s) + QHR_logLobs_norm + 
     $        log10(5.1918e-2 / (QHR_lambda(iY) * QHR_frecomb(iY)))

c     Compute the corresponding contribution to chi2
c     Added if iFIR_IsFIRcOn.EQ.1 to cope with -1 (predict-only) option!
         chi2_Y = 0.
         if (iQHR_IsQHRcOn.EQ.1) then
            if (QHR_ModlogY(iY).GE.QHR_logY_upp(iY)) then
               z = (QHR_ModlogY(iY) - QHR_logY_upp(iY)) /QHR_ErrlogY(iY)
               chi2_Y = z**2
            else if (QHR_ModlogY(iY).LE.QHR_logY_low(iY)) then
               z = (QHR_ModlogY(iY) - QHR_logY_low(iY)) /QHR_ErrlogY(iY)
               chi2_Y = z**2
            end if
         end if
         QHR_chi2_Y(iY) = chi2_Y * QHR_Chi2ScaleFactor(iY) *chi2_OPTICAL
         s_chi2_Ys = s_chi2_Ys + QHR_chi2_Y(iY) 

      end do

c     For ELR-fits, compute the predicted A/B emission line ratio and
c     the associated chi2
      if (IsELROn.EQ.1) then
         ELR_ModlogR = QHR_ModlogY(iY_ELR_A) - QHR_ModlogY(iY_ELR_B)

         chi2_ELR = 0.
         if (ELR_ModlogR.GE.ELR_logR_upp) then
            z = (ELR_ModlogR - ELR_logR_upp) / ELR_ErrlogR
            chi2_ELR = z**2
         else if (ELR_ModlogR.LE.ELR_logR_low) then
            z = (ELR_ModlogR - ELR_logR_low) / ELR_ErrlogR
            chi2_ELR = z**2
         end if

         ELR_chi2_logR = chi2_ELR * ELR_Chi2ScaleFactor * chi2_OPTICAL
         s_chi2_Ys = s_chi2_Ys + ELR_chi2_logR
      end if

c     Here's the result of this function
      QHR_CalcQHRchi2 = s_chi2_Ys
c ***************************************************************************


c ***************************************************************************
c     FINAL OUTPUT: Calling this function with WhereFrom = 'OUT' means
c     we are done, and all that is left to do is compute
c     output-purpose-only QHR_* things: 
c          QHR_log_QH0  & QHR_log_QHeff  - [scalars]
c          QHR_x_QH0(j) & QHR_x_QHeff(j) - [arrays]
c          QHR_x_Y(j,iY)                 - [matrix]
c     These are all passed back to the calling routine via the common
c     QHR_OutputStuff.
      if (WhereFrom.EQ.'OUT') then

c     QHR_log_QH0   = log of the total rate of ionizing photons
c     QHR_log_QHeff = log of the rate of ionizing photons actually used up in
c                 photoionization (the rest is eaten by dust).
         s_QH0   = 0.
         s_QHeff = 0.
         do j=1,NQHR_base
            aux_AV_tot = AV_tot(j)
            if (AV_tot(j).LT.0.0) aux_AV_tot = 0.

            z = 0.4 * aux_AV_tot * QHR_q_norm
            s_QH0 = s_QH0 + QHR_Q2Lnorm40(j) * x(j) * 10.**z

            z = 0.4 * aux_AV_tot * (QHR_q_norm - 0.785 * QHRbeta_I)
            s_QHeff = s_QHeff + QHR_Q2Lnorm40(j) * x(j) * 10.**z
         end do

         if ((s_QH0.LE.0.).OR.(s_QHeff.LE.0)) then
            write (*,*) '[QHR_CalcQHRchi2] HELP  s_QH0=', s_QH0 , 
     $           's_QHeff=',s_QHeff
            write (*,*) '                  Run aborted :('
            stop
         end if

         QHR_log_QH0   = log10(s_QH0)   + QHR_logLobs_norm + 40.
         QHR_log_QHeff = log10(s_QHeff) + QHR_logLobs_norm + 40.

c     QHR_x_QH0(j)   = fraction of QH0 associated w/population j
c     QHR_x_QHeff(j) = fraction of QHeff associated w/population j
         do j=1,NQHR_base
            aux_AV_tot = AV_tot(j)
            if (AV_tot(j).LT.0.0) aux_AV_tot = 0.

            z = 0.4 * aux_AV_tot * QHR_q_norm
            QHR_x_QH0(j) = (QHR_Q2Lnorm40(j) * x(j) * 10.**z) /s_QH0

            z = 0.4 * aux_AV_tot * (QHR_q_norm - 0.785 * QHRbeta_I)
            QHR_x_QHeff(j) = (QHR_Q2Lnorm40(j) * x(j) * 10.**z) /s_QHeff
         end do

c     QHR_x_Y(j,iY) = fractional contribution to line iY from population j
         do iY=1,NQHR_Ys
            s = 0.
            do j=1,NQHR_base
               aux_AV_tot = AV_tot(j)
               if (AV_tot(j).LT.0.0) aux_AV_tot = 0.
               z = 0.4 * aux_AV_tot * (QHR_q_norm - QHR_qlambda(iY) - 
     $              0.785 * QHRbeta_I)
               s = s + QHR_Q2Lnorm40(j) * x(j) * 10.**z
            end do

            do j=1,NQHR_base
               aux_AV_tot = AV_tot(j)
               if (AV_tot(j).LT.0.0) aux_AV_tot = 0.
               z = 0.4 * aux_AV_tot * (QHR_q_norm - QHR_qlambda(iY) - 
     $              0.785 * QHRbeta_I)
               QHR_x_Y(j,iY) = (QHR_Q2Lnorm40(j) * x(j) * 10.**z) / s
            end do
         end do
      end if
c ***************************************************************************


      end
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR


c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
c QHR                                                                     QHR
c QHR                      End of QHR-Routines Block                      QHR
c QHR                                                                     QHR
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR






c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO
c PHO                                                                     PHO
c PHO                     Start of PHO-Routines Block                     PHO
c PHO                                                                     PHO
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO


c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO
c     Initialize PHO-related things

c     Here I (1) read filter curves, pad and resample to dl = 1 Angs and
c     normalize them; (2) compute transmission-weighted mean lambdas &
c     the filter's curve std dev; (3) compute the q-factor applicable to
c     the filter's mean lambda; (4) resample each base spectrum (after
c     redshifting it to the source's redshift!!) to the same lambdas as
c     the padded + resampled filters & compute the 3D-array
c     PHO_3Mat(j,iY,i_AV), which gives the filter luminosity wrt to that
c     at the normalization wavelength (l_norm) for each population j,
c     filter iY and extinction AV_tot(i_AV) in PHOAV_grid(i_AV) - the
c     AV-effect on PHO_3Mat is only a within-the-filter thing (and thus
c     ~ negligible to 1st order, except when AV grows crazily).

c     The main output of this routine are:
c
c       PHO_MeanLambda(iY)   = integral of lambda * T(lambda) * dlambda
c       PHO_StdDevLambda(iY) = sqrt[<l**2> - <l>**2] = a measure of the
c                              filter's width (just for output purposes)
c       PHO_q_MeanLambda(iY) = q(l) at the MeanLambda of the filter
c       PHO_q_norm = q_norm for PHO-routines
c       NPHO_base  = N_base for PHO-routines
c       PHOdAV & PHOAV_grid(i_AV) = the AV-grid for PHO-work

c       PHO_3Mat(j=1..N_base;,iY=1..NPHO_Ys;i_AV=1...NPHO_AV), defined
c     as the integral, over the filter range, of
c
c         t(l) * b_j(l) * exp[-tau_V * [q(l) - q(l_MeanLambda)] * dl
c
c     where t(l) is the normalized filter transmission curve, b_j(l) is
c     the j-th base-spectrum normalized to l_norm, tauV is the optical
c     depth and q(l) is the usual A(l)/A(V).

c     OBS2: some of the steps, like normalizing base spectra to their
c     flux at l_norm, and reddening-curve definitions are also done
c     later-on in the main routine, but I need to do them here and now
c     for PHO-purposes.
c
c     OBS2: PHO_redshift & NPHO_Ys belong to PHP_DataStuff-common, but,
c     for compactness, are passed as arguments to this routine. (They
c     are not changed, anyway.)

c     Cid@Lagoa - 11/Jan/2011

c     AKI! falta acertar redshifting ....
c???? ! NAO TENHO CERTZEZA DE EM QUE PONTO NORMAZLIAR?? ANTES OU DEPOIS DE APLICAR z???

!AKI! cosmetic change Nl,l <==> NOl_base,Ol_base for similarity with QHRInit???

      SUBROUTINE PHOInitialization(PHO_redshift,NPHO_Ys,i_verbose,N_base
     $     ,NOl_base,Ol_base,l_norm,red_law_option,etc_dir)

c ***************************************************************************
c     Defs
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
      parameter (NPHO_Ys_max=100)
      parameter (NPHO_Nl_max=5000)
!IIB! Changed Nl_max to NAVgrid_max!
      parameter (NAVgrid_max = 1001)

c     Original base things
      real Ol_base(Nl_max)
      real Of_base(Nb_max,Nl_max)
      common/Obase/ Of_base

      character*3 red_law_option
      character*100 arq , etc_dir

c     Auxiliar lambda & transmission arrays for filters ...
      real l0(NPHO_Nl_max) , l1(NPHO_Nl_max) , l2(NPHO_Nl_max)
      real T0(NPHO_Nl_max) , T1(NPHO_Nl_max) , T2(NPHO_Nl_max)
      real l0_aux(NPHO_Nl_max)
      
c     Auxiliar lambda & base-flux arrays for base-spectra....
      real lb0(Nl_max) , lb1(Nl_max)
      real fb0(Nl_max) , fb1(Nl_max) , fb2(Nl_max)
      
c     These arrays/matrices are used only here - don't need to be passed
c     back nowhere!!!
      integer NPHO_Nl_Filter(NPHO_Ys_max)
      real q_lambda(NPHO_Nl_max) ,
     $     PHO_l_Filter(NPHO_Ys_max,NPHO_Nl_max) ,
     $     PHO_T_Filter(NPHO_Ys_max,NPHO_Nl_max)

*PHO* filter name & file name common for PHO-work
      character*10  PHO_name(NPHO_Ys_max)
      character*100 PHO_Filter_file(NPHO_Ys_max)
      common/PHO_NamesStuff/ PHO_name , PHO_Filter_file

*PHO* Base-info-common. 
!IIB! Changed Nl_max to NAVgrid_max!
      common/PHO_BaseStuff/ NPHO_AV , PHOAV_grid(NAVgrid_max) ,
     $     PHOAV_low , PHOAV_upp , PHOdAV , 
     $     PHO_MeanLambda(NPHO_Ys_max) , PHO_StdDevLambda(NPHO_Ys_max) ,
     $     PHO_q_MeanLambda(NPHO_Ys_max) ,
     $     PHO_3Mat(Nb_max,NPHO_Ys_max,NAVgrid_max) ,
     $     PHO_q_norm , NPHO_base
c ***************************************************************************


c ***************************************************************************
c     Resets. OBS: PHOAV_grid(), PHO_MeanLambda(), PHO_StdDevLambda() &
c     PHO_3Mat(,,) already zero'ed in the MAIN routine!
      do i=1,NPHO_Nl_max
         l0(i) = 0.
         l1(i) = 0.
         l2(i) = 0.
         T0(i) = 0.
         T1(i) = 0.
         T2(i) = 0.
         q_lambda(i) = 0.
         l0_aux = 0.
      end do
      do i=1,Nl_max
         lb0(i) = 0.
         lb1(i) = 0.
         fb0(i) = 0.
         fb1(i) = 0.
         fb2(i) = 0.
      end do
      do iY=1,NPHO_Ys_max
         NPHO_Nl_Filter(iY) = 0
         do i=1,NPHO_Nl_max
            PHO_l_Filter(iY,i) = 0.
            PHO_T_Filter(iY,i) = 0.
         end do
      end do

c     Defines PHO_q_norm, passed on to other routines via the
c     PHO_BaseStuff common.
      call Calc_red_term(red_law_option,l_norm,NOl_base,Ol_base,
     $      fb0,PHO_q_norm)

c     Store N_base in NPHO_base, such that it is available through the
c     PHO_BaseStuff common above.
      NPHO_base = N_base

c     Sanity check on NPHO_AV
      if (NPHO_AV.GT.NAVgrid_max) then
         write (*,'(a,2(a,1x,i5))')  '[PHOInitialization] ATT! ' ,
     &        '  NPHO_AV = ', NPHO_AV , ' too large! > ', NAVgrid_max
         stop
      end if

c     Setup PHOAV_grid for PHO_3Mat(j=1..N_base;,iY=1..NPHO_Ys;i_AV=1...NPHO_AV)
c     OBS: NPHO_AV, PHOAV_low & PHOAV_upp were read from arq_config and
c          passed on via common/PHO_BaseStuff.
      PHOdAV = (PHOAV_upp - PHOAV_low) / (NPHO_AV - 1)
      do i_AV=1,NPHO_AV
         PHOAV_grid(i_AV) = PHOAV_low + (i_AV - 1) * PHOdAV
      end do

c     Screen output
      if (i_verbose.GT.0) write (*,'(a,i4,4(2x,a,f6.3))') 
     &     '[PHOInit..] AVgrid: NPHO_AV=',NPHO_AV,
     &     'PHOAV_low=',PHOAV_low , 'PHOAV_upp=',PHOAV_upp ,
     &     'PHOdAV=',PHOdAV , ' PHO_q_norm=',PHO_q_norm
c ***************************************************************************


c ***************************************************************************
c                 Read, pad, resample & normalize filters
c ***************************************************************************
c     Read filter curves, pad them, resample to dl = 1 Angs, and
c     normalize. Also, compute transmission-weighted mean lambdas and
c     the std dev of the filter curve. The output of this block is:
c
c       NPHO_Nl_Filter(iY)    = array with # of points (lambdas) in filter iY
c       PHO_l_Filter(iY,il)   = matrix of lambdas in filter iY
c       PHO_T_Filter(iY,il)   = matrix of transmission coef. in filter iY
c       PHO_MeanLambda(iY)    = integral of lambda * T(lambda) * dlambda
c       PHO_StdDevLambda(iY)  = sqrt[<l**2> - <l>**2] = a measure of the
c                              filter's width (just for output purposes)
c
c     It is assumed that the transmission curves in files
c     PHO_Filter_file(iY) are given in a lambda,T(lambda) format,
c     possibly with **header** lines starting with a # in the 1st
c     column.

      i_last_non_blank_char = 0
      do i=1,99
         if ((etc_dir(i:i+1).EQ.' ').AND.(i_last_non_blank_char.EQ.0))
     &        i_last_non_blank_char = i - 1
      end do

      do iY=1,NPHO_Ys

c     Read filter file PHO_Filter_file(iY), assumed to be in etc_dir.
c     ==> T0(l0)
         arq = etc_dir(1:i_last_non_blank_char) // PHO_Filter_file(iY)
         call ReadSpecWithHeader(arq,NPHO_Nl_max,N0,l0,T0)
         if (i_verbose.GT.0) write (*,9) iY , PHO_Filter_file(iY) ,
     $        PHO_name(iY) , N0 , l0(1) , l0(N0)
 9       format ('[PHOInit...]',i3,2x,a20,2x,a10,2x,
     $        ' N=',i5,2x,f9.2,'==>',f9.2,' Angs')

         
c     Test whether lambda-range is too large to fit in the allocated arrays!
c     (after padding & resampling to 1 Angs below...)
c     ATT: filters may get very wide in the IR!
         ntest = l0(N0) - l0(1) + 10
         if (ntest.GT.NPHO_Nl_max) then
            write (*,*) '[PHOInit] Filter too wide! ' , iY ,
     $           l0(1) , ' ==> ', l0(N0) , NPHO_Nl_max
            write (*,*) '          Run aborted:-('
            stop
         end if

c     Pad filter-curve with a T=zero @ 2 Angs before and 2 Angs after the
c     filter lambda-range. T0(l0) ===> T1(l1)
         N1     = N0 + 2
         l1(1)  = l0(1) - 2.0
         l1(N1) = l0(N0) + 2.0
         T1(1)  = 0.0
         T1(N1) = 0.0
         do i1=2,N1-1
            l1(i1) = l0(i1-1)
            T1(i1) = T0(i1-1)
         end do

c     Resample padded filter to dl = 1 Angs - use Interpol_Spec. T1(l1) ==> T2(l2)
         l2_ini = int(l1(1))
         N2     = int(l1(N1) - l1(1)) + 1
         do i2=1,N2
            l2(i2) = l2_ini + (i2 - 1) * 1.0
         end do
         call Interpol_Spec(N1,l1,T1,N2,l2,T2)

c     Normalize padded+resampled filter & compute weighted-mean-lambda.
c     For checking purposes, computes the std dev of the filter curve,
c     ie., PHO_StdDevLambda(iY) = sqrt[<l**2> - <l>**2], a measure of
c     its width.
         sum_T = 0.
         do i2=1,N2
            sum_T = sum_T  + T2(i2)
         end do
         sum_lT  = 0.
         sum_l2T = 0.
         do i2=1,N2
            T2(i2) = T2(i2) / sum_T
            sum_lT  = sum_lT + l2(i2) * T2(i2)
            sum_l2T = sum_l2T + l2(i2)**2 * T2(i2)
         end do

c     Store things into proper names
         PHO_MeanLambda(iY)   = sum_lT
         PHO_StdDevLambda(iY) = sqrt(abs(sum_l2T - sum_LT**2))

         NPHO_Nl_Filter(iY) = N2
         do i2=1,N2
            PHO_l_Filter(iY,i2) = l2(i2)
            PHO_T_Filter(iY,i2) = T2(i2)
         end do
!         write (*,*) 'test ', iY , PHO_MeanLambda(iY) , l2(1) , l2(N2) ,
!     $        T2(1) , T2(N2) , ' >> ', l2(int(N2/2)) , T2(int(N2/2)) , PHO_StdDevLambda(iY)

c     Sanity check: Mean lambda must be within filter's range!
         if ((PHO_MeanLambda(iY).LT.l2(1)).OR.
     $       (PHO_MeanLambda(iY).GT.l2(N2))) then
            write (*,*) '[PHOInit] Mean filter lambda = ',
     $           PHO_MeanLambda(iY) , ' is outside filter boundaries: '
     $           , l2(1) , ' ==> ' , l2(N2) , ' iY=' , iY
            write (*,*) '          Run aborted:-('
            stop
         end if
      end do
c ***************************************************************************


c ***************************************************************************
c                    Compute PHO_q_MeanLambda(iY)
c ***************************************************************************
c     A few lines above we called Calc_red_term with fb0 as an auxiliar
c     array for red_term. Here We use it to define the corresponding
c     q_lambda array (stored in fb1), from which we interpolate (using
c     calc_Fnorm) to obtain the q-value at the mean lambda of each
c     filter. 
c     The output of this block is PHO_q_MeanLambda(iY).
      do i=1,NOl_base
         fb1(i) = PHO_q_norm - 2.5 * fb0(i)
      end do

      do iY=1,NPHO_Ys
         l_MeanLambda = PHO_MeanLambda(iY)
         call calc_Fnorm(NOl_base,Ol_base,fb1,l_MeanLambda,q_MeanLambda)
         PHO_q_MeanLambda(iY) = q_MeanLambda
      end do
c ***************************************************************************


c ***************************************************************************
c                  Compute PHO_3Mat(j,iY,i_AV) "3D-array"
c ***************************************************************************
c     Normalize base specra, redshift them, & then, for each filter iY,
c     compute the PHO_3Mat(j,iY,i_AV) 3D-array... used in PHO-calculations.

c     Start (outer) loop in j=1...N_base
      do j=1,N_base

c     Normalize j-th base spectrum ==> fb0(lb0)
         Nb0 = NOl_base
         do ib0=1,Nb0
            lb0(ib0) = Ol_base(ib0)
            fb0(ib0) = Of_base(j,ib0)
         end do
         call calc_Fnorm(Nb0,lb0,fb0,l_norm,fbase_norm)
         do ib0=1,Nb0
            fb0(ib0) = fb0(ib0) / fbase_norm
         end do

c     Redshifting normalized base spectrum to PHO_redshift ==>fb1(lb1)
c     ATT! THIS IS NOT RIGHT ..???.. need (1+z)**3 factors & etc... <** AKI FIX IT!!
         write (*,*) '[PHOInit ...] ATT: Fix redshifting <-- AKI!!!'
         Nb1 = Nb0
         do ib1=1,Nb1
!RAFA! Now we are going to do the redshift corrections on the base spectra
            lb1(ib1) = lb0(ib1) * (1. + PHO_redshift) 
            fb1(ib1) = fb0(ib1) / (1. + PHO_redshift)
c            fb1(ib1) = fb0(ib1)
         end do

c     -----------------------------------------------------------------------
c     Start (middle) loop in filters iY=1...NPHO_Ys
         do iY=1,NPHO_Ys

c     Store iY-filter things into aux variables N0, l0 & T0
            N0 = NPHO_Nl_Filter(iY)
            do i0=1,N0
               l0(i0) = PHO_l_Filter(iY,i0)
               T0(i0) = PHO_T_Filter(iY,i0)
            end do

c     Sanity check: Test whether the redshifted base spectrum fb1(lb1)
c     covers the range of filter-iY, otherwise interpolation will fail!
            if ((l0(1).LT.lb1(1)).OR.(l0(N0).GT.lb1(Nb1))) then
               write (*,*) '[PHOInit]' , 
     $              ' Filter overshoots redshifted base spectrum!'
               write (*,*) ' Filter ' , iY , ' : ' , l0(1) , ' ==> ',
     $              l0(N0) ,' |  Base ' , j , ' : ' , lb1(1) , ' ==> ',
     $              lb1(Nb1) , ' |  z = ' , PHO_redshift
               write (*,*) '          Run aborted:-('
               stop            
            end if

c     Resample redshifted j-th base spectrum to 1 Angs in the range of
c     the (already resampled) iY-filter curve: fb1(lb1) ==> fb2(l0)
            call Interpol_Spec(Nb1,lb1,fb1,N0,l0,fb2)

c     Compute array q_lambda = q(l0) = AV(l0) / AV(l_norm) in l0-grid.
c     ==> q_lambda(l0)
c     OBS: red_term(l0) = -0.4 * (q(l0) - q_norm) is stored in aux-variable T1.
c Rafa Al mover lambda de la base hay que tener en cuenta donde se aplica la extincion
c          Hay que hacer la misma correccion.
            do i0=1,N0
               l0_aux(i0) = l0(i0)/(1. + PHO_redshift)
            end do
            call Calc_red_term(red_law_option,l_norm,N0,l0_aux,
     $      T1,q_norm)
            do i0=1,N0
               q_lambda(i0) = q_norm - 2.5 * T1(i0)
            end do

c     .......................................................................
c     Compute PHO_3Mat(j,iY,i_AV), defined as the integral, over the
c     filter range, of
c
c         t(l) * b_j(l) * exp[-tau_V * [q(l) - q(l_MeanLambda)] * dl
c
c     where t(l) is the nomarlized filter transmission curve, b_j(l) is
c     the j-th base-spectrum normalized to l_norm, tauV is the optical
c     depth and q(l) is the usual A(l)/A(V).
c     OBS: Notice that dl = 1 Angs is implicit below.

c     Start (inner) loop in predefined PHOAV_grid. ==> PHO_3Mat(j,iY,i_AV)
            do i_AV=1,NPHO_AV
               tauV = PHOAV_grid(i_AV) * 0.4 / log10(exp(1.))
               s  = 0.
               s_filter = 0.
               do i0=1,N0
c                  s = s + T0(i0) * fb2(i0) *
!RAFA¡ We are going to change the integral to compute AB_magnitude
                  s = s + T0(i0) * fb2(i0) * l0(i0) *
     $                exp(-tauV * (q_lambda(i0) - PHO_q_MeanLambda(iY)))
!RAFA¡
                  s_filter = s_filter + T0(i0) / l0(i0)
               end do
c               PHO_3Mat(j,iY,i_AV) = s
!RAFA¡             
               PHO_3Mat(j,iY,i_AV) = s/s_filter
            end do
c     -----------------------------------------------------------------------

c     close middle loop: goto next filter iY, same base element j
         end do

c     close outer loop: goto next base element
      end do
c ***************************************************************************


c ***************************************************************************
c     Screen output
      if (i_verbose.GT.0) then
         istep = float(NPHO_AV) / float(14)
         ilast = 14 * istep
 11      format ('[PHOInit..]',i3,1x,i3,1x,a10,2(1x,i5),1x,f6.4,
     $        ' | ',a,201(1x,f6.2))
 12      format (157('-'))
 13      format ('[PHOInit..]  j  iY Filter      <l>  sig_l q(<l>) ', 
     &        '| AV_grid=',14(1x,f6.2))
         write (*,12)
         write (*,13) (PHOAV_grid(i_AV),i_AV=1,ilast,istep)
         write (*,12)
         do j=1,N_base
            do iY=1,NPHO_Ys
               write (*,11) j , iY , PHO_name(iY) ,
     $              int(PHO_MeanLambda(iY) +  0.5) , 
     $              int(PHO_StdDevLambda(iY) + 0.5) ,
     $              PHO_q_MeanLambda(iY) , 
     $              'log3Mat=' ,
     $              (log10(PHO_3Mat(j,iY,i_AV)),i_AV=1,ilast,istep)
            end do
            if ((j.NE.N_base).AND.(NPHO_Ys.GT.1)) write (*,*) 
         end do
         write (*,12)
         write (*,*)
      end if
c ***************************************************************************


      end
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO



c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO
!     Computes the PHO-part of chi2...

!     Given x(j) and AV_tot(j) population and extinction vectors from my
!     MEx-STARLIGHT fits, and the stuff initialized in routine
!     PHOInitialization (passed though the common/PHO_BaseStuff), this
!     function computes the predicted filter-luminosities (Y) for the
!     NPHO_Ys photometric bands informed in the arq_ETCinfo and passed
!     on via common/PHO_DataStuff/. The main output is the chi2 of the
!     comparison of these predicted Ys with the observed ones.

!     As in the FIRc-routines, the value of chi2 is = 0 between lower &
!     upper limits, rising parabolically beyond these limits, such that
!     the likelyhood is a flat-top with gaussian wings. The chi2 of each
!     predicted band is scaled such that it can compete with the optical
!     fit (same scheme as for the FIR work).

!     WhereFrom is just a string which tells from where this function
!     was called. When it is = 'OUT', then some quantities are computed
!     for output purposes only - see common/PHO_OutputStuff/.
      
!     A NOTE ABOUT AV < 0: In principle, PHOAV_grid does NOT contain AV
!     < 0 (although that, in principle, could be set in arq_config
!     ...). Hence, there is no entry in PHO_3Mat(j,iY,i_AV)
!     corresponding to AV < 0. However, the effect of AV < 0 (as long as
!     not too large) upon PHO_3Mat(j,iY,i_AV) should be small (the
!     reddening effect in this matrix is small by definition, since it
!     applies to within-the-filter only!). AV < 0 will still affect the
!     10**[-0.4 * AV * (q_norm - q_MeanLambda)] term in the predicted
!     band luminosity! This is **good**, because a AV<0 should lead to
!     large effects in UV-luminosities....

!     Cid@Lagoa - 11/Jan/2011


      FUNCTION PHO_CalcPHOchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)

c ***************************************************************************
c     Definitions
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nl_max = 14000,Nb_max = 300)
      parameter (NPHO_Ys_max=100)
!IIB! Changed Nl_max to NAVgrid_max!
      parameter (NAVgrid_max = 1001)
      real x(Nb_max) , AV_tot(Nb_max)
      character*3 WhereFrom

*PHO* filter name & file name common for PHO-work
      character*10  PHO_name(NPHO_Ys_max)
      character*100 PHO_Filter_file(NPHO_Ys_max)
      common/PHO_NamesStuff/ PHO_name , PHO_Filter_file

*PHO* Base-info-common. 
!IIB! Changed Nl_max to NAVgrid_max!
      common/PHO_BaseStuff/ NPHO_AV , PHOAV_grid(NAVgrid_max) ,
     $     PHOAV_low , PHOAV_upp , PHOdAV , 
     $     PHO_MeanLambda(NPHO_Ys_max) , PHO_StdDevLambda(NPHO_Ys_max) ,
     $     PHO_q_MeanLambda(NPHO_Ys_max) ,
     $     PHO_3Mat(Nb_max,NPHO_Ys_max,NAVgrid_max) ,
     $     PHO_q_norm , NPHO_base

*PHO* Data-info-common for PHO-work
      common/PHO_DataStuff/ PHO_logY_TOT(NPHO_Ys_max) ,
     $     PHO_YFrac2Model(NPHO_Ys_max) , PHO_ErrlogY(NPHO_Ys_max)  ,
     $     PHO_RangelogY(NPHO_Ys_max) , 
     $     PHO_Chi2ScaleFactor(NPHO_Ys_max) ,
     $     PHO_logY_obs(NPHO_Ys_max) ,
     $     PHO_logY_low(NPHO_Ys_max) , PHO_logY_upp(NPHO_Ys_max) ,
     $     PHO_GlobalChi2ScaleFactor , 
     $     PHO_redshift , PHO_logLobs_norm , NPHO_Ys , iPHO_IsPHOcOn

*PHO* Common for predicted PHO-things
      common/PHO_OutputStuff/ PHO_ModlogY(NPHO_Ys_max) , 
     $     PHO_chi2_Y(NPHO_Ys_max) , PHO_x_Y(Nb_max,NPHO_Ys_max)


c     Reset output chi2
      PHO_CalcPHOchi2 = 0.

c     If PHO is not turned on, just set the PHO-chi2 to zero & go back!
c     OBS: In non-PHOc-fits, I may temporarely "turn PHO On" just to
c     call this routine to obtain the model PHO-related quantities...
      if (iPHO_IsPHOcOn.EQ.0) return
c ***************************************************************************


c ***************************************************************************
c     For each filter, compute the predicted line luminosity, which is
c     the sum of luminosities coming from each population, each one with
c     its own total extinction. The log of each of these NPHO_Ys
c     band-luminosities is used to compute the chi2, in the fashion
c     described above.
      s_chi2_Ys = 0.
      do iY=1,NPHO_Ys

c     Add upp contributions of all base elements to filter-iY.
         s = 0.
         do j=1,NPHO_base

c     Find index of entry in PHOAV_grid which best matches AV_tot(j)
c     ATT: AV < 0 will lead to AV_grid = 0 = 1st element in PHOAV_grid.
c          Cannot bluen PHO-models here, but see below!)
            ind_PHOAV = int(0.5 + (AV_tot(j) - PHOAV_low) / PHOdAV) + 1
            if (ind_PHOAV.LT.1) ind_PHOAV = 1
            if (ind_PHOAV.GT.NPHO_AV) ind_PHOAV = NPHO_AV

c     Computes this base-element's addition to the integral s.
c     ATT: Note that here AV<0 will work!! 
c     !!!Rafa: He cambiado z1 = -0.4 etc. .....por z1 = 0.4 etc.
            z1 = 0.4 * AV_tot(j) * (PHO_q_norm - PHO_q_MeanLambda(iY))
            s  = s + x(j) * (10.**z1) * PHO_3Mat(j,iY,ind_PHOAV)
         end do

c     Sanity check
         if (s.LE.0.) then
            write (*,*) '[PHO_CalcPHOchi2] HELP!  s = ',s, ' iY=',iY
            write (*,*) '                  Run aborted :('
            stop
         end if

c     Model logY for filter iY
c         PHO_ModlogY(iY) = log10(s) + PHO_logLobs_norm
c !RAFA¡ We want to compute AB_magnitude 
         PHO_ModlogY(iY) = -2.5*log10(s) -2.41 
     &     - 2.5*PHO_logLobs_norm
     
c     Compute the corresponding contribution to total chi2
c     Added if iPHO_IsPHOcOn.EQ.1 to cope with -1 (predict-only) option!
         chi2_Y = 0.
         if (iPHO_IsPHOcOn.EQ.1) then
            if (PHO_ModlogY(iY).GE.PHO_logY_upp(iY)) then
               z = (PHO_ModlogY(iY) - PHO_logY_upp(iY)) /PHO_ErrlogY(iY)
               chi2_Y = z**2
            else if (PHO_ModlogY(iY).LE.PHO_logY_low(iY)) then
               z = (PHO_ModlogY(iY) - PHO_logY_low(iY)) /PHO_ErrlogY(iY)
               chi2_Y = z**2
            end if
         end if
         PHO_chi2_Y(iY) = chi2_Y * PHO_Chi2ScaleFactor(iY) *chi2_OPTICAL
         s_chi2_Ys = s_chi2_Ys + PHO_chi2_Y(iY) 

!         write (*,*) 'test-AKI PHO_CalcPHOchi2- ', iY , PHO_name(iY) ,
!     $        ' ModlogY = ', PHO_ModlogY(iY) , chi2_Y , PHO_chi2_Y(iY) ,
!     $        PHO_Chi2ScaleFactor(iY) , chi2_OPTICAL , s_chi2_Ys
      end do

c     Here's the result of this function
      PHO_CalcPHOchi2 = s_chi2_Ys
c ***************************************************************************


c ***************************************************************************
c     FINAL OUTPUT: Calling this function with WhereFrom = 'OUT' means
c     we are done, and all that is left to do is compute (for output
c     purposes only) PHO_x_Y(j,iY): the filter-by-filter pop vector
c     matrix, defined as the fractional contribution to filter iY from
c     population j, and passed back to the calling routine via the
c     common PHO_OutputStuff.
      if (WhereFrom.EQ.'OUT') then
         do iY=1,NPHO_Ys
            s = 0.
            do j=1,NPHO_base
               ind_PHOAV = int(0.5 + (AV_tot(j) - PHOAV_low)/PHOdAV) + 1
               if (ind_PHOAV.LT.1) ind_PHOAV = 1
               if (ind_PHOAV.GT.NPHO_AV) ind_PHOAV = NPHO_AV
c Rafa; he vuelto a cambiar z1 = -0.4 etc..... por z1 = 0.4etc....
               z1 = 0.4*AV_tot(j) * (PHO_q_norm - PHO_q_MeanLambda(iY))
               PHO_x_Y(j,iY) = x(j) *(10.**z1)* PHO_3Mat(j,iY,ind_PHOAV)
               s  = s + PHO_x_Y(j,iY)
            end do

            do j=1,NPHO_base
               PHO_x_Y(j,iY) = PHO_x_Y(j,iY) / s
            end do
         end do
      end if
c ***************************************************************************


      end
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO


c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO
c PHO                                                                     PHO
c PHO                      End of PHO-Routines Block                      PHO
c PHO                                                                     PHO
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO






c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c ETC                                                                     ETC
c ETC                     Start of ETC-Routines Block                     ETC
c ETC                                                                     ETC
c ETC                                                                     ETC
c ETC Read_ETC_Info(i_verbose,arq)                                        ETC
c ETC ETC_CalcETCchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)                    ETC
c ETC                                                                     ETC
c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC


c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c     Read general ETC-info routine. Reads a file designed to contain
c     FIR and/or QHR and/or PHO data to be used as extra constraints in
c     the optical fits.
c
c     This routine works like this: First search arq to find the line
c     numbers correspondig to key-string: FIR, QHR and PHO. Then reads
c     the blocks available. Not that the FIR/QHR/PHO info can come in
c     any order. All that maters is that the corresponding key-string is
c     found, in which case the following line contain data in the
c     expected format. If a key string is not found as the 1st entry in
c     any line of arq, no data corresponding data is read. A quite
c     flexible (if not excedingly elegant) coding:-)

c     Except for FIRlogLobs_norm & iFIR_IsFIRcOn, all
c     common/FIR_DataStuff/ variables are defined here.
c
c     Except for QHR_qlambda(), QHR_logLobs_norm & iQHR_IsQHRcOn, all
c     common/QHR_DataStuff/ variables are defined here.


c     Cid@Lagoa - 07/Jan/2011

      SUBROUTINE Read_ETC_Info(i_verbose,arq)

c ***************************************************************************
c     Defs
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (NQHR_Ys_max=10)
      parameter (NPHO_Ys_max=100)
      character*3 char3
      character*100 arq
 
*FIR* Data-info-common for FIR-work
      common/FIR_DataStuff/ FIRlogLobs_norm , 
     $     FIR_logLFIR_TOT , FIR_LFIRFrac2Model ,
     $     FIR_logLFIR_obs , FIR_ErrlogLFIR ,
     $     FIR_logLFIR_low , FIR_logLFIR_upp , FIR_RangelogLFIR ,
     $     FIRChi2ScaleFactor , iFIR_IsFIRcOn

*QHR* Data-info-common for QHR-work
      common/QHR_DataStuff/ QHR_lambda(NQHR_Ys_max) ,
     $     QHR_frecomb(NQHR_Ys_max)  ,
     $     QHR_logY_TOT(NQHR_Ys_max) , QHR_YFrac2Model(NQHR_Ys_max) ,
     $     QHR_ErrlogY(NQHR_Ys_max)  , QHR_RangelogY(NQHR_Ys_max) ,
     $     QHR_Chi2ScaleFactor(NQHR_Ys_max) , 
     $     QHR_logY_obs(NQHR_Ys_max) ,
     $     QHR_logY_upp(NQHR_Ys_max) , QHR_logY_low(NQHR_Ys_max)  ,
     $     QHR_qlambda(NQHR_Ys_max) ,
     $     QHR_GlobalChi2ScaleFactor ,
     $     QHR_logLobs_norm , 
     $     NQHR_Ys , iQHR_IsQHRcOn

*ELR* Data-info for ELR-work - goes together w/QHR-things!
      common/ELR_DataStuff/ ELR_lambdaA , ELR_lambdaB , 
     $     ELR_logRint , ELR_logRobs , 
     $     ELR_ErrlogR , ELR_logR_low , ELR_logR_upp , 
     $     ELR_RangelogR , ELR_Chi2ScaleFactor , 
     $     iY_ELR_A , iY_ELR_B , IsELROn , 
     $     IsELRInputInconsistentWarningFlag

*PHO* filter name & file name common for PHO-work
      character*10  PHO_name(NPHO_Ys_max)
      character*100 PHO_Filter_file(NPHO_Ys_max)
      common/PHO_NamesStuff/ PHO_name , PHO_Filter_file

*PHO* Data-info-common for PHO-work
      common/PHO_DataStuff/ PHO_logY_TOT(NPHO_Ys_max) ,
     $     PHO_YFrac2Model(NPHO_Ys_max) , PHO_ErrlogY(NPHO_Ys_max)  ,
     $     PHO_RangelogY(NPHO_Ys_max) , 
     $     PHO_Chi2ScaleFactor(NPHO_Ys_max) ,
     $     PHO_logY_obs(NPHO_Ys_max) ,
     $     PHO_logY_low(NPHO_Ys_max) , PHO_logY_upp(NPHO_Ys_max) ,
     $     PHO_GlobalChi2ScaleFactor , 
     $     PHO_redshift , PHO_logLobs_norm , NPHO_Ys , iPHO_IsPHOcOn

c     OBS: All resets were made in the MAIN routine:-)
c ***************************************************************************


c ***************************************************************************
c     Finds where (and if) FIR, QHR & PHO information are located in
c     arq, finding the line numbers for the 3 key strings FIR/QHR/PHO.
c     OBS1: The strings must be the 1st entries in their respective lines,
c           as for all info read realted to each FIR/QHR/PHO block.
c     OBS2: Assumes arq does not have > 1000 lines!
c     OBS3: The file MUST NOT CONTAIN BLANCK LINES!!!
      iFIR = 0
      iQHR = 0
      iPHO = 0
      open (unit=1,file=arq,status='old')
      do i=1,1000
         read (1,*,end=401) char3
         if (char3.EQ.'FIR') iFIR = i
         if (char3.EQ.'QHR') iQHR = i
         if (char3.EQ.'PHO') iPHO = i
      end do
 401  close(1)

c     Screen check
      if (i_verbose.GT.0) then
         write (*,'(3(a,i3))') '[Read_ETC_Info] iFIR = ' , iFIR , 
     $        '  iQHR  = ', iQHR , '  iPHO = ', iPHO
      end if

c     A series of sanity checks...
      if ((iFIR.LE.0).AND.(iFIR_IsFIRcOn.EQ.1)) then
         write (*,'(a,a,a)') '[Read_ETC_Info] You asked for a FIRc fit',
     $        ' but provided no FIR-data in ', arq 
         write (*,'(a)') '        Run aborted:-('
         stop
      end if

      if ((iQHR.LE.0).AND.(iQHR_IsQHRcOn.EQ.1)) then
         write (*,'(a,a,a)') '[Read_ETC_Info] You asked for a QHRc fit',
     $        ' but provided no QHR-data in ', arq 
         write (*,'(a)') '        Run aborted:-('
         stop
      end if

      if ((iPHO.LE.0).AND.(iPHO_IsPHOcOn.EQ.1)) then
         write (*,'(a,a,a)') '[Read_ETC_Info] You asked for a PHOc fit',
     $        ' but provided no PHO-data in ', arq 
         write (*,'(a)') '        Run aborted:-('
         stop
      end if

      if ((iFIR.GT.0).AND.(iFIR_IsFIRcOn.NE.1)) then
         write (*,'(a,a)') '[Read_ETC_Info] ***WARNING*** ', 
     $    'FIR-data will be IGNORED since you asked for a non-FIRc-fit!'
         write (*,'(a,a)') '[Read_ETC_Info] ... arq = ', arq
      end if

      if ((iQHR.GT.0).AND.(iQHR_IsQHRcOn.EQ.0)) then
         write (*,'(a,a)') '[Read_ETC_Info] ***WARNING*** ', 
     $    'QHR-data will be IGNORED since you asked for a non-QHRc-fit!'
         write (*,'(a,a)') '[Read_ETC_Info] ... arq = ', arq
      end if
      
      if ((iPHO.GT.0).AND.(iPHO_IsPHOcOn.EQ.0)) then
         write (*,'(a,a)') '[Read_ETC_Info] ***WARNING*** ', 
     $    'PHO-data will be IGNORED since you asked for a non-PHOc-fit!'
         write (*,'(a,a)') '[Read_ETC_Info] ... arq = ', arq
      end if

      if ((iQHR.GT.0).AND.(iQHR_IsQHRcOn.EQ.-1)) then
         write (*,'(a,a)') '[Read_ETC_Info] ATT: ', 
     $    'QHR-info will only be used to PREDICT QHR-stuff!'
         write (*,'(a,a)') '[Read_ETC_Info]      arq = ', arq
      end if

      if ((iPHO.GT.0).AND.(iPHO_IsPHOcOn.EQ.-1)) then
         write (*,'(a,a)') '[Read_ETC_Info] ATT: ', 
     $    'PHO-info will only be used to PREDICT PHO-stuff!'
         write (*,'(a,a)') '[Read_ETC_Info]      arq = ', arq
      end if
c ***************************************************************************


c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR
c     Read FIR-input - the lines following the "FIR" key string in arq
      if ((iFIR.GT.0).AND.(iFIR_IsFIRcOn.EQ.1)) then

         write (*,'(a,a)') '[Read_ETC_Info] ... Reading FIR-info from '
     $        ,arq
         open (unit=1,file=arq,status='old')
         do i=1,iFIR
            read (1,*)
         end do 
         read (1,*) FIR_logLFIR_TOT
         read (1,*) FIR_LFIRFrac2Model
         read (1,*) FIR_ErrlogLFIR
         read (1,*) FIR_RangelogLFIR
         read (1,*) FIRChi2ScaleFactor
         close(1)

c     As the name implies, FIR_logLFIR_TOT is the total (say, IRAS)
c     FIR-luminosity, but the user may want to use only a FRACTION of it
c     as a constraint (say, because the IRAS data covers a much wider
c     region than the optical spectrum, or Spitzer images show that most
c     of the FIR comes from a companion not included in the optical
c     slit...). This is what the entry FIR_LFIRFrac2Model is for. The
c     fit is carried out considering an "observed" FIR luminosity
c     10**FIR_logLFIR_obs which is a fraction FIR_LFIRFrac2Model of
c     FIR_logLFIR_TOT.
         FIR_logLFIR_obs = FIR_logLFIR_TOT + log10(FIR_LFIRFrac2Model)

c     Define upper & lower soft-limits for the FIR luminosity
c===> Natalia's suggestion: Range >/< 0 makes obs = upper/lower limit!
         if (FIR_RangelogLFIR.LT.0.) then
            FIR_logLFIR_low = FIR_logLFIR_obs
            FIR_logLFIR_upp = FIR_logLFIR_obs + abs(FIR_RangelogLFIR)
         else
            FIR_logLFIR_upp = FIR_logLFIR_obs
            FIR_logLFIR_low = FIR_logLFIR_obs - FIR_RangelogLFIR
         end if

c    Screen check ...
         if (i_verbose.GT.0) then
            write (*,'(a,f7.3)') '[FIR]>> log L_FIR_TOT = ' , 
     &           FIR_logLFIR_TOT 
            write (*,'(a,f7.3)') '[FIR]>> Frac to Model = ' , 
     &           FIR_LFIRFrac2Model
            write (*,'(a,f7.3)') '[FIR]>> log L_FIR_obs = ' , 
     &           FIR_logLFIR_obs 
            write (*,'(a,f7.3)') '[FIR]>> Error in log L_FIR = +/-' , 
     &           FIR_ErrlogLFIR
            write (*,'(a,f8.2)') '[FIR]>> Range in log L_FIR = ' , 
     &           FIR_RangelogLFIR
            write (*,'(2(a,f7.3))') '[FIR]>> log L_FIR range from ' , 
     &           FIR_logLFIR_low , '   --> ' , FIR_logLFIR_upp
            write (*,'(a,f7.3,a)')'[FIR]>> Scale factor for chi2_FIR = '
     $           ,FIRChi2ScaleFactor , ' times chi2_OPTICAL'
         end if

c    Basic sanity checks on input FIR-info...
         if ((FIR_logLFIR_TOT.LT.0.0).OR.
     $       (FIR_logLFIR_obs.LT.0.0).OR.
     $       (FIR_LFIRFrac2Model.LE.0.0).OR.
     $       (FIR_LFIRFrac2Model.GT.1.0).OR.
     $       (FIR_ErrlogLFIR.LE.0.0).OR.
     $       (FIRChi2ScaleFactor.LT.0.0)) then
            write (*,'(a)') 
     $           '[Read_ETC_Info] Something wrong with FIR the input!!'
            stop
         end if
      end if
c FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR-FIR


c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR
c     Read QHR-input - the lines following the "QHR" key string in arq

c     QHR-info:

c     NQHR_s = number os lines used
c     QHR_GlobalChi2ScaleFactor = each line has a weight specified by
c         QHR_Chi2ScaleFactor(iY), read below. The sum of the chi2's of
c         all lines is expected to have a global weight of
c         QHR_GlobalChi2ScaleFactor wrt to the optical fit. This feature
c         was introduced quite late in the development, and is expected
c         to be useful when more than 1 line is used.

c     ELR-info: IsELROn , ELR_lambdaA & ELR_lambdaB , ELR_ErrlogR ,
c         ELR_RangelogR , ELR_Chi2ScaleFactor.
c
c     IsELROn = 0/1 = N/Y decide if 2 of the lines will be used to
c         perform a EmissionLineRatio-constraint (eg Hbeta/Hgamma,
c         Halpha/Hbeta).
c     ELR_lambdaA , ELR_lambdaB = the lambdas of the lines. There
c         recobination factors and fluxes MUST be given below in the
c         QHR-entries!
c     ELR_ErrlogR , ELR_RangelogR , ELR_Chi2ScaleFactor error, range &
c         chi2-scaling for the logR fit. 
c
c         This info is used to predict L_A / L_B and compare it to the
c         observed value.

c     QHR_lambda(iY) is the lambda in Angs of the emission line.
c         QHR_frecomb(iY) is the fraction of ionizing photons which end
c         up producing an iY-line-photon, such that the iY-line
c         luminosity can be computed as f_Inoue * h nu(iY) /
c         f_recomb(iY)] * QH, where f_Inoue = accounts for the fact that
c         some of the ionizing radiation is eaten up by dust [and thus
c         is not available for photoionization; see below].

c     The next input info follows the same logic as for FIR-work (except
c     that now we may have more than one constraint if NQHR_Ys > 1):
c     QHR_logY_TOT(iY) = total Y-line-luminosity, in Lsun.
c     QHR_YFrac2Model(iY) = fraction of Y-luminosity to be modeled.
c     QHR_ErrlogY(iY) = error in logY
c     QHR_RangelogY(iY) = Range in logY
c     QHR_Chi2ScaleFactor(iY) = weight of this logY-value compared to
c     the optical fit & scaled by QHR_GlobalChi2ScaleFactor 
c ==> QHR_qlambda(iY) = A(lambda(iY)) / AV - this is defined elsewhere!


c     ATT: The "abs" in the if below is to allow for reading QHR-info
c     also when iQHR_IsQHRcOn = -1, ie, in the predict-only option. (We
c     need the lines info in order to predict them!)
      if ((iQHR.GT.0).AND.(abs(iQHR_IsQHRcOn).EQ.1)) then

         write (*,'(a,a)') '[Read_ETC_Info] ... Reading QHR-info from '
     $        ,arq
         open (unit=1,file=arq,status='old')
         do i=1,iQHR
            read (1,*)
         end do 

c     Read & check number of lines to be used
         read (1,*) NQHR_Ys
         if (NQHR_Ys.GT.NQHR_Ys_max) then
            write (*,*)  
     $           '[Read_ETC_Info] OOPS! Too many QHR-Ys: ' , NQHR_Ys , 
     $           '! Maximum is = ', NQHR_Ys_max
            write (*,*) ' Run aborted :-('
            stop
         end if

c     Read ELR info: On/Off, the 2 lines, error, range & chi2-scaling
         read (1,*) IsELROn , ELR_lambdaA , ELR_lambdaB , ELR_ErrlogR ,
     $        ELR_RangelogR , ELR_Chi2ScaleFactor

c     Sanity check: ELR requires >= 2 lines!  Also, for the sake of a
c     logical output, reset to 0 ELR things if they are not used (ie.,
c     when IsELROn = 0). Finaly, make sure ELR will have chi2 = 0 in
c     the predict-QHR-only-mode!
         if ((IsELROn.EQ.1).AND.(NQHR_Ys.LT.2).AND.
     $        (ELR_Chi2ScaleFactor.GT.0.)) then
            write (*,*)  
     $           '[Read_ETC_Info] OOPS! ELR needs >= 2 lines,' , 
     $           'but you have only ' , NQHR_Ys , ' !!'
            write (*,*) ' Run aborted :-('
            stop
         end if
         if (IsELROn.EQ.0) then
            ELR_lambdaA = 0.
            ELR_lambdaB = 0.
            ELR_ErrlogR = 0.
            ELR_RangelogR = 0.
            ELR_Chi2ScaleFactor = 0.
         end if
         if (iQHR_IsQHRcOn.EQ.-1) ELR_Chi2ScaleFactor = 0.

c     Global Chi2 Scaling factor for QHR
         read (1,*) QHR_GlobalChi2ScaleFactor

c     Read data for each QHRelated emission line
         do iY=1,NQHR_Ys
            read (1,*) QHR_lambda(iY) , QHR_frecomb(iY) ,
     $           QHR_logY_TOT(iY) , QHR_YFrac2Model(iY) ,
     $           QHR_ErrlogY(iY) , QHR_RangelogY(iY) ,
     $           QHR_Chi2ScaleFactor(iY) 

c    As the name implies, QHR_logY() is the total line-luminosity, but
c    the user may want to use only a FRACTION of it as a constraint
c    (say, because of a potential AGN-contribution or due to aperture
c    issues). This is what the entry QHR_YFrac2Model() is for. The fit
c    is carried out considering "observed" line luminosities of
c    10**QHR_logY_obs which are a fraction QHR_YFrac2Model of
c    10**QHR_logY_TOT.
            QHR_logY_obs(iY) = QHR_logY_TOT(iY) +
     $           log10(QHR_YFrac2Model(iY))

c    Define upper & lower soft-limits for the Y (log) luminosity
c===> Natalia's suggestion: Range >/< 0 makes obs = upper/lower limit!
         if (QHR_RangelogY(iY).LT.0.) then
            QHR_logY_low(iY) = QHR_logY_obs(iY)
            QHR_logY_upp(iY) = QHR_logY_obs(iY) + abs(QHR_RangelogY(iY))
         else
            QHR_logY_upp(iY) = QHR_logY_obs(iY)
            QHR_logY_low(iY) = QHR_logY_obs(iY) - QHR_RangelogY(iY)
         end if

c     For the sake of a "rational output", reset Chi2ScaleFactor to zero
c     when QHR-stuff is only going to be predicted (not fitted)!
            if (iQHR_IsQHRcOn.EQ.-1) QHR_Chi2ScaleFactor(iY) = 0.
         end do


c     .......................................................................
c     Find ELR lines lambda_A & lambda_B in QHR-line-list, define
c     intrinsic & observed line ratios, lower & upper values for the
c     ratio R. Also, do some tests on the way.
         if (IsELROn.EQ.1) then
            iY_ELR_A = 0
            iY_ELR_B = 0
            do iY=1,NQHR_Ys
               if (int(QHR_lambda(iY)).EQ.int(ELR_lambdaA)) iY_ELR_A=iY
               if (int(QHR_lambda(iY)).EQ.int(ELR_lambdaB)) iY_ELR_B=iY
            end do

c     Check if lambda_A & lambda_B were found among the QHR_lambda's...
            if ((iY_ELR_A.EQ.0).OR.(iY_ELR_B.EQ.0)) then
               write (*,'(a)') 
     $              '[Read_ETC_Info] Error with the ELR+QHR input!!'
               write (*,*) ' Could not find ELR lines ' , ELR_lambdaA ,
     $              ' and ' , ELR_lambdaB , ' in the QHR list...'
               write (*,*) ' Run aborted :-('
               stop
            end if

c     Check if there are actual line luminosities for lines A & B ...
            if ((ELR_Chi2ScaleFactor.GT.0.).AND.
     $          ((QHR_logY_TOT(iY_ELR_A).LE.0.).OR.
     $           (QHR_logY_TOT(iY_ELR_B).LT.0.))) then
               write (*,'(a)') 
     $              '[Read_ETC_Info] Error with the ELR+QHR input!!'
               write (*,*) ' No decent luminosities for  ELR lines ' ,
     $              ELR_lambdaA ,' and ' , ELR_lambdaB , ' !'
               write (*,*) ' Run aborted :-('
               stop
            end if

c     Intrinsic (dust-free) L_A/L_B line ratio (computed from lambda's &
c     frecomb's given in QHR-line-list) & its log
            ELR_Rint = (QHR_lambda(iY_ELR_B) * QHR_frecomb(iY_ELR_B)) /
     $                 (QHR_lambda(iY_ELR_A) * QHR_frecomb(iY_ELR_A))
            ELR_logRint = log10(ELR_Rint)

c     Observed L_A/L_B line ratio & its log
            ELR_logRobs = QHR_logY_TOT(iY_ELR_A) -QHR_logY_TOT(iY_ELR_B)
            ELR_Robs = 10.**ELR_logRobs

!ATT! In my 1st trial runs (@IAA/Jan/2011) I faced the problem that my
!     ETCinfo files occasionally implied an AVneb(A,B) < 0! (For
!     instance, Haalpha/Hbeta < 2.97 or Hbeta/Hgamma < 2.146). Even
!     though the associated errors were large, the runs were aborted in
!     the previous version of the code (a tmp elr4.for). Here I,
!     instead, set Robs to Rint plus an epsilon (for numerical
!     safety/paranoia) and turn on the flag IsELRInputInconsistent, to
!     remind the reader, on output, that he screwed up with the input
!     info! BUT notice that this is only done if ELR_Chi2ScaleFactor >
!     0, ie, if the ELR is in fact being used seriously!

c     ELR-sanity check: Robs must be > Rint if lambda_A > lambda_B & vice-versa!
c     If not, MODIFY Robs & WARN!
            IsELRInputInconsistentWarningFlag = 0
            if (ELR_Chi2ScaleFactor.GT.0.) then
               if (((ELR_lambdaA.GT.ELR_lambdaB).AND.
     $              (ELR_Robs.LT.ELR_Rint)).OR.
     $              ((ELR_lambdaA.LT.ELR_lambdaB).AND.
     $              (ELR_Robs.GT.ELR_Rint))) then
                  write (*,'(a,a)') 
     $                 '[Read_ETC_Info] Error with the ELR+QHR input:' , 
     $                 ' Inconsistent line ratios !!! '
                  write (*,*) ELR_lambdaA , ELR_lambdaB , ' Rint= ',
     $                 ELR_Rint ,' but Robs=', ELR_Robs
                  write (*,*) ' I am turning Robs to slightly > than ' ,
     $                 'Rint ...and turning on a FIX-IT flag!'

                  ELR_Robs    = ELR_Rint + 1.e-3
                  ELR_logRobs = ELR_logRint + 1.e-3
                  IsELRInputInconsistentWarningFlag = 1
               end if
            end if

c     Define upper & lower soft-limits for logR
c===> Natalia's suggestion: Range >/< 0 makes obs = upper/lower limit!
            if (ELR_RangelogR.LT.0.) then
               ELR_logR_low = ELR_logRobs
               ELR_logR_upp = ELR_logRobs + abs(ELR_RangelogR)
            else
               ELR_logR_upp = ELR_logRobs
               ELR_logR_low = ELR_logRobs - ELR_RangelogR
            end if
         end if
c     .......................................................................


c     Re-scale QHR_Chi2ScaleFactor(iY) & ELR_Chi2ScaleFactor such that
c     their sum equals the desired QHR_GlobalChi2ScaleFactor.
c     Notice protection for sum = s = 0! [Cid@IAA - 28/Jan/2011]
         if (iQHR_IsQHRcOn.EQ.1)  then
            s = 0.
            if (IsELROn.EQ.1) s = ELR_Chi2ScaleFactor 
            do iY=1,NQHR_Ys
               s = s + QHR_Chi2ScaleFactor(iY) 
            end do
            if (s.GT.0.) then
               if (IsELROn.EQ.1) ELR_Chi2ScaleFactor =
     $              ELR_Chi2ScaleFactor * QHR_GlobalChi2ScaleFactor / s
               do iY=1,NQHR_Ys
                  QHR_Chi2ScaleFactor(iY) = QHR_Chi2ScaleFactor(iY) *
     $                 QHR_GlobalChi2ScaleFactor / s
               end do
            end if
         end if 


c     Screen check ... (ELH & QHR)
         if (i_verbose.GT.0) then
            write (*,'(a,i3)') '[Read_ETC_Info]>> NQHR_Ys = ', NQHR_Ys
            write (*,'(a,f9.3)') 
     $           '[Read_ETC_Info]>> QHR_GlobalChi2ScaleFactor = ',
     $           QHR_GlobalChi2ScaleFactor
            if ((IsELROn.EQ.1).AND.(ELR_Chi2ScaleFactor.GT.0.0)) then
               write (*,'(2(a,i5),2(a,f9.4))') '      ==> [ELR]>> log' ,
     $              int(ELR_lambdaA) , ' /', int(ELR_lambdaB) , ' =' ,
     $              ELR_Robs , ' while intrinsic ratio is =' ,
     $              ELR_Rint
               write (*,'(3(a,f9.4))') '      ==> [ELR]>> ErrlogR =',
     $              ELR_ErrlogR , '   RangelogR=', RangelogR ,
     $              '    Chi2ScaleFactor=',ELR_Chi2ScaleFactor
            end if
            do iY=1,NQHR_Ys
               write (*,'(a,i3,12(1x,a,f9.3))') 
     $              '      ==> [QHR]>> iY=', iY ,
     $              'lambda='     , QHR_lambda(iY) , 
     $              'frecomb= '   , QHR_frecomb(iY) ,
     $              'logY='       , QHR_logY_TOT(iY) , 
     $              'Frac2Model=' , QHR_YFrac2Model(iY) ,
     $              'ErrlogY='    , QHR_ErrlogY(iY) , 
     $              'Range='      , QHR_RangelogY(iY) ,
     $              'Chi2Scale='  , QHR_Chi2ScaleFactor(iY) 
            end do
         end if

c     Basic sanity checks on input QHR-info...
         do iY=1,NQHR_Ys
            if ((QHR_lambda(iY).LE.0.).OR.
     $           (QHR_frecomb(iY).LE.0.).OR.
     $           ((QHR_logY_TOT(iY).LE.0.).AND.
     $            (QHR_Chi2ScaleFactor(iY).GT.0.)).OR.
     $           ((QHR_YFrac2Model(iY).LE.0.).AND.
     $            (QHR_Chi2ScaleFactor(iY).GT.0.)).OR.
     $           ((QHR_ErrlogY(iY).LE.0.).AND.
     $            (QHR_Chi2ScaleFactor(iY).GT.0.)).OR.
     $           (QHR_Chi2ScaleFactor(iY).LT.0.)) then
               write (*,'(a)')  
     $           '[Read_ETC_Info] Something wrong with QHR the input!!'
               stop
            end if
         end do

         close(1)
      end if
c QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR-QHR


c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO
c     Read PHO-input - the lines following the "PHO" key string in arq
c     ATT: The "abs" in the if below is to allow for reading PHO-info
c     also when iPHO_IsPHOcOn = -1, ie, in the predict-only option. (We
c     need the filters info in order to predict them!)
      if ((iPHO.GT.0).AND.(abs(iPHO_IsPHOcOn).EQ.1)) then

         write (*,'(a,a)') '[Read_ETC_Info] ... Reading PHO-info from '
     $        ,arq
         open (unit=1,file=arq,status='old')
         do i=1,iPHO
            read (1,*)
         end do 
         read (1,*) PHO_redshift

c     Read & check number of bands to be used
         read (1,*) NPHO_Ys
         if (NPHO_Ys.GT.NPHO_Ys_max) then
            write (*,*)  
     $           '[Read_ETC_Info] OOPS! Too many PHO-Ys: ' , NPHO_Ys , 
     $           '! Maximum is = ', NPHO_Ys_max
            stop
         end if

c     Global Chi2 Scaling factor for PHO
         read (1,*) PHO_GlobalChi2ScaleFactor

c     Read data for each PHOtometric band
         do iY=1,NPHO_Ys
            read (1,*) PHO_name(iY) , PHO_Filter_file(iY) , 
     $           PHO_logY_TOT(iY) , PHO_YFrac2Model(iY) ,
     $           PHO_ErrlogY(iY) , PHO_RangelogY(iY) ,
     $           PHO_Chi2ScaleFactor(iY) 

c    As the name implies, PHO_logY() is the total band-luminosity, but
c    the user may want to use only a FRACTION of it as a constraint
c    (say, due to aperture issues). This is what the entry
c    PHO_YFrac2Model() is for. The fit is carried out considering
c    "observed" band luminosities of 10**PHO_logY_obs which are a
c    fraction PHO_YFrac2Model of 10**PHO_logY_TOT.
            PHO_logY_obs(iY) = PHO_logY_TOT(iY) +
     $           log10(PHO_YFrac2Model(iY))

c    Define upper & lower soft-limits for the Y (log) luminosity
c===> Natalia's suggestion: Range >/< 0 makes obs = upper/lower limit!
         if (PHO_RangelogY(iY).LT.0.) then
            PHO_logY_low(iY) = PHO_logY_obs(iY)
            PHO_logY_upp(iY) = PHO_logY_obs(iY) + abs(PHO_RangelogY(iY))
         else
            PHO_logY_upp(iY) = PHO_logY_obs(iY)
            PHO_logY_low(iY) = PHO_logY_obs(iY) - PHO_RangelogY(iY)
         end if

c     For the sake of a "rational output", reset Chi2ScaleFactor to zero
c     when PHO-stuff is only going to be predicted (not fitted)!
            if (iPHO_IsPHOcOn.EQ.-1) PHO_Chi2ScaleFactor(iY) = 0.
         end do

c     Re-scale PHO_Chi2ScaleFactor(iY) such that their sum equals the
c     desired PHO_GlobalChi2ScaleFactor
c     Notice protection for sum = s = 0! [Cid@IAA - 28/Jan/2011]
         if (iPHO_IsPHOcOn.EQ.1)  then
            s  = 0.
            do iY=1,NPHO_Ys
               s = s + PHO_Chi2ScaleFactor(iY) 
            end do
            if (s.GT.0.) then
               do iY=1,NPHO_Ys
                  PHO_Chi2ScaleFactor(iY) = PHO_Chi2ScaleFactor(iY) *
     $                 PHO_GlobalChi2ScaleFactor / s
               end do
            end if
         end if

c     Screen check ...
         if (i_verbose.GT.0) then
            write (*,'(a,i3)')'[Read_ETC_Info]>> NPHO_Ys = ', NPHO_Ys
            write (*,'(a,f9.3)') 
     $           '[Read_ETC_Info]>> PHO_GlobalChi2ScaleFactor = ',
     $           PHO_GlobalChi2ScaleFactor
            do iY=1,NPHO_Ys
               write (*,'(a,i3,2(1x,a,1x,a10),12(1x,a,f8.3))') 
     $              '[Read_ETC_Info]>> iY=', iY ,
     $              'name='       , PHO_name(iY) , 
     $              'Filter_file=', PHO_Filter_file(iY) , 
     $              'logY='       , PHO_logY_TOT(iY) , 
     $              'Frac2Model=' , PHO_YFrac2Model(iY) ,
     $              'ErrlogY='    , PHO_ErrlogY(iY) , 
     $              'Range='      , PHO_RangelogY(iY) ,
     $              'Chi2Scale='  , PHO_Chi2ScaleFactor(iY) 
            end do
         end if

c     Basic sanity checks on input PHO-info...
c         do iY=1,NPHO_Ys
c            if ((PHO_logY_TOT(iY).LE.0.).OR.
c     $           (PHO_YFrac2Model(iY).LE.0.).OR.
c     $           (PHO_ErrlogY(iY).LE.0.).OR.
c     $           (PHO_Chi2ScaleFactor(iY).LT.0.)) then
c               write (*,'(a)')  
c     $           '[Read_ETC_Info] Something wrong with PHO the input!!'
c               stop
c            end if
c         end do
c
         close(1)
      end if

c     Another basic sanity checks on input PHO-info...
      if ((NPHO_Ys.LE.0).AND.(abs(iPHO_IsPHOcOn).EQ.1)) then
         write (*,'(a,a,i4,a,a)')  
     $        '[Read_ETC_Info] You wanna use/predict PHO-things,',
     $        ' but says NPHO_Ys = ',  NPHO_Ys, 
     $        ' !! ==> Check info in' , arq 
         write (*,'(a)') '        Run aborted:-('
         stop
      end if
c PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO-PHO


      end
c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC



c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c     "Intermediary-function" which calls others to get FIR, QHR & PHO
c     chi2's and adds them as a result. Also, store them in
c     common/ETC_Chi2s/ for output purposes.
c
c     Cid@Lagoa - 07/Jan/2011

      FUNCTION ETC_CalcETCchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)

c     Definitions
      implicit real (a-h,k-m,o-z)
      implicit integer (i,j,n)
      parameter (Nb_max = 300, NQHR_Ys_max=10)
      real x(Nb_max) , AV_tot(Nb_max)
      character*3 WhereFrom
      common/ETC_Chi2s/ chi2_FIR , chi2_QHR , chi2_PHO

      chi2_FIR = FIR_CalcFIRchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)
      chi2_QHR = QHR_CalcQHRchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)
      chi2_PHO = PHO_CalcPHOchi2(x,AV_tot,chi2_OPTICAL,WhereFrom)
      ETC_CalcETCchi2 = chi2_FIR + chi2_QHR + chi2_PHO

      return
      end
c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC


c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c ETC                                                                     ETC
c ETC                      End of ETC-Routines Block                      ETC
c ETC                                                                     ETC
c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC
c ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC-ETC


