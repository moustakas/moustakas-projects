pro calibrations_sfr, atlasdust, atlasnodust, nfgsdust, nfgsnodust, agn=agn, write=write
; jm05mar25uofa

; collect all the various star-formation rate calibrations here     
    
    path = atlas_path(/projects)+'sfrs/'

    if keyword_set(agn) then suffix = '_agn' else suffix = ''
    if (n_elements(atlasdust) eq 0L) then atlasdust = read_atlas_sfrs_sample(atlasnodust=atlasnodust,agn=agn)
    if (n_elements(nfgsdust) eq 0L) then nfgsdust = read_nfgs_sfrs_sample(nfgsnodust=nfgsnodust,agn=agn)
    
    if n_elements(snrcut) eq 0L then snrcut = 3.0

    Uinfo = im_filterspecs(filterlist='bessell_U.par')

    lsun = 3.826D33             ; [erg/s]
    haconst = 7.9D-42           ; K98 conversion L(IR) --> SFR(IR)
    loghaconst = alog10(haconst)
    hgsnrcut = 7.0
    hghbewcut = 10.0
    hbewcut = 1.0
    hahb_k04 = 2.85
    hahb = 2.86
;   hahb = return_tbalmer(/hahb)
    kluband = k_lambda(Uinfo.weff,/charlot)
    klhb = k_lambda(4861,/odonnell)
    kloii = k_lambda(3727,/odonnell)
    haoii = 1.0                 ; assumed intrinsic ratio!
    mbolsun = 4.74         ; bolometric absolute solar magnitude [mag]

    lusun = lsun*10^(-0.4*(uinfo.solarmags-mbolsun))

    dust_fraction = 0.5
    
    hbsfr_lumcut = 8.7
    hbsfr_coeff_lolum = 2.46D-41
    hbsfr_coeff_hilum = [2.45D-43,0.23]
;   hbsfr_coeff = [-0.32872866,0.21257953]

    Usfr_obs_coeff     = 1.4D-43 ; 3.6D-44
    Usfr_obs_coeff_err = 1.1D-43 ; 2.7D-44

    oiisfr_lumcut = 8.7
    oiisfr_coeff     = [4.2D-45,0.37]
    oiisfr_coeff_err = [0.0,0.0]
    
; initialize the output data structures
    
    atlas = {$
      galaxy:                  '', $
      hgsnrcut:          hgsnrcut, $ ; applied to Hb/Hg calibrations
      hbewcut:            hbewcut, $ ; applied to SFR(Hb) calibrations
      hghbewcut:        hghbewcut, $ ; applied to Hb/Hg calibrations
      sfr_ha:              -999.0, $ ; SFR(Ha) reference SFR
      sfr_ha_err:          -999.0, $
      sfr_hb_uncor:        -999.0, $ ; SFR(Hb) uncorrected
      sfr_hb_uncor_err:    -999.0, $
      sfr_hb_hbhg:         -999.0, $ ; SFR(Hb) extinction-corrected using Hb/Hg
      sfr_hb_hbhg_err:     -999.0, $
      sfr_hb_best:         -999.0, $ ; SFR(Hb) best (this paper)
      sfr_hb_best_err:     -999.0, $
      sfr_oii_uncor:       -999.0, $ ; SFR([O II]) uncorrected
      sfr_oii_uncor_err:   -999.0, $
      sfr_oii_hbhg:        -999.0, $ ; SFR([O II]) extinction-corrected using Hb/Hg
      sfr_oii_hbhg_err:    -999.0, $
      sfr_oii_hahb:        -999.0, $ ; SFR([O II]) extinction-corrected using Ha/Hb
      sfr_oii_hahb_err:    -999.0, $
      sfr_oii_k98:         -999.0, $ ; SFR([O II]) using Kennicutt 1998
      sfr_oii_k98_err:     -999.0, $
      sfr_oii_r02:         -999.0, $ ; SFR([O II]) using Rosa-Gonzalez et al. 2002
      sfr_oii_r02_err:     -999.0, $
      sfr_oii_k04:         -999.0, $ ; SFR([O II]) using Kewley et al. 2004
      sfr_oii_k04_err:     -999.0, $
      sfr_oii_k04_theory:  -999.0, $
      sfr_oii_k04_theory_err: -999.0, $
      sfr_oii_best:        -999.0, $ ; SFR([O II]) best (this paper)
      sfr_oii_best_err:    -999.0, $
      sfr_uband_uncor:     -999.0, $
      sfr_uband_uncor_err: -999.0, $
      sfr_uband_cram:      -999.0, $ ; Cram et al. (1998)
      sfr_uband_cram_err:  -999.0, $
      sfr_uband_best:      -999.0, $ ; SFR(U) best (this paper)
      sfr_uband_best_err:  -999.0}
    nfgs = atlas

    atlas = replicate(atlas,n_elements(atlasdust))
    nfgs = replicate(nfgs,n_elements(nfgsdust))

    atlas.galaxy = atlasnodust.galaxy
    atlas.sfr_ha = atlasnodust.sfr_h_alpha
    atlas.sfr_ha_err = atlasnodust.sfr_h_alpha_err

    nfgs.galaxy = nfgsnodust.galaxy
    nfgs.sfr_ha = nfgsnodust.sfr_h_alpha
    nfgs.sfr_ha_err = nfgsnodust.sfr_h_alpha_err

; --------------------------------------------------
; SFR(U) uncorrected
; --------------------------------------------------

; Atlas    

    indx = where((atlasdust.u_lum_obs gt -900.0),nindx)

    atlas[indx].sfr_uband_uncor = atlasdust[indx].u_lum_obs + alog10(lusun) + alog10(Usfr_obs_coeff)
    atlas[indx].sfr_uband_uncor_err = sqrt(atlasdust[indx].u_lum_obs_err^2 + $
      (Usfr_obs_coeff_err/Usfr_obs_coeff/alog(10.0))^2)

; NFGS    
    
    indxnfgs = where((nfgsdust.u_lum_obs gt -900.0),nindxnfgs)

    nfgs[indxnfgs].sfr_uband_uncor = nfgsdust[indxnfgs].u_lum_obs + alog10(lusun) + alog10(Usfr_obs_coeff)
    nfgs[indxnfgs].sfr_uband_uncor_err = sqrt(nfgsdust[indxnfgs].u_lum_obs_err^2 + $
      (Usfr_obs_coeff_err/Usfr_obs_coeff/alog(10.0))^2)

; --------------------------------------------------
; SFR(U) based on Cram et al. (1998)
; --------------------------------------------------

; Atlas    

    indx = where((atlasdust.u_lum_obs gt -900.0),nindx)

    atlas[indx].sfr_uband_cram = alog10(9.51D-45) + atlasdust[indx].u_lum_obs + alog10(lusun)
    atlas[indx].sfr_uband_cram_err = 0.0

; NFGS    
    
    indxnfgs = where((nfgsdust.u_lum_obs gt -900.0),nindxnfgs)

    nfgs[indxnfgs].sfr_uband_cram = alog10(9.51D-45) + nfgsdust[indxnfgs].u_lum_obs + alog10(lusun)
    nfgs[indxnfgs].sfr_uband_cram_err = 0.0

; --------------------------------------------------
; SFR(U) best (this paper's calibration)
; --------------------------------------------------    
    
; Atlas    

    indx = where((atlasdust.u_lum_obs gt -900.0) and (atlasdust.b_lum_obs gt -900),nindx)

    lubandobs = atlasdust[indx].u_lum_obs[0] + alog10(lusun)
    loglb = atlasdust[indx].b_lum_obs

    atlas[indx].sfr_uband_best = alog10(uband_sfr(loglb,lubandobs))    
    atlas[indx].sfr_uband_best_err = 0.0
    
; NFGS    
    
    indxnfgs = where((nfgsdust.u_lum_obs gt -900) and (nfgsdust.b_lum_obs gt -900),nindxnfgs)

    lubandobs_nfgs = nfgsdust[indxnfgs].u_lum_obs[0] + alog10(lusun)
    loglb_nfgs = nfgsdust[indxnfgs].b_lum_obs

    nfgs[indxnfgs].sfr_uband_best = alog10(uband_sfr(loglb_nfgs,lubandobs_nfgs))    
    nfgs[indxnfgs].sfr_uband_best_err = 0.0
    
; --------------------------------------------------
; SFR(Hb) uncorrected
; --------------------------------------------------

; Atlas    

    indx = where((atlasdust.h_beta_ew_uncor[0] ge hbewcut),nindx)

    atlas[indx].sfr_hb_uncor = atlasdust[indx].h_beta_lum[0] + alog10(lsun) + loghaconst + alog10(hahb)
    atlas[indx].sfr_hb_uncor_err = atlasdust[indx].h_beta_lum[1]

; NFGS    
    
    indxnfgs = where((nfgsdust.h_beta_ew_uncor[0] ge hbewcut),nindxnfgs)

    nfgs[indxnfgs].sfr_hb_uncor = nfgsdust[indxnfgs].h_beta_lum[0] + alog10(lsun) + loghaconst + alog10(hahb)
    nfgs[indxnfgs].sfr_hb_uncor_err = nfgsdust[indxnfgs].h_beta_lum[1]

; --------------------------------------------------
; SFR(Hb) extinction-corrected using Hb/Hg      
; --------------------------------------------------

; Atlas    

    indx = where((atlasnodust.ebv_hbhg_err gt 0.0) and (atlasdust.h_beta_ew[0] ge hghbewcut) and $
      (atlasdust.h_gamma[0]/atlasdust.h_gamma[1] gt hgsnrcut),nindx)

    AHb = atlasnodust[indx].ebv_hbhg*klhb
    AHb_err = atlasnodust[indx].ebv_hbhg_err*klhb

    atlas[indx].sfr_hb_hbhg = atlasdust[indx].h_beta_lum[0] + 0.4*AHb + alog10(lsun) + loghaconst + alog10(hahb)
    atlas[indx].sfr_hb_hbhg_err = sqrt(atlasdust[indx].h_beta_lum[1]^2 + (0.4*AHb_err)^2)

; NFGS    
    
    indxnfgs = where((nfgsnodust.ebv_hbhg_err gt 0.0) and (nfgsdust.h_beta_ew[0] ge hghbewcut) and $
      (nfgsdust.h_gamma[0]/nfgsdust.h_gamma[1] gt hgsnrcut),nindxnfgs)

    if (nindxnfgs ne 0L) then begin

       AHbnfgs = nfgsnodust[indxnfgs].ebv_hbhg*klhb
       AHbnfgs_err = nfgsnodust[indxnfgs].ebv_hbhg_err*klhb

       nfgs[indxnfgs].sfr_hb_hbhg = nfgsdust[indxnfgs].h_beta_lum[0] + 0.4*AHbnfgs + $
         alog10(lsun) + loghaconst + alog10(hahb)
       nfgs[indxnfgs].sfr_hb_hbhg_err = sqrt(nfgsdust[indxnfgs].h_beta_lum[1]^2 + (0.4*AHbnfgs_err)^2)

    endif

; --------------------------------------------------
; SFR(Hb) best (this paper's calibration)
; --------------------------------------------------    
    
; Atlas    

    indx = where((atlasdust.h_beta_ew_uncor[0] ge hbewcut) and (atlasdust.b_lum_obs gt -900),nindx)

    lhbobs = atlasdust[indx].h_beta_lum[0] + alog10(lsun)
    loglb = atlasdust[indx].b_lum_obs

    atlas[indx].sfr_hb_best = alog10(hb_sfr(loglb,lhbobs))    
    atlas[indx].sfr_hb_best_err = 0.0
    
;   lhbobs = lsun*10^atlasdust[indx].h_beta_lum[0]
;   Blum = 10^atlasdust[indx].b_lum_obs
;   
;   lolum = where(alog10(Blum) lt hbsfr_lumcut,nlolum,comp=hilum,ncomp=nhilum)
;   hbsfr = lhbobs*0.0
;   
;   if (nlolum ne 0L) then hbsfr[lolum] = hbsfr_coeff_lolum*lhbobs[lolum]
;   if (nhilum ne 0L) then hbsfr[hilum] = hbsfr_coeff_hilum[0] * $
;     lhbobs[hilum]*Blum[hilum]^hbsfr_coeff_hilum[1]
;
;   atlas[indx].sfr_hb_best     = alog10(hbsfr)
;   atlas[indx].sfr_hb_best_err = 0.0

; NFGS    
    
    indxnfgs = where((nfgsdust.h_beta_ew_uncor[0] ge hbewcut) and (nfgsdust.b_lum_obs gt -900),nindxnfgs)

    lhbobs_nfgs = nfgsdust[indxnfgs].h_beta_lum[0] + alog10(lsun)
    loglb_nfgs = nfgsdust[indxnfgs].b_lum_obs

    nfgs[indxnfgs].sfr_hb_best = alog10(hb_sfr(loglb_nfgs,lhbobs_nfgs))    
    nfgs[indxnfgs].sfr_hb_best_err = 0.0
    
;   lhbobs_nfgs = lsun*10^nfgsdust[indxnfgs].h_beta_lum[0]
;   Blum_nfgs = 10^nfgsdust[indxnfgs].b_lum_obs
;
;   lolum_nfgs = where(alog10(Blum_nfgs) lt hbsfr_lumcut,nlolum_nfgs,comp=hilum_nfgs,ncomp=nhilum_nfgs)
;   hbsfr_nfgs = lhbobs_nfgs*0.0
;   
;   if (nlolum_nfgs ne 0L) then hbsfr_nfgs[lolum_nfgs] = hbsfr_coeff_lolum*lhbobs_nfgs[lolum_nfgs]
;   if (nhilum_nfgs ne 0L) then hbsfr_nfgs[hilum_nfgs] = hbsfr_coeff_hilum[0] * $
;     lhbobs_nfgs[hilum_nfgs]*Blum_nfgs[hilum_nfgs]^hbsfr_coeff_hilum[1]
;
;   nfgs[indxnfgs].sfr_hb_best     = alog10(hbsfr_nfgs)
;   nfgs[indxnfgs].sfr_hb_best_err = 0.0
    
; --------------------------------------------------
; SFR([O II]) uncorrected
; --------------------------------------------------

; Atlas    

    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut),nindx)

    atlas[indx].sfr_oii_uncor = atlasdust[indx].oii_3727_lum[0] + alog10(lsun) + loghaconst + alog10(haoii)
    atlas[indx].sfr_oii_uncor_err = atlasdust[indx].oii_3727_lum[1]

; NFGS    
    
    indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut),nindxnfgs)

    nfgs[indxnfgs].sfr_oii_uncor = nfgsdust[indxnfgs].oii_3727_lum[0] + alog10(lsun) + loghaconst + alog10(haoii)
    nfgs[indxnfgs].sfr_oii_uncor_err = nfgsdust[indxnfgs].oii_3727_lum[1]

; --------------------------------------------------
; SFR([O II]) extinction-corrected using Hb/Hg      
; --------------------------------------------------

; Atlas    

    indx = where((atlasnodust.ebv_hbhg_err gt 0.0) and (atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasdust.h_gamma[0]/atlasdust.h_gamma[1] gt hgsnrcut) and (atlasdust.h_beta_ew[0] gt hghbewcut),nindx)

    Aoii = atlasnodust[indx].ebv_hbhg*kloii
    Aoii_err = atlasnodust[indx].ebv_hbhg_err*kloii

    atlas[indx].sfr_oii_hbhg = atlasdust[indx].oii_3727_lum[0] + 0.4*Aoii + alog10(lsun) + loghaconst + alog10(haoii)
    atlas[indx].sfr_oii_hbhg_err = sqrt(atlasdust[indx].oii_3727_lum[1]^2 + (0.4*Aoii_err)^2)

; NFGS    
    
    indxnfgs = where((nfgsnodust.ebv_hbhg_err gt 0.0) and (nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
      (nfgsdust.h_gamma[0]/nfgsdust.h_gamma[1] gt hgsnrcut) and (nfgsdust.h_beta_ew[0] gt hghbewcut),nindxnfgs)

    if (nindxnfgs ne 0L) then begin
       
       Aoiinfgs = nfgsnodust[indxnfgs].ebv_hbhg*kloii
       Aoiinfgs_err = nfgsnodust[indxnfgs].ebv_hbhg_err*kloii

       nfgs[indxnfgs].sfr_oii_hbhg = nfgsdust[indxnfgs].oii_3727_lum[0] + 0.4*Aoiinfgs + $
         alog10(lsun) + loghaconst + alog10(haoii)
       nfgs[indxnfgs].sfr_oii_hbhg_err = sqrt(nfgsdust[indxnfgs].oii_3727_lum[1]^2 + (0.4*Aoiinfgs_err)^2)

    endif
; --------------------------------------------------
; SFR([O II]) extinction-corrected using Ha/Hb
; --------------------------------------------------

; Atlas    

    indx = where((atlasnodust.ebv_hahb_err gt 0.0) and (atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut),nindx)

    atlas[indx].sfr_oii_hahb = atlasnodust[indx].oii_3727_lum[0] + alog10(lsun) + loghaconst + alog10(haoii)
    atlas[indx].sfr_oii_hahb_err = atlasnodust[indx].oii_3727_lum[1]

; NFGS    
    
    indxnfgs = where((nfgsnodust.ebv_hahb_err gt 0.0) and (nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut),nindx)

    nfgs[indxnfgs].sfr_oii_hahb = nfgsnodust[indxnfgs].oii_3727_lum[0] + alog10(lsun) + loghaconst + alog10(haoii)
    nfgs[indxnfgs].sfr_oii_hahb_err = nfgsdust[indxnfgs].oii_3727_lum[1]

; --------------------------------------------------
; SFR([O II]) from Rosa-Gonzalez et al. 2002
; --------------------------------------------------

    oiiconst = 8.4D-41
    oiiconst_err = 0.2*oiiconst ; NOTE!!!!

; Atlas    
    
    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut),nindx)

    atlas[indx].sfr_oii_r02 = atlasdust[indx].oii_3727_lum[0] + alog10(lsun) + alog10(oiiconst)
    atlas[indx].sfr_oii_r02_err = sqrt(atlasdust[indx].oii_3727_lum[1]^2 + (oiiconst_err/oiiconst/alog(10.0))^2)

; NFGS    
    
    indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut),nindxnfgs)

    nfgs[indxnfgs].sfr_oii_r02 = nfgsdust[indxnfgs].oii_3727_lum[0] + alog10(lsun) + alog10(oiiconst)
    nfgs[indxnfgs].sfr_oii_r02_err = sqrt(nfgsdust[indxnfgs].oii_3727_lum[1]^2 + (oiiconst_err/oiiconst/alog(10.0))^2)

; --------------------------------------------------
; SFR([O II]) from Kennicutt 1998
; --------------------------------------------------

    oiiconst = 1.4D-41
    oiiconst_err = 0.4D-41

; Atlas    
    
    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut),nindx)

    atlas[indx].sfr_oii_k98 = atlasdust[indx].oii_3727_lum[0] + alog10(lsun) + alog10(oiiconst)
    atlas[indx].sfr_oii_k98_err = sqrt(atlasdust[indx].oii_3727_lum[1]^2 + (oiiconst_err/oiiconst/alog(10.0))^2)

; NFGS    
    
    indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut),nindxnfgs)

    nfgs[indxnfgs].sfr_oii_k98 = nfgsdust[indxnfgs].oii_3727_lum[0] + alog10(lsun) + alog10(oiiconst)
    nfgs[indxnfgs].sfr_oii_k98_err = sqrt(nfgsdust[indxnfgs].oii_3727_lum[1]^2 + (oiiconst_err/oiiconst/alog(10.0))^2)

; --------------------------------------------------
; SFR([O II]) from Kewley et al. 2004
; --------------------------------------------------

; Atlas
    
    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasdust.oiii_5007[0]/atlasdust.oiii_5007[1] gt snrcut) and $
      (atlasdust.h_beta[0]/atlasdust.h_beta[1] gt snrcut),nindx)

    oii3727  = atlasdust[indx].oii_3727
    oiii4959 = atlasdust[indx].oiii_4959
    oiii5007 = atlasdust[indx].oiii_5007
    hbeta  = atlasdust[indx].h_beta

    lumoii_o = lsun*10^atlasdust[indx].oii_3727_lum[0]

; [1] compute the intrinsic luminosity from eq. (18)

    lumoii_i = 3.11D-20*lumoii_o^1.495

; [2] estimate E(B-V) from eq. (16)

    ebv = (0.174 * alog10(lumoii_i) - 6.84) > 0.0

; [3] correct R23 for the E(B-V) estimated in step [2]; we have to
; trick IUNRED_LINEDUST() into using the right E(B-V) value

    line = {$
      linename:       ['OII_3727','H_BETA','OIII_4959','OIII_5007','H_ALPHA'], $
      oii_3727:       [0.0,0.0], $
      oii_3727_wave:  3727.0, $
      h_beta:         [0.0,0.0], $
      h_beta_wave:    4861.0, $
      oiii_4959:      [0.0,0.0], $
      oiii_4959_wave: 4959.0, $
      oiii_5007:      [0.0,0.0], $
      oiii_5007_wave: 5007.0, $
      h_alpha:        [0.0,0.0], $
      h_alpha_wave:   6563.0}
    line = replicate(line,nindx)

    factor = hahb_k04*10^(0.4*(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))*ebv)

    halpha = hbeta
    halpha[0,*] = hbeta[0,*]*factor

    line.oii_3727  = oii3727
    line.h_beta    = hbeta
    line.oiii_4959 = oiii4959
    line.oiii_5007 = oiii5007
    line.h_alpha   = halpha

    linenodust = iunred_linedust(line,snrcut=snrcut,/silent,/nopropagate)

; [4] use eq. (11) to estimate the abundance from the Z94 diagnostic;
; note that Kewley's definition of R23 does not include [O III] 4959! 

    r23 = alog10((linenodust.oii_3727[0]+linenodust.oiii_4959[0]+linenodust.oiii_5007[0])/linenodust.h_beta[0])
;   r23 = alog10((linenodust.oii_3727[0]+linenodust.oiii_5007[0])/linenodust.h_beta[0])
    logoh = 9.265 - 0.33*r23 - 0.202*r23^2 - 0.207*r23^3 - 0.333*r23^4

    abund = im_abundance(linenodust,snrcut=snrcut)
    logoh_me = abund.zstrong_12oh_zkh94

; [5] use eq. (10) to compute the reddening and abundance corrected
; SFR([O II],Z) estimate

    sfr_oii_theory = logoh*0.0-999.0

    numer = 7.9D-42*lumoii_i
    denom = -1857.24+612.693*logoh-67.0264*logoh^2+2.43209*logoh^3
    pos = where(denom gt 0.0,npos)
    if (npos ne 0L) then begin
       atlas[indx[pos]].sfr_oii_k04_theory = alog10(numer[pos]/denom[pos])
       atlas[indx[pos]].sfr_oii_k04_theory_err = 0.1 ; NOTE!
    endif

    atlas[indx].sfr_oii_k04 = alog10(7.9D-42*lumoii_i/(-1.75*logoh+16.73)) ; 
    atlas[indx].sfr_oii_k04_err = 0.1 ; NOTE!!!

; NFGS

    indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
      (nfgsdust.oiii_5007[0]/nfgsdust.oiii_5007[1] gt snrcut) and $
      (nfgsdust.h_beta[0]/nfgsdust.h_beta[1] gt snrcut),nindx)

    oii3727  = nfgsdust[indxnfgs].oii_3727
    oiii4959 = nfgsdust[indxnfgs].oiii_4959
    oiii5007 = nfgsdust[indxnfgs].oiii_5007
    hbeta  = nfgsdust[indxnfgs].h_beta

    lumoii_o = lsun*10^nfgsdust[indxnfgs].oii_3727_lum[0]

; [1] compute the intrinsic luminosity from eq. (18)

    lumoii_i = 3.11D-20*lumoii_o^1.495
    lumoii_true = lsun*10^nfgsnodust[indxnfgs].oii_3727_lum[0]

; [2] estimate E(B-V) from eq. (16)

    ebv = (0.174 * alog10(lumoii_i) - 6.84) > 0.02

; [3] correct R23 for the E(B-V) estimated in step [2]; we have to
; trick IUNRED_LINEDUST() into using the right E(B-V) value

    line = {$
      linename:       ['OII_3727','H_BETA','OIII_4959','OIII_5007','H_ALPHA'], $
      oii_3727:       [0.0D,0.0D], $
      oii_3727_wave:  3727.0, $
      h_beta:         [0.0D,0.0D], $
      h_beta_wave:    4861.0, $
      oiii_4959:      [0.0D,0.0D], $
      oiii_4959_wave: 4959.0, $
      oiii_5007:      [0.0D,0.0D], $
      oiii_5007_wave: 5007.0, $
      h_alpha:        [0.0D,0.0D], $
      h_alpha_wave:   6563.0}
    line = replicate(line,nindx)

    factor = hahb_k04*10^(0.4*(k_lambda(4861,/odonnell)-k_lambda(6563,/odonnell))*ebv)

    halpha = hbeta
    halpha[0,*] = hbeta[0,*]*factor

    line.oii_3727  = oii3727
    line.h_beta    = hbeta
    line.oiii_4959 = oiii4959
    line.oiii_5007 = oiii5007
    line.h_alpha   = halpha

    linenodust = iunred_linedust(line,snrcut=snrcut,/silent,/nopropagate)

; [4] use eq. (11) to estimate the abundance from the Z94 diagnostic;
; note that Kewley's definition of R23 does not include [O III] 4959! 

    r23 = alog10((linenodust.oii_3727[0]+linenodust.oiii_4959[0]+linenodust.oiii_5007[0])/linenodust.h_beta[0])
;   r23 = alog10((linenodust.oii_3727[0]+linenodust.oiii_5007[0])/linenodust.h_beta[0])
    logoh = 9.265 - 0.33*r23 - 0.202*r23^2 - 0.207*r23^3 - 0.333*r23^4

    abund = im_abundance(linenodust,snrcut=snrcut)
    logoh_me = abund.zstrong_12oh_zkh94

; [5] use eq. (10) to compute the reddening and abundance corrected
; SFR([O II],Z) estimate

    sfr_oii_theory = logoh*0.0-999.0

    numer = 7.9D-42*lumoii_i
    denom = -1857.24+612.693*logoh-67.0264*logoh^2+2.43209*logoh^3
    pos = where(denom gt 0.0,npos)
    if (npos ne 0L) then begin
       nfgs[indxnfgs[pos]].sfr_oii_k04_theory = alog10(numer[pos]/denom[pos])
       nfgs[indxnfgs[pos]].sfr_oii_k04_theory_err = 0.1 ; NOTE!
    endif

    nfgs[indxnfgs].sfr_oii_k04 = alog10(7.9D-42*lumoii_i/(-1.75*logoh+16.73))
    nfgs[indxnfgs].sfr_oii_k04_err = 0.1 ; NOTE!!!

;;; ---------------------------------------------------------------------------    
;;; try to reproduce the Kewley et al. (2004) results
;;; ---------------------------------------------------------------------------    
;;
;;    k04 = rsex('04kewley_table1.txt')
;;    match, long(strmid(nfgsdust[indxnfgs].drift_file,0,3)), k04.id, indx1, indx2
;;
;;    my_ebv_k04 = get_ebv(nfgsnodust[indxnfgs[indx1]].hahb,temp=10771.0,/hahb,/ccm)
;;    
;;    good = where(nfgsnodust[indxnfgs[indx1]].ebv_hahb_err gt 0.0)
;;    plot, nfgsnodust[indxnfgs[indx1[good]]].ebv_hahb, k04[indx2[good]].e_bv_real-$
;;      nfgsnodust[indxnfgs[indx1[good]]].ebv_hahb, ps=4, xsty=3, ysty=3
;;    jj = im_stats(k04[indx2[good]].e_bv_real-nfgsnodust[indxnfgs[indx1[good]]].ebv_hahb,/verbose)
;;
;;    plot, my_ebv_k04[good], k04[indx2[good]].e_bv_real-my_ebv_k04, ps=4, xsty=3, ysty=3
;;    jj = im_stats(k04[indx2[good]].e_bv_real-my_ebv_k04[good],/verbose)
;;    niceprint, my_ebv_k04[good], k04[indx2[good]].e_bv_real
;;
;;; --
;;    
;;    plot, k04[indx2].e_bv_est, k04[indx2].e_bv_est-ebv[indx1], ps=4
;;    jj = im_stats(k04[indx2].e_bv_est-ebv[indx1],/verbose)
;;
;;    good = where(nfgsnodust[indxnfgs[indx1]].oii_3727_lum[1] gt 0.0)
;;    plot, nfgsnodust[indxnfgs[indx1[good]]].oii_3727_lum[0], k04[indx2[good]].l_o2 - $
;;      nfgsnodust[indxnfgs[indx1[good]]].oii_3727_lum[0], ps=4, xsty=3, ysty=3
;;    jj = im_stats(k04[indx2[good]].l_o2-nfgsnodust[indxnfgs[indx1[good]]].oii_3727_lum[0],/verbose)
;;
;;stop
;;    
;;; first check that her SFR(Ha) indeed comes from the
;;; reddening-corrected L(Ha) --> kinda
;;
;;    plot, 7.9D-42*lsun*10^k04.l_ha, k04.sfr_ha_k98, ps=4, /xlog, /ylog, $
;;      xsty=3, ysty=3, xrange=[0.01,100], yrange=[0.01,100]
;;    plot, k04.sfr_ha_k98, alog10(7.9D-42*lsun*10^k04.l_ha/k04.sfr_ha_k98), ps=4, /xlog, $
;;      xsty=3, ysty=3, xrange=[0.01,100]
;;
;;; try to reproduce Figure 16, left; i can get lisa's claimed scatter
;;; of 0.21 dex only if i reject all 2-sigma outliers
;;    
;;    good = where((k04.sfr_ha_k98 gt 0.0) and (k04.sfr_o2_eq18 gt 0.0))
;;    plot, k04[good].sfr_ha_k98, k04[good].sfr_o2_eq18, ps=4, /xlog, /ylog, $
;;      xsty=3, ysty=3, xrange=[0.01,100], yrange=[0.01,100]
;;    djs_oplot, 10^!x.crange, 10^!y.crange
;;    plot, k04[good].sfr_ha_k98, alog10(k04[good].sfr_ha_k98/k04[good].sfr_o2_eq18), ps=4, /xlog, $
;;      xsty=3, ysty=3, xrange=[0.01,100], yrange=[-1.0,1.0]
;;    djs_oplot, 10^!x.crange, [0,0]
;;    jj = im_stats(alog10(k04[good].sfr_ha_k98/k04[good].sfr_o2_eq18),/verbose)
;;
;;stop    
;;    
;;; try to reproduce Figure 16, right
;;    
;;    good = where((k04.sfr_ha_k98 gt 0.0) and (k04.sfr_o2_eq15 gt 0.0))
;;    plot, k04[good].sfr_ha_k98, k04[good].sfr_o2_eq15, ps=4, /xlog, /ylog, $
;;      xsty=3, ysty=3, xrange=[0.01,100], yrange=[0.01,100]
;;    djs_oplot, 10^!x.crange, 10^!y.crange
;;    plot, k04[good].sfr_ha_k98, alog10(k04[good].sfr_ha_k98/k04[good].sfr_o2_eq15), ps=4, /xlog, $
;;      xsty=3, ysty=3, xrange=[0.01,100], yrange=[-1.0,1.0]
;;    djs_oplot, 10^!x.crange, [0,0]
;;    jj = im_stats(alog10(k04[good].sfr_ha_k98/k04[good].sfr_o2_eq15),/verbose)
;;    
;;; try to reproduce Figure 15
;;
;;    plot, k04.l_o2+alog10(lsun), k04.e_bv_real, ps=4, $
;;      xsty=3, ysty=3, xrange=[38,44], yrange=[0,1.0]
;;    xaxis = findgen((!x.crange[1]-!x.crange[0])/0.01+1)*0.01+!x.crange[0]
;;    yaxis = 0.174*xaxis - 6.84
;;    chop = where(yaxis gt 0.0)
;;    djs_oplot, xaxis[chop], yaxis[chop], line=2
;;; ---------------------------------------------------------------------------    
    
; --------------------------------------------------
; SFR([O II]) best (this paper's calibration)
; --------------------------------------------------    
    
; Atlas    

    indx = where((atlasdust.oii_3727[0]/atlasdust.oii_3727[1] gt snrcut) and $
      (atlasdust.b_lum_obs gt -900),nindx)

    loiiobs = atlasdust[indx].oii_3727_lum[0] + alog10(lsun)
    loglb = atlasdust[indx].b_lum_obs

    atlas[indx].sfr_oii_best = alog10(oii_sfr(loglb,loiiobs))    
    atlas[indx].sfr_oii_best_err = 0.0
    
;   loiiobs = lsun*10^atlasdust[indx].oii_3727_lum[0]
;   Blum = 10^atlasdust[indx].b_lum_obs
;   
;   lolum = where(alog10(Blum) lt oiisfr_lumcut,nlolum,comp=hilum,ncomp=nhilum)
;   oiisfr = loiiobs*0.0
;   
;   if (nhilum ne 0L) then begin
;      oiisfr[hilum] = oiisfr_coeff[0]*loiiobs[hilum]*Blum[hilum]^oiisfr_coeff[1]
;      atlas[indx[hilum]].sfr_oii_best     = alog10(oiisfr[hilum])
;      atlas[indx[hilum]].sfr_oii_best_err = 0.0
;   endif
    
; NFGS    
    
    indxnfgs = where((nfgsdust.oii_3727[0]/nfgsdust.oii_3727[1] gt snrcut) and $
      (nfgsdust.b_lum_obs gt -900),nindxnfgs)

    loiiobs_nfgs = nfgsdust[indxnfgs].oii_3727_lum[0] + alog10(lsun)
    loglb_nfgs = nfgsdust[indxnfgs].b_lum_obs

    nfgs[indxnfgs].sfr_oii_best = alog10(oii_sfr(loglb_nfgs,loiiobs_nfgs))    
    nfgs[indxnfgs].sfr_oii_best_err = 0.0
    
;   loiiobs_nfgs = lsun*10^nfgsdust[indxnfgs].oii_3727_lum[0]
;   Blum_nfgs = 10^nfgsdust[indxnfgs].b_lum_obs
;   
;   lolum_nfgs = where(alog10(Blum_nfgs) lt oiisfr_lumcut,nlolum_nfgs,comp=hilum_nfgs,ncomp=nhilum_nfgs)
;   oiisfr_nfgs = loiiobs_nfgs*0.0
;   
;   if (nhilum_nfgs ne 0L) then begin
;      oiisfr_nfgs[hilum_nfgs] = oiisfr_coeff[0]*loiiobs_nfgs[hilum_nfgs]*Blum_nfgs[hilum_nfgs]^oiisfr_coeff[1]
;      nfgs[indxnfgs[hilum_nfgs]].sfr_oii_best     = alog10(oiisfr_nfgs[hilum_nfgs])
;      nfgs[indxnfgs[hilum_nfgs]].sfr_oii_best_err = 0.0
;   endif

; --------------------------------------------------    
; Write out
; --------------------------------------------------    

    if keyword_set(write) then begin

       atlasfile = 'atlas_sfrs'+suffix+'.fits'
       splog, 'Writing '+path+atlasfile+'.'
       mwrfits, atlas, path+atlasfile, /create

       nfgsfile = 'nfgs_sfrs'+suffix+'.fits'
       splog, 'Writing '+path+nfgsfile+'.'
       mwrfits, nfgs, path+nfgsfile, /create

    endif

return
end
