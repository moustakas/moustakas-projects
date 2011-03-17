pro write_ages_agn_catalog
; jm08jul10nyu - generate a catalog of galaxy classifications based on
;   a variety of spectrophotometric and multiwavelength criteria; the
;   output file (both ASCII and FITS) contain 39943 and have been
;   line-matched to all the other AGES catalogs; this routine should
;   be run after every new version of ISPEC
; jm08aug29nyu - use v2.0 ispec catalog
; jm08nov01nyu - use v2.1 ispec catalog and v2.1 unfluxed catalog and
;   new ICLASSIFICATION code; remove chi2cut requirement; increase
;   broad_sigmacut from 250 to 300 km/s

    ewhbcor = 2.5 ; from mzsampleplots
    ewhacor = 0.8

    chi2cut1 = 1000.0 ; 5.0
    sigmacut1 = 250.0

    bpt_snrcut1 = 3.0
    broad_snrcut1 = 10.0
    broad_sigmacut1 = 300.0 ; 250.0
    max_broad_sigma1 = 5000.0 ; maximum value (usually means not a good fit)

    nev_ewcut1 = 2.0
    nev_snrcut1 = 3.0
    
;   vv = 'v1.0' ; uses v1.1 emission-line catalog
;   vv = 'v1.1' ; uses v2.0 emission-line catalog
    vv = 'v1.2' ; see nov01 notes, above

    version = ages_version(/ispec)
    apath = ages_path(/analysis)
    outpath = ages_path(/catalogs)
    
    aa = read_ages(/ancillary)
    ii = read_ages(/ispec);,/notweak)
    uu = read_ages(/unfluxed)
    cc = mrdfits(apath+'catalog.cat.noguidestars.fits.gz',1)
    nobj = n_elements(cc)
    
; classification codes:
;    -1 = no classification attempted
;     0 = classification not possible
;     2 = AGN
;     4 = SF
;     8 = SF/AGN

    class = {$
      ra:  -999.0D,$
      dec: -999.0D,$
      z:   -999.0, $
      pass:     -1,$
      aper:     -1,$
      bpt:      -1,$
;     nev:      -1,$
      broad:    -1}
;     xray:     -1}
    class = replicate(class,nobj)

    format = [$
      ['RA','F12.7'], $ ; for the output file
      ['DEC','F12.7'],$
      ['Z','F12.5'],  $           
      ['PASS','I5'],  $
      ['APER','I5'],  $
      ['BPT','I5'],   $
;     ['NEV','I5'],   $
      ['BROAD','I5']]
;     ['XRAY','I5']]

;   class = replicate({ra: -999.0D, dec: -999.0D, z: -999.0, $
;     pass: -1, aper: -1, bpt: -1, nev: -1, broad: -1, xray: -1},nobj)
;   format = [['RA','F12.7'],['DEC','F12.7'],['Z','F12.5'],$ ; for the output file
;     ['PASS','I5'],['APER','I5'],['BPT','I5'],['NEV','I5'],$
;     ['BROAD','I5'],['XRAY','I5']]

    spherematch, aa.ra, aa.dec, 15.0D*cc.cat_ra, $
      cc.cat_dec, 1.5/3600.0, m1, m2
;   niceprint, aa[m1[0:10]].ra, cc[m2[0:10]].cat_ra*15.0D

    class.ra = cc.cat_ra
    class.dec = cc.cat_dec
    class[m2].pass = aa[m1].pass
    class[m2].aper = aa[m1].aper
    class[m2].z = aa[m1].z

; #########################
; BPT    
; #########################

; fluxed plates    
    these = where((aa[m1].z gt 0.01) and (aa[m1].z lt 1.0) and $
;     (strtrim(aa[m1].class,2) eq 'GALAXY') and $
      (aa[m1].pass ne 106) and (aa[m1].pass ne 110) and $
      (aa[m1].pass ne 209) and (aa[m1].pass ne 310) and $
      (aa[m1].pass ne 311))
    bpt = iclassification(ii[m1[these]],snrcut_class=bpt_snrcut1,$
      niihacut=-0.3,chi2cut=chi2cut1,ratios=iratios)

; this code uses upper limits    
;   unk = where(strtrim(bpt.bpt_pure_nii_class,2) eq 'Unknown')
;   sf = where(strtrim(bpt.bpt_pure_nii_class,2) eq 'HII')
;   agn = where(strtrim(bpt.bpt_pure_nii_class,2) eq 'AGN')
;   help, sf, agn, these
; this code defines mixed classes and does not use upper limits
    unk =   where(strtrim(iratios.final_class,2) eq 'UNKNOWN')
    sf =    where(strtrim(iratios.final_class,2) eq 'SF')
    agn =   where(strtrim(iratios.final_class,2) eq 'AGN')
    sfagn = where(strtrim(iratios.final_class,2) eq 'SF/AGN')
    help, sf, agn, sfagn, these

    class[m2[these[unk]]].bpt = 0
    class[m2[these[agn]]].bpt = 2
    class[m2[these[sf]]].bpt = 4
    class[m2[these[sfagn]]].bpt = 8

need to fix this -- we need to distinguish the unknowns between no classification attempted and inconclusive classification    
    
    
stop    
    
; unfluxed plates; statistically correct H-alpha and H-beta for
; absorption 
    uthese = where((aa[m1].z gt 0.01) and (aa[m1].z lt 1.0) and $
;     (strtrim(aa[m1].class,2) eq 'GALAXY') and $
      ((aa[m1].pass eq 106) or (aa[m1].pass eq 110) or $
      (aa[m1].pass eq 209) or (aa[m1].pass eq 310) or $
      (aa[m1].pass eq 311)),nuthese)

    uucor = uu
;;    hb = where((uu.h_beta[0] gt 0.0) and (uu.h_beta[1] gt 0.0) and $
;;      (uu.h_beta_ew[0] gt 0.0) and (uu.h_beta_ew[1] gt 0.0),nhb)
;;    if (nhb ne 0L) then begin
;;       hbfactor = uu[hb].h_beta_ew[0]/ewhbcor
;;       uucor[hb].h_beta_ew[0] = uu[hb].h_beta_ew[0] + ewhbcor
;;       uucor[hb].h_beta[0] = uu[hb].h_beta[0]*hbfactor
;;       uucor[hb].h_beta[1] = uu[hb].h_beta[1]*hbfactor
;;;      im_plothist, uu[hb].h_beta_ew[0]/ewhbcor, xr=[0,10], bin=0.1
;;    endif
;;    ha = where((uu.h_alpha[0] gt 0.0) and (uu.h_alpha[1] gt 0.0) and $
;;      (uu.h_alpha_ew[0] gt 0.0) and (uu.h_alpha_ew[1] gt 0.0),nha)
;;    if (nha ne 0L) then begin
;;       hafactor = uu[ha].h_alpha_ew[0]/ewhacor
;;       uucor[ha].h_alpha_ew[0] = uucor[ha].h_alpha_ew[0] + ewhacor
;;       uucor[ha].h_alpha[0] = uucor[ha].h_alpha[0]*hafactor
;;       uucor[ha].h_alpha[1] = uucor[ha].h_alpha[1]*hafactor
;;;      im_plothist, uu[ha].h_alpha_ew[0]/ewhacor, xr=[0,10], bin=0.1
;;    endif
    
    ubpt = iclassification(uucor[m1[uthese]],snrcut_class=bpt_snrcut1,$
      niihacut=-0.3,chi2cut=chi2cut1,ratios=uratios,doplot=0)

    uunk =   where(strtrim(uratios.final_class,2) eq 'UNKNOWN')
    usf =    where(strtrim(uratios.final_class,2) eq 'SF')
    uagn =   where(strtrim(uratios.final_class,2) eq 'AGN')
    usfagn = where(strtrim(uratios.final_class,2) eq 'SF/AGN')
    help, usf, uagn, usfagn, uthese

    class[m2[uthese[uunk]]].bpt = 0
    class[m2[uthese[uagn]]].bpt = 2
    class[m2[uthese[usf]]].bpt = 4
    class[m2[uthese[usfagn]]].bpt = 8

; #########################
; [Ne V]
; #########################

    these = where(ii[m1].nev_3426[1] gt 0.0,nthese)
;   djs_plot, ii[m1[these]].nev_3426_ew[0], ii[m1[these]].nev_3426[0]/ii[m1[these]].nev_3426[1], $
;     ps=6, /xlog, /ylog, xr=[0.001,1000], yr=[0.001,1000], sym=0.2
;   djs_oplot, nev_ewcut1*[1,1], 10^!y.crange, line=0, color='green'
;   djs_oplot, 10^!x.crange, nev_snrcut1*[1,1], line=0, color='green'

    nev_agn = where((ii[m1[these]].nev_3426[0]/ii[m1[these]].nev_3426[1] gt nev_snrcut1) and $
      (ii[m1[these]].nev_3426_ew[0] gt nev_ewcut1) and (ii[m1[these]].nev_3426_chi2 le chi2cut1) and $
      (ii[m1[these]].nev_3426_sigma[0] lt sigmacut1))
;   ages_display_spectrum, aa[m1[these[agn[0:10]]]], labeltype=5, xr=[3300,3600], plottype=3

;   class[m2[these]].nev = 0
;   class[m2[these[nev_agn]]].nev = 2

; #########################
; broad
; #########################

;;    jj = where((aa[m1].pass ne 311) and (aa[m1].pass ne 106) and (aa[m1].pass ne 110) and $
;;      (aa[m1].pass ne 209) and (aa[m1].pass ne 310) and $
;;      (uu[m1].h_alpha[0]/uu[m1].h_alpha[1] gt broad_snrcut1) and (uu[m1].h_alpha_sigma[0] ge broad_sigmacut1) and $
;;      (uu[m1].h_beta[0]/uu[m1].h_beta[1] gt broad_snrcut1) and (uu[m1].h_beta_sigma[0] ge broad_sigmacut1) and $
;;;     (uu[m1].oiii_4959[0]/uu[m1].oiii_4959[1] gt broad_snrcut1) and $
;;      (uu[m1].oiii_5007[0]/uu[m1].oiii_5007[1] gt broad_snrcut1))

;     (aa[m1].pass ne 311) and $
;     (aa[m1].pass ne 106) and (aa[m1].pass ne 110) and $
;     (aa[m1].pass ne 209) and (aa[m1].pass ne 310) and $

; identify all the broad-line objects; since we're using the
; unfluxed spectra here, we can also classify plates 106, 110, 209,
; 310, and 311    
    
    broad = where($
;     (aa[m1].main_flag eq 1) and $
; H-alpha + H-beta
      (((uu[m1].h_alpha[0]/uu[m1].h_alpha[1] gt broad_snrcut1) and (uu[m1].h_alpha_sigma[0] ge broad_sigmacut1) and $
      (uu[m1].h_beta[0]/uu[m1].h_beta[1] gt broad_snrcut1) and (uu[m1].h_beta_sigma[0] ge broad_sigmacut1) and $
      (uu[m1].h_alpha_sigma[0] lt max_broad_sigma1) and (uu[m1].h_beta_sigma[0] lt max_broad_sigma1)) or $
; H-beta + Mg II
      ((uu[m1].h_beta[0]/uu[m1].h_beta[1] gt broad_snrcut1) and (uu[m1].h_beta_sigma[0] ge broad_sigmacut1) and $
      (uu[m1].mgii_2800[0]/uu[m1].mgii_2800[1] gt broad_snrcut1) and (uu[m1].mgii_2800_sigma[0] ge broad_sigmacut1) and $
      (uu[m1].h_beta_sigma[0] lt max_broad_sigma1) and (uu[m1].mgii_2800_sigma[0] lt max_broad_sigma1)) or $
; Mg II alone
      ((uu[m1].mgii_2800[0]/uu[m1].mgii_2800[1] gt broad_snrcut1) and (uu[m1].mgii_2800_sigma[0] ge broad_sigmacut1) and $
      (uu[m1].mgii_2800_sigma[0] lt max_broad_sigma1)) or $
; Mg II + C III
      ((uu[m1].mgii_2800[0]/uu[m1].mgii_2800[1] gt broad_snrcut1) and (uu[m1].mgii_2800_sigma[0] ge broad_sigmacut1) and $
      (uu[m1].ciii_1908[0]/uu[m1].ciii_1908[1] gt broad_snrcut1) and (uu[m1].ciii_1908_sigma[0] ge broad_sigmacut1) and $
      (uu[m1].mgii_2800_sigma[0] lt max_broad_sigma1) and (uu[m1].ciii_1908_sigma[0] lt max_broad_sigma1))))
    class[m2[broad]].broad = 16

;   ages_display_spectrum, uu[m1[broad]], labeltype=2, plottype=2, /unfluxed, $
;     /postscript, psname=outpath+'ages_agn_broad.ps', /pdf
    
;;; for all z, classified as QSO by Schlegel, and non broad-line AGN by me
;;    w1 = where((class[m2].broad gt 0) and (strtrim(aa[m1].class,2) eq 'QSO'))
;;; for z>1, classified as QSO by Schlegel, and non broad-line AGN by
;;; me; these are all quasars, so classify them as such
;;    w2 = where((class[m2].broad lt 0) and (strtrim(aa[m1].class,2) eq 'QSO') and (aa[m1].z gt 1.0))
;;    ages_display_spectrum, uu[m1[w2]], labeltype=5, plottype=2, /unfluxed, $
;;      /postscript, psname='ages_broad0.schlegel.qso.z.gt.1.ps', /pdf
;;; for z<1, classified as QSO by Schlegel, and non broad-line AGN by
;;; me; visual-inspect all of these
;;    w3 = where((class[m2].broad lt 0) and (strtrim(aa[m1].class,2) eq 'QSO') and (aa[m1].z lt 1.0))
;;    ages_display_spectrum, uu[m1[w3]], labeltype=5, plottype=2, /unfluxed, $
;;      /postscript, psname='ages_broad0.schlegel.qso.z.lt.1.ps', /pdf
;;; for z<1, classified as non-QSO by Schlegel, and broad-line AGN by
;;; me; visual-inspect all of these
;;    w4 = where((class[m2].broad gt 0) and (strtrim(aa[m1].class,2) ne 'QSO') and (aa[m1].z lt 1.0))
;;    ages_display_spectrum, uu[m1[w4]], labeltype=5, plottype=2, /unfluxed, $
;;      /postscript, psname='ages_broad1.schlegel.galaxy.z.lt.1.ps', /pdf
;;
;;; visual inspected: ages_broad0.schlegel.qso.z.lt.1.pdf
;;    visual_broad = ['ages_'+$
;;      '417/247',$
;;      '309/112',$
;;      '606/086',$
;;      '710/093',$
;;      '508/077',$
;;      '613/238',$
;;      '112/232',$
;;      '423/055',$
;;      '102/229',$
;;      '210/248',$
;;      '312/227',$
;;      '213/207',$
;;      '211/249',$
;;      '114/004',$
;;      '301/121',$
;;      '109/158',$
;;      '309/254',$
;;      '611/234',$
;;      '107/282']
;;    class[where_array(class.galaxy,visual_broad)].broad = 2
    
;   ages_display_spectrum, uu[m1[w3[0:5]]], labeltype=5, plottype=2, /unfluxed
;;  ages_display_spectrum, uu[m1[w3]], labeltype=5, plottype=2, /unfluxed, /postscript, psname='junk.ps', /pdf
    
; how many broadline AGN at z<1 were classified as 'GALAXY'? - only 4!
;   w1 = where((strtrim(aa[m1[broad]].class,2) eq 'GALAXY') and (aa[m1[broad]].z lt 1.0))
;   ages_display_spectrum, uu[m1[broad[w1]]], labeltype=5, plottype=2, /unfluxed

;;; how well do the QSO and broadline classifications agree - very well!
;;    junk = where((aa[m1].pass ne 311) and $
;;      (aa[m1].pass ne 106) and (aa[m1].pass ne 110) and $
;;      (aa[m1].pass ne 209) and (aa[m1].pass ne 310) and $
;;      (strtrim(aa[m1].class,2) eq 'QSO') and (aa[m1].z lt 1.0))
;;    w2 = cmset_op(broad,'and',junk,/not2) ; - 5 
;;    w2 = cmset_op(junk,'and',broad,/not2) ; - 95
;;    ages_display_spectrum, uu[m1[w2[0:5]]], labeltype=5, plottype=2, /unfluxed

; write out

    outfile = outpath+'catalog.moustakas.agn.'+vv+'.fits'
    splog, 'Writing '+outfile
    mwrfits, class, outfile, /create

    outfile = outpath+'catalog.moustakas.agn.'+vv
    splog, 'Writing '+outfile

    openw, lun, outfile, /get_lun
    printf, lun, '# J. Moustakas, NYU'
    printf, lun, '# Contact john.moustakas@gmail.com with questions'
    printf, lun, '# '+im_today()
    printf, lun, '# '
    printf, lun, '# Spectrophotometric and multiwavelength classification of'
    printf, lun, '# AGN and star-forming galaxies based on the '+version+' AGES'
    printf, lun, '# emission-line catalog.  Plates 106, 110, 209, 310, and 311 '
    printf, lun, '# have not been flux-calibrated; therefore, the '+version+' catalog'
    printf, lun, '# of unfluxed plates was used.  The first five columns identify the AGES'
    printf, lun, '# spectrum used and the adopted redshift; the subsequent'
    printf, lun, '# columns contain the results of various classification schemes:'
    printf, lun, '# '
    printf, lun, '#   +BPT: [NII]/Ha vs [OIII]/Hb diagnostic diagram; objects were '
    printf, lun, '# classified if all four lines were detected with S/N>3 and line-width '
    printf, lun, '# sigma<500 km/s as follows: star-forming (SF): below and to the left '
    printf, lun, '# of the Kauffmann et al. (2003) boundary; AGN: above and to the right '
    printf, lun, '# of the Kewley et al. (2001); between the two curves: composite, SF/AGN.'
    printf, lun, '# '
    printf, lun, '#   +BROAD: Broad-line (type 1) AGN were identified using the following '
    printf, lun, '# criteria applied to lines with S/N>10 and line-width 300<sigma<5000 km/s: '
    printf, lun, '# (1) H-alpha AND H-beta; OR (2) H-beta AND Mg II 2800; OR; (3) Mg II 2800; '
    printf, lun, '# OR (4) Mg II 2800 AND CIII 1908'
    printf, lun, '# '
;   printf, lun, '#   +BPT: [NII]/Ha vs [OIII]/Hb diagnostic diagram; objects were'
;   printf, lun, '# classified using the BPT diagram and the Kauffmann et al.'
;   printf, lun, '# (2003) boundary if all four lines were detected with S/N>3; '
;   printf, lun, '# alternatively, if [NII] and Ha were detected with S/N>3 and'
;   printf, lun, '# log([NII]/Ha)>-0.3 then the object was classified as an AGN.'
;   printf, lun, '# '
;   printf, lun, '# '
;   printf, lun, '#   +NEV: If [Ne V] 3426 was detected with S/N>5 and EW([NeV])>1 A'
;   printf, lun, '# (rest) then the object was classified as an AGN.'
;   printf, lun, '# '
    printf, lun, '# Notes:'
    printf, lun, '# '
    printf, lun, '#   + A variety of tests have been performed to ensure that these '
    printf, lun, '# results are as reliable as possible; however, with a sample as '
    printf, lun, '# large as this it is very difficult to ensure 100% fidelity.  Thus, '
    printf, lun, '# visual inspection of the spectra is highly recommended.  In the '
    printf, lun, '# future, additional techniques for identifying AGN (using, e.g., '
    printf, lun, '# X-ray and mid-infrared colors) will be implemented.'
    printf, lun, '# '
    printf, lun, '# The classes are stored as integers with the following meaning:'
    printf, lun, '# '
    printf, lun, '#    -1: no classification attempted (insufficient data)'
    printf, lun, '#     0: inconclusive classification'
    printf, lun, '#     2: AGN'
    printf, lun, '#     4: Star-Forming (SF)'
    printf, lun, '#     8: SF/AGN Composite'
    printf, lun, '#    16: Broad-line AGN'
    printf, lun, '# '
    printf, lun, '# 1 ra'
    printf, lun, '# 2 dec'
    printf, lun, '# 3 z'
    printf, lun, '# 4 pass'
    printf, lun, '# 5 aper'
    printf, lun, '# 6 bpt'
    printf, lun, '# 7 broad'
;   printf, lun, '# 7 nev'
    struct_print, class, format=format, lun=lun, /no_head
    free_lun, lun
;   for jj = 0L, nobj-1L do printf, lun, 
;   README

stop    
    

return
end
    
