;+
; NAME:
;       SINGS_ANCILLARY_DATA
;
; PURPOSE:
;       Compile ancillary data for SINGS.
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jul 25, U of A
;-

pro sings_ancillary_data, table, write=write

    version = sings_version(/ancillary)
    datapath = sings_path(/analysis)
    outpath = datapath

    outname = 'sings_ancillary_data_'+version+'.fits'
    
    basicname = 'sings_ned.fits.gz'
    photoname = 'sings_ned_photo.fits.gz'
    distname = 'sings_distances.fits.gz'
    diamname = 'sings_diameters.fits.gz'

;   leda = sings_read_leda()
    
    write_ancillary_data, table, datapath=datapath, outpath=outpath, $
      basicname=basicname, photoname=photoname, distname=distname, $
      diamname=diamname, leda=leda, outname=outname;, write=write
    ngalaxy = n_elements(table)

    indx = speclinefit_locate(table,'ngc1266') ; no RC3 photometry
    table[indx].rc3_ubv[1] = 14.00
    table[indx].rc3_ubv_err[1] = 0.36
    table[indx].rc3_ubv_ref[1] = 'LEDA'

    dmod = im_d2dmod(table[indx].distance,err_dmod=err_dmod,$
      err_dist=table[indx].distance_err)
    table[indx].rc3_ubv_absmag[1] = table[indx].rc3_ubv[1] - dmod
    table[indx].rc3_ubv_absmag_err[1] = sqrt(table[indx].rc3_ubv_err[1]^2.0 + err_dmod^2.0)

; add gas mass tags and UBVRIJHKs photometry tags  

    moretags = {$
      log_mhi:                   -999.0, $
      log_mh2:                   -999.0, $
      log_mdust:                 -999.0, $
      log_ldust:                 -999.0, $
      qpah:                      -999.0, $
      umean:                     -999.0, $
      mhi_ref:                    '...', $
      mh2_ref:                    '...', $
      bvri_ref:                   '...', $
      bvri:             fltarr(4)-999.0, $
      bvri_err:         fltarr(4)-999.0, $
      bvri_absmag:      fltarr(4)-999.0, $
      bvri_absmag_err:  fltarr(4)-999.0, $
      twomass_ref:                   '...', $
      twomass:             fltarr(3)-999.0, $
      twomass_err:         fltarr(3)-999.0, $
      twomass_absmag:      fltarr(3)-999.0, $
      twomass_absmag_err:  fltarr(3)-999.0}
    table = struct_addtags(table,replicate(moretags,ngalaxy))

; start with the Leroy+08 THINGS measurements

    leroy = rsex(datapath+'08leroy_table4.sex')
    match, strtrim(strupcase(table.galaxy),2), strtrim(strupcase(leroy.galaxy),2), m1, m2

; M(HI); no upper limits
    these = where((table[m1].log_mhi lt -900.0) and (leroy[m2].log_mhi gt -900.0),nthese)
    struct_print, leroy[m2[these]]
    dist_corr = 2.0*alog10(table[m1[these]].distance/leroy[m2[these]].distance) ; distance correction
    table[m1[these]].log_mhi = leroy[m2[these]].log_mhi + dist_corr
    table[m1[these]].mhi_ref = 'leroy08a'
; M(H2); deal with 5-sigma upper limits
    these = where((table[m1].log_mh2 lt -900.0) and (leroy[m2].log_mh2 gt -900.0),nthese)
    struct_print, leroy[m2[these]]
    pos = where(leroy[m2[these]].log_mh2 gt 0.0,comp=neg)
    table[m1[these[pos]]].log_mh2 =  leroy[m2[these[pos]]].log_mh2 + dist_corr
    table[m1[these[neg]]].log_mh2 = -leroy[m2[these[neg]]].log_mh2 + dist_corr
    table[m1[these]].mh2_ref = 'leroy08a'

;   ww = where(table.log_mhi lt -900.0)
;   struct_print, struct_trimtags(table,select=['galaxy',$
;     'log_mhi','mhi_ref','log_mh2','mh2_ref'])
    
; next use the Draine+07 compilation; don't copy over objects
; that already have masses from Leroy+08

    draine = rsex(datapath+'07draine.sex')
    match, strtrim(strupcase(table.galaxy),2), strtrim(strupcase(draine.galaxy),2), m1, m2
    niceprint, table[m1].galaxy, draine[m2].galaxy

;   draine_copy = struct_trimtags(draine[m2],except=['GALAXY','DISTANCE','TYPE'])
;   indx1 = speclinefit_locate(table,draine.galaxy)
;   these = where(table[indx1].log_mhi lt -900.0)
;   indx1 = indx1[these]
;   junk = table[m1] & struct_assign, draine_copy[m2], junk, /nozero & table[m1] = junk

    dist_corr = 2.0*alog10(table[m1].distance/draine[m2].distance) ; distance correction

; M(HI)    
    neg = where((table[m1].log_mhi lt -900.0) and (draine[m2].log_mhi gt -900.0) and (draine[m2].log_mhi lt 0.0),nneg)
    pos = where((table[m1].log_mhi lt -900.0) and (draine[m2].log_mhi gt 0.0),npos)
    table[m1[pos]].log_mhi = draine[m2[pos]].log_mhi + dist_corr[pos]
    table[m1[neg]].log_mhi = -(abs(draine[m2[neg]].log_mhi) + dist_corr[neg])
    if (npos gt 0) then table[m1[pos]].mhi_ref = 'draine07a'
    if (nneg gt 0) then table[m1[neg]].mhi_ref = 'draine07a'

; M(H2)    
    neg = where((table[m1].log_mh2 lt -900.0) and (draine[m2].log_mh2 gt -900.0) and (draine[m2].log_mh2 lt 0.0),nneg)
    pos = where((table[m1].log_mh2 lt -900.0) and (draine[m2].log_mh2 gt 0.0),npos)
    table[m1[pos]].log_mh2 = draine[m2[pos]].log_mh2 + dist_corr[pos]
    table[m1[neg]].log_mh2 = -(abs(draine[m2[neg]].log_mh2) + dist_corr[neg])
    if (npos gt 0) then table[m1[pos]].mh2_ref = 'draine07a'
    if (nneg gt 0) then table[m1[neg]].mh2_ref = 'draine07a'

; some other fun quantites    
    good = where((draine[m2].log_mdust gt -900.0),ngood)
    table[m1[good]].log_mdust = draine[m2[good]].log_mdust + dist_corr[good]
    good = where((draine[m2].log_ldust gt -900.0),ngood)
    table[m1[good]].log_ldust = draine[m2[good]].log_ldust + dist_corr[good]

;   ww = where(table.log_mhi lt -900.0)
;   struct_print, struct_trimtags(table,select=['galaxy',$
;     'log_mhi','mhi_ref','log_mh2','mh2_ref'])

;; then add the Walter+07 (private communication) gas measurements from
;; THINGS (superseded by Leroy+08, above)
;
;    walt = rsex(datapath+'07walter_things.sex')
;    indx = speclinefit_locate(table,walt.galaxy)
;
;    these = where((table[indx].log_mhi lt -900.0) and (walt.mhi gt -900.0),nthese)
;    if (nthese ne 0L) then begin
;       struct_print, walt[these]
;       dist_corr = 2.0*alog10(table[indx[these]].distance/walt[these].distance) ; distance correction
;       table[indx[these]].log_mhi = alog10(walt[these].mhi*1D8) + dist_corr
;       table[indx[these]].mhi_ref = 'walter07_priv'
;    endif

; then add the Walter+07 gas measurements

    walt_gal = ['IC2574','HOLMBERGII','HOLMBERGI','DDO165','DDO053','M81DWB','M81DWA']
    walt_dist = [4.02,3.39,3.84,4.57,3.56,5.3,3.55]
    walt_mhi = [14.75,5.95,1.40,1.72,0.60,0.25,0.12]
    walt_indx = speclinefit_locate(table,walt_gal)

    walt_these = where(table[walt_indx].log_mhi lt -900.0,nthese)
    if (nthese ne 0L) then begin
       dist_corr = 2.0*alog10(table[walt_indx[walt_these]].distance/walt_dist[walt_these]) ; distance correction
       table[walt_indx[walt_these]].log_mhi = alog10(walt_mhi[walt_these]*1D8) + dist_corr
       table[walt_indx[walt_these]].mhi_ref = 'walter07a'
    endif
    
; finally supplement with the Kennicutt+03 gas masses

    kenn = rsex(datapath+'2003_kennicutt_table2.sex')
    kenn = struct_trimtags(kenn,select=['GALAXY','DISTANCE','LOG_MHI','LOG_MH2'])
;   struct_print, kenn

    these = where((table.log_mhi lt -900.0) and ((kenn.log_mhi gt -900.0) or $
      (kenn.log_mhi gt -900.0)),nthese)
    if (nthese ne 0L) then begin
;      struct_print, kenn[these]
       dist_corr = 2.0*alog10(table[these].distance/kenn[these].distance) ; distance correction
; M(HI)
       pos = where((kenn[these].log_mhi gt -900.0) and (kenn[these].log_mhi gt 0.0),npos)
       if (npos ne 0L) then table[these[pos]].log_mhi = kenn[these[pos]].log_mhi + dist_corr[pos]
       neg = where((kenn[these].log_mhi gt -900.0) and (kenn[these].log_mhi lt 0.0),nneg)
       if (nneg ne 0L) then table[these[neg]].log_mhi = -(abs(kenn[these[neg]].log_mhi) + dist_corr[neg])
       if (npos gt 0) then table[these[pos]].mhi_ref = 'kenn03b'
       if (nneg gt 0) then table[these[neg]].mhi_ref = 'kenn03b'

; M(H2)
       pos = where((kenn[these].log_mh2 gt -900.0) and (kenn[these].log_mh2 gt 0.0),npos)
       if (npos ne 0L) then table[these[pos]].log_mh2 = kenn[these[pos]].log_mh2 + dist_corr[pos]
       neg = where((kenn[these].log_mh2 gt -900.0) and (kenn[these].log_mh2 lt 0.0),nneg)
       if (nneg ne 0L) then table[these[neg]].log_mh2 = -(abs(kenn[these[neg]].log_mh2) + dist_corr[neg])
       if (npos gt 0) then table[these[pos]].mh2_ref = 'kenn03b'
       if (nneg gt 0) then table[these[neg]].mh2_ref = 'kenn03b'
    endif

; Hyperleda values; m21 = -2.5 log f + 17.40 
;   hyper_gal = ['NGC0584',
;   hyper_himag = [16.94,
;   hyper_himag_err = [0.39
;   hyper_mhi = [
;   hyper_indx = speclinefit_locate(table,hyper_gal)
;
;   hyper_these = where(table[hyper_indx].log_mhi lt -900.0,nthese)
;   if (nthese ne 0L) then begin
;      dist_corr = 2.0*alog10(table[hyper_indx[hyper_these]].distance/hyper_dist[hyper_these]) ; distance correction
;      table[hyper_indx[hyper_these]].log_mhi = alog10(hyper_mhi[hyper_these]*1D8) + dist_corr
;      table[hyper_indx[hyper_these]].mhi_ref = 'hyperleda'
;   endif
    
;   ww = where(table.log_mhi lt -900.0)
;   struct_print, struct_trimtags(table,select=['galaxy',$
;     'log_mhi','mhi_ref','log_mh2','mh2_ref'])

; read and parse the Dale et al. 2006 photometry    

    flux = rsex(datapath+'sings_photometry_2006dale.sex')
    ferr = rsex(datapath+'sings_photometry_2006dale.uncertainty.sex')

    flux.uv_1528    = flux.uv_1528*1D-3    ; [mJy-->Jy]
    flux.uv_2271    = flux.uv_2271*1D-3    ; [mJy-->Jy]
    flux.radio_20cm = flux.radio_20cm*1D-3 ; [mJy-->Jy]
    
    ferr.uv_1528    = ferr.uv_1528*1D-3    ; [mJy-->Jy]
    ferr.uv_2271    = ferr.uv_2271*1D-3    ; [mJy-->Jy]
    ferr.radio_20cm = ferr.radio_20cm*1D-3 ; [mJy-->Jy]

; IRAC photometry

    irac_Jy = [[flux.irac_3_6],[flux.irac_4_5],[flux.irac_5_8],[flux.irac_8]]
    irac_Jy_err = [[ferr.irac_3_6],[ferr.irac_4_5],[ferr.irac_5_8],[ferr.irac_8]]

    good = where(irac_Jy_err gt 0.0,comp=crap)
    irac = irac_Jy*0.0
    irac_err = irac_Jy_err*0.0

    irac[good] = -2.5*alog10(irac_Jy[good]*1D-23)-48.6              ; [AB]
    irac_err[good] = 2.5*irac_Jy_err[good]/irac_Jy[good]/alog(10.0) ; [AB]
    irac[crap] = -999.0 & irac_err[crap] = -999.0
    
    table.irac_ref = 'dale07a'
    table.irac = transpose(irac)
    table.irac_err = transpose(irac_err)

    for iband = 0L, 3L do begin
       good = where((table.distance gt -900.0) and (table.irac[iband] gt -900.0),ngood)
       if (ngood ne 0L) then begin
          dmod = im_d2dmod(table[good].distance,err_dmod=err_dmod,$
            err_dist=table[good].distance_err)
          table[good].irac_absmag[iband] = table[good].irac[iband] - dmod
          table[good].irac_absmag_err[iband] = sqrt(table[good].irac_err[iband]^2.0 + err_dmod^2.0)
       endif
    endfor

; MIPS photometry
    
    mips_Jy = [[flux.mips_24],[flux.mips_70],[flux.mips_160]]
    mips_Jy_err = [[ferr.mips_24],[ferr.mips_70],[ferr.mips_160]]

    good = where(mips_Jy_err gt 0.0,comp=crap)
    mips = mips_Jy*0.0
    mips_err = mips_Jy_err*0.0

    mips[good] = -2.5*alog10(mips_Jy[good]*1D-23)-48.6              ; [AB]
    mips_err[good] = 2.5*mips_Jy_err[good]/mips_Jy[good]/alog(10.0) ; [AB]
    mips[crap] = -999.0 & mips_err[crap] = -999.0
    
    table.mips_ref = 'dale07a'
    table.mips = transpose(mips)
    table.mips_err = transpose(mips_err)

    for iband = 0L, 2L do begin
       good = where((table.distance gt -900.0) and (table.mips[iband] gt -900.0),ngood)
       if (ngood ne 0L) then begin
          dmod = im_d2dmod(table[good].distance,err_dmod=err_dmod,$
            err_dist=table[good].distance_err)
          table[good].mips_absmag[iband] = table[good].mips[iband] - dmod
          table[good].mips_absmag_err[iband] = sqrt(table[good].mips_err[iband]^2.0 + err_dmod^2.0)
       endif
    endfor
    
; GALEX photometry

    galex_Jy = [[flux.uv_1528],[flux.uv_2271]]
    galex_Jy_err = [[ferr.uv_1528],[ferr.uv_2271]]

    good = where(galex_Jy_err gt 0.0,comp=crap)
    galex = galex_Jy*0.0
    galex_err = galex_Jy_err*0.0

    galex[good] = -2.5*alog10(galex_Jy[good]*1D-23)-48.6              ; [AB]
    galex_err[good] = 2.5*galex_Jy_err[good]/galex_Jy[good]/alog(10.0) ; [AB]
    galex[crap] = -999.0 & galex_err[crap] = -999.0
    
    table.galex_ref = 'dale07a'
    table.galex = transpose(galex)
    table.galex_err = transpose(galex_err)

    for iband = 0L, 1L do begin
       good = where((table.distance gt -900.0) and (table.galex[iband] gt -900.0),ngood)
       if (ngood ne 0L) then begin
          dmod = im_d2dmod(table[good].distance,err_dmod=err_dmod,$
            err_dist=table[good].distance_err)
          table[good].galex_absmag[iband] = table[good].galex[iband] - dmod
          table[good].galex_absmag_err[iband] = sqrt(table[good].galex_err[iband]^2.0 + err_dmod^2.0)
       endif
    endfor
    
;   filterlist = ['bessell_'+['B','V','R','I'],'twomass_'+['J','H','Ks']]+'.par'

    filterlist = ['bessell_'+['B','V','R','I']]+'.par'
    weff = k_lambda_eff(filterlist=filterlist)

    vega_Jy_dale = [4270.0,3670.0,2840.0,2550.0] ; [Vega flux densities, Jy]
    vega_Jy_ned = [4260.0,3640.0,2890.0,2280.0]  ; [http://nedwww.ipac.caltech.edu/help/photoband.lst]

    vega2ab_dale = -2.5*alog10(vega_Jy_dale*1D-23)-48.6 ; (see below)
    vega2ab_ned = -2.5*alog10(vega_Jy_ned*1D-23)-48.6   ; (see below)
    vega2ab = k_vega2ab(filterlist=filterlist,/kurucz)

;   vega2ab_dale = vega2ab_dale*0.0
;   vega2ab_ned = vega2ab_ned*0.0
;   vega2ab = vega2ab*0.0
    
    amcor = [0.25,0.14,0.08,0.053] ; air mass correction terms (see below)

; ###########################################################################
; BVRI photometry; as of 2007-Sep-31 the airmass correction is
; obsolete, as Danny decided to publish an erratum to his optical
; fluxes (see email archive)
; ###########################################################################
;
; Danny's extinction correction was made using the Li & Draine (2001)
; extinction curve; I use O'Donnell (1986) here, but note that this is
; an insignificant correction


; The Vega-->AB conversion stuff is a mess.  The STScI data analyst
; used the Vega_Jy_Dale values, above, but these lead to *very*
; different vega2ab values compared to K-correct (see above).
;
; IDL> niceprint, vega2ab_dale, vega2ab_ned, vega2ab
;    -0.17606816       -0.17352247     -0.0923588
;   -0.011663635     -0.0027519332      0.0184169
;     0.26670568        0.24775692       0.201340
;     0.38365107        0.50516441       0.451995
;
; This matters because using either vega2ab_dale or vega2ab_ned leads
; to (B-V)_Vega colors that are way too red (e.g., on the red
; sequence).  However, they used these values to convert the RC3
; BV_Vega magnitudes to AB.  So I have to convert the RC3 magnitudes
; from AB-->Vega using vega2ab_dale, but for the rest of the
; sample I'm going to use vega2ab
;
; One final note: here are two equivalent ways of getting the Vega
; magnitude: 
;     
; IDL> bmag_ab = -2.5*alog10(flux[good].b*1D-23) - 48.6
; IDL> bmag_vega_1 = bmag_ab - vega2ab_dale[0]
; IDL> bmag_vega_2 = -2.5*alog10(flux[good].b/vega_Jy_dale[0])
    
    table.bvri_ref = 'dale07a'
    
    good = where((flux.b gt 0.0) and (ferr.b gt 0.0),ngood)
    if (ngood ne 0L) then begin
       bmag_vega_dale = -2.5*alog10(flux[good].b*1D-23) - 48.6 - vega2ab_dale[0]
       bmag_vega = -2.5*alog10(flux[good].b*1D-23) - 48.6 - vega2ab[0]
       table[good].bvri[0] = bmag_vega + table[good].ebv_mw * $
         (k_lambda(4400.0,/li)-k_lambda(weff[0],/odonnell))
       table[good].bvri_err[0] = 2.5*ferr[good].b/flux[good].b/alog(10.0)
; now use the Dale Vega magnitude for the RC3 galaxies so that it
; matches the published RC3 values
       bindx_rc3 = speclinefit_locate(table[good],'NGC'+['0855',$
         '2915','3031','3621','4569','5055','5408','6946'])
       table[good[bindx_rc3]].bvri[0] = bmag_vega_dale[bindx_rc3] + $
         table[good[bindx_rc3]].ebv_mw * (k_lambda(4400.0,/li)-k_lambda(weff[0],/odonnell))
    endif
;   niceprint, table.galaxy, bmag_vega_dale, table.bvri[0], table.rc3_ubv[1]

    good = where((flux.v gt 0.0) and (ferr.v gt 0.0),ngood)
    if (ngood ne 0L) then begin
       vmag_vega_dale = -2.5*alog10(flux[good].v*1D-23) - 48.6 - vega2ab_dale[1]
       vmag_vega = -2.5*alog10(flux[good].v*1D-23) - 48.6 - vega2ab[1]
       table[good].bvri[1] = vmag_vega + table[good].ebv_mw * $
         (k_lambda(5500.0,/li)-k_lambda(weff[1],/odonnell))
       table[good].bvri_err[1] = 2.5*ferr[good].v/flux[good].v/alog(10.0)
; now use the Dale Vega magnitude for the RC3 galaxies so that it
; matches the published RC3 values
       vindx_rc3 = speclinefit_locate(table[good],'NGC'+['0855',$
         '3031','3034','4569','5055','5408'])
       table[good[vindx_rc3]].bvri[1] = vmag_vega_dale[vindx_rc3] + $
         table[good[vindx_rc3]].ebv_mw * (k_lambda(5500.0,/li)-k_lambda(weff[1],/odonnell))
    endif
;   niceprint, table.galaxy, vmag_vega_dale, table.bvri[1], table.rc3_ubv[2]

; interlude to deal with a problem; for some galaxies mixing the Dale
; and RC3 photometry leads to funky results; for example, the B-V
; color of NGC3621 (B=RC3, V=Dale) is ~0.93, but if I use RC3 for both
; B and V then the color is ~0.54, which is *much* more consistent
; with its morphology; the other problem galaxies are: NGC2915,
; NGC6946, and NGC3034

;   w1 = speclinefit_locate(table,'NGC'+['0855',$
;     '2915','3031','3621','4569','5055','5408','6946'])
;   w2 = speclinefit_locate(table,'NGC'+['0855','3031','3034','4569','5055','5408'])
;   niceprint, table[w1].galaxy, $
;;    table[w1].bvri[0], table[w1].rc3_ubv[1], table[w1].bvri[1], table[w1].rc3_ubv[2], $
;     table[w1].rc3_ubv[1]-table[w1].rc3_ubv[2], table[w1].bvri[0]-table[w1].bvri[1]
; NGC0855       0.637841       0.630604
; NGC2915       0.290405       0.191121
; NGC3031       0.868407       0.864376
; NGC3621       0.538005       0.931025
; NGC4569       0.672617       0.667243
; NGC5055       0.702169       0.693172
; NGC5408       0.490538       0.484311
; NGC6946       0.451898       0.713748
;
;   niceprint, table[w2].galaxy, $
;;    table[w2].bvri[0], table[w2].rc3_ubv[1], table[w2].bvri[1], table[w2].rc3_ubv[2], $
;     table[w2].rc3_ubv[1]-table[w2].rc3_ubv[2], table[w2].bvri[0]-table[w2].bvri[1]
; NGC0855       0.637841       0.630604
; NGC3031       0.868407       0.864376
; NGC3034       0.730925     -0.0625954
; NGC4569       0.672617       0.667243
; NGC5055       0.702169       0.693172
; NGC5408       0.490538       0.484311
;
; so, for the four galaxies mention above I'm going to use the
; RC3 photometry for both B and V; I'm going to cheat, though,
; and keep the Dale *uncertainties*

    fixme = speclinefit_locate(table,'NGC'+['2915','3621','3034','6946'])
    table[fixme].bvri[0] = table[fixme].rc3_ubv[1]
    table[fixme].bvri[1] = table[fixme].rc3_ubv[2]
    niceprint, table[fixme].galaxy, table[fixme].bvri[0], table[fixme].rc3_ubv[1], $
      table[fixme].bvri[1], table[fixme].rc3_ubv[2]

; now deal with R and I    
    
    good = where((flux.r gt 0.0) and (ferr.r gt 0.0),ngood)
    if (ngood ne 0L) then begin
       table[good].bvri[2] = -2.5*alog10(flux[good].r*1D-23) - 48.6 - vega2ab[2] + $
         table[good].ebv_mw*(k_lambda(6600.0,/li)-k_lambda(weff[2],/odonnell))
       table[good].bvri_err[2] = 2.5*ferr[good].r/flux[good].r/alog(10.0)
    endif
    
    good = where((flux.i gt 0.0) and (ferr.i gt 0.0),ngood)
    if (ngood ne 0L) then begin
       table[good].bvri[3] = -2.5*alog10(flux[good].i*1D-23) - 48.6 - vega2ab[3] + $
         table[good].ebv_mw*(k_lambda(8100.0,/li)-k_lambda(weff[3],/odonnell))
       table[good].bvri_err[3] = 2.5*ferr[good].i/flux[good].i/alog(10.0)
    endif

; compute absolute magnitudes
    for iband = 0L, 3L do begin
       good = where((table.distance gt -900.0) and (table.bvri[iband] gt -900.0),ngood)
       if (ngood ne 0L) then begin
          dmod = im_d2dmod(table[good].distance,err_dmod=err_dmod,$
            err_dist=table[good].distance_err)
          table[good].bvri_absmag[iband] = table[good].bvri[iband] - dmod
          table[good].bvri_absmag_err[iband] = sqrt(table[good].bvri_err[iband]^2.0 + err_dmod^2.0)
       endif
    endfor

; ###########################################################################    
    
; make a plot    

    plotsym, 0, 1, /fill
    im_openclose, datapath+'dale_vs_rc3.ps', xsize=8.5, ysize=8.0, /postscript
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal
;   dfpsplot, /square, /color
    bw = speclinefit_locate(table,'NGC'+['0855','2915','3031','3621','4569','5055','5408','6946'])
    vw = speclinefit_locate(table,'NGC'+['0855','3031','3034','4569','5055','5408'])
    splog, 'B-band'
    niceprint, table[bw].galaxy, table[bw].rc3_ubv[1]+table[bw].ebv_mw*k_lambda(weff[0],/odonnell), $
      table[bw].rc3_ubv[1], table[bw].bvri[0], table[bw].rc3_ubv[1]-table[bw].bvri[0]
    splog, 'V-band'
    niceprint, table[vw].galaxy, table[vw].rc3_ubv[2]+table[vw].ebv_mw*k_lambda(weff[1],/odonnell), $
      table[vw].rc3_ubv[2], table[vw].bvri[1], table[vw].rc3_ubv[2]-table[vw].bvri[1]
; --------------------------------------------------
    good = where((table.bvri[0] gt -900.0) and (table.rc3_ubv[1] gt -900.0))
    ploterror, table[good].bvri[0], table[good].rc3_ubv[1], table[good].bvri_err[0], table[good].rc3_ubv_err[1], $
      ps=8, xthick=2.0, ythick=2.0, charthick=2.0, charsize=1.8, xtitle='B [Dale]', $
      ytitle='B [RC3]', xrange=[6,16], yrange=[6,16], position=pos
    oploterror, table[bw].bvri[0], table[bw].rc3_ubv[1], table[bw].bvri_err[0], table[bw].rc3_ubv_err[1], $
      color=djs_icolor('red'), errcolor=djs_icolor('red'), ps=8, thick=2.0
    oplot, !x.crange, !y.crange, line=0, thick=2.0
    jj = im_stats(table[good].bvri[0]-table[good].rc3_ubv[1],/verbose,/baremin)
    jj = im_stats(table[bw].bvri[0]-table[bw].rc3_ubv[1],/verbose,/baremin,/no_head)
; --------------------------------------------------
    good = where((table.bvri[1] gt -900.0) and (table.rc3_ubv[2] gt -900.0))
    ploterror, table[good].bvri[1], table[good].rc3_ubv[2], table[good].bvri_err[1], table[good].rc3_ubv_err[2], $
      ps=8, xthick=2.0, ythick=2.0, charthick=2.0, charsize=1.8, xtitle='V [Dale]', $
      ytitle='V [RC3]', xrange=[6,16], yrange=[6,16], position=pos
    oploterror, table[vw].bvri[1], table[vw].rc3_ubv[2], table[vw].bvri_err[1], table[vw].rc3_ubv_err[2], $
      color=djs_icolor('red'), errcolor=djs_icolor('red'), ps=8
    oplot, !x.crange, !y.crange, line=0, thick=2.0
    jj = im_stats(table[good].bvri[1]-table[good].rc3_ubv[2],/verbose,/baremin,/no_head)
    jj = im_stats(table[vw].bvri[1]-table[vw].rc3_ubv[2],/verbose,/baremin,/no_head)
    dfpsclose

; ###########################################################################
; 2MASS photometry; use Danny's photometry because he went to the raw
; data and measured the fluxes for the dwarfs, which are not in the
; official 2MASS or Jarrett catalogs; Danny's fluxes have been
; corrected to zero airmass; here is some correspondence with him
; dated 2007-Feb-10:
;
;     2MASS provides calibration to 1.0 airmasses, so I corrected to 0
;     airmasses.  I can double check with Jarrett, but the equations
;     and text at the below links confirm my thinking.
;     http://spider.ipac.caltech.edu/staff/roc/2mass/calib/calib.plan.v4.2.html
;     http://www.ipac.caltech.edu/2mass/releases/second/doc/sec4_8.html
;     Also, look at the 2mass large galaxy atlas numbers for n6946 on
;     ned.  you'll see that the listed flux is 6.27 Jy, which means
;     that's the value corrected for milky way extinction.  so my
;     2mass numbers additionally correct for mw extinction and airmass
;     1->0.  I will double-check things, though, on monday.
;
; for the Vega2AB zero-point he used:
;    http://nedwww.ipac.caltech.edu/help/photoband.lst
;
;   System           Central  Freq.    Band       Flux       Reference
;                     wavel.    Hz      width   conversion*  
;                    A(unless          dnu/nu  W m^-2 Hz^-1         
;                    (noted)                  (F=k10^(-0.4m))
;    ----------      -------  --------  -----    --------     -------------------
;
;    J   (2MASS)       1.25um  2.40E+14  0.21     1.592E-23    20002MASX.2.......:
;    H      "          1.65um  1.82E+14  0.18     1.024E-23    20002MASX.2.......: 
;    K_s    "          2.17um  1.38E+14  0.15     6.668E-24    20002MASX.2.......:  
;
; compare with Blanton (the factor of 1D3 converts W/m2 to erg/s/cm2):
;
; IDL> niceprint, -2.5*alog10([1.592E-23,1.024E-23,6.668E-24]*1D3)-48.6, $
;   k_vega2ab(filterlist='twomass_'+['J','H','Ks']+'.par',/kurucz)
;     0.89514384       0.909857
;      1.3742516        1.38826
;      1.8400126        1.85430
;
; ###########################################################################

; 2MASS photometry

    twomass_Jy = [[flux.j],[flux.h],[flux.ks]]
    twomass_Jy_err = [[ferr.j],[ferr.h],[ferr.ks]]

    good = where(twomass_Jy_err gt 0.0,comp=crap)
    twomass = twomass_Jy*0.0
    twomass_err = twomass_Jy_err*0.0

    twomass[good] = -2.5*alog10(twomass_Jy[good]*1D-23)-48.6                 ; [AB]
    twomass_err[good] = 2.5*twomass_Jy_err[good]/twomass_Jy[good]/alog(10.0) ; [AB]
;   twomass[crap] = -999.0 & twomass_err[crap] = -999.0
    
    table.twomass_ref = 'dale07a'
    table.twomass = transpose(twomass)
    table.twomass_err = transpose(twomass_err)

    for iband = 0L, 2L do begin
       good = where((table.distance gt -900.0) and (table.twomass[iband] gt -900.0),ngood)
       if (ngood ne 0L) then begin
          dmod = im_d2dmod(table[good].distance,err_dmod=err_dmod,$
            err_dist=table[good].distance_err)
          table[good].twomass_absmag[iband] = table[good].twomass[iband] - dmod
          table[good].twomass_absmag_err[iband] = sqrt(table[good].twomass_err[iband]^2.0 + err_dmod^2.0)
       endif
    endfor
    
; write out

    if keyword_set(write) then begin
       splog, 'Writing '+outpath+outname+'.'
       mwrfits, table, outpath+outname, /create
       spawn, ['gzip -f '+outpath+outname], /sh
    endif

return
end

;;  splog, 'Computing absolute magnitudes and luminosities.'
;;
;;  band = ['B','V','R','I']
;;  tags = band
;;  tags_err = tags+'_err'
;;  abstags = 'M_'+band
;;  abstags_err = abstags+'_err'
;;  lumtags = band+'_lum'
;;  lumtags_err = lumtags+'_err'
;;
;;  for iband = 0L, n_elements(tags)-1L do begin
;;
;;     true = tag_exist(table,tags[iband],index=tagsindx)
;;     true = tag_exist(table,tags_err[iband],index=tagsindx_err)
;;     true = tag_exist(table,abstags[iband],index=abstagsindx)
;;     true = tag_exist(table,abstags_err[iband],index=abstagsindx_err)
;;     true = tag_exist(table,lumtags[iband],index=lumtagsindx)
;;     true = tag_exist(table,lumtags_err[iband],index=lumtagsindx_err)
;;
;;     good = where((table.distance gt -900.0) and (table.(tagsindx) gt -900.0),ngood)
;;     if (ngood ne 0L) then begin
;;
;;        mags = im_absolute_magnitudes(band[iband],table[good].(tagsindx),$
;;          table[good].distance,mag_err=table[good].(tagsindx_err),$
;;          distance_err=table[good].distance_err)
;;
;;        table[good].(abstagsindx)     = mags.absmag
;;        table[good].(abstagsindx_err) = mags.absmag_err
;;        
;;        table[good].(lumtagsindx)     = mags.lum
;;        table[good].(lumtagsindx_err) = mags.lum_err
;;        
;;     endif
;;     
;;  endfor

