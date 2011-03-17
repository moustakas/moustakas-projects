pro kcorrect_03lilly, phot, result, write=write
; jm06nov20nyu - based on MASS_05SAVAGLIO
; jm08apr23nyu - totally re-written
; jm08sep20nyu - another rewrite to use IM_KCORRECT()
    
; define the cosmology and some constants

    omega0 = 0.3 & omegal0 = 0.7 & h100 = 0.7
    
    ugriz_band_shift = 0.1 ; z=0.1
    ubvri_band_shift = 0.0 ; z=0.0

    vname = 'default.nolines' ; 'default'

; read the photometry    
    
    path = getenv('CATALOGS_DIR')+'/03lilly/'

    if (n_elements(phot) eq 0L) then begin

       phot1 = rsex(path+'cfrs_phot_from_savaglio.dat')
       ngalaxy = n_elements(phot1)

       phot = {galaxy: '', z: 0.0, $
         phot_b: -999.0, phot_b_err: -999.0, $
         phot_v: -999.0, phot_v_err: -999.0, $
         phot_i: -999.0, phot_i_err: -999.0, $
         phot_k: -999.0, phot_k_err: -999.0}
       phot = replicate(phot,ngalaxy)

       phot.galaxy = phot1.galaxy
       phot.z      = phot1.z

; the cut on one magnitude of error excludes the K-band magnitude for
; CFRS_03.1375 (z=0.635)       

       apcor = 0                ; use the aperture photometry
       
       if keyword_set(apcor) then begin
          
          good = where((phot1.i_iso lt 99.99) and (phot1.b3 lt 99.99) and (phot1.b3_err lt 1.00),ngood)
          if (ngood ne 0L) then begin
             phot[good].phot_b = phot1[good].i_iso + (phot1[good].b3-phot1[good].i3)
             phot[good].phot_b_err = phot1[good].b3_err
          endif
          
          good = where((phot1.i_iso lt 99.99) and (phot1.v3 lt 99.99) and (phot1.v3_err lt 1.00),ngood)
          if (ngood ne 0L) then begin
             phot[good].phot_v = phot1[good].i_iso + (phot1[good].v3-phot1[good].i3)
             phot[good].phot_v_err = phot1[good].v3_err
          endif
          
          good = where((phot1.i_iso lt 99.99) and (phot1.i3 lt 99.99) and (phot1.i3_err lt 1.00),ngood)
          if (ngood ne 0L) then begin
             phot[good].phot_i = phot1[good].i_iso + (phot1[good].i3-phot1[good].i3)
             phot[good].phot_i_err = phot1[good].i3_err
          endif
          
          good = where((phot1.i_iso lt 99.99) and (phot1.k3 lt 99.99) and (phot1.k3_err lt 1.00),ngood)
          if (ngood ne 0L) then begin
             phot[good].phot_k = phot1[good].i_iso + (phot1[good].k3-phot1[good].i3)
             phot[good].phot_k_err = phot1[good].k3_err
          endif

       endif else begin

          good = where((phot1.b3 lt 99.99) and (phot1.b3_err lt 1.00),ngood)
          if (ngood ne 0L) then begin
             phot[good].phot_b = phot1[good].b3
             phot[good].phot_b_err = phot1[good].b3_err
          endif
          
          good = where((phot1.v3 lt 99.99) and (phot1.v3_err lt 1.00),ngood)
          if (ngood ne 0L) then begin
             phot[good].phot_v = phot1[good].v3
             phot[good].phot_v_err = phot1[good].v3_err
          endif
          
          good = where((phot1.i3 lt 99.99) and (phot1.i3_err lt 1.00),ngood)
          if (ngood ne 0L) then begin
             phot[good].phot_i = phot1[good].i3
             phot[good].phot_i_err = phot1[good].i3_err
          endif
          
          good = where((phot1.k3 lt 99.99) and (phot1.k3_err lt 1.00),ngood)
          if (ngood ne 0L) then begin
             phot[good].phot_k = phot1[good].k3
             phot[good].phot_k_err = phot1[good].k3_err
          endif

       endelse

    endif else ngalaxy = n_elements(phot)
       
; convert the observed photometry to maggies and store

    obsfilters = ['bessell_'+['B','V','I'],'twomass_Ks']+'.par'
    obsbands = 'phot_'+['b','v','i','k']
    obsvega2ab = k_vega2ab(filterlist=obsfilters,/kurucz)*0.0 ; already in AB
    nobsfilter = n_elements(obsfilters)

    splog, 'Converting observed photometry to maggies.'
    obsmaggies = convert_to_maggies(phot,obsbands,obsvega2ab,$
      maggies_invvar=obsmaggies_ivar,mags=mobs_ab,err_mags=mobs_ab_err,$
      ivar_mags=mobs_ab_ivar)
    setzero = where(obsmaggies_ivar eq 1D-32,nsetzero)
    if (nsetzero ne 0L) then obsmaggies_ivar[setzero] = 0.0

; compute k-corrections    
    
    splog, 'Computing ugriz k-corrections'
    ugriz_kcorrect = im_kcorrect(phot.z,obsmaggies,obsmaggies_ivar,$
      obsfilters,band_shift=ugriz_band_shift,chi2=chi2,/sdss,$ ; note /SDSS
      coeffs=coeffs,rmaggies=rmaggies,vname=vname,mass=mass,$
      absmag=ugriz_absmag,ivarabsmag=ugriz_absmag_ivar,intsfh=intsfh,$
      clineflux=cflux,omega0=omega0,omegal0=omegal0,/silent)

    splog, 'Computing UBVRI k-corrections'
    ubvri_kcorrect = im_kcorrect(phot.z,obsmaggies,obsmaggies_ivar,$
      obsfilters,band_shift=ubvri_band_shift,/vega,$ ; note /VEGA
      vname=vname,absmag=ubvri_absmag,ivarabsmag=ubvri_absmag_ivar,$
      omega0=omega0,omegal0=omegal0,/silent)

    result_template = {$
;     z:                          -999.0, $
      abmaggies:      fltarr(nobsfilter), $
      abmaggies_ivar: fltarr(nobsfilter), $
      mass:                       -999.0, $
      intsfh:                     -999.0, $
      coeffs:                  fltarr(5), $
      chi2:                       -999.0, $
      cflux_3727:                 -999.0, $
      cflux_4861:                 -999.0, $
      cflux_4959:                 -999.0, $
      cflux_5007:                 -999.0, $
      cflux_6563:                 -999.0, $
      ugriz_absmag:      fltarr(5)-999.0, $
      ugriz_absmag_ivar: fltarr(5)-999.0, $
      ugriz_kcorrect:    fltarr(5)-999.0, $
      ubvri_absmag:      fltarr(5)-999.0, $
      ubvri_absmag_ivar: fltarr(5)-999.0, $
      ubvri_kcorrect:    fltarr(5)-999.0}
    result = replicate(result_template,ngalaxy)
    result = struct_addtags(phot,temporary(result))

;   result.z                 = phot.z
    result.abmaggies         = obsmaggies
    result.abmaggies_ivar    = obsmaggies_ivar
    result.mass              = alog10(mass)   ; Chabrier IMF
    result.intsfh            = alog10(intsfh) ; Chabrier IMF
    result.coeffs            = coeffs
    result.chi2              = chi2
    result.ugriz_absmag      = ugriz_absmag
    result.ugriz_absmag_ivar = ugriz_absmag_ivar
    result.ugriz_kcorrect    = ugriz_kcorrect
    result.ubvri_absmag      = ubvri_absmag
    result.ubvri_absmag_ivar = ubvri_absmag_ivar
    result.ubvri_kcorrect    = ubvri_kcorrect

    result.cflux_3727 = reform(cflux[0,*])
    result.cflux_4861 = reform(cflux[1,*])
    result.cflux_4959 = reform(cflux[2,*])
    result.cflux_5007 = reform(cflux[3,*])
    result.cflux_6563 = reform(cflux[4,*])

    if keyword_set(write) then begin
       outfile_kcorr = path+'03lilly_kcorr.fits'
       splog, 'Writing '+outfile_kcorr
       mwrfits, result, outfile_kcorr, /create
       spawn, 'gzip -f '+outfile_kcorr, /sh
    endif
       
return
end
