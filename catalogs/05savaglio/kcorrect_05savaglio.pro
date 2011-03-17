pro kcorrect_05savaglio, phot, result, write=write
; jm06oct24nyu - derive stellar masses for the GDDS
; jm06nov19nyu - rewritten based on AGES_KCORRECT
; jm08apr09nyu - another major re-write
; jm08sep18nyu - another rewrite to use IM_KCORRECT()
    
; define the cosmology and some constants

    omega0 = 0.3 & omegal0 = 0.7 & h100 = 0.7
    
    ugriz_band_shift = 0.1 ; z=0.1
    ubvri_band_shift = 0.0 ; z=0.0

    vname = 'default.nolines' ; 'default'

; read the photometry    
    
    path = getenv('CATALOGS_DIR')+'/05savaglio/'
;   sava = rsex(path+'05savaglio.dat')

    if (n_elements(phot) eq 0L) then begin
       
       phot1 = im_read_fmr(path+'04abraham_table4.dat')
       goodphot1 = where(phot1.z gt 0.0,ngalaxy)

       phot = {galaxy: '', z: 0.0, $
         phot_b: -999.0, phot_b_err: -999.0, $
         phot_v: -999.0, phot_v_err: -999.0, $
         phot_r: -999.0, phot_r_err: -999.0, $
         phot_i: -999.0, phot_i_err: -999.0, $
         phot_z: -999.0, phot_z_err: -999.0, $
         phot_h: -999.0, phot_h_err: -999.0, $
         phot_k: -999.0, phot_k_err: -999.0}
       phot = replicate(phot,ngalaxy)

       phot.galaxy = phot1[goodphot1].id
       phot.z      = phot1[goodphot1].z

       phot.phot_b = phot1[goodphot1].bmag & phot.phot_b_err = phot1[goodphot1].ebmag
       phot.phot_v = phot1[goodphot1].vmag & phot.phot_v_err = phot1[goodphot1].evmag 
       phot.phot_r = phot1[goodphot1].rmag & phot.phot_r_err = phot1[goodphot1].ermag
       phot.phot_i = phot1[goodphot1].imag & phot.phot_i_err = phot1[goodphot1].eimag
       phot.phot_z = phot1[goodphot1].zmag & phot.phot_z_err = phot1[goodphot1].ezmag
       phot.phot_h = phot1[goodphot1].hmag & phot.phot_h_err = phot1[goodphot1].ehmag
       phot.phot_k = phot1[goodphot1].kmag & phot.phot_k_err = phot1[goodphot1].ekmag

; replace upper limits with "no detection"; two objects have zero
; R-band magnitude error: set this equal to 0.01
       
       for i = 0L, ngalaxy-1L do begin
          for j = 0L, n_tags(phot[0])-1L do begin
             if (size(phot[i].(j),/type) ne 7L) then begin
                if (phot[i].(j) eq -9.99) then begin
                   phot[i].(j) = -999.0 ; magnitude error
                   phot[i].(j-1) = -999.0 ; magnitude
                endif
             endif
          endfor
       endfor

       zero = where(phot.phot_r_err eq 0.0)
       phot[zero].phot_r_err = 0.01

    endif else ngalaxy = n_elements(phot)
       
; convert the observed photometry to maggies and store

    obsfilters = 'gdds_'+['B','V','R','I','z','H','K']+'.par'
    obsbands = 'phot_'+['b','v','r','i','z','h','k']
    obsvega2ab = k_vega2ab(filterlist=obsfilters,/kurucz)
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
       outfile_kcorr = path+'05savaglio_kcorr.fits'
       splog, 'Writing '+outfile_kcorr
       mwrfits, result, outfile_kcorr, /create
       spawn, 'gzip -f '+outfile_kcorr, /sh
    endif
       
return
end
