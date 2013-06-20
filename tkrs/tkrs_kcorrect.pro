pro tkrs_kcorrect, zcat, result, write=write
; jm08apr23nyu - stolen from DEEP2_KCORRECT_ALL
; jm08sep20nyu - now a wrapper on IM_KCORRECT()
    
; define the cosmology and some constants

    omega0 = 0.3 & omegal0 = 0.7 & h100 = 0.7
    
    ugriz_band_shift = 0.1 ; z=0.1
    ubvri_band_shift = 0.0 ; z=0.0

    vname = 'default.nolines' ; 'default'

; read the TKRS redshift catalog

    if (n_elements(zcat) eq 0L) then zcat = read_tkrs(/zcat)
    ngalaxy = n_elements(zcat)

; convert the observed photometry to maggies and store

    splog, 'Converting observed photometry to maggies.'
    obsfilters = 'goods_'+['acs_'+['f435w','f606w','f775w','f850lp'],$
      ['J','H','Ks']+'_isaac_etc']+'.par'
    nobsfilter = n_elements(obsfilters)
    tkrs_to_maggies, zcat, obsmaggies, obsmaggies_ivar

; compute k-corrections    
    
    splog, 'Computing ugriz k-corrections'
    ugriz_kcorrect = im_kcorrect(zcat.z,obsmaggies,obsmaggies_ivar,$
      obsfilters,sdss_filterlist(),band_shift=ugriz_band_shift,$
      chi2=chi2,coeffs=coeffs,rmaggies=rmaggies,vname=vname,mass=mass,$
      absmag=ugriz_absmag,ivarabsmag=ugriz_absmag_ivar,$
      clineflux=cflux,omega0=omega0,omegal0=omegal0,/silent)

    splog, 'Computing UBVRI k-corrections'
    ubvri_kcorrect = im_kcorrect(zcat.z,obsmaggies,obsmaggies_ivar,$
      obsfilters,bessell_filterlist(),band_shift=ubvri_band_shift,$ 
      vname=vname,absmag=ubvri_absmag,ivarabsmag=ubvri_absmag_ivar,$
      omega0=omega0,omegal0=omegal0,/silent)

    result_template = {$
      zobj:                       -999.0, $
      abmaggies:      fltarr(nobsfilter), $
      abmaggies_ivar: fltarr(nobsfilter), $
      mass:                       -999.0, $
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
    result = struct_addtags(zcat,temporary(result))

    result.zobj              = zcat.z
    result.abmaggies         = obsmaggies
    result.abmaggies_ivar    = obsmaggies_ivar
    result.mass              = alog10(mass)   ; Chabrier IMF
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
       outpath = tkrs_path()
       outfile_kcorr = outpath+'tkrs_kcorr.fits'
       splog, 'Writing '+outfile_kcorr
       mwrfits, result, outfile_kcorr, /create
       spawn, 'gzip -f '+outfile_kcorr, /sh
    endif
       
return
end
