;+
; NAME:
;   BUILD_AGES_PHOTOMETRY
;
; PURPOSE:
;   Build the parent AGES photometric catalogs.
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2006 Sep 19, NYU - written
;   jm09may12nyu - overhaul
;   jm09aug27ucsd - another overhaul
;-

pro build_ages_photometry, phot, clobber=clobber

    common ndwfs_phot, cat1, zmerge1, kcorr1, ndwfsb1, ndwfsr1, ndwfsi1, $
      codes1, weight1, ebv1, bootes1, mips1, sdss1, sdssstars1, tmass1, $
      galex1, wsrt1, xray1

    catpath = ages_path(/catalogs)
    mycatpath = ages_path(/mycatalogs)

    vv = ages_version(/phot)
    outfile = mycatpath+'ages_photometry_'+vv+'.fits'
    if file_test(outfile+'.gz',/regular) and $
      (keyword_set(clobber) eq 0) then begin
       splog, 'FITS file '+outfile+'.gz exists; use /CLOBBER'
       return
    endif
    
; basic AGES catalogs
    if (n_elements(cat1) eq 0L) then begin
       splog, 'Reading '+catpath+'catalog.cat.noguidestars.fits.gz'
       cat1 = mrdfits(catpath+'catalog.cat.noguidestars.fits.gz',1,silent=0)
    endif
    ngal = n_elements(cat1)

; read the ZMERGE catalog to get the PASS and APER    
    if (n_elements(zmerge1) eq 0L) then begin
       splog, 'Reading '+catpath+'catalog.zmerge.fits.gz'
       zmerge1 = mrdfits(catpath+'catalog.zmerge.fits.gz',1,silent=0)
    endif
    
; store Eisenstein's corrected photometry
    if (n_elements(kcorr1) eq 0L) then begin
       splog, 'Reading '+catpath+'catalog.kcorr.v3.fits.gz'
       kcorr1 = mrdfits(catpath+'catalog.kcorr.v3.fits.gz',1,silent=0)
    endif

; store the original NDWFS BwRI photometry
    if (n_elements(ndwfsb1) eq 0L) then begin
       splog, 'Reading '+catpath+'catalog.ndwfsb.fits.gz'
       ndwfsb1 = mrdfits(catpath+'catalog.ndwfsb.fits.gz',1,silent=0)
    endif
    if (n_elements(ndwfsr1) eq 0L) then begin
       splog, 'Reading '+catpath+'catalog.ndwfsr.fits.gz'
       ndwfsr1 = mrdfits(catpath+'catalog.ndwfsr.fits.gz',1,silent=0)
    endif
    if (n_elements(ndwfsi1) eq 0L) then begin
       splog, 'Reading '+catpath+'catalog.ndwfsi.fits.gz'
       ndwfsi1 = mrdfits(catpath+'catalog.ndwfsi.fits.gz',1,silent=0)
    endif

;; store the original NDWFS photometry
;    if (n_elements(ndwfs1) eq 0L) then begin
;       splog, 'Reading '+catpath+'catalog.ndwfs.short.fits.gz'
;       ndwfs1 = mrdfits(catpath+'catalog.ndwfs.short.fits.gz',1,silent=0)
;    endif

;; usno catalog
;    if (n_elements(usno1) eq 0L) then begin
;       splog, 'Reading '+catpath+'catalog.usno.fits.gz'
;       usno1 = mrdfits(catpath+'catalog.usno.fits.gz',1,silent=0)
;    endif

; targeting codes and selection weights    
    if (n_elements(codes1) eq 0L) then begin
       splog, 'Reading '+catpath+'catalog.codes.fits.gz'
       codes1 = mrdfits(catpath+'catalog.codes.fits.gz',1,silent=0)
    endif
    if (n_elements(weight1) eq 0L) then begin
       splog, 'Reading '+catpath+'catalog.spectweight.fits.gz'
       weight1 = mrdfits(catpath+'catalog.spectweight.fits.gz',1,silent=0)
    endif

; Galactic reddening    
    if (n_elements(ebv1) eq 0L) then begin
       splog, 'Reading '+catpath+'catalog.ebv.fits.gz'
       ebv1 = mrdfits(catpath+'catalog.ebv.fits.gz',1,silent=0)
    endif

; BOOTES - see AGES_MATCH_BOOTES
    if (n_elements(bootes1) eq 0L) then begin
       suffix = '2010b' ; '2009b'
       splog, 'Reading '+mycatpath+'ages_bootes_'+suffix+'.fits.gz'
       bootes1 = mrdfits(mycatpath+'ages_bootes_'+suffix+'.fits.gz',1,silent=0)
    endif    
    
;; zBootes - see AGES_MATCH_ZBOOTES
;    if (n_elements(zband1) eq 0L) then begin
;       splog, 'Reading '+mycatpath+'ages_zbootes.fits.gz'
;       zband1 = mrdfits(mycatpath+'ages_zbootes.fits.gz',1,silent=0)
;    endif

;; UBootes - see AGES_MATCH_UBOOTES
;    if (n_elements(uband1) eq 0L) then begin
;       splog, 'Reading '+mycatpath+'ages_ubootes.fits.gz'
;       uband1 = mrdfits(mycatpath+'ages_ubootes.fits.gz',1,silent=0)
;    endif

; Spitzer/MIPS
    if (n_elements(mips1) eq 0L) then begin
       splog, 'Reading '+catpath+'catalog.mips.fits.gz'
       mips1 = mrdfits(catpath+'catalog.mips.fits.gz',1,silent=0)
    endif

; SDSS galaxies - see AGES_MATCH_SDSS
    if (n_elements(sdss1) eq 0L) then begin
       splog, 'Reading '+mycatpath+'ages.sdss.phot.dr72.fits.gz'
       sdss1 = mrdfits(mycatpath+'ages.sdss.phot.dr72.fits.gz',1,silent=0)
    endif    
    
; SDSS stars - see AGES_MATCH_SDSS
    if (n_elements(sdssstars1) eq 0L) then begin
       splog, 'Reading '+mycatpath+'ages.sdss.stars.phot.dr72.fits.gz'
       sdssstars1 = mrdfits(mycatpath+'ages.sdss.stars.phot.dr72.fits.gz',1,silent=0)
    endif    
    
; 2MASS - see AGES_MATCH_TWOMASS
    if (n_elements(tmass1) eq 0L) then begin
       splog, 'Reading '+mycatpath+'ages.twomass.phot.fits.gz'
       tmass1 = mrdfits(mycatpath+'ages.twomass.phot.fits.gz',1,silent=0)
    endif    

; GALEX
    if (n_elements(galex1) eq 0L) then begin
       splog, 'Reading '+mycatpath+'ages_galex_gr6.fits.gz'
       galex1 = mrdfits(mycatpath+'ages_galex_gr6.fits.gz',1,silent=0)
    endif

; WSRT
    if (n_elements(wsrt1) eq 0L) then begin
       splog, 'Reading '+mycatpath+'ages_wsrt_02devries.fits.gz'
       wsrt1 = mrdfits(mycatpath+'ages_wsrt_02devries.fits.gz',1,silent=0)
;      splog, 'Reading '+catpath+'catalog.wsrt.fits.gz'
;      wsrt1 = mrdfits(catpath+'catalog.wsrt.fits.gz',1,silent=0)
    endif

; X-ray
    if (n_elements(xray1) eq 0L) then begin
       splog, 'Reading '+mycatpath+'ages_xbootes.fits.gz'
       xray1 = mrdfits(mycatpath+'ages_xbootes.fits.gz',1,silent=0)
    endif

; initialize the output data structure
    phot = {$
      ages_id:             0L,$ ; zero-indexed identification number
      pass:                -1,$
      aper:                -1,$
      ra:                0.0D,$ ; from SPECINFO/ZMERGE
      dec:               0.0D,$ ; from SPECINFO/ZMERGE
      z:                  0.0,$ ; ZMERGE redshift
;     ndwfs_galaxy:        '',$ ; NDWFS galaxy name
      ebv_mw:             0.0,$ ; Milky Way E(B-V)
      infiber_i:       -999.0,$ ; I-band light fraction
      infiber_r:       -999.0}  ; I-band light fraction
    phot = replicate(phot,ngal)

    phot.ages_id  = lindgen(ngal)
    phot.ra       = cat1.ra*15.0D
    phot.dec      = cat1.dec
    phot.z        = weight1.redshift
    phot.ebv_mw   = ebv1.ebv_mw

    nospec = where(zmerge1.pass ne 0)
    phot[nospec].pass = zmerge1[nospec].pass
    phot[nospec].aper = zmerge1[nospec].aper

; AGES weights and targeting codes; redefine the structure tag names
; here so that we can recast LONG-->INT, DOUBLE-->FLOAT, etc., where
; the additional bits are not needed
    moretags = replicate({$
      newcode: 0L, $
      oldcode:  0, $
      rcode:    0, $
      field:    0, $
      bstar:    0, $
      qshort:   0, $
      gshort:   0, $
      gbright:  0, $
      grand:    0, $
      galaxy:   0, $
      quasar:   0, $
      point:    0, $
      igood:    0, $
      rgood:    0, $
      bgood:    0, $
      spec_weight:    0.0,$ ; spectroscopic completeness (~2.1%)
      target_weight:  0.0,$ ; sparse-sampling correction 
      fiber_weight:   0.0,$ ; fiber incompleteness (~4.3%)
      main_weight:      0,$
      spec_yesno:       0,$
      z_yesno:          0,$
      main_flag:        0,$ ; Eisenstein's definition of MAIN (union of MAIN_WEIGHT, SPEC_YESNO, and Z_YESNO)
      imain:            0,$ ; my definition of MAIN, except the magnitude cut
      bw_auto:     -999.0,$ ; old photometry
      r_auto:      -999.0,$
      i_auto:      -999.0,$
      bw_obs:      -999.0,$ ; corrected photometry
      r_obs:       -999.0,$
      i_obs:       -999.0,$
      i_tot:       -999.0},ngal) ; Eisenstein-like total I-band magnitude
    phot = struct_addtags(temporary(phot),moretags)
    struct_assign, codes1, phot, /nozero
    struct_assign, weight1, phot, /nozero

    phot.main_flag = (weight1.main_weight ge 1) and $
      (weight1.spec_yesno ge 1) and (weight1.z_yesno ge 1)

; add the BOOTES photometry; do not use M. Brown's FUV/NUV
; photometry  
    splog, 'Adding BOOTES UBwRIzJHKs[ch1-4]'

stop
    
    phot = struct_addtags(temporary(phot),bootes1)
;   phot = struct_addtags(temporary(phot),struct_trimtags(bootes1,$
;     except=['fuv*','nuv*','field']))
;   phot = struct_addtags(temporary(phot),struct_trimtags(bootes1,except=['field']))

; store the old observed photometry; also, until M. Brown fixes a
; problem with the *new* mag_auto values, use the old ones
    phot.bw_auto = ndwfsb1.bw_mag_auto
    phot.r_auto = ndwfsr1.r_mag_auto
    phot.i_auto = ndwfsi1.i_mag_auto

;   splog, 'Kludge!!!!!!!'
;   phot.i_mag_auto = phot.i_auto

; also store Eisenstein's corrected photometry, for comparison
    phot.bw_obs = kcorr1.bw_obs
    phot.r_obs = kcorr1.r_obs
    phot.i_obs = kcorr1.i_obs

;   kk = mrdfits('catalog.kcorr.v3.fits.gz',1)
;   djs_plot, kk.i_obs, ndwfs.itot, ps=4

; compute the total I-band magnitude using Eisenstein's kludge; WEIGHT
; is a statistic that allows us to linearly combine the two I-band
; magnitudes (I_avg and I_fainter); when I_Kron and I_R,Kron differ by
; a lot, use I_fainter, otherwise use I_avg

; the code below shows that I get the same answer as Eisenstein if I
; use the non PSF-matched photometry
;   good = where($ ; non PSF-matched aperture colors
;     (phot.i_auto gt 0.0) and (phot.i_auto lt 90.0) and $
;     (phot.r_auto gt 0.0) and (phot.r_auto lt 90.0) and $
;     (ndwfs1.r_aper_06 gt 0.0) and (ndwfs1.r_aper_06 lt 90.0) and $
;     (ndwfs1.i_aper_06 gt 0.0) and (ndwfs1.i_aper_06 lt 90.0),ngood)
;   ircolor = ndwfs1[good].i_aper_06-ndwfs1[good].r_aper_06
    
    good = where($ ; PSF-matched aperture colors
      (phot.i_auto gt 0.0) and (phot.i_auto lt 90.0) and $
      (phot.r_auto gt 0.0) and (phot.r_auto lt 90.0) and $
      (phot.r_mag_aper_06 gt 0.0) and (phot.r_mag_aper_06 lt 90.0) and $
      (phot.i_mag_aper_06 gt 0.0) and (phot.i_mag_aper_06 lt 90.0),ngood)
    ircolor = phot[good].i_mag_aper_06-phot[good].r_mag_aper_06

    ikron = phot[good].i_auto            ; I_Kron
    irband = phot[good].r_auto + ircolor ; I_R = R_Kron+(I-R)_6"
    iavg = (ikron + irband)/2.0   ; average the two magnitudes
    ifainter = ikron>irband       ; choose the fainter of the two
    weight = exp(-((ikron-irband)/0.2)^2.0)

    itot_flux = (1.0-weight)*10^(-0.4*ifainter)+weight*10^(-0.4*iavg)
    itot = -2.5*alog10(itot_flux)

    phot.i_tot = phot.i_mag_auto ; use I-auto by default
    phot[good].i_tot = itot
;   plot, phot.i_tot, phot.i_tot-phot.i_obs, ps=4, ysty=3, yr=[-0.5,0.5], xr=[18,20]

;   ww = where(phot.main_weight gt 0 and phot.spec_yesno gt 0 and phot.z_yesno gt 0 and $
;     phot.i_auto gt 0 and phot.i_auto lt 90 and phot.i_tot gt 0 and phot.i_tot lt 90,nww)
;   ss = im_stats(phot[ww].i_tot-phot[ww].i_auto,/ver,sigrej=5)
;   print, total(abs(phot[ww].i_tot-phot[ww].i_auto) gt 0.1)/float(nww)
;   print, total(abs(phot[ww].i_tot-phot[ww].i_auto) gt 0.5)/float(nww)
    
; compute the approximate R- and I-band light fractions (i.e., the
; fraction of light within the 1.5" diameter AGES fiber); use
; the *unconvolved* catalogs
;
; NOTE!  _APER = 1", _APER1 = 2" (diameter!)

;    good = where((phot.i_mag_aper_01 gt 0.0) and (phot.i_mag_aper_01 lt 90.0) and $
;      (phot.i_mag_aper_02 gt 0.0) and (phot.i_mag_aper_02 lt 90.0) and $
;      (phot.i_mag_auto gt 0.0) and (phot.i_mag_auto lt 90.0),ngood)
;    for jj = 0L, ngood-1L do begin
;       mag_1_5 = interpol([phot[good[jj]].i_mag_aper_01,$ ; 1.5" aperture magnitude
;         phot[good[jj]].i_mag_aper_02],[1.0,2.0],1.5)
;       phot[good[jj]].infiber_i = 10.0^(-0.4*(mag_1_5-phot[good[jj]].i_mag_auto))
;    endfor
;    old = phot.infiber_i
    
    good = where($
      (ndwfsi1.i_mag_aper gt 0.0) and (ndwfsi1.i_mag_aper lt 90.0) and $
      (ndwfsi1.i_mag_aper1 gt 0.0) and (ndwfsi1.i_mag_aper1 lt 90.0) and $
      (ndwfsi1.i_mag_auto gt 0.0) and (ndwfsi1.i_mag_auto lt 90.0),ngood)
    for jj = 0L, ngood-1L do begin
       mag_1_5 = interpol([ndwfsi1[good[jj]].i_mag_aper,$ ; 1.5" aperture magnitude
         ndwfsi1[good[jj]].i_mag_aper1],[1.0,2.0],1.5)
       phot[good[jj]].infiber_i = 10.0^(-0.4*(mag_1_5-ndwfsi1[good[jj]].i_mag_auto))
    endfor

    good = where($
      (ndwfsr1.r_mag_aper gt 0.0) and (ndwfsr1.r_mag_aper lt 90.0) and $
      (ndwfsr1.r_mag_aper1 gt 0.0) and (ndwfsr1.r_mag_aper1 lt 90.0) and $
      (ndwfsr1.r_mag_auto gt 0.0) and (ndwfsr1.r_mag_auto lt 90.0),ngood)
    for jj = 0L, ngood-1L do begin
       mag_1_5 = interpol([ndwfsr1[good[jj]].r_mag_aper,$ ; 1.5" aperture magnitude
         ndwfsr1[good[jj]].r_mag_aper1],[1.0,2.0],1.5)
       phot[good[jj]].infiber_r = 10.0^(-0.4*(mag_1_5-ndwfsr1[good[jj]].r_mag_auto))
    endfor

;; zBootes; see the Cool+07 paper for the quality cuts used
;    splog, 'Adding zBootes'
;    select = ['match','object_position','alpha_j2000','delta_j2000',$
;      'x_image','y_image','mag_auto','magerr_auto','nobs','photflag',$
;      'fieldname','duplicate']
;    phot = struct_addtags(temporary(phot),im_struct_trimtags(zband1,$
;      select=select,newtags='z_'+select))

;; UBootes - temporary catalog from Fuyan Bian
;    splog, 'Adding UBootes'
;    select = ['match','object_position','ra','dec',$
;      'mag','mag_err','field']
;    phot = struct_addtags(temporary(phot),im_struct_trimtags(uband1,$
;      select=select,newtags='u_'+select))
    
; MIPS
    splog, 'Adding MIPS'
    moretags = replicate({$
      phot_mips24:     -999.0, $
      phot_mips24_err: -999.0},ngal)
    phot = struct_addtags(temporary(phot),moretags)
      
    good = where((mips1.f24 gt 0.0) and (mips1.f24_err gt 0.0),ngood)
    phot[good].phot_mips24     = mips1[good].f24 ; mJy
    phot[good].phot_mips24_err = mips1[good].f24_err

; SDSS...
    splog, 'Adding SDSS'
    sdss = {$
      sdss_galaxy: 0, sdss_star: 0, sdss_ra: -999.0D, sdss_dec: -999.0D, sdss_type: -1, $
      psfflux: fltarr(5), psfflux_ivar: fltarr(5), $
      petroflux: fltarr(5), petroflux_ivar: fltarr(5), $
      modelflux: fltarr(5), modelflux_ivar: fltarr(5), $
      extinction: fltarr(5)}
    sdss = replicate(sdss,ngal)

; ...galaxies    
    isgal = where(sdss1.ra gt -900.0)
    sdss[isgal].sdss_galaxy = 1
    sdss[isgal].sdss_ra = sdss1[isgal].ra
    sdss[isgal].sdss_dec = sdss1[isgal].dec
    sdss[isgal].sdss_type = sdss1[isgal].objc_type ;=3
    sdss[isgal].petroflux = sdss1[isgal].petroflux
    sdss[isgal].petroflux_ivar = sdss1[isgal].petroflux_ivar
    sdss[isgal].modelflux = sdss1[isgal].modelflux
    sdss[isgal].modelflux_ivar = sdss1[isgal].modelflux_ivar
    sdss[isgal].extinction = sdss1[isgal].extinction

; ...stars
    isstar = where(sdssstars1.ra gt -900.0)
    sdss[isstar].sdss_star = 1
    sdss[isstar].sdss_ra = sdssstars1[isstar].ra
    sdss[isstar].sdss_dec = sdssstars1[isstar].dec
    sdss[isstar].sdss_type = sdssstars1[isstar].objc_type ;=3
    sdss[isstar].psfflux = sdssstars1[isstar].psfflux
    sdss[isstar].psfflux_ivar = sdssstars1[isstar].psfflux_ivar
    sdss[isstar].extinction = sdssstars1[isstar].extinction

    phot = struct_addtags(phot,sdss)

; 2MASS
    splog, 'Adding 2MASS'
    tmass = struct_addtags(replicate({twomass_match: 0},ngal),tmass1)
    tmass[where(tmass1.ra gt -999.0)].twomass_match = 1
    select = ['TWOMASS_MATCH','*_EXT*']
    phot = struct_addtags(phot,struct_trimtags(tmass,select=select))

; WSRT radio catalog 
    splog, 'Adding WSRT'
    moretags = replicate({$
      wsrt_match:          0, $ ; match to WSRT radio catalog
      wsrt_type:          '', $
      wsrt_flux:      -999.0, $
      wsrt_flux_err:  -999.0},ngal)
    phot = struct_addtags(temporary(phot),moretags) 
    
    phot.wsrt_match = wsrt1.match
    phot.wsrt_type = wsrt1.stype ; S=point source, M=resolved, E=complex
    phot.wsrt_flux = wsrt1.si1_4ghz       ; 1.4GHz flux [mJy]
    phot.wsrt_flux_err = wsrt1.e_si1_4ghz ; error [mJy]
; old code
;   phot.wsrt_match = wsrt1.match
;   phot.wsrt_flag = wsrt1.flag
;   phot.wsrt_flux = wsrt1.flux
;   phot.wsrt_flux_err = wsrt1.flux_err

; XBootes
    splog, 'Adding XBootes'
    select = ['match','object_position','optprob','bayprob',$
      'flux','fluxerr','hr','hrp','hrm']
    phot = struct_addtags(temporary(phot),im_struct_trimtags(xray1,$
      select=select,newtags='x_'+select))
    
; GALEX
    splog, 'Adding GALEX'
    phot = struct_addtags(temporary(phot),$
      struct_trimtags(galex1,select=['nuv_*','fuv_*']))

; define the MAIN sample, minus the magnitude cut, which should be
; project specific
    windowfile = ages_path(/window)+'ages_fields.ply'
    win = im_is_in_window(windowfile,ra=phot.ra,dec=phot.dec)
    istar = (phot.psfflux[2]*1E-9 ge 10^(-0.4*19.0)) and (phot.sdss_star eq 1)
    phot.imain = (phot.gshort gt 0) and (phot.gshort ne 2048) and $
      (phot.bgood eq 1) and (phot.rgood eq 1) and (win eq 1) and $
      (phot.field ge 1) and (phot.field le 15) and (istar eq 0)
    
; write out    
    im_mwrfits, phot, outfile, clobber=clobber

return
end    
