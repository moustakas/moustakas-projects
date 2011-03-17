;
; NAME:
;   SG1120_KCORRECT
;
; PURPOSE:
;   Compute K-corrections and rest-frame luminosities (and stellar
;   masses) for SG1120 using K-correct.
;
; INPUTS
;   use_zpt - see SG1120_TO_MAGGIES
;
; KEYWORD PARAMETERS:
;   sdss - use the SDSS photometry 
;   nowrite - do not write out the results
;   tweak - see SG1120_TWEAK_ZEROPOINTS 
;
; OUTPUTS:
;   allsg1120 - photometry data structure
;   allresult - stellar masses and K-corrections
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Jun 15, NYU - rewritten
;   jm09jul07nyu - added the V-band test stuff; note that /SDSS should
;     only be used for testing, as the photometry is *very* noisy 
;   jm09aug05ucsd - rederive K-corrections using the revised
;     zeropoints (see VIMOS_CALIBRATE and LDSS3_CALIBRATE); also now
;     include the FLAMINGOS Ks-band photometry
;-

pro sg1120_kcorrect, sg1120, result, use_zpt=use_zpt, sdss=sdss, $
  nowrite=nowrite, tweak=tweak

    if keyword_set(sdss) then begin
       splog, 'It is not recommended that you use /SDSS!'
       return
    endif
    
; define the cosmology and some constants
    version = sg1120_version(/kcorrect)
    analysis_path = sg1120_path(/analysis)

    omega0 = 0.3 & omegal0 = 0.7 & h100 = 0.7
    
    ugriz_band_shift = 0.1 ; z=0.1
    ubvri_band_shift = 0.0 ; z=0.0
    vname = 'default.nolines' ; 'default'

; read the merged photometric catalog; only compute k-corrections for
; objects with well-measured redshifts; see
; BUILD_SG1120_PARENT_CATALOG for the meaning of SUFFIX 
    suffix = 'refR2.0chi2'
    pversion = sg1120_version(/parent)
    photfile = analysis_path+'sg1120_parent_'+suffix+'_'+pversion+'.fits.gz'
    splog, 'Reading '+photfile
    allsg1120 = mrdfits(photfile,1,/silent)
    nallgalaxy = n_elements(allsg1120)

; define a subsample of objects for calibrating the V-band, otherwise
; fit everything with a redshift and photometry in at least three
; bandpasses     
    if keyword_set(tweak) then begin
       keepindx = where((allsg1120.z gt 0.3) and (allsg1120.z lt 0.4) and $
         (allsg1120.phot_b gt 0.0) and (allsg1120.phot_r gt 0.0) and $
         (allsg1120.phot_gprime gt 0.0) and (allsg1120.phot_rprime gt 0.0) and $
         (allsg1120.phot_ks gt 0.0))
       outsuffix = suffix+'_'+version+'.tweak'
    endif else begin
       keepindx = where((allsg1120.phot_b gt 0.0) and (allsg1120.phot_v gt 0.0) and $
         (allsg1120.phot_r gt 0.0) and (allsg1120.z gt 0.0),ngal)
       outsuffix = suffix+'_'+version
    endelse
    sg1120 = allsg1120[keepindx]

; convert the observed photometry to maggies
    sg1120_to_maggies, sg1120, obsmaggies, obsmaggies_ivar, $
      sdss=sdss, filterlist=obsfilters, use_zpt=use_zpt
    nobsfilter = n_elements(obsfilters)
    ubvri_filterlist = ['bessell_U.par','bessell_B.par',$
      'bessell_V.par','bessell_R.par','bessell_I.par']
    ugriz_filterlist = ['sdss_u0.par','sdss_g0.par',$
      'sdss_r0.par','sdss_i0.par','sdss_z0.par']

;   obsmaggies_ivar[5:8,*] = 0.0 ; SDSS

; compute k-corrections    
;   splog, 'Computing ugriz k-corrections'
;   ugriz_kcorrect = im_kcorrect(sg1120.z,obsmaggies,obsmaggies_ivar,$
;     obsfilters,band_shift=ugriz_band_shift,chi2=chi2,/sdss,$ ; note /SDSS
;     coeffs=coeffs,rmaggies=rmaggies,vname=vname,mass=mass,$
;     absmag=ugriz_absmag,ivarabsmag=ugriz_absmag_ivar,intsfh=intsfh,$
;     /silent,out_filterlist=ugriz_filterlist)

    splog, 'Computing UBVRI k-corrections'
    ubvri_kcorrect = im_kcorrect(sg1120.z,obsmaggies,obsmaggies_ivar,$
      obsfilters,band_shift=ubvri_band_shift,out_filterlist=ubvri_filterlist,$
      synth_inmaggies=bestmaggies,vname=vname,absmag=ubvri_absmag,$
      ivarabsmag=ubvri_absmag_ivar,clineflux=cflux,mass=mass,coeffs=coeffs,$
      chi2=chi2,/reset,/silent,/vega)

; pack into a structure and write out
    result_template = {$
      id:                             0L, $
      z:                          -999.0, $
      maggies:        fltarr(nobsfilter)-999.0, $
      ivarmaggies:    fltarr(nobsfilter)-999.0, $
      bestmaggies:    fltarr(nobsfilter)-999.0, $
      mass:                       -999.0, $
;     intsfh:                     -999.0, $
      coeffs:                  fltarr(5)-999.0, $
      chi2:                       -999.0, $
;     cflux_3727:                 -999.0, $
;     cflux_4861:                 -999.0, $
;     cflux_4959:                 -999.0, $
;     cflux_5007:                 -999.0, $
;     cflux_6563:                 -999.0, $
;     ugriz_absmag:      fltarr(5)-999.0, $
;     ugriz_absmag_ivar: fltarr(5)-999.0, $
;     ugriz_kcorrect:    fltarr(5)-999.0, $
      ubvri_absmag:      fltarr(5)-999.0, $
      ubvri_absmag_ivar: fltarr(5)-999.0, $
      ubvri_kcorrect:    fltarr(5)-999.0}
    allresult = replicate(result_template,nallgalaxy)

    allresult.id = allsg1120.id
    result = allresult[keepindx]

    result.z                 = sg1120.z
    result.maggies           = obsmaggies
    result.ivarmaggies       = obsmaggies_ivar
    result.bestmaggies       = bestmaggies
    result.mass              = alog10(mass)   ; Chabrier IMF, h=0.7
;   result.intsfh            = alog10(intsfh) ; Chabrier IMF, h=0.7
    result.coeffs            = coeffs
    result.chi2              = chi2
;   result.ugriz_absmag      = ugriz_absmag
;   result.ugriz_absmag_ivar = ugriz_absmag_ivar
;   result.ugriz_kcorrect    = ugriz_kcorrect
    result.ubvri_absmag      = ubvri_absmag
    result.ubvri_absmag_ivar = ubvri_absmag_ivar
    result.ubvri_kcorrect    = ubvri_kcorrect

;   result.cflux_3727 = reform(cflux[0,*])
;   result.cflux_4861 = reform(cflux[1,*])
;   result.cflux_4959 = reform(cflux[2,*])
;   result.cflux_5007 = reform(cflux[3,*])
;   result.cflux_6563 = reform(cflux[4,*])
    
; just keep the objects used to calibrate the V-band, otherwise, write
; out a line-matched k-correction catalog; also make a QAplot
    
    if (keyword_set(nowrite) eq 0) then begin
       kcorrfile = analysis_path+'sg1120_kcorr_'+outsuffix+'.fits'
       psfile = repstr(kcorrfile,'.fits','.ps')
       im_mwrfits, result, kcorrfile
       if (keyword_set(tweak) eq 0) then $
         wsex, result, outfile=repstr(kcorrfile,'.fits','.dat')
;      im_mwrfits, allresult, kcorrfile
;      wsex, allresult, outfile=repstr(kcorrfile,'.fits','.dat')
       kcorrect_qaplot, result, obsfilters, psfile=psfile, /clobber
    endif

; test code when running /CALIBRATE_VBAND
;  analysis_path = sg1120_path(/analysis)
;  suffix = 'refR2.0chi2.Vcalib'
;  kcorrfile = analysis_path+'sg1120_kcorr_'+suffix+'.fits.gz'
;  rr = mrdfits(kcorrfile,1)
;
;  diff = -2.5*alog10(rr.maggies[1]/rr.bestmaggies[1])
;  plot, rr.z, diff, ps=4
;  ss = im_stats(diff,/ver,sigrej=3.0)
;
;  scale[jj] = total(rr.ivarmaggies[*,jj]*maggies[*,jj]*rmaggies[*,jj])/$
;    total(ivarmaggies[*,jj]*rmaggies[*,jj]^2.0)
   
return
end
