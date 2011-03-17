function calibsfr_kcorrect, in_redshift, in_maggies, $
  in_ivarmaggies, filterlist=filterlist
; compute additional K-corrections for both the AGES and SDSS samples 

    ngal = n_elements(in_redshift)
    nfilt = n_elements(filterlist)
    
    kcorr = {$
      k_bestmaggies:            fltarr(nfilt), $
      k_mass:                          -999.0, $
      k_coeffs:                     fltarr(5), $
      k_chi2:                          -999.0, $
      k_cflux:                fltarr(5)-999.0, $
      k_uvflux:               fltarr(2)-999.0, $
      k_galex_absmag:         fltarr(2)-999.0, $
      k_galex_absmag_ivar:    fltarr(2)-999.0, $
      k_galex_kcorrect:       fltarr(2)-999.0, $
      k_ugriz_absmag:         fltarr(5)-999.0, $
      k_ugriz_absmag_ivar:    fltarr(5)-999.0, $
      k_ugriz_kcorrect:       fltarr(5)-999.0, $
      k_ubvrijhk_absmag:      fltarr(8)-999.0, $
      k_ubvrijhk_absmag_ivar: fltarr(8)-999.0, $
      k_ubvrijhk_kcorrect:    fltarr(8)-999.0}
    kcorr = replicate(kcorr,ngal)

    ugriz_kcorrect = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,sdss_filterlist(),band_shift=0.1,chi2=chi2,mass=mass,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=ugriz_absmag,$
      ivarabsmag=ugriz_absmag_ivar,clineflux=cflux,uvflux=uvflux,$
      /silent,vname=vname)

    galex_kcorrect = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,galex_filterlist(),band_shift=0.0,absmag=galex_absmag,$
      ivarabsmag=galex_absmag_ivar,/silent,vname=vname) ; AB

    kcorrect = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,[bessell_filterlist(),twomass_filterlist()],$
      band_shift=0.0,absmag=absmag,ivarabsmag=absmag_ivar,$
      /vega,/silent,vname=vname) ; Vega           
    
    kcorr.k_bestmaggies = bestmaggies
    kcorr.k_mass = alog10(mass) ; h=0.7, Chabrier
    kcorr.k_coeffs = coeffs
    kcorr.k_chi2 = chi2
    kcorr.k_cflux = cflux
    kcorr.k_uvflux = uvflux
    kcorr.k_ugriz_absmag = ugriz_absmag
    kcorr.k_ugriz_absmag_ivar = ugriz_absmag_ivar
    kcorr.k_ugriz_kcorrect = ugriz_kcorrect
    kcorr.k_galex_absmag = galex_absmag
    kcorr.k_galex_absmag_ivar = galex_absmag_ivar
    kcorr.k_galex_kcorrect = galex_kcorrect
    kcorr.k_ubvrijhk_absmag = absmag
    kcorr.k_ubvrijhk_absmag_ivar = absmag_ivar
    kcorr.k_ubvrijhk_kcorrect = kcorrect

return, kcorr
end

pro build_calibsfr_sample, sample
; jm10feb21ucsd - build the AGES/CALIBSFR parent samples

    common calibsfr_sample, allphot, lowz, galex, twomass

; select the parent sample    
    agespath = ages_path(/analysis)
    calibsfrpath = ages_path(/projects)+'calibsfr/'
    h100 = 0.7

    splog, '#########################'
    splog, 'Building the AGES sample'

    if (n_elements(allphot) eq 0) then begin
       phot_version = ages_version(/photometry)
       allphot = mrdfits(agespath+'ages_photometry_'+$
         phot_version+'.fits.gz',1)
    endif

    allispec = read_ages_gandalf(/ppxf)
    phot = allphot[allispec.ages_id]
    
; basic parameters of the selection    
    sample_zmin = 0.05
    sample_zmax = 0.35

    vname = 'default.nolines'
    ifilt = 'ndwfs_I.par'
    ibright = 15.0              ; [Vega]
    ifaint = 19.95              ; [Vega]

    ewcut1 = 1.0
    snrcut1 = 1.0
    
; define the parent sample
    splog, 'Converting to maggies'
    ages_to_maggies, phot, maggies, ivarmaggies, $
      filterlist=filterlist, nozpoffset=nozpoffset
    nfilt = n_elements(filterlist)
    bwindx = where(strmatch(filterlist,'*_Bw*'))
    rindx = where(strmatch(filterlist,'*_R*'))
    iindx = where(strmatch(filterlist,'*_I*'))

    splog, 'Applying the window function'
    istar = (phot.psfflux[2]*1E-9 gt 10^(-0.4*19.0)) and $ ; SDSS type=star and r<19
      (phot.sdss_star eq 1)

    parent = where((phot.z ge sample_zmin) and (phot.z le sample_zmax) and $             ; redshift cuts
      (phot.gshort gt 0) and (phot.i_tot le ifaint) and (phot.i_tot ge ibright) and $    ; galaxy, 15<I<19.95
      (istar eq 0) and (phot.phot_mips24 gt 0.0) and (ivarmaggies[bwindx,*] gt 0.0) and $
      (ivarmaggies[rindx,*] gt 0.0) and (ivarmaggies[iindx,*] gt 0.0) and $
      (allispec.h_alpha_ew[0] gt ewcut1) and (allispec.h_beta_ew[0] gt ewcut1) and $
      (allispec.h_alpha[0]/allispec.h_alpha[1] gt snrcut1) and $
      (allispec.h_beta[0]/allispec.h_beta[1] gt snrcut1),ngal)
    splog, 'Ngal = ', ngal

    ispec = allispec[parent]
    sample = struct_addtags(phot[parent],struct_trimtags(ispec,$
      except=['ages_id','pass','aper','z']))

; compute K-corrections
    splog, 'Computing K-corrections'
    kcorr = calibsfr_kcorrect(sample.z,maggies[*,parent],$
      ivarmaggies[*,parent],filterlist=filterlist)
    sample = struct_addtags(temporary(sample),kcorr)

; convert the measured emission-line EWs to observed and
; reddening-corrected luminosities
    dust = iunred_linedust(sample,snrcut=0.0)
    sample = struct_addtags(temporary(sample),$
      struct_trimtags(dust,select=['*hahb*','*hbha*']))
    
    lum = {$
      lha:         0.0,$
      lha_err:     0.0,$
      lha_cor:     0.0,$
      lha_cor_err: 0.0}
    lum = replicate(lum,ngal)

    dlum = dluminosity(sample.z,/cm)
    factor = sample.k_cflux[4]*4.0*!dpi*dlum^2
    lum.lha = alog10(factor*sample.h_alpha_ew[0]) ; [log, erg/s]
    lum.lha_err = sample.h_alpha_ew[1]/sample.h_alpha_ew[0]/alog(10.0)

    dustfactor = 10.0D^(0.4*sample.ebv_hahb*k_lambda(6562.8,/odon))
    lum.lha_cor = lum.lha + alog10(dustfactor)
    lum.lha_cor_err = lum.lha_err + alog10(dustfactor)
    
    sample = struct_addtags(temporary(sample),lum)

; compute L(IR) and L(24) using the Chary & Elbaz (2001) and Dale &
; Helou (2002) models
    ir = {$
      lir_ce01:     0.0,$ ; Chary & Elbaz (2001) L(IR)
      lir_ce01_err: 0.0,$
      lir_dh02:     0.0,$ ; Dale & Helou (2002) L(IR)
      lir_dh02_err: 0.0,$
      l24_ce01:     0.0,$
      l24_ce01_err: 0.0,$
      l24_dh02:     0.0,$
      l24_dh02_err: 0.0}
    ir = replicate(ir,ngal)
    
    lir_ce01 = im_f24tolir(sample.z,sample.phot_mips24,$
      f24_err=sample.phot_mips24_err,/chary,err_lir=lir_ce01_err,$
      l24=l24_ce01,err_l24=l24_ce01_err)
    ir.lir_ce01 = alog10(lir_ce01*3.826D33) ; [log, erg/s]
    ir.l24_ce01 = alog10(l24_ce01*3.826D33)
    ir.lir_ce01_err = lir_ce01_err/lir_ce01/alog(10.0)
    ir.l24_ce01_err = l24_ce01_err/l24_ce01/alog(10.0)
    
    lir_dh02 = im_f24tolir(sample.z,sample.phot_mips24,$
      f24_err=sample.phot_mips24_err,/dale,err_lir=lir_dh02_err,$
      l24=l24_dh02,err_l24=l24_dh02_err)
    ir.lir_dh02 = alog10(lir_dh02*3.826D33) ; [log, erg/s]
    ir.l24_dh02 = alog10(l24_dh02*3.826D33)
    ir.lir_dh02_err = lir_dh02_err/lir_dh02/alog(10.0)
    ir.l24_dh02_err = l24_dh02_err/l24_dh02/alog(10.0)

    sample = struct_addtags(temporary(sample),ir)
    
; write out       
    im_mwrfits, sample, calibsfrpath+'ages_calibsfr.fits', /clobber
       
return
end
