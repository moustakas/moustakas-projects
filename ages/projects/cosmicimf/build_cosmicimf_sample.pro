function cosmicimf_kcorrect, in_redshift, in_maggies, $
  in_ivarmaggies, filterlist=filterlist
; compute additional K-corrections for both the AGES and SDSS samples 

    ngal = n_elements(in_redshift)
    nfilt = n_elements(filterlist)
    
    kcorr = {$
      k_filterlist:             strarr(nfilt),$
      k_maggies:                fltarr(nfilt),$
      k_ivarmaggies:            fltarr(nfilt),$
      k_bestmaggies:            fltarr(nfilt),$
      k_mass:                          -999.0,$
      k_coeffs:                     fltarr(5),$
      k_chi2:                          -999.0,$
      k_cflux:                fltarr(5)-999.0,$
      k_uvflux:               fltarr(2)-999.0,$
      k_galex_absmag:         fltarr(2)-999.0,$
      k_galex_absmag_ivar:    fltarr(2)-999.0,$
      k_galex_kcorrect:       fltarr(2)-999.0,$
      k_ugriz_absmag:         fltarr(5)-999.0,$
      k_ugriz_absmag_ivar:    fltarr(5)-999.0,$
      k_ugriz_kcorrect:       fltarr(5)-999.0,$
      k_ubvrijhk_absmag:      fltarr(8)-999.0,$
      k_ubvrijhk_absmag_ivar: fltarr(8)-999.0,$
      k_ubvrijhk_kcorrect:    fltarr(8)-999.0}
    kcorr = replicate(kcorr,ngal)

    ugriz_kcorrect = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,sdss_filterlist(),band_shift=0.1,chi2=chi2,mass=mass,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=ugriz_absmag,$
      ivarabsmag=ugriz_absmag_ivar,clineflux=cflux,uvflux=uvflux,$
      /silent,vname=vname)

    galex_kcorrect = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,galex_filterlist(),band_shift=0.0,absmag=galex_absmag,$
      ivarabsmag=galex_absmag_ivar,/vega,/silent,vname=vname) ; Vega

    kcorrect = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,[bessell_filterlist(),twomass_filterlist()],$
      band_shift=0.0,absmag=absmag,ivarabsmag=absmag_ivar,$r=[
      /vega,/silent,vname=vname) ; Vega           
    
    kcorr.k_filterlist = filterlist
    kcorr.k_maggies = in_maggies
    kcorr.k_ivarmaggies = in_ivarmaggies
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

pro build_cosmicimf_sample, sample
; jm10mar16ucsd - build the AGES/COSMICIMF parent sample

; select the parent sample    
    agespath = ages_path(/analysis)
    cosmicimfpath = ages_path(/projects)+'cosmicimf/'

; read the parent AGES photometry
    splog, '#########################'
    splog, 'Building the AGES sample'

    phot_version = ages_version(/photometry)
    allphot = mrdfits(agespath+'ages_photometry_'+$
      phot_version+'.fits.gz',1)

; basic parameters of the selection    
    sample_zmin = 0.05
    sample_zmax = 0.75

    vname = 'default.nolines'
    ifilt = 'ndwfs_I.par'
    ibright = 15.0              ; [Vega]
    ifaint = 19.95              ; [Vega]

; define the parent sample
    splog, 'Converting to maggies'
    ages_to_maggies, allphot, maggies, ivarmaggies, $
      filterlist=filterlist, nozpoffset=nozpoffset
    nfilt = n_elements(filterlist)
    bwindx = where(strmatch(filterlist,'*_Bw*'))
    rindx = where(strmatch(filterlist,'*_R*'))
    iindx = where(strmatch(filterlist,'*_I*'))

    splog, 'Applying the window function'
    window = ages_isin_window(allphot.ra,allphot.dec)            ; within the 2004 survey region
    istar = (allphot.psfflux[2]*1E-9 gt 10^(-0.4*19.0)) and $    ; SDSS type=star and r<19
      (allphot.sdss_star eq 1)

    parent = where((allphot.z ge sample_zmin) and (allphot.z le sample_zmax) and $ ; redshift cuts
      (allphot.gshort gt 0) and (allphot.gshort ne 2048) and $                     ; main galaxy
      (allphot.field ge 1) and (allphot.field le 15) and $                         ; field assignment cut
      (allphot.i_tot le ifaint) and (allphot.i_tot ge ibright) and $ ; 15<I<19.95
      (istar eq 0) and (window eq 1) and (ivarmaggies[bwindx,*] gt 0.0) and $
      (ivarmaggies[rindx,*] gt 0.0) and (ivarmaggies[iindx,*] gt 0.0),ngal)
;     (allphot.phot_mips24 gt 0.27),ngal) ; require f(24)>0.27 mJy
    splog, 'Ngal = ', ngal

; trim to just the tags we need...
    sample = allphot[parent]
    sample = struct_trimtags(sample,select=['ages_id','ra','dec',$
      'z','field','*_weight','i_tot','i_mag_auto','*_mag_aper_04',$
      '*mips*','x_*','galex_detect'])

; compute K-corrections; don't use IRAC/ch2-4
    splog, 'Computing additional K-corrections'
    in_maggies = maggies[*,parent]
    in_ivarmaggies = ivarmaggies[*,parent]
    toss = where(strmatch(filterlist,'*_ch[2-4].par*',/fold))
    in_ivarmaggies[toss,*] = 0.0

    kcorr = cosmicimf_kcorrect(sample.z,in_maggies,$
      in_ivarmaggies,filterlist=filterlist)
    sample = struct_addtags(temporary(sample),kcorr)

; check the largest outliers    
;   check = where(kcorr.k_chi2 gt im_quantile(sample.k_chi2,quant=0.99) and $
;     sample.phot_mips24 gt 0.27,ncheck)
;   tt = tag_names(kcorr)
;   junk = im_struct_trimtags(kcorr[check],select=tt,newtags=repstr(tt,'K_',''))
;   junk = struct_addtags(junk,replicate({z: 0.0},ncheck))
;   junk.z = sample[check].z
;   kcorrect_qaplot, junk, psfile='test.ps', in_filterlist=filterlist, $
;     vname=vname, /clobber

; compute additional ancillary quantities that we are going to need:
; ZMIN and ZMAX with and without evolution, and the final galaxy
; weight as the product of all the various selection terms
    splog, 'Computing zmin/zmax'
    moretags = replicate({zmin_noevol: 0.0, zmax_noevol: 0.0, $
      zmin_evol: 0.0, zmax_evol: 0.0, final_weight: 0.0},ngal)
    sample = struct_addtags(temporary(sample),moretags)
    sample.final_weight = sample.spec_weight*sample.target_weight*$
      sample.fiber_weight*sample.catalog_weight

    q0 = 1.6 & q1 = 0.0 & qz0 = 0.1 ; following Eisenstein+09
    zminzmax_evol = im_zminzmax(sample.z,sample.i_tot,$
      sample.k_coeffs,bright=ibright,faint=ifaint,$
      filter=ifilt,vname=vname,q0=q0,q1=q1,qz0=qz0)

    q0 = 0.0 & q1 = 0.0 & qz0 = 0.0 ; no evolution
    zminzmax_noevol = im_zminzmax(sample.z,sample.i_tot,$
      sample.k_coeffs,bright=ibright,faint=ifaint,$
      filter=ifilt,vname=vname,q0=q0,q1=q1,qz0=qz0)

    sample.zmin_noevol = zminzmax_noevol.zmin
    sample.zmax_noevol = zminzmax_noevol.zmax
    sample.zmin_evol = zminzmax_evol.zmin
    sample.zmax_evol = zminzmax_evol.zmax

; write out       
    im_mwrfits, sample, cosmicimfpath+'ages_cosmicimf.fits', /clobber

return
end
