pro bcgsfhs_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber
; jm13dec29siena - do SED-fitting

    prefix = 'bcgsfhs'
    
    sersicpath = bcgsfhs_path(/sersic)
    isedfit_dir = bcgsfhs_path(/isedfit)
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; read the sample
    sample = read_bcgsfhs_sample(/zsort)
    struct_print, sample
    ncl = n_elements(sample)

    filterlist = clash_filterlist(/dropbluest,short=short)
    nfilt = n_elements(filterlist)

; gather the photometry
;    cat = rsex(isedfit_dir+'flx_iso.dec26')

; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfile) then begin
       spsmodels = 'fsps_v2.4_miles'
       imf = 'salp'
       nmodel = 5000L
; age and metallicity are free
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, use_redshift=sample.z, $
         spsmodels=spsmodels, imf=imf, igm=0, redcurve='none', AV=[0.0,0.0], $
         sfhgrid=1, nmodel=nmodel, age=[6.2,11.2], tau=[0.0,2.0], $
         Zmetal=[0.01,0.03], /delayed, clobber=clobber
; metallicity is fixed at super-solar
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, use_redshift=sample.z, $
         spsmodels=spsmodels, imf=imf, igm=0, redcurve='none', AV=[0.0,0.0], $
         sfhgrid=2, nmodel=nmodel, age=[6.2,11.2], tau=[0.0,2.0], $
         Zmetal=[0.03,0.03], /delayed, clobber=clobber, /append
    endif

; --------------------------------------------------
; build the Monte Carlo grids    
    if keyword_set(build_grids) then begin
       isedfit_montegrids, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, clobber=clobber
    endif

; --------------------------------------------------
; calculate the model photometry 
    if keyword_set(model_photometry) then begin
       isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber
    endif

; --------------------------------------------------
; do the fitting!
    if keyword_set(isedfit) then begin
       for ic = 0, ncl-1 do begin
          cluster = strtrim(sample[ic].shortname,2)
          outprefix = prefix+'_'+cluster

          phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
          naper = n_elements(phot[0].photradius_kpc)+1 ; radial bins + integrated
          nband = n_elements(phot)

          maggies = fltarr(nfilt,naper)
          ivarmaggies = fltarr(nfilt,naper)
          for ii = 0, nfilt-1 do begin
             this = where(short[ii] eq strtrim(phot.band,2))
             maggies[ii,*] = [phot[this].maggies_int,phot[this].maggies]
             ivarmaggies[ii,*] = [phot[this].ivarmaggies_int,phot[this].ivarmaggies]
          endfor
          
          isedfit, isedfit_paramfile, maggies, ivarmaggies, replicate(sample[ic].z,naper), $
            thissfhgrid=thissfhgrid, isedfit_dir=isedfit_dir, isedfit_results=ised, $
            isedfit_post=isedpost, clobber=clobber, outprefix=outprefix
       endfor
    endif 

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         absmag_filterlist=bessell_filterlist(), band_shift=0.0, $
         clobber=clobber
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       for ic = 0, ncl-1 do begin
          cluster = strtrim(sample[ic].shortname,2)
          outprefix = prefix+'_'+cluster
          isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
            montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
            clobber=clobber, /xlog, galaxy=galaxy, index=index, $
            outprefix=outprefix
       endfor
    endif
    
return
end
