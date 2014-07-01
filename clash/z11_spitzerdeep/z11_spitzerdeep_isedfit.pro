pro z11_spitzerdeep_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber, photoz=photoz, lowz=lowz
; jm14jun23siena

    prefix = 'z11_spitzerdeep'
    if keyword_set(photoz) then prefix = 'z11_spitzerdeep_photoz'
    if keyword_set(lowz) then prefix = 'z11_spitzerdeep_lowz'

    isedfit_dir = getenv('CLASH_PROJECTS')+'/z11_spitzerdeep/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; gather the photometry
    cat = read_z11_spitzerdeep()

    filterlist = z11_filterlist()
    nfilt = n_elements(filterlist)

; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfile) then begin
       spsmodels = 'fsps_v2.4_miles'
       imf = 'chab'
       if keyword_set(photoz) then begin
          nmodel = 20000L
          write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
            prefix=prefix, filterlist=filterlist, zminmax=[0.1,12.0], nzz=100, $
            spsmodels=spsmodels, imf=imf, redcurve='charlot', /igm, pburst=0.2, $
            sfhgrid=1, nmodel=nmodel, age=[0.01,12.0], AV=[0.0,3.0], tau=[0.0,5.0], $
            Zmetal=[0.0008,0.03], oiiihb=[-0.3,1.0], /nebular, /delayed, /flatav, $
            modelchunksize=1000L, clobber=clobber
;         write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
;           prefix=prefix, filterlist=filterlist, zminmax=[8.0,12.0], nzz=50, $
;           spsmodels=spsmodels, imf=imf, redcurve='charlot', /igm, pburst=0.2, $
;           sfhgrid=1, nmodel=nmodel, age=[0.01,0.65], AV=[0.0,2.0], tau=[0.0,1], $
;           Zmetal=[0.0008,0.03], oiiihb=[-0.5,1.1], /nebular, /delayed, /flatav, $
;           modelchunksize=5000L, clobber=clobber
       endif 
       if keyword_set(lowz) then begin
          nmodel = 5000L
          write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
            prefix=prefix, filterlist=filterlist, zminmax=[2.5,2.7], nzz=3, $
            spsmodels=spsmodels, imf=imf, redcurve='charlot', /igm, $
            sfhgrid=1, nmodel=nmodel, age=[0.01,2.6], AV=[0.0,3.0], tau=[0.0,5.0], $
            Zmetal=[0.0008,0.03], oiiihb=[0.0,1.0], /nebular, /delayed, /flatav, $
            clobber=clobber
          write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
            prefix=prefix, filterlist=filterlist, zminmax=[2.5,2.7], nzz=3, $
            spsmodels=spsmodels, imf=imf, redcurve='calzetti', /igm, $
            sfhgrid=2, nmodel=nmodel, age=[0.01,2.6], AV=[0.0,3.0], tau=[0.0,5.0], $
            Zmetal=[0.0008,0.03], oiiihb=[0.0,1.0], /nebular, /delayed, /flatav, $
            /append
       endif
       if keyword_set(photoz) eq 0 and keyword_set(lowz) eq 0 then begin
          nmodel = 5000L
          write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
            prefix=prefix, filterlist=filterlist, zminmax=[10.5,10.9], nzz=5, $
            spsmodels=spsmodels, imf=imf, redcurve='none', /igm, $
            sfhgrid=1, nmodel=nmodel, age=[0.01,0.45], AV=[0.0,0.0], tau=[0.0,0.5], $
            Zmetal=[0.0008,0.019], oiiihb=[0.0,1.0], /nebular, /delayed, $
            clobber=clobber
       endif
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
       z11_to_maggies, cat, maggies, ivarmaggies, filterlist=filt
       if keyword_set(lowz) then z = cat.z*0+2.6 ; force the low-z solution
       if keyword_set(photoz) eq 0 and keyword_set(lowz) eq 0 then z = cat.z ; cat.bpz
       isedfit, isedfit_paramfile, maggies, ivarmaggies, z, thissfhgrid=thissfhgrid, $
         isedfit_dir=isedfit_dir, outprefix=outprefix, isedfit_results=ised, $
         isedfit_post=isedpost, clobber=clobber, photoz=photoz, index=index
    endif 

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         absmag_filterlist=galex_filterlist(), band_shift=0.0, $
         clobber=clobber
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, galaxy=cat.galaxy, yrange=[30.5,23], $
         xrange=[0.3,7.0]*1D4, nsigma=2.0, index=index, outprefix=outprefix
    endif
    
return
end
