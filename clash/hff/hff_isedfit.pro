pro hff_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber, photoz=photoz
; jm12sep17siena

    usehawki = 0

    if keyword_set(photoz) then prefix = 'hff_photoz' else prefix = 'hff'
    isedfit_dir = getenv('CLASH_PROJECTS')+'/hff/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; gather the photometry
    cat = rsex(isedfit_dir+'flx_iso.jan03')
    cat.bpz = abs(cat.bpz)

;   oldcat = rsex(isedfit_dir+'a2744_13dec03.lis')
;   cat = rsex(isedfit_dir+'a2744_13dec08.lis')
;   cat.bpz = oldcat.bpz
;   cat = rsex(isedfit_dir+'flx_iso.lis')

    if keyword_set(photoz) eq 0 then begin
       cat = cat[where(cat.bpz ge 7.0 and cat.bpz lt 10.01)]
       cat = cat[reverse(sort(cat.bpz))]
;      cat[where(cat.id eq 1580)].bpz = 5.45
    endif

;   if keyword_set(photoz) eq 0 then cat = cat[where(cat.bpz ge 7.0)]
    galaxy = 'ID'+strtrim(cat.id,2)

;   cat = cat[multisort(cat.bpz,cat.name)]
;   use_redshift = cat.bpz       ; custom redshift array
;   cat = cat[multisort(cat.zb,cat.name)]
;   use_redshift = cat.zb       ; custom redshift array

    filterlist = hff_filterlist(/useirac,usehawki=usehawki)
    nfilt = n_elements(filterlist)

; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfile) then begin
       spsmodels = 'fsps_v2.4_miles'
       imf = 'chab'
       if keyword_set(photoz) then begin
          nmodel = 5000 
          redcurve = 'calzetti'
          write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
            prefix=prefix, filterlist=filterlist, zminmax=[0.5,11.0], nzz=50, $
            spsmodels=spsmodels, imf=imf, redcurve=redcurve, /igm, $
            sfhgrid=1, nmodel=nmodel, age=[0.01,8.0], AV=[0.0,3.0], tau=[0.01,6.0], $
            Zmetal=[0.0008,0.03], oiiihb=[0.0,1.0], /nebular, /delayed, /flatav, $
            modelchunksize=1000L, clobber=clobber
       endif else begin
          nmodel = 10000L
          write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
            prefix=prefix, filterlist=filterlist, zminmax=[7.0,11], nzz=20, $
            spsmodels=spsmodels, imf=imf, redcurve=redcurve, /igm, $
            sfhgrid=1, nmodel=nmodel, age=[0.01,0.75], AV=[0.0,0.0], tau=[0.01,1.0], $
            Zmetal=[0.0008,0.019], oiiihb=[0.0,1.0], /nebular, /delayed, $
            clobber=clobber
       endelse
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
       hff_to_maggies, cat, maggies, ivarmaggies, /nJy, usehawki=usehawki, filterlist=filt
       if keyword_set(photoz) eq 0 then z = cat.bpz
;      if keyword_set(photoz) eq 0 then z = cat.zb
;      if keyword_set(photoz) then index = where(cat.irac_ch1_flux gt 0.0 or cat.irac_ch2_flux gt 0.0)
       isedfit, isedfit_paramfile, maggies, ivarmaggies, z, thissfhgrid=thissfhgrid, $
         isedfit_dir=isedfit_dir, outprefix=outprefix, isedfit_results=ised, $
         isedfit_post=isedpost, clobber=clobber, photoz=photoz, index=index
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
;      index = where(cat.irac_ch1_flux gt 0.0 or cat.irac_ch2_flux gt 0.0)

;      cat = cat[where(cat.bpz ge 7.0 and cat.bpz lt 10.01)]
;      cat = cat[reverse(sort(cat.bpz))]
;      index = where(cat.id ne 1580 and cat.id ne 490 and $
;        (cat.irac_ch1_flux gt 0.0 or cat.irac_ch2_flux gt 0.0))
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, galaxy=galaxy, yrange=[30.5,23], $
         xrange=[0.3,7.0]*1D4, nsigma=2.0, index=index
    endif
    
return
end
