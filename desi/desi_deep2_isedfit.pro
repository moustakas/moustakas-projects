pro desi_deep2_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, clobber=clobber
; jm13dec18siena - fit the parent sample of DEEP2 galaxies for the
; DESI project
    
    prefix = 'desi_deep2'
    isedfit_dir = getenv('IM_PROJECTS_DIR')+'/desi/deep2/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    filterlist = deep2_filterlist()

    cat = mrdfits(isedfit_dir+'deep2_zcat.fits.gz',1)
    
; --------------------------------------------------
; write the parameter file
    if keyword_set(write_paramfile) then begin
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels='fsps_v2.4_miles', $
         imf='chab', redcurve='charlot', /igm, zminmax=[0.7,1.6], nzz=30.0, $
         nmodel=1000L, age=[0.1,7.2], tau=[0.01,7], Zmetal=[0.004,0.03], $
         pburst=0.2, interval_pburst=2.0, clobber=clobber
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
; generate the model photometry QAplots
    if keyword_set(qaplot_models) then begin
       thesefilters = ['deep_B','deep_R','deep_I']
       isedfit_qaplot_models, isedfit_paramfile, cat.maggies, $
         cat.ivarmaggies, cat.z, isedfit_dir=isedfit_dir, $
         thissfhgrid=thissfhgrid, thesefilters=thesefilters, clobber=clobber
    endif
    
; --------------------------------------------------
; fit!
    if keyword_set(isedfit) then begin
       isedfit, isedfit_paramfile, cat.maggies, cat.ivarmaggies, $
         cat.zbest, ra=cat.ra, dec=cat.dec, isedfit_dir=isedfit_dir, $
         thissfhgrid=thissfhgrid, clobber=clobber
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
       galaxy = strtrim(cat.objno,2)+'/'+strtrim(cat.source,2)
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, nrandom=50, galaxy=galaxy
    endif

return
end
