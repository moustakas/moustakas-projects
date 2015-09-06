pro clash_redgals_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber
; jm14sep07siena

; echo "decals_dr1_isedfit, /ised, /qaplot_sed, /cl" | /usr/bin/nohup idl > & logall.log & 
; echo "decals_dr1_isedfit, /ised, /qaplot_sed, /nosdss, /cl" | /usr/bin/nohup idl > & lognosdss.log & 
; echo "decals_dr1_isedfit, /ised, /qaplot_sed, /nodecals, /cl" | /usr/bin/nohup idl > & lognodecals.log & 

    prefix = 'clash_redgals'
    isedfit_dir = getenv('CLASH_PROJECTS')+'/redgals/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    clash = clash[sort(clash.z)]

    cat = rsex(isedfit_dir+'redgals.cat')

    zobj = cat.zml
    clash_redgals_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist

stop

; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfile) then begin
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels='fsps_v2.4_miles', $
         imf='salp', redcurve='charlot', /igm, use_redshift=clash.z, $
         nmodel=10000L, age=[1.0,11.5], tau=[0.1,10.0], Zmetal=[0.005,0.03], $
         clobber=clobber
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
       isedfit, isedfit_paramfile, maggies, ivarmaggies, zobj, $
         isedfit_dir=isedfit_dir, outprefix=outprefix, isedfit_results=ised, $
         isedfit_post=isedpost, clobber=clobber, thissfhgrid=thissfhgrid, $
         photoz=photoz, index=index;, ra=cat.ra, dec=cat.dec
    endif 

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         absmag_filterlist=galex_filterlist(), band_shift=0.0, $
         clobber=clobber, outprefix=outprefix
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
;      index = lindgen(50)
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, nrandom=40, galaxy=galaxy, $
         index=index, outprefix=outprefix
    endif
    
return
end
