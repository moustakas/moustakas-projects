pro clash_spitzer_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, clobber=clobber, $
  thissfhgrid=thissfhgrid
; jm13aug22siena - 

    isedfit_dir = getenv('IM_RESEARCH_DIR')+'/projects/clash/spitzer/'
    montegrids_dir = isedfit_dir+'montegrids/'

    prefix = 'spitzer'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

;   cc = rsex(isedfit_dir+'abell_209_final.cat')
;   keep = where(cc.f160w_mag lt 90 and cc.redshift gt 0.0 and cc.redshift lt 100.0)
;   wsex, cc[keep], outfile='abell_209_final_hstirac.cat'
    cat = rsex(isedfit_dir+'abell_209_final_hstirac.cat')
    clash_to_maggies, cat, maggies, ivarmaggies, /useirac, $
      filterlist=filterlist
    zz = cat.redshift

    zminmax = minmax(zz)
    zbin = 0.01 ; fix this!

; --------------------------------------------------
; (mandatory) choose your priors: write the iSEDfit parameter file 
    if keyword_set(write_paramfile) then begin
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels='fsps_v2.4_miles', $
         imf='salp', redcurve='charlot', /igm, zminmax=zminmax, zbin=zbin, $
         nmodel=1000L, age=[0.1,11.0], tau=[0.1,10.0], Zmetal=[0.004,0.03], $
         pburst=0.1, interval_pburst=1.5, tburst=[0.1,11.0], /delayed, $
         /nebular, clobber=clobber
    endif

; --------------------------------------------------
; (mandatory) build the Monte Carlo grids    
    if keyword_set(build_grids) then begin
       isedfit_montegrids, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, clobber=clobber
    endif

; --------------------------------------------------
; (mandatory) calculate the model photometry 
    if keyword_set(model_photometry) then begin
       isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber
    endif

; --------------------------------------------------
; (optional) generate the model photometry QAplots
    if keyword_set(qaplot_models) then begin
       thesefilters = ['clash_wfc3_f390w.par','clash_acs_f814w',$
         'clash_wfc3_f110w','clash_wfc3_f160w']
       isedfit_qaplot_models, isedfit_paramfile, maggies, ivarmaggies, $
         zz, isedfit_dir=isedfit_dir, thissfhgrid=thissfhgrid, $
         thesefilters=thesefilters, clobber=clobber
    endif
    
; --------------------------------------------------
; (mandatory) fit!
    if keyword_set(isedfit) then begin
       isedfit, isedfit_paramfile, maggies, ivarmaggies, zz, $
         ra=cat.ra, dec=cat.dec, isedfit_dir=isedfit_dir, $
         thissfhgrid=thissfhgrid, clobber=clobber
    endif 

; --------------------------------------------------
; (optional) compute K-corrections
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber
    endif 

; --------------------------------------------------
; (optional) generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       isedfit_qaplot_sed, isedfit_paramfile, $ ; nrandom=50, $
         isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, $
         thissfhgrid=thissfhgrid, clobber=clobber, /xlog
    endif

return
end
