pro hff_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, clobber=clobber
; jm12sep17siena

    usehawki = 0
    
    prefix = 'hizgalaxies'
    isedfit_dir = getenv('CLASH_PROJECTS')+'/hizgalaxies/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; gather the photometry
    cat = rsex(isedfit_dir+'sed_flx_moustakas.inp')
;   cat = rsex(isedfit_dir+'sed_flx.inp')
    cat = cat[multisort(cat.zb,cat.name)]
    use_redshift = cat.zb       ; custom redshift array

    filterlist = clash_filterlist(/useirac,usehawki=usehawki)
    nfilt = n_elements(filterlist)

; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfile) then begin
       spsmodels = 'fsps_v2.4_miles'
       imf = 'chab'
       redcurve = 'none'
       nmodel = 100
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, use_redshift=use_redshift, $
         spsmodels=spsmodels, imf=imf, redcurve=redcurve, /igm, $
         sfhgrid=1, nmodel=nmodel, age=[0.005,0.6], AV=[0.0,0.0], tau=[0.01,1.0], $
         Zmetal=[0.0008,0.03], oiiihb=[0.0,1.0], /nebular, /delayed, $
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
; generate the model photometry QAplots
    if keyword_set(qaplot_models) then begin
       thesefilters = ['clash_wfc3_f275w','clash_wfc3_f390w','clash_acs_f475w',$
         'clash_acs_f814w.par','clash_wfc3_f160w.par']
       hff_to_maggies, cat, maggies, ivarmaggies, /nJy, usehawki=usehawki
       isedfit_qaplot_models, isedfit_paramfile, maggies, $
         ivarmaggies, cat.z, isedfit_dir=isedfit_dir, $
         thissfhgrid=thissfhgrid, thesefilters=thesefilters, $
         zcolor_pdffile='', clobber=clobber
    endif
    
; --------------------------------------------------
; do the fitting!
    if keyword_set(isedfit) then begin
       hff_to_maggies, cat, maggies, ivarmaggies, /nJy, usehawki=usehawki
       isedfit, isedfit_paramfile, maggies, ivarmaggies, cat.zb, thissfhgrid=thissfhgrid, $
         isedfit_dir=isedfit_dir, outprefix=outprefix, isedfit_results=ised, $
         isedfit_post=isedpost, clobber=clobber
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
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, galaxy=cat.name, yrange=[27,21]
    endif
    
return
end
