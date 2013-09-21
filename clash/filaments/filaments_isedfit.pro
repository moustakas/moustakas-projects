pro filaments_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, clobber=clobber
; jm12sep17siena

    prefix = 'filaments'
    isedfit_dir = getenv('CLASH_PROJECTS')+'/filaments/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; gather the photometry for RXJ1532
    cat = rsex(isedfit_dir+'Filaments_0_6_SED_SextractorFMT.txt')
    cat = struct_addtags(replicate({name: '', z: 0.0},n_elements(cat)),cat)
    cat.name = 'Filament '+string(cat.id,format='(I2.2)')
    
    cat = rsex(isedfit_dir+'Redenning_Map_Cells_0_118_SED_SextractorFMT.txt')
    cat = struct_addtags(replicate({name: '', z: 0.0},n_elements(cat)),cat)
    cat.name = 'Grid '+string(cat.id,format='(I3.3)')
    cat.z = 0.363
    zminmax = [0.363,0.363]
    nzz = 1
    
;   use_redshift = cat.z ; custom redshift array

    filterlist = clash_filterlist()
    nfilt = n_elements(filterlist)

; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfile) then begin
       spsmodels = 'fsps_v2.4_miles'
       imf = 'chab'
       redcurve = 'calzetti'
       nmodel = 1000
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, zminmax=zminmax, nzz=nzz, $
         spsmodels=spsmodels, imf=imf, redcurve=redcurve, /igm, $
         sfhgrid=1, nmodel=nmodel, age=[0.005,0.5], tau=[0.1,5.0], $
         Zmetal=[0.0008,0.03], AV=[0.0,5.0], pburst=0.0, interval_pburst=1.0, $
         oiiihb=[0.0,1.0], /nebular, /flatAV, /delayed, clobber=clobber

;      write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
;        prefix=prefix, filterlist=filterlist, use_redshift=use_redshift, $
;        spsmodels=spsmodels, imf=imf, redcurve=redcurve, /igm, $
;        sfhgrid=2, nmodel=nmodel, age=[0.05,6.0], tau=[0.01,6.0], $
;        Zmetal=[0.0008,0.03], AV=[0.0,3.0], pburst=0.1, interval_pburst=1.0, $
;        nebular=0, /flatAV, /delayed, /append, clobber=clobber
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

;; --------------------------------------------------
;; generate the model photometry QAplots
;    if keyword_set(qaplot_models) then begin
;       thesefilters = ['galex_NUV','sdss_g0','sdss_r0','sdss_i0','wise_w1']
;       isedfit_qaplot_models, isedfit_paramfile, cat.maggies, $
;         cat.ivarmaggies, cat.z, isedfit_dir=isedfit_dir, $
;         thissfhgrid=thissfhgrid, thesefilters=thesefilters, clobber=clobber
;    endif
    
; --------------------------------------------------
; do the fitting!
    if keyword_set(isedfit) then begin
       clash_to_maggies, cat, maggies, ivarmaggies, /nodustcorr
       isedfit, isedfit_paramfile, maggies, ivarmaggies, cat.z, thissfhgrid=thissfhgrid, $
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
         clobber=clobber, /xlog, nrandom=20, galaxy=cat.name
    endif
    
return
end
