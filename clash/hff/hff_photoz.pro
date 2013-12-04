pro hizgalaxies_photoz, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, isedfit_photoz=isedfit_photoz, qaplot_photoz=qaplot_photoz, $
  thissfhgrid=thissfhgrid, clobber=clobber
; jm13dec04siena 

    usehawki = 0
    
    prefix = 'hff_photoz'
    isedfit_dir = getenv('CLASH_PROJECTS')+'/hizgalaxies/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; gather the photometry
    cat = rsex(isedfit_dir+'flx_iso.lis')
;   cat = rsex(isedfit_dir+'sed_flx_moustakas.inp')
;   cat = rsex(isedfit_dir+'sed_flx.inp')
;   cat = cat[multisort(cat.zb,cat.name)]
;   use_redshift = cat.zb       ; custom redshift array

    filterlist = clash_filterlist(/useirac,usehawki=usehawki)
    nfilt = n_elements(filterlist)

; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfile) then begin
       spsmodels = 'fsps_v2.4_miles'
       imf = 'chab'
       redcurve = 'calzetti'
       nmodel = 1000
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, zminmax=[1.0,13.0], nzz=100, $
         spsmodels=spsmodels, imf=imf, redcurve=redcurve, /igm, $
         sfhgrid=1, nmodel=nmodel, age=[0.01,6.0], AV=[0.0,3.0], tau=[0.01,10.0], $
         Zmetal=[0.0008,0.03], oiiihb=[-0.5,1.0], /nebular, /delayed, /flatav, $
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

;; --------------------------------------------------
;; generate the model photometry QAplots
;    if keyword_set(qaplot_models) then begin
;       thesefilters = ['clash_wfc3_f275w','clash_wfc3_f390w','clash_acs_f475w',$
;         'clash_acs_f814w.par','clash_wfc3_f160w.par']
;       hizgalaxies_to_maggies, cat, maggies, ivarmaggies, /nJy
;       isedfit_qaplot_models, isedfit_paramfile, maggies, $
;         ivarmaggies, cat.z, isedfit_dir=isedfit_dir, $
;         thissfhgrid=thissfhgrid, thesefilters=thesefilters, $
;         zcolor_pdffile='', clobber=clobber
;    endif
    
; --------------------------------------------------
; get the redshifts!
    if keyword_set(isedfit_photoz) then begin
       hizgalaxies_to_maggies, cat, maggies, ivarmaggies, /nJy, usehawki=usehawki
       isedfit_photoz, isedfit_paramfile, maggies, ivarmaggies, thissfhgrid=thissfhgrid, $
         isedfit_dir=isedfit_dir, outprefix=outprefix, isedfit_photoz_results=ised_photoz, $
         isedfit_photoz_post=ised_photoz_post, clobber=clobber
    endif 

;; --------------------------------------------------
;; compute K-corrections
;    if keyword_set(kcorrect) then begin
;       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
;         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
;         absmag_filterlist=bessell_filterlist(), band_shift=0.0, $
;         clobber=clobber
;    endif 

; --------------------------------------------------
; generate photoz QAplots
    if keyword_set(qaplot_photoz) then begin
       isedfit_qaplot_photoz, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, galaxy=cat.name, yrange=[27,21]
    endif
    
return
end
