pro desi_deep2_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, clobber=clobber
; jm13dec18siena - fit the parent sample of DEEP2 galaxies for the
; DESI project
    
    version = 'v1.1'

    prefix = 'desi_deep2'
    isedfit_dir = getenv('IM_PROJECTS_DIR')+'/desi/templates/'+version+'/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

;   filterlist = deep2_filterlist()
;   cat = mrdfits(isedfit_dir+'deep2_zcat.fits.gz',1)
    
; fit everything in DR4 so that I can use the iSEDfit results for both
; targeting tests and template simulations
    cat = read_deep2_zcat(photo=phot)
    deep2_to_maggies, phot, maggies, ivarmaggies, /unwise, $
      filterlist=filterlist

    zminmax = [0.1,2.0]
    index = where(cat.zbest ge zminmax[0] and cat.zbest le zminmax[1])
    ngal = n_elements(cat)

; --------------------------------------------------
; write the parameter file
    if keyword_set(write_paramfile) then begin
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels='fsps_v2.4_miles', $
         imf='chab', redcurve='charlot', /igm, zminmax=zminmax, nzz=40.0, $
         nmodel=10000L, age=[0.1,7.2], tau=[0.01,7], Zmetal=[0.004,0.03], $
         pburst=0.2, interval_pburst=1.0, clobber=clobber
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels='fsps_v2.4_miles', $
         imf='chab', redcurve='charlot', /igm, zminmax=zminmax, nzz=40.0, $
         nmodel=10000L, age=[0.1,7.2], tau=[0.01,7], Zmetal=[0.004,0.03], $
         pburst=0.2, interval_pburst=1.0, /nebular, /append
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
       isedfit_qaplot_models, isedfit_paramfile, maggies, $
         ivarmaggies, cat.zbest, isedfit_dir=isedfit_dir, $
         thissfhgrid=thissfhgrid, thesefilters=thesefilters, clobber=clobber
    endif
    
; --------------------------------------------------
; fit!
    if keyword_set(isedfit) then begin
;      outprefix = 'unwise'
;      index = where(phot.w1_nanomaggies_ivar ne 0 and cat.zbest ge zminmax[0] and $
;        cat.zbest le zminmax[1])
       isedfit, isedfit_paramfile, maggies, ivarmaggies, $
         cat.zbest, ra=cat.ra, dec=cat.dec, isedfit_dir=isedfit_dir, $
         thissfhgrid=thissfhgrid, clobber=clobber, index=index, $
         outprefix=outprefix
    endif 

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
;      index = [9,23,29]
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         absmag_filterlist=sdss_filterlist(), band_shift=0.0, $
         clobber=clobber, index=index, outprefix=outprefix
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
;      outprefix = 'unwise'
;      index = (where(phot.w1_nanomaggies_ivar ne 0 and cat.zbest ge zminmax[0] and $
;        cat.zbest le zminmax[1]))[0:30]
       galaxy = 'DEEP2/'+strtrim(cat.objno,2);+'/'+strtrim(cat.source,2)
       these = shuffle_indx(ngal,num=25)
;      these = where(cat[index].objno eq 12024524) & yrange = [24,20]
;      these = where(cat[index].objno eq 12024078)
;      these = where(cat[index].objno eq 12101118)
;      these = where(cat[index].objno eq 12015944)
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, galaxy=galaxy, index=index[these];, yrange=yrange
;        yrange=[26,15] ;, outprefix=outprefix
    endif

return
end
