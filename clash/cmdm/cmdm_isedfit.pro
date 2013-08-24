pro cmdm_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  models=models, isedfit=isedfit, qaplot=qaplot, clobber=clobber, $
  thissfhgrid=thissfhgrid, fordoron=fordoron
; jm13jul18siena

    prefix = 'cmdm'
    datapath = getenv('IM_PROJECTS_DIR')+'/clash/cmdm/'
    isedfit_dir = datapath+'isedfit/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; read the catalog and then specify the redshifts and filters
    cat = mrdfits(datapath+'cmdm_cat.fits.gz',1)
    ngal = n_elements(cat)

    filterlist = clash_subaru_filterlist()
    nfilt = n_elements(filterlist)
    
    nmodel = 500L ; 30000L

; --------------------------------------------------
; do the preliminaries: build the parameter files and the Monte Carlo
; grids
    if keyword_set(write_paramfile) then begin
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, zminmax=[0.40,0.48], nzz=10, $
         spsmodels='fsps_v2.4_miles', imf='salp', redcurve='calzetti', /igm, $
         sfhgrid=1, nmodel=nmodel, age=[0.1,9.2], tau=[0.01,5.0], $
         Zmetal=[0.003,0.03], pburst=0.1, interval_pburst=1.5, $
         oiiihb=[-1.0,1.0], /nebular, /delayed, clobber=clobber
    endif

    if keyword_set(build_grids) then begin
       isedfit_montegrids, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, clobber=clobber
    endif

; --------------------------------------------------
; build the models
    if keyword_set(models) then begin
       isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, clobber=clobber, thissfhgrid=thissfhgrid
    endif

; --------------------------------------------------
; do the fitting!
    if keyword_set(isedfit) then begin
       clash_subaru_to_maggies, cat, maggies, ivarmaggies
       isedfit, isedfit_paramfile, maggies, ivarmaggies, cat.zobj, ised, $
         isedfit_dir=isedfit_dir, $
         outprefix=outprefix, thissfhgrid=thissfhgrid, clobber=clobber
    endif 

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
       indx = shuffle_indx(ngal,num=20)
       isedfit_qaplot, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         outprefix=outprefix, clobber=clobber, $
         yrange=yrange, /xlog, result=result, index=indx
    endif

; --------------------------------------------------
; package the results for Doron
    if keyword_set(fordoron) then begin
stop       

;       ised = read_isedfit(isedfit_paramfile,isedfit_dir=isedfit_dir,$
;         supergrid_paramfile=supergrid_paramfile)
;       ngal = n_elements(ised)
;       mstar = isedfit_reconstruct_posterior(isedfit_paramfile,$
;         supergrid_paramfile=supergrid_paramfile,thissupergrid=thissupergrid,$
;         isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir)
;
;       
;       ndraw = isedfit_ndraw()
;       out = struct_addtags(cat,struct_trimtags(ised,select=['mass_50','mass_err']))
;       out = struct_addtags(out,replicate({pofm: fltarr(ndraw)},ngal))
;       out.pofm = mstar
;
;       wsex, struct_trimtags(out,except='pofm'), outfile=datapath+'cmdm_isedfit.txt'
;       openw, lun, datapath+'cmdm_isedfit_pofm.txt', /get_lun
;       for ii = 0L, ngal-1 do printf, lun, 
;       free_lun, lun
;       im_mwrfits, out, datapath+'cmdm_isedfit.fits', clobber=clobber
    endif
    
return
end
