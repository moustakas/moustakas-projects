pro irclusters_write_paramfile, paramfile, igm=igm, super=super

    sfhgridstring = strtrim(super.sfhgrid,2)
    redcurvestring = redcurve2string(super.redcurve)
    synthmodels = strtrim(super.synthmodels,2)
    imf = strtrim(super.imf,2)
    
    splog, 'Writing '+paramfile
    zrange = string(zminmax[0],format='(F4.2)')+','+string(zminmax[1],$
      format='(F4.2)')+','+nzz+' # [minz,maxz,dz]'
    openw, lun, paramfile, /get_lun
    printf, lun, 'synthmodels          bc03'
    printf, lun, 'imf                  chab'
    printf, lun, 'sfhgrid              '+sfhgridstring
    printf, lun, 'redcurve             calzetti'
    printf, lun, 'prefix               '+prefix
    printf, lun, 'redshift             '+zrange
    printf, lun, 'igm                  '+igm+' # [0=no, 1=yes]'
    printf, lun, 'maxold               0 # [0=no, 1=yes]'
    printf, lun, 'filterlist           ndwfs_Bw.par,ndwfs_R.par,ndwfs_I.par,newfirm_J.par,newfirm_H.par,newfirm_Ks.par,spitzer_irac_ch1.par,spitzer_irac_ch2.par,spitzer_irac_ch3.par,spitzer_irac_ch4.par'
    free_lun, lun 

return
end

pro irclusters_isedfit, make_montegrids=make_montegrids, models=models, $
  isedfit=isedfit, qaplot=qaplot, chi2grid_qaplot=chi2grid_qaplot, clobber=clobber, $
  debug=debug
; jm10oct26ucsd - measure stellar masses for galaxies in
; Brodwin's sample of IR-selected clusters
; jm11aug30ucsd - everything updated

    iopath = ages_path(/projects)+'irclusters/'
    isedfit_sfhgrid_dir = iopath+'montegrids/'
    sfhgrid_paramfile = getenv('IDL_PROJECTS_DIR')+'/ages/projects/irclusters/irclusters_sfhgrid.par'

    paramfile = iopath+'irclusters_isedfit.par'
    param = read_isedfit_paramfile(paramfile)

    cat = read_irclusters(maggies=maggies,ivarmaggies=ivarmaggies)
;   test = lindgen(200)+500
;   cat = cat[test]
;   maggies = maggies[*,test]
;   ivarmaggies = ivarmaggies[*,test]
    ngal = n_elements(cat)

; --------------------------------------------------
; build the Monte Carlo grids
    if keyword_set(make_montegrids) then begin
; bursts, calzetti
       build_isedfit_sfhgrid, 1, synthmodels='bc03', imf='chab', redcurve=0, $
         /make_montegrid, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
         sfhgrid_paramfile=sfhgrid_paramfile
; no bursts, calzetti
       build_isedfit_sfhgrid, 1, synthmodels='bc03', imf='chab', redcurve=0, $
         /make_montegrid, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
         sfhgrid_paramfile=sfhgrid_paramfile
    endif

; --------------------------------------------------
; build the models
    if keyword_set(models) then begin
       irclusters_write_paramfile, paramfile, field=field[ii], zminmax=zminmax, $
         nzz=nzz, filters=filters, igm=igm, super=super[gg]
       isedfit_models, paramfile, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
         iopath=iopath, clobber=clobber
    endif

; --------------------------------------------------
; do the fitting!  do not use the NEWFIRM photometry nor IRAC/ch3-4
    if keyword_set(isedfit) then begin
       toss = where(strmatch(param.filterlist,'*ch3*') or $
         strmatch(param.filterlist,'*ch4*'))
       ivarmaggies[toss,*] = 0.0
;; fit at both the photometric redshift...
;       index = where((cat.z ge min(param.redshift)) and $
;         (cat.z le max(param.redshift)))
;       isedfit, paramfile, maggies, ivarmaggies, cat.z, $
;         iopath=iopath, clobber=clobber, index=index
;; ...and at the cluster redshift...
;;      res = mrdfits('BwRIJHKsirac_zclust_bc03_chab_calzetti_sfhgrid02.fits.gz',1)
;;      index = (where((res.mass_err lt 0.01) and (res.zobj gt 1.2) and (res.zobj lt 1.3) ))[0:49]
;       isedfit, paramfile, maggies, ivarmaggies, cat.zclust, $
;         iopath=iopath, outprefix=param.prefix+'_zclust', $
;         clobber=clobber;, index=index
;; ...and then refit at the photometric redshift excluding IRAC 
       index = where((cat.z ge min(param.redshift)) and $
         (cat.z le max(param.redshift)))
;      index = [200,239]
       alltoss = where(strmatch(param.filterlist,'*spitzer*'))
       ivarmaggies[alltoss,*] = 0.0
       isedfit, paramfile, maggies, ivarmaggies, cat.z, $
         iopath=iopath, outprefix=repstr(param.prefix,'irac',''), $
         clobber=clobber, index=index
    endif 

; --------------------------------------------------
; make some QAplots of a random subset of objects
    if keyword_set(qaplot) then begin
       qaindx = long(randomu(seed,100)*ngal)
       qaindx = qaindx[uniq(qaindx,sort(qaindx))]
       galaxy = 'Galaxy '+string(qaindx,format='(I4.4)')
       
; includes irac, photoz
       isedfit_qaplot, paramfile, iopath=iopath, index=qaindx, $
         galaxy=galaxy, clobber=clobber, outprefix=outprefix
; includes irac, zclust
       isedfit_qaplot, paramfile, iopath=iopath, index=qaindx, $
         galaxy=galaxy, clobber=clobber, outprefix=param.prefix+'_zclust'
; no irac, photoz
       isedfit_qaplot, paramfile, iopath=iopath, index=qaindx, $
         galaxy=galaxy, clobber=clobber, outprefix=repstr(param.prefix,'irac','')
    endif

    if keyword_set(chi2grid_qaplot) then begin
       isedfit_chi2grid_qaplot, paramfile, isedfit, iopath=iopath, galaxy=galaxy, $
         index=index, clobber=clobber, outprefix=outprefix
    endif

return
end
