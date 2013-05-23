pro redmapper_isedfit, prelim=prelim, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber, bcgs=bcgs, firstchunk=firstchunk, $
  lastchunk=lastchunk
; jm13mar28siena

    prefix = 'redmapper'
    catalogs_dir = redmapper_path(/catalogs,version=ver)
    isedfit_dir = redmapper_path(/isedfit)
    
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    sfhgrid_paramfile = isedfit_dir+prefix+'_sfhgrid.par'
    supergrid_paramfile = isedfit_dir+prefix+'_supergrid.par'

; global parameters    
    filterlist = redmapper_filterlist()
    zmin = 0.05
    zmax = 0.6
    nzz = 50

; SFHgrid priors    
    nage = 20
    nmonte = 1000 ; 1000
    Z = [0.004,0.04]
    minage = 3.0
    maxage = 13.0
;   AV = [0.0,3.0]
    tau = [0.01,5.0]
    delayed = 1
    
; supergrid parameters    
    sfhgrid = 1
    supergrid = 1
    synthmodels = 'fsps'

    imf = 'chab'
    redcurve = 1 ; Charlot & Fall

; --------------------------------------------------
; do the preliminaries: build the parameter files and the Monte Carlo
; grids
    if keyword_set(prelim) then begin
       write_isedfit_paramfile, filterlist, prefix=prefix, minz=zmin, $
         maxz=zmax, nzz=nzz, igm=1, isedfit_dir=isedfit_dir, clobber=clobber

; delayed
       write_sfhgrid_paramfile, sfhgrid_paramfile, params, clobber=clobber, $
         sfhgrid=sfhgrid, nage=nage, nmonte=nmonte, Z=Z, minage=minage, $
         maxage=maxage, AV=AV, tau=tau, flatAV=flatAV, delayed=delayed

       write_supergrid_paramfile, supergrid_paramfile, supergrid=supergrid, $
         sfhgrid=sfhgrid, synthmodels=synthmodels, imf=imf, redcurve=redcurve, $
         clobber=clobber

       build_montegrids, sfhgrid_paramfile, supergrid_paramfile=supergrid_paramfile, $
         isedfit_dir=isedfit_dir, clobber=clobber
    endif

; --------------------------------------------------
; build the models
    if keyword_set(models) then begin
       isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         supergrid_paramfile=supergrid_paramfile, clobber=clobber
    endif

; --------------------------------------------------
; do the fitting
    if keyword_set(isedfit) then begin
       phot = mrdfits(catalogs_dir+'redmapper_'+ver+'_phot.fits.gz',1)
       ngal = n_elements(phot)

; split the members catalog into chunks since it's so big
       if keyword_set(bcgs) then begin
          outprefix = 'bcgs'
          phot = phot[where(phot.isbcg)]

          isedfit, isedfit_paramfile, phot.maggies, phot.ivarmaggies, $
            phot.z, result, isedfit_dir=isedfit_dir, clobber=clobber, $
            supergrid_paramfile=supergrid_paramfile, outprefix=outprefix, $
            sfhgrid_paramfile=sfhgrid_paramfile, galchunksize=2000L
       endif else begin
          nchunk = 8
          chunksize = ceil(ngal/float(nchunk))
          if n_elements(firstchunk) eq 0 then firstchunk = 0
          if n_elements(lastchunk) eq 0 then lastchunk = nchunk-1

          for ii = firstchunk, lastchunk do begin
             splog, 'Working on CHUNK '+strtrim(ii,2)+'/'+strtrim(lastchunk+1,2)
             splog, im_today()
             t0 = systime(1)
             outprefix = 'redmapper_chunk'+string(ii,format='(I0)')
             these = lindgen(chunksize)+ii*chunksize
             these = these[where(these lt ngal)]
;            these = these[0:99] ; test!
             
             isedfit, isedfit_paramfile, phot[these].maggies, phot[these].ivarmaggies, $
               phot[these].z, result, isedfit_dir=isedfit_dir, clobber=clobber, $
               supergrid_paramfile=supergrid_paramfile, outprefix=outprefix, $
               sfhgrid_paramfile=sfhgrid_paramfile, galchunksize=1000L
             splog, 'Total time (min) = '+strtrim((systime(1)-t0)/60.0,2)
          endfor
       endelse 
    endif

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
       phot = mrdfits(catalogs_dir+'redmapper_'+ver+'_phot.fits.gz',1)
       allcl = phot.mem_match_id
       cl = allcl[uniq(allcl,sort(allcl))]

;      jj = mrdfits('testcl1_fsps_chab_charlot_sfhgrid01.fits.gz',1)
;      index = where(phot.mem_match_id eq 12)
;      outprefix = 'testcl1'
       
       jj = mrdfits('redmapper_fsps_chab_charlot_sfhgrid01.fits.gz',1)
       phot = phot[where(phot.isbcg,ngal)]
       these = where(jj.mass_50 gt 11.5 and jj.mass_50 lt 12.5,nthese)
;      these = where(jj.chi2 lt 5,nthese)
       index = these[shuffle_indx(nthese,num=50)]
;;       index = shuffle_indx(ngal,num=50)
;;      index = where(jj.mass_50 gt 12.0)
;;      index = lindgen(50)
       galaxy = 'BCG '+string(phot[index].mem_match_id,format='(I6.6)')

       isedfit_qaplot, isedfit_paramfile, supergrid_paramfile=supergrid_paramfile, $
         isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, index=index, $
         galaxy=galaxy, outprefix=outprefix, clobber=clobber, /xlog
    endif
    
    
return
end
    
