pro irclusters_isedfit, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber, use_zcluster=use_zcluster
; jm10oct26ucsd - measure stellar masses for galaxies in Brodwin's sample of IR-selected clusters
; jm11aug30ucsd - everything updated
; jm11dec16ucsd - and again, a major rewrite

; echo "clash_bcg_isedfit, /model, /ised, /clob" | idl > & ~/clash.bcg.log & 

    irpath = ages_path(/projects)+'irclusters/'
    isedpath = irpath+'isedfit/'

    isedfit_sfhgrid_dir = irpath+'montegrids/'
    sfhgrid_paramfile = getenv('IDL_PROJECTS_DIR')+'/ages/projects/irclusters/irclusters_sfhgrid.par'

; read the supergrid parameter file
    super = get_irclusters_supergrid(supergrid,nsuper=nsuper)
    struct_print, super

    igm = '0'
    nzz = '100'
    zlog = '0'
    zminmax = [0.01D,2D]
    prefix = 'irclusters'

; loop on each supergrid
    for gg = 0, nsuper-1 do begin
       splog, 'working on grid '+strtrim(super[gg].supergrid,2)

       paramfile = isedpath+prefix+'_supergrid'+string(super[gg].supergrid,$
         format='(i2.2)')+'_isedfit.par'
       irclusters_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
         nzz=nzz, zlog=zlog, igm=igm, super=super[gg], filters=filters
       
; build the models
       if keyword_set(models) then begin
          isedfit_models, paramfile, iopath=isedpath, clobber=clobber, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       endif
       
; optionally fit at the cluster redshift  
       if keyword_set(isedfit) then begin
          cat = read_irclusters(maggies=maggies,ivarmaggies=ivarmaggies)

          toss = where(strmatch(filters,'*ch3*') or strmatch(filters,'*ch4*'))
          ivarmaggies[toss,*] = 0.0

          delvarx, index, outprefix
          if keyword_set(use_zcluster) then begin
             zobj = cat.zclust
             outprefix = prefix+'_zclust'
          endif else begin
; only fit objects at 0<z<2
;            index = where(total(maggies gt 0,1) lt 3) ; test!!
;            index = [index,index+1]
             zobj = cat.z
             index = where((cat.z ge zminmax[0]) and (cat.z le zminmax[1]),nindex)
          endelse
          isedfit, paramfile, maggies, ivarmaggies, zobj, iopath=isedpath, $
            clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, index=index, outprefix=outprefix
       endif       

; make some QAplots of a random subset of objects
       if keyword_set(qaplot) then begin
          cat = read_irclusters(maggies=maggies,ivarmaggies=ivarmaggies)
          if keyword_set(use_zcluster) then begin
             outprefix = prefix+'_zclust'
             qaindx = shuffle_indx(n_elements(cat),num=30)
          endif else begin
             index = where((cat.z ge zminmax[0]) and (cat.z le zminmax[1]),nindex)
             qaindx = index[shuffle_indx(nindex,num=30)]
          endelse
          galaxy = 'Galaxy '+string(qaindx,format='(I5.5)')
          isedfit_qaplot, paramfile, result, iopath=isedpath, galaxy=galaxy, $
            index=qaindx, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
            outprefix=outprefix, /xlog
       endif
    endfor

return
end

;; --------------------------------------------------
;; do the fitting!  do not use the NEWFIRM photometry nor IRAC/ch3-4
;    if keyword_set(isedfit) then begin
;       toss = where(strmatch(param.filterlist,'*ch3*') or $
;         strmatch(param.filterlist,'*ch4*'))
;       ivarmaggies[toss,*] = 0.0
;;; fit at both the photometric redshift...
;;       index = where((cat.z ge min(param.redshift)) and $
;;         (cat.z le max(param.redshift)))
;;       isedfit, paramfile, maggies, ivarmaggies, cat.z, $
;;         iopath=iopath, clobber=clobber, index=index
;;; ...and at the cluster redshift...
;;;      res = mrdfits('BwRIJHKsirac_zclust_bc03_chab_calzetti_sfhgrid02.fits.gz',1)
;;;      index = (where((res.mass_err lt 0.01) and (res.zobj gt 1.2) and (res.zobj lt 1.3) ))[0:49]
;;       isedfit, paramfile, maggies, ivarmaggies, cat.zclust, $
;;         iopath=iopath, outprefix=param.prefix+'_zclust', $
;;         clobber=clobber;, index=index
;;; ...and then refit at the photometric redshift excluding IRAC 
;       index = where((cat.z ge min(param.redshift)) and $
;         (cat.z le max(param.redshift)))
;;      index = [200,239]
;       alltoss = where(strmatch(param.filterlist,'*spitzer*'))
;       ivarmaggies[alltoss,*] = 0.0
;       isedfit, paramfile, maggies, ivarmaggies, cat.z, $
;         iopath=iopath, outprefix=repstr(param.prefix,'irac',''), $
;         clobber=clobber, index=index
;    endif 

;; --------------------------------------------------
;    if keyword_set(qaplot) then begin
;       cat = read_irclusters(maggies=maggies,ivarmaggies=ivarmaggies)
;       ngal = n_elements(cat)
;       qaindx = long(randomu(seed,100)*ngal)
;       qaindx = qaindx[uniq(qaindx,sort(qaindx))]
;       galaxy = 'Galaxy '+string(qaindx,format='(I4.4)')
;       
;; includes irac, photoz
;       isedfit_qaplot, paramfile, iopath=iopath, index=qaindx, $
;         galaxy=galaxy, clobber=clobber, outprefix=outprefix
;; includes irac, zclust
;       isedfit_qaplot, paramfile, iopath=iopath, index=qaindx, $
;         galaxy=galaxy, clobber=clobber, outprefix=param.prefix+'_zclust'
;; no irac, photoz
;       isedfit_qaplot, paramfile, iopath=iopath, index=qaindx, $
;         galaxy=galaxy, clobber=clobber, outprefix=repstr(param.prefix,'irac','')
;    endif
