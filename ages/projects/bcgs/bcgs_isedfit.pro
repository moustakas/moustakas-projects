pro bcgs_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, nzz=nzz, $
  filters=filters, igm=igm, super=super

    sfhgridstring = strtrim(super.sfhgrid,2)
    redcurvestring = redcurve2string(super.redcurve)
    synthmodels = strtrim(super.synthmodels,2)
    imf = strtrim(super.imf,2)
    
    splog, 'writing '+paramfile
    zrange = string(zminmax[0],format='(f4.2)')+','+string(zminmax[1],$
      format='(f4.2)')+','+nzz+' # [minz,maxz,dz]'
    openw, lun, paramfile, /get_lun
    printf, lun, 'synthmodels          '+synthmodels
    printf, lun, 'imf                  '+imf
    printf, lun, 'sfhgrid              '+sfhgridstring
    printf, lun, 'redcurve             '+redcurvestring
    printf, lun, 'prefix               '+prefix
    printf, lun, 'redshift             '+zrange
    printf, lun, 'igm                  '+igm+' # [0=no, 1=yes]'
    printf, lun, 'maxold               0 # [0=no, 1=yes]'
    printf, lun, 'filterlist           '+strjoin(filters,',')
    free_lun, lun 

return
end

pro bcgs_isedfit, sdss=sdss, supergrid=supergrid, models=models, $
  isedfit=isedfit, clobber=clobber
; jm11apr06ucsd -

; bcgs_isedfit, /model, /ised, /clob, supergrid=[1,2,3]    
; bcgs_isedfit, /model, /ised, /clob, supergrid=[4,5,6]
    
    isedpath = ages_path(/proj)+'bcgs/isedfit/'
    sfhgrid_basedir = ages_path(/projects)+'bcgs/isedfit/montegrids/'
    sfhgrid_paramfile = ages_path(/projects)+'bcgs/isedfit/bcgs_sfhgrid.par'

; read the supergrid parameter file    
    supergrid_paramfile = ages_path(/projects)+'bcgs/isedfit/bcgs_supergrid.par'
    super = yanny_readone(supergrid_paramfile)
    if (n_elements(supergrid) ne 0) then begin
       match2, super.supergrid, supergrid, m1, m2
       if (total(m2 eq -1) ne 0) then message, 'unknown supergrid!'
       match, super.supergrid, supergrid, m1, m2
       srt = sort(m2) & m1 = m1[srt] & m2 = m2[srt]
       super = super[m1]
    endif
    struct_print, super
    nsuper = n_elements(super)

; isedfit parameters
    igm = '0'
    if keyword_set(sdss) then begin
       nzz = '25'
       prefix = 'sdss_bcgs'
       zminmax = [0.02,0.26]
    endif else begin
       nzz = '90'
       prefix = 'bcgs'
       zminmax = [0.1,1.8]
    endelse

; loop on each supergrid
    for gg = 0, nsuper-1 do begin
       splog, 'working on grid '+strtrim(super[gg].supergrid,2)

       filters = bcgs_filterlist(sdss=sdss)
       nfilt = n_elements(filters)

       paramfile = isedpath+prefix+'_supergrid'+string(super[gg].supergrid,$
         format='(i2.2)')+'_isedfit.par'
       if im_file_test(paramfile,clobber=clobber) then return
       
; --------------------------------------------------
; build the models
       if keyword_set(models) then begin
          bcgs_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
            nzz=nzz, filters=filters, igm=igm, super=super[gg]
          isedfit_models, paramfile, iopath=isedpath, clobber=clobber, $
            sfhgrid_basedir=sfhgrid_basedir
       endif

; --------------------------------------------------
; do the fitting!
       if keyword_set(isedfit) then begin
          bcgs_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
            nzz=nzz, filters=filters, igm=igm, super=super[gg]
          maggies = read_bcgs_photometry(sdss=sdss,ivarmaggies=ivarmaggies, $
            zobj=zobj,galaxy=galaxy,filterlist=filters)
          isedfit, paramfile, maggies, ivarmaggies, zobj, iopath=isedpath, $
            clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
            sfhgrid_basedir=sfhgrid_basedir, debug=debug, outprefix=outprefix, $
            galchunksize=1000
       endif       

; --------------------------------------------------
; make some QAplots
       if keyword_set(qaplot) then begin
          isedfit_qaplot, paramfile, isedfit, iopath=iopath, galaxy=galaxy, $
            index=index, clobber=clobber, outprefix=outprefix
       endif
       
       if keyword_set(chi2grid_qaplot) then begin
          isedfit_chi2grid_qaplot, paramfile, isedfit, iopath=iopath, galaxy=galaxy, $
            index=index, clobber=clobber, outprefix=outprefix
       endif
       
; --------------------------------------------------
; measure rest-frame quantities
       if keyword_set(measure) then begin
          isedfit_measure, paramfile, measure, isedfit, iopath=iopath, $
         clobber=clobber, outprefix=outprefix
       endif
    endfor                      ; close SFHGRID loop

return
end
