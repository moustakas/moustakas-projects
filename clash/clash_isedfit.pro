pro clash_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, nzz=nzz, $
  filters=filters, igm=igm, super=super

    sfhgridstring = strtrim(super.sfhgrid,2)
    redcurvestring = redcurve2string(super.redcurve)
    synthmodels = strtrim(super.synthmodels,2)
    imf = strtrim(super.imf,2)
    
    splog, 'Writing '+paramfile
    zrange = string(zminmax[0],format='(f4.2)')+','+string(zminmax[1],$
      format='(f4.2)')+','+nzz+' # [minz,maxz,nz]'
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

pro clash_isedfit, supergrid=supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber
; jm11may05ucsd -

; clash_isedfit, /model, /ised, /clob, supergrid=[1,2,3]    
    
    isedpath = clash_path(/ised)
    catpath = clash_path(/cat)
    isedfit_sfhgrid_dir = clash_path(/monte)
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/clash_sfhgrid.par'

; read the supergrid parameter file    
    supergrid_paramfile = getenv('CLASH_DIR')+'/clash_supergrid.par'
    super = yanny_readone(supergrid_paramfile)
    if (n_elements(supergrid) ne 0) then begin
       match2, super.supergrid, supergrid, m1, m2
       if (total(m2 eq -1) ne 0) then message, 'Unknown supergrid!'
       match, super.supergrid, supergrid, m1, m2
       srt = sort(m2) & m1 = m1[srt] & m2 = m2[srt]
       super = super[m1]
    endif
    struct_print, super
    nsuper = n_elements(super)

; isedfit parameters
    igm = '1'
    nzz = '3'
    prefix = 'a383'
    zminmax = [1.0,1.1]

; loop on each supergrid
    for gg = 0, nsuper-1 do begin
       splog, 'working on grid '+strtrim(super[gg].supergrid,2)

       filters = clash_filterlist()
       nfilt = n_elements(filters)

       paramfile = isedpath+prefix+'_supergrid'+string(super[gg].supergrid,$
         format='(i2.2)')+'_isedfit.par'
       if im_file_test(paramfile,clobber=clobber) then return
       
; --------------------------------------------------
; build the models
       if keyword_set(models) then begin
          clash_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
            nzz=nzz, filters=filters, igm=igm, super=super[gg]
          isedfit_models, paramfile, iopath=isedpath, clobber=clobber, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       endif

; --------------------------------------------------
; do the fitting!
       if keyword_set(isedfit) then begin
          clash_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
            nzz=nzz, filters=filters, igm=igm, super=super[gg]
; read the catalog, and fit both with and without the first two bands
          cat = mrdfits(catpath+'a383.fits.gz',1)
          cat = [cat,cat]
          clash_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist
          ivarmaggies[0:1,1] = 0 ; do not use!

          isedfit, paramfile, maggies, ivarmaggies, cat.z, iopath=isedpath, $
            clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, debug=debug, outprefix=outprefix
       endif       

; --------------------------------------------------
; make some QAplots
       if keyword_set(qaplot) then begin
          cat = mrdfits(catpath+'a383.fits.gz',1)
          cat = [cat,cat]
          isedfit_qaplot, paramfile, result, iopath=iopath, galaxy=cat.galaxy, $
            index=index, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
            outprefix=outprefix
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
