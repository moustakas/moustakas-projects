pro mz_isedfit, sfhgrid=sfhgrid, imf=imf, synthmodels=synthmodels, $
  redcurve=redcurve, models=models, isedfit=isedfit, qaplot=qaplot, $
  measure=measure, clobber=clobber, sdss=sdss, debug=debug
; jm10jul29ucsd - derive stellar masses for the MZ AGES and SDSS
; samples 

    isedpath = ages_path(/projects)+'mz/isedfit/'
    sfhgrid_basedir = isedpath+'montegrids/'
    
; defaults    
    if (n_elements(imf) eq 0) then imf = 'chab'
    if (n_elements(synthmodels) eq 0) then synthmodels = 'bc03'
    if (n_elements(redcurve) eq 0) then redcurve = 1 ; charlot

    sfhgrid_paramfile = mz_sfhgrid_parfile()
    mzgrid = read_sfhgrid_paramfile(sfhgrid,sfhgrid_paramfile=sfhgrid_paramfile)
    ngrid = n_elements(mzgrid)

; specify some parameters
    if keyword_set(sdss) then begin
       field = 'sdss'
       igm = '0'
       nzz = '25'
       vagc = mz_get_vagc(zminmax=zminmax,sample=sample,$
         letter=letter,poststr=poststr)
    endif else begin
       field = 'ages'
       igm = '0'
       nzz = '60'
       zbins = mz_zbins(zmin=zmin,zmax=zmax)
       zminmax = [zmin,zmax]
    endelse
    
    filters = mz_filterlist(sdss=sdss)
       
; loop on each grid
    for gg = 0, ngrid-1 do begin
       redcurvestring = redcurve2string(redcurve,params=mzgrid[gg])
       sfhgridstring = strtrim(mzgrid[gg].sfhgrid,2)

; --------------------------------------------------
; always make the parameter file 
       paramfile = isedpath+'mz_'+field+'_isedfit.par'
       if im_file_test(paramfile,clobber=clobber) then return
       splog, 'Writing '+paramfile
       zrange = string(zminmax[0],format='(F4.2)')+','+string(zminmax[1],$
         format='(F4.2)')+','+nzz+' # [minz,maxz,dz]'
       openw, lun, paramfile, /get_lun
       printf, lun, 'synthmodels          '+synthmodels
       printf, lun, 'imf                  '+imf
       printf, lun, 'sfhgrid              '+sfhgridstring
       printf, lun, 'redcurve             '+redcurvestring
       printf, lun, 'prefix               '+field
       printf, lun, 'redshift             '+zrange
       printf, lun, 'igm                  '+igm+' # [0=no, 1=yes]'
       printf, lun, 'maxold               0 # [0=no, 1=yes]'
       printf, lun, 'filterlist           '+strjoin(filters,',')
       free_lun, lun 

; --------------------------------------------------
; build the models
       if keyword_set(models) then isedfit_models, paramfile, $
         iopath=isedpath, sfhgrid_basedir=sfhgrid_basedir, $
         clobber=clobber

; --------------------------------------------------
; do the fitting!  
       if keyword_set(isedfit) then begin
          if keyword_set(sdss) then begin
             post = read_vagc_garching(sample=sample,$
               letter=letter,poststr=poststr,/postlss)
             maggies = mz_get_maggies(post,/sdss,ivarmaggies=ivarmaggies)
             zobj = post.z
          endif else begin
             phot = read_ages(/photo)
             index = where((phot.imain eq 1) and (phot.z ge zmin) and (phot.z le zmax))
             maggies = mz_get_maggies(phot,ivarmaggies=ivarmaggies)
             zobj = phot.z
          endelse

          isedfit, paramfile, maggies, ivarmaggies, zobj, iopath=isedpath, $
            clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
            sfhgrid_basedir=sfhgrid_basedir, debug=debug, outprefix=outprefix, $
            galchunksize=750, index=index
       endif 
    endfor ; close SFHGRID loop
    
return
end
