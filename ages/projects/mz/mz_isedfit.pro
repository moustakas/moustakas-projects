pro mz_isedfit, sdss=sdss, sfhgrid=sfhgrid, imf=imf, synthmodels=synthmodels, $
  redcurve=redcurve, models=models, isedfit=isedfit, clobber=clobber
; jm10jul29ucsd - derive stellar masses for the MZ AGES and SDSS
; samples 

    isedpath = mz_path(/ised)
    isedfit_sfhgrid_dir = mz_path(/montegrids)
    sfhgrid_paramfile = getenv('IDL_PROJECTS_DIR')+'/ages/projects/mz/mz_sfhgrid.par'
    
; defaults    
    if (n_elements(imf) eq 0) then imf = 'chab'
    if (n_elements(synthmodels) eq 0) then synthmodels = 'bc03'
    if (n_elements(redcurve) eq 0) then redcurve = 1 ; charlot

    mzgrid = read_sfhgrid_paramfile(sfhgrid,sfhgrid_paramfile=sfhgrid_paramfile)
    ngrid = n_elements(mzgrid)

    h100 = mz_h100(omega0=omega0,omegal=omegal)
    
; specify some parameters
    if keyword_set(sdss) then begin
       field = 'sdss'
       igm = '0'
       nzz = '40'
       zlog = '1'
       vagc = mz_get_vagc(zminmax=zminmax,sample=sample,$
         letter=letter,poststr=poststr)
    endif else begin
       field = 'ages'
       igm = '1'
       nzz = '50'
       zlog = '0'
       zbins = mz_zbins(zmin=zmin,zmax=zmax)
       zminmax = [zmin,zmax]
    endelse
    
    filters = mz_filterlist(sdss=sdss)
       
; loop on each grid
    for gg = 0, ngrid-1 do begin
       redcurvestring = redcurve2string(redcurve)
       sfhgridstring = string(mzgrid[gg].sfhgrid,format='(I2.2)')
;      sfhgridstring = strtrim(mzgrid[gg].sfhgrid,2)

; --------------------------------------------------
; always make the parameter file 
       paramfile = isedpath+'mz_'+field+'_sfhgrid'+sfhgridstring+'_isedfit.par'
       if im_file_test(paramfile,clobber=clobber) then return
       splog, 'Writing '+paramfile
       zrange = string(zminmax[0],format='(F4.2)')+','+string(zminmax[1],$
         format='(F4.2)')+','+nzz+','+zlog+' # [minz,maxz,dz,log?]'
       openw, lun, paramfile, /get_lun
       printf, lun, 'h100                 '+string(h100,format='(F4.2)')
       printf, lun, 'omega0               '+string(omega0,format='(F4.2)')
       printf, lun, 'omegal               '+string(omegal,format='(F4.2)')
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
       if keyword_set(models) then begin
          isedfit_models, paramfile, iopath=isedpath, clobber=clobber, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       endif

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
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, outprefix=outprefix, $
            galchunksize=500, index=index
       endif 
    endfor ; close SFHGRID loop
    
return
end
