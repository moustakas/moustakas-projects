function read_bcgs_isedfit, sdss=sdss, supergrid=supergrid, rows=rows, $
  post=post
; jm11apr07ucsd - read the desired output from ISEDFIT

    bcgspath = ages_path(/projects)+'bcgs/isedfit/'

    if keyword_set(sdss) then begin
       prefix = 'sdss_bcgs'
    endif else begin
       prefix = 'bcgs'
    endelse
    
    supergrid_paramfile = ages_path(/projects)+'bcgs/isedfit/bcgs_supergrid.par'
    super = yanny_readone(supergrid_paramfile)
    if (n_elements(supergrid) ne 0) then begin
       this = where(super.supergrid eq supergrid)
       if (this[0] eq -1) then message, 'Unknown supergrid!'
       super = super[this]
    endif else begin
       splog, 'SUPERGRID input required'
       return, -1
    endelse
    nsuper = n_elements(super)

; build the filename and then read
    superstring = string(super.supergrid,format='(I2.2)')
    sfhgridstring = string(super.sfhgrid,format='(I2.2)')
    redcurvestring = redcurve2string(super.redcurve)
    synthmodels = strtrim(super.synthmodels,2)
    imf = strtrim(super.imf,2)
    
    isedfile = bcgspath+prefix+'_'+synthmodels+'_'+imf+'_'+$
      redcurvestring+'_sfhgrid'+sfhgridstring+'.fits.gz'
    isedfile_post = bcgspath+prefix+'_'+synthmodels+'_'+imf+'_'+$
      redcurvestring+'_sfhgrid'+sfhgridstring+'_post.fits.gz'

    if (file_test(isedfile) eq 0) then begin
       splog, 'iSEDFIT file '+isedfile+' not found!'
       return, -1
    endif
    splog, 'Reading '+isedfile
    ised = mrdfits(isedfile,1,rows=rows)

; optionally read the posterior distribution (just the rows we care
; about) 
    if arg_present(post) then begin
       if (file_test(isedfile_post) eq 0) then begin
          splog, 'iSEDFIT file '+isedfile_post+' not found!'
          return, -1
       endif
       splog, 'Reading '+isedfile_post
       post = mrdfits(isedfile_post,1,rows=ised.isedfit_id)
    endif
    
return, ised
end
