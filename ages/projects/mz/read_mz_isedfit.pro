function read_mz_isedfit, sdss=sdss, sfhgrid=sfhgrid, imf=imf, $
  synthmodels=synthmodels, redcurve=redcurve, isedfile=isedfile, $
  post=post
; jm10aug23ucsd - read the desired output from ISEDFIT

    mzpath = ages_path(/projects)+'mz/isedfit/'
    if keyword_set(sdss) then field = 'sdss' else field = 'ages'
    
; read the parameter file but optionally overwrite SFHGRID    
    if (n_elements(sfhgrid) eq 0) then begin
       splog, 'SFHGRID input required'
       return, -1
    endif

    if (n_elements(imf) eq 0) then imf = 'chab'
    if (n_elements(synthmodels) eq 0) then synthmodels = 'bc03'
    if (n_elements(redcurve) eq 0) then redcurve = 1 ; charlot
    redcurvestring = redcurve2string(redcurve)

    isedfile = mzpath+field+'_'+synthmodels+'_'+imf+'_'+$
      redcurvestring+'_sfhgrid'+string(sfhgrid,format='(I2.2)')+'.fits.gz'
    isedfile_post = mzpath+field+'_'+synthmodels+'_'+imf+'_'+$
      redcurvestring+'_sfhgrid'+string(sfhgrid,format='(I2.2)')+'_post.fits.gz'

    if (file_test(isedfile) eq 0) then begin
       splog, 'iSEDFIT file '+isedfile+' not found!'
       return, -1
    endif
    splog, 'Reading '+isedfile
    ised = mrdfits(isedfile,1)

; also read the posterior distributions
    if arg_present(post) then begin
       if (file_test(isedfile_post) eq 0) then begin
          splog, 'iSEDFIT file '+isedfile_post+' not found!'
          return, -1
       endif
       splog, 'Reading '+isedfile_post
       post = mrdfits(isedfile_post,1)
    endif
    
return, ised
end
