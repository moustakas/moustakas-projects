function read_tkrs_rotcurve, rotfile
; jm13jul03siena - read a TKRS-style rotation curve file

    pixscale = 0.119 ; [arcsec/pixel] (see Weiner+06)
    
; Files for each individual galaxy follow the naming:
; <tkrs id>__<rest wavelength>.fitdata, etc.  The rest wavelength
; is for the line chosen for ROTCURVE analysis.
; 
; Formats: 
; 
; 0000428__5007.fitdata - output from ROTCURVE fitting.  Fields are
;   ID
;   intensity along slit
;   sigma along slit (spatial), arcsec
;   sigma along slit (spatial) corrected for seeing, arcsec
;   rot curve systemic velocity offset, km/sec
;   rot curve scale radius (may be fixed value), arcsec
;   error on rot curve scale radius (may be fixed)
;   rot curve fit velocity Vrot, km/sec
;   error on Vrot (not very statistically well-defined)
;   rot curve fit dispersion sigma_2d, km/sec
;   error on sigma_2d (not very statistically well-defined)
;   number of good velocity data points used in the fit
;   chi-squared of the fit
; 
; 0000428__5007.rotcurv - data and models generated by ROTCURVE as a
; function of row in the 2-d CCD data.  Columns are:
;   number of pixel row
;   y = displacement along slit, arcsec
;   total flux along row (incl continuum)
;   intensity of emission line fit
;   velocity of em line
;   error on velocity
;   dispersion of em line 
;   error on dispersion
;   rot curve unblurred model intensity
;   rot curve unblurred model velocity
;   rot curve unblurred model dispersion
;   rot curve blurred model intensity 
;   rot curve blurred model velocity
;   rot curve blurred model dispersion
;   ifgood = did the program decide that this was a good datapoint to use
;     in fitting? (1-yes, 0-no)
; 
; list.all.id.rcnotes - quality and comments on the rotation curve fits
; Qualities 3 and 4 are fairly reliable.  Qualities 1 and 2 are not
; reliable and the ROTCURVE output values must not be trusted!!!
; columns are:
;   ID
;   RC quality code
;   RC problem code
;   comments
; 
;   quality codes:
;   4 - good
;   3 - ok, minor problems
;   2 - semi dubious
;   1 - bad
;   problem codes:
;   0 - none
;   1 - poor fit
;   2 - bad outlier data points included
;   3 - good data points excluded
;   4 - bad center
;   5 - small spatial extent
;   6 - bad data
;   7 - oversmoothing/small intensity profile width
; 
; list.all.id.rcqual.yrange - summary of RC qualities and extent of 
; data along slit (y-range).  columns:
;   ID
;   RC quality code
;   RC problem code
;   ID
;   minimum row number with good data
;   maximum row number with good data
;   extent of good data in rows
;   extent of good data in arcsec
; 
; list.goodrc.id - list of IDs with decent quality rotation curve fits.
; This is everything with a quality code 3 or 4.  
; *** DO NOT USE *** the fit values from objects with qualities 1 and 2.
; 
; Readme.plotprof.macro - supermongo macros for plotting the .rotcurv files
; 
; sm.plotall.in2 - redirect this file into sm to make plots of
; rotation and dispersion profile for each galaxy.
; 
; rcfit.*.ps - plots produced by supermongo.
; 
; # make archive
; tar zchvf ktrs_rc_data.tgz 00*.fitdata 00*.rotcurv Readme.rc_data \
;   list.all.id.rcnotes list.all.id.rcqual.yrange list.goodrc.id \
;   Readme.plotprof.macro sm.plotall.in2 rcfit.00*.ps
    
    readcol, rotfile, row, yoff, totflux, flux, vel, $
      errvel, sigma, errsigma, model_unblur_flux, $
      model_unblur_vel, model_unblur_sigma, $
      model_blur_flux, model_blur_vel, model_blur_sigma, $
      good, /silent, format='I,F,F,F,F,F,F,F,F,F,F,F,F,F,I'
    nrow = n_elements(row)
    curve = replicate({row: 0, $
      yoff: 0.0, $
      totflux: 0.0,$
      flux: 0.0, $
      vel: 0.0, $
      errvel: 0.0, $
      sigma: 0.0, $
      errsigma: 0.0, $
      model_unblur_flux: 0.0, $
      model_unblur_vel: 0.0, $
      model_unblur_sigma: 0.0, $
      model_blur_flux: 0.0, $
      model_blur_vel: 0.0, $
      model_blur_sigma: 0.0, $
      good: 0},nrow)

    curve.row = row
    curve.yoff = yoff
    curve.totflux = totflux
    curve.flux = flux
    curve.vel = vel
    curve.errvel = errvel
    curve.sigma = sigma
    curve.errsigma = errsigma
    curve.model_unblur_flux = model_unblur_flux
    curve.model_unblur_vel = model_unblur_vel
    curve.model_unblur_sigma = model_unblur_sigma
    curve.model_blur_flux = model_blur_flux
    curve.model_blur_vel = model_blur_vel
    curve.model_blur_sigma = model_blur_sigma
    curve.good = good
return, curve
end
