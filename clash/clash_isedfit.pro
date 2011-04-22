pro hizea_isedfit, models=models, write_paramfile=write_paramfile, $
  isedfit=isedfit, qaplot=qaplot, clobber=clobber, debug=debug, $
  noirac=noirac
; jm10dec20ucsd - derive stellar masses for the HIZEA sample
; jm11apr06ucsd - major updates

    iopath = hizea_path(/isedfit)
    paramfile = iopath+'hizea_isedfit.par'

    sfhgrid_basedir = hizea_path(/monte)
    sfhgrid_paramfile = hizea_path(/mass)+'hizea_sfhgrid.par'

; --------------------------------------------------
; write the parameter file
    if keyword_set(write_paramfile) then begin
       splog, 'Writing '+paramfile
       openw, lun, paramfile, /get_lun
       printf, lun, 'synthmodels          bc03'
       printf, lun, 'imf                  chab'
       printf, lun, 'sfhgrid              1'
       printf, lun, 'redcurve             charlot'
       printf, lun, 'prefix               hizea'
       printf, lun, 'redshift             0.35,0.91,30'
       printf, lun, 'igm                  1'
       printf, lun, 'maxold               0 # [0=no, 1=yes]'
       printf, lun, 'filterlist           '+strjoin(hizea_filterlist(),',')
       free_lun, lun 
    endif

; --------------------------------------------------
; build the models
    if keyword_set(models) then isedfit_models, paramfile, $
      iopath=iopath, sfhgrid_basedir=sfhgrid_basedir, clobber=clobber

; --------------------------------------------------
; do the fitting!  
    if keyword_set(isedfit) then begin
       phot = mrdfits(hizea_path(/sdss)+'hizea_galex_sdss_spitzer.fits.gz',1)
       maggies = phot.maggies
       ivarmaggies = phot.ivarmaggies
       zobj = phot.z
       filt = hizea_filterlist()

       if keyword_set(noirac) then begin
          toss = where(strmatch(filt,'*irac*'))
          ivarmaggies[toss,*] = 0.0
          outprefix = 'hizea_noirac'
       endif
       isedfit, paramfile, maggies, ivarmaggies, zobj, result, $
         iopath=iopath, outprefix=outprefix, galchunksize=galchunk, $
         sfhgrid_paramfile=sfhgrid_paramfile, sfhgrid_basedir=sfhgrid_basedir, $
         clobber=clobber, debug=debug, index=index
    endif 

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
       if keyword_set(noirac) then outprefix = 'hizea_noirac'

       sample = mrdfits(hizea_path(/sdss)+'hizea_galex_sdss_spitzer.fits.gz',1)
;      index = where(sample.maggies[8] gt 0) ; has ch1 flux
;      galaxy = sample[index].galaxy

       keck = rsex(hizea_path(/sdss)+'keckao_sample.sex')

       ra = 15D*im_hms2dec(keck.ra) 
       dec = im_hms2dec(keck.dec)
       spherematch, sample.ra, sample.dec, ra, dec, 3.0/3600.0, m1, m2
       srt = sort(m2) & m1 = m1[srt] & m2 = m2[srt]
       index = m1
       galaxy = sample[index].galaxy
       
       isedfit_qaplot, paramfile, isedfit, iopath=iopath, galaxy=galaxy, $
         index=index, sfhgrid_basedir=sfhgrid_basedir, clobber=clobber, $
         outprefix=outprefix
    endif

return
end
