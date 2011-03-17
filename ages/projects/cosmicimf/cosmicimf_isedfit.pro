pro cosmicimf_isedfit, imf, make_parfiles=make_parfiles, models=models, $
  isedfit=isedfit, qaplot=qaplot, measure=measure, clobber=clobber, $
  debug=debug
; jm10mar19ucsd - measure stellar masses for the AGES/COSMICIMF
; project 

    rootpath = ages_path(/projects)+'cosmicimf/isedfit/'

; initialize the IMF parameters
    alpha2 = [2.35,cosmicimf_slope()]
    nimf = n_elements(alpha2)
    imf = ['Salpeter','cosmicimf_'+string(alpha2[1:nimf-1],format='(F4.2)')]
    imfstr = repstr(repstr(imf,'Salpeter','salp'),'Kroupa02','kroupa')

; --------------------------------------------------
; build the parameter files
    if keyword_set(make_parfiles) then begin
       for ii = 0, nimf-1 do begin
          imfpath = rootpath+imf[ii]+'/'
          case imf[ii] of 
             'Salpeter': imfstr = 'salp'
             'Kroupa02': imfstr = 'kroupa'
             else: if strmatch(imf[ii],'*cosmicimf*') then $
               imfstr = imf[ii] else message, 'Fix me'
          endcase
          if (file_test(rootpath+imf[ii]+'/',/dir) eq 0) then $
            spawn, 'mkdir -p '+rootpath+imf[ii]
          paramfile = rootpath+imf[ii]+'/'+imf[ii]+'_isedfit.par'
          openw, lun, paramfile, /get_lun
          printf, lun, '# Input parameters for ISEDFIT_MODELS'
          printf, lun, 'synthmodels          pegase'
          printf, lun, 'imf                  '+imfstr
          printf, lun, 'sfhgrid              2' ; 3
          printf, lun, 'redcurve             calzetti' ; none
          printf, lun, 'prefix               BwRIJHKsirac'
          printf, lun, 'redshift             0.01,0.8,0.02 # [minz,maxz,dz]'
          printf, lun, 'igm                  1 # [0=no, 1=yes]'
          printf, lun, 'filterlist           galex_FUV.par,galex_NUV.par,ndwfs_Bw.par,ndwfs_R.par,ndwfs_I.par,newfirm_J.par,newfirm_H.par,newfirm_Ks.par,spitzer_irac_ch1.par,spitzer_irac_ch2.par,spitzer_irac_ch3.par,spitzer_irac_ch4.par'
          free_lun, lun
       endfor
    endif

; build the SFH grids    
;   for ii = 1, 15 do build_isedfit_sfhgrid, synthmodels='pegase', imf=imf[ii], sfhgrid=3, /clobber
;   for ii = 16, 16 do build_isedfit_sfhgrid, synthmodels='pegase', imf=imf[ii], sfhgrid=3, /clobber
;   for ii = 0, nimf-1 do build_isedfit_sfhgrid, synthmodels='pegase', imf=imfstr[ii], sfhgrid=2, /clobber

; loop on each IMF
;   for ii = 0, nimf/2-1 do begin
;   for ii = nimf/2, nimf-1 do begin
    for ii = 0, nimf-1 do begin
       iopath = rootpath+imf[ii]+'/'
       paramfile = imf[ii]+'_isedfit.par'

; build the models
       if keyword_set(models) then isedfit_models, $
         paramfile, iopath=iopath, clobber=clobber

; read the AGES photometry
       maggies = read_ages_photometry(ivarmaggies=ivarmaggies,$
         index=index,filterlist=filterlist,redshift=redshift)

; do the fitting!
       if keyword_set(isedfit) then begin
          isedfit, paramfile, maggies, ivarmaggies, redshift, result, $
            iopath=iopath, outprefix=outprefix, nminphot=nminphot, $
            clobber=clobber, debug=debug, index=index
       endif 

; make a QAplot
       if keyword_set(qaplot) then begin
;         index = where(total 
;         index = index[0:10]
          isedfit_qaplot, paramfile, isedfit, iopath=iopath, galaxy=galaxy, $
            index=index, clobber=clobber, outprefix=outprefix
stop
       endif

; measure rest-frame quantities
       if keyword_set(measure) then begin
          isedfit_measure, paramfile, measure, isedfit, iopath=iopath, $
            clobber=clobber, outprefix=outprefix
       endif
    endfor
       
return
end
