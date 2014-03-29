;+
; NAME:
;   BUILD_DESI_DEEP2_TEMPLATES
;
; PURPOSE:
;   Generate emission-line galaxy spectra for DESI using the DEEP2
;   spectroscopy and photometry.
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   minwave - minimum output wavelength (default 3000 A)
;   maxwave - approximate output wavelength (default 10,000 A) 
;   velpixsize - spectrum pixel size [default 20 km/s]
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   I have verified that my BUILD_EMLINE script yields the same
;   emission-line spectrum as GANDALF, which was used to fit the DEEP2
;   spectra. 
; 
;   All line fluxes are defined with respect to [OII] 3729.
; 
; ToDo:
;   1. Add more nebular emission lines.
;   2. Allow the nebular emission-line ratios to vary.
;   3. Build stacks of spectra.
;   4. Add some jitter to the emission-line redshifts compared to the
;      continuum redshifts.  Could also vary the line-widths somewhat.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2014 Jan 09, Siena
;   jm14mar11siena - major update
;
; Copyright (C) 2014, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function build_emline, emspectrum, logwave=logwave, lineflux=lineflux, $
  linesigma=linesigma, linewave=linewave, zshift=zshift
; internal support routine: build a single emission line
;   logwave - desired (rest-frame) wavelength vector [A]
;   lineflux - integrated (rest-frame) line-flux [erg/s/cm2]
;   linesigma - kinematic velocity width [km/s]
;   linewave - rest-frame emission-line wavelength [A]
;   zshift - radial velocity shift relative to the systemic redshift 

; allow the line to have a radial velocity shift with respect to the
; systemic redshift/velocity
    if n_elements(zshift) eq 0 then zshift = 0.0 

    light = 2.99792458D5
    sigma = linesigma/light/alog(10)       ; line-width [log-10 Angstrom]
    amplitude = lineflux/alog(10)/linewave ; line-amplitude [erg/s/cm2/A]
    linewave1 = alog10(linewave*(1+zshift))
    emspectrum += amplitude*exp(-0.5*(logwave-linewave1)^2/$ ; [erg/s/cm2/A, rest]
      sigma^2)/(sqrt(2.0*!pi)*sigma)

return, emspectrum
end

pro build_desi_deep2_templates, minwave=minwave, maxwave=maxwave, $
  velpixsize=velpixsize, debug=debug
; jm13dec18siena - take the output of BUILD_DESI_DEEP2_TEMPLATE_SAMPLE
;   and generate the full-resolution spectra

    version = 'v1.0'
    templatepath = getenv('IM_PROJECTS_DIR')+'/desi/templates/'
    
    light = im_light(/kms)
    light_ang = im_light(/Ang)

; simulation parameter defaults
    if n_elements(minwave) eq 0 then minwave = 3000D  ; [A]
    if n_elements(maxwave) eq 0 then maxwave = 1D4  ; [A]
    if n_elements(velpixsize) eq 0 then velpixsize = 20D ; [km/s]

    crval1 = alog10(minwave)           ; starting wavelength [log-10 A]
    cdelt1 = velpixsize/light/alog(10) ; [pixel size in log-10 A]
    npix = round((alog10(maxwave)-crval1)/cdelt1+1)

;   crval1 = alog(minwave)      ; starting wavelength [log-10 A]
;   cdelt1 = velpixsize/light   ; [pixel size in log-10 A]
;   npix = round((alog(maxwave)-crval1)/cdelt1+1)

    obswave = crval1+dindgen(npix)*cdelt1 ; log-10 spacing
    obswave_edges = k_lambda_to_edges(10D^obswave)

    weff = k_lambda_eff(filterlist=deep2_filterlist())
    nband = n_elements(weff)

; specify the lines (other than [OII]) we care about; adjust the [OII]
; doublet ratio according to some prior; this will need to go inside
; the FOR loop when I do this with Monte Carlo
    linewave = [3726.032D,4861.33D,4958.91D,5006.84D,6562.8D]
    lineratio = [0.73,0.4577,0.094,0.2821,2.002]
    nline = n_elements(linewave)
    
; read the fitted sample and select a fiducial subset 
    info = mrdfits(templatepath+'desi_deep2_template_sample.fits.gz',1)
;   ngal = n_elements(info)

    keep = where(info.oii_3727[0] gt 8D-17 and info.ugriz[2] gt 23.5,ngal)
    outinfo = struct_addtags(struct_trimtags(info[keep],select=['OBJNO','RA','DEC',$
      'Z','SIGMA_KMS']),replicate({oii_3726: 0.0, $
      oii_3729: 0.0, oii: 0.0, oii_ew: 0.0, $ ; lineratio: fltarr(nline), $
      mag_u: 0.0, mag_g: 0.0, mag_r: 0.0, mag_i: 0.0, mag_z: 0.0, $
;     ugriz: fltarr(5), 
      logmstar: 0.0, logsfr: 0.0},ngal))
    outinfo.oii_3726 = info[keep].oii_3727_1[0]
    outinfo.oii_3729 = info[keep].oii_3727_2[0]
    outinfo.oii = info[keep].oii_3727[0]
    outinfo.oii_ew = info[keep].oii_3727_ew[0]
    
; divide the sample into chunks
;   chunksize = 5L
;   chunksize = 500L
    chunksize = ngal
    nchunk = ceil(ngal/float(chunksize))
    
;   for ichunk = 0, 0 do begin
    for ichunk = 0, nchunk-1 do begin
       i1 = ichunk*chunksize
       i2 = ((ichunk*chunksize+chunksize)<ngal)-1L
       nthese = i2-i1+1L
       these = lindgen(nthese)+i1

       ised = read_isedfit(templatepath+'desi_deep2_paramfile.par',/getmodels,$
         isedfit_dir=templatepath,montegrids_dir=templatepath+'montegrids/',$
         index=info[these].isedfit_id,/flam)

; initialize the output stack of spectra
       outflux = fltarr(nthese,npix)
       
       for igal = 0, nthese-1 do begin
; restore the best-fit iSEDfit model and resample to constant
; log10-lambda; note that we're ignoring the change in spectral 
; resolution of the underlying continuum spectrum 
          ised1 = ised[igal]
          restwave = obswave-alog10(1+info[these[igal]].z)

; observed [erg/s/cm2/A]                    
          continuum = im_log_rebin(im_double(ised1.wave),ised1.flux,/log10,$  
            minwave=minwave,maxwave=maxwave,vsc=velpixsize);,outwave=obswave1)
          if n_elements(continuum) ne npix then message, 'Problem here!'
          
;; get the *rest-frame* continuum flux around [OII]
;          cindx = where(obswave gt 3727.0*(1+info[these[igal]].z)*0.995 and $
;            obswave lt 3727.0*(1+info[these[igal]].z)*1.005)
;          cflux = djs_median(continuum[cindx]*(1+info[these[igal]].z))
;;         djs_plot, wave[these], cflux[these]

; rebuild the measured emission-line spectrum at the new velocity
; resolution, all relative to [OII] 3729; this step assumes that our 
; emission-line fitting is spot-on
          linesigma = info[these[igal]].sigma_kms
          reflineflux = info[these[igal]].oii_3727_2[0] ; [erg/s/cm2]
          emspectrum = build_emline(continuum*0,logwave=restwave,zshift=zshift,$
            lineflux=reflineflux,linesigma=linesigma,linewave=3728.8149D)
       
;         lineflux = reflineflux*ispec[these[igal]].oii_3727_1[0]/$
;           ispec[these[igal]].oii_3727_2[0]
;         emspectrum = build_emline(emspectrum,logwave=restwave,zshift=zshift,$
;           lineflux=lineflux,linesigma=linesigma,$
;           linewave=ispec[these[igal]].oii_3727_1_wave)

;         djs_plot, ised1.wave, ised1.flux, xr=[6000,8000]              
;         djs_oplot, [1,1]*3727*(1+info[igal].z), !y.crange, color='red'
;         djs_oplot, !x.crange, info[igal].cflux_3727*[1,1], color='red'
          
; now do the other nebular lines
          for ll = 0, nline-1 do emspectrum = build_emline(emspectrum,$
            logwave=restwave,zshift=zshift,lineflux=lineratio[ll]*reflineflux,$
            linesigma=linesigma,linewave=linewave[ll])

; shift to the observed frame and build the final spectrum
          emspectrum = emspectrum/(1+info[these[igal]].z)
          flux = continuum + emspectrum ; [erg/s/cm2/A]
          outflux[igal,*] = flux

; synthesize photometry and pack in the output structure
;         ugrizch12 = k_project_filters(obswave_edges,flux,$
;           filterlist=[sdss_filterlist(),wise_filterlist(/short)])
          ugriz = -2.5*alog10(reform(k_project_filters(obswave_edges,flux,$
            filterlist=[sdss_filterlist()])))
          
;         outinfo[these[igal]].lineratio = lineratio
          outinfo[these[igal]].logmstar = ised[igal].mstar_50
          outinfo[these[igal]].logsfr = ised[igal].sfr_50
;         outinfo[these[igal]].ugriz = ugriz
          outinfo[these[igal]].mag_u = reform(ugriz[0,*])
          outinfo[these[igal]].mag_g = reform(ugriz[1,*])
          outinfo[these[igal]].mag_r = reform(ugriz[2,*])
          outinfo[these[igal]].mag_i = reform(ugriz[3,*])
          outinfo[these[igal]].mag_z = reform(ugriz[4,*])
          
          if keyword_set(debug) then begin
             mag = maggies2mag(info[these[igal]].maggies,magerr=magerr,$
               ivarmaggies=info[these[igal]].ivarmaggies)
             abflux = -2.5*alog10((flux*(10D^obswave)^2.0/light_ang)>1D-50)-48.6
             yrange = [max(abflux),min(abflux)]
             djs_plot, 10D^obswave/1D4, abflux, psym=10, xsty=1, $
               ysty=3, xrange=xrange, yrange=yrange, $
               xtitle='Wavelength (\mu'+'m)', ytitle='Flux'
             oploterror, weff/1D4, mag, magerr, psym=6, symsize=2.0, $
               color=cgcolor('green'), errcolor=cgcolor('green')
             cc = get_kbrd(1)
          endif
       endfor ; close galaxy loop

       if nchunk eq 1 then begin
          outfile = templatepath+'desi_deep2_templates_'+version+'.fits'
       endif else begin
          outfile = templatepath+'desi_deep2_templates_'+version+'_chunk'+$
            string(ichunk+1,format='(I2.2)')+'.fits'
       endelse

; build the spectrum header
       mkhdr, fluxhdr, outflux, /extend
       sxdelpar, fluxhdr, 'DATE'
       sxdelpar, fluxhdr, 'COMMENT'
       sxaddpar, fluxhdr, 'OBJTYPE', 'ELG-DEEP2', 'object type'
       sxaddpar, fluxhdr, 'FLUXUNIT', 'erg/s/cm2/A', ' spectrum flux units'
       sxaddpar, fluxhdr, 'CRVAL1', crval1, ' reference log10(Angstrom)'
       sxaddpar, fluxhdr, 'CDELT1', cdelt1, ' delta log10(Angstrom)'
       sxaddpar, fluxhdr, 'LOGLAM', 1, ' log10 spaced wavelengths'
       sxaddpar, fluxhdr, 'VELSCALE', velpixsize, ' pixel size in km/s'
       
; we need to do some MWRFITS jujitsu to get the metadata header right 
       mwrfits, 0, outfile, /create
       mwrfits, outinfo, outfile, /silent

       nheader = [$
;        'DESI/DEEP2 galaxy templates',$
;        'Generated by J. Moustakas (jmoustakas@siena.edu) on '+im_today(),$
         'EXTNAME'+" = '"+string('METADATA','(A-8)')+"'           /"]

       units = [$
         'OBJNO,,DEEP2/DR4 unique object number',$
         'RA,degrees,right ascension (J2000)',$
         'DEC,degrees,declination (J2000)',$
         'Z,,DEEP2 heliocentric redshift',$
         'SIGMA_KMS,km/s,emission line velocity width',$
         'OII_3726,erg/s/cm2,[OII]\lambda3726 emission-line flux',$
         'OII_3729,erg/s/cm2,[OII]\lambda3729 emission-line flux',$
         'OII,erg/s/cm2,total [OII]\lambda\lambda3726,29 emission-line flux',$
         'OII_EW,Angstrom,rest-frame [OII]\lambda\lambda3726,29 emission-line equivalent width',$
;        'LINERATIO,,adopted nebular emission-line ratios relative to [OII] \lambda3729',$
;        'UGRIZ,mag,synthesized SDSS ugriz AB photometry',$
         'MAG_U,mag,synthesized SDSS u-band AB mag',$
         'MAG_G,mag,synthesized SDSS g-band AB mag',$
         'MAG_R,mag,synthesized SDSS r-band AB mag',$
         'MAG_I,mag,synthesized SDSS i-band AB mag',$
         'MAG_Z,mag,synthesized SDSS z-band AB mag',$
         'LOGMSTAR,Msun,log10(stellar mass) (Chabrier, h=0.7)',$
         'LOGSFR,Msun/yr,log10(star formation rate) (Chabrier, h=0.7)']
       metahdr =  im_update_header(outfile,units,nheader)

; finally write everything out       
       im_mwrfits, transpose(outflux), outfile, fluxhdr, /clobber, /nogzip
       im_mwrfits, outinfo, outfile, metahdr, /append, /gzip
    endfor                      ; close chunk loop

return
end

;; read the observed DEEP2 spectrum       
;         ss = mrdfits('testspec.fits.gz',1)
;         good = where(ss.wave gt 0)
;;        datawave = exp(ss.wave[good])*(1+ispec[these[igal]].zline_forbidden)
;         datawave = exp(ss.wave[good])*(1+ss.z)
;         dataflux = ss.flux[good]/(1+ss.z)
;         dataflux = ss.flux[good]/ss.continuum[good]*$
;           interpolate(cflux,findex(wave,datawave))
