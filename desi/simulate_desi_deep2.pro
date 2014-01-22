;+
; NAME:
;   SIMULATE_DESI_DEEP2
;
; PURPOSE:
;   Generate emission-line galaxy spectra for DESI using DEEP2.
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
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2014 Jan 09, Siena
;
; Copyright (C) 2013, John Moustakas
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

pro simulate_desi_deep2, minwave=minwave, maxwave=maxwave, $
  velpixsize=velpixsize, debug=debug
; jm13dec18siena - simulate DESI spectra using DEEP2

    desipath = getenv('IM_PROJECTS_DIR')+'/desi/deep2/'

    light = im_light(/kms)
    light_ang = im_light(/Ang)

; simulation parameter defaults
    if n_elements(minwave) eq 0 then minwave = 3000D  ; [A]
    if n_elements(maxwave) eq 0 then maxwave = 1D4  ; [A]
    if n_elements(velpixsize) eq 0 then velpixsize = 20D ; [km/s]

    coeff0 = alog10(minwave)             ; starting wavelength [log-10 A]
    coeff1 = alog10(velpixsize/light+1D) ; [pixel size in log-10 A]
;   npix = round((alog10(maxwave)-coeff0)/coeff1)+1
;   outwave = coeff0+dindgen(npix)*coeff1

; specify the lines (other than [OII]) we care about 
    linewave = [4861.33D,4958.91D,5006.84D,6562.8D]
    lineratio = [0.4577,0.094,0.2821,2.002]
    nline = n_elements(linewave)
    
; read the fitted sample
    if keyword_set(debug) then zcat = mrdfits(desipath+'deep2_zcat.fits.gz',1)
    ispec = mrdfits(desipath+'deep2_ispec.fits.gz',1)
    ngal = n_elements(ispec)

    weff = k_lambda_eff(filterlist=deep2_filterlist())
    nband = n_elements(weff)

; divide the sample into chunks
    chunksize = 500L
    nchunk = ceil(ngal/float(chunksize))
    
    for ichunk = 0, 0 do begin
;   for ichunk = 0, nchunk-1 do begin
       i1 = ichunk*chunksize
       i2 = ((ichunk*chunksize+chunksize)<ngal)-1L
       nthese = i2-i1+1L
       these = lindgen(nthese)+i1

       ised = read_isedfit(desipath+'desi_deep2_paramfile.par',/getmodels,$
         isedfit_dir=desipath,montegrids_dir=desipath+'montegrids/',$
         index=these,/flam)
       
       for igal = 0, nthese-1 do begin
;; read the observed DEEP2 spectrum       
;         ss = mrdfits('testspec.fits.gz',1)
;         good = where(ss.wave gt 0)
;;        datawave = exp(ss.wave[good])*(1+ispec[these[igal]].zline_forbidden)
;         datawave = exp(ss.wave[good])*(1+ss.z)
;         dataflux = ss.flux[good]/(1+ss.z)
;         dataflux = ss.flux[good]/ss.continuum[good]*$
;           interpolate(cflux,findex(wave,datawave))

; restore the best-fit iSEDfit model and resample to constant
; log-lambda; note that we're ignoring the change in spectral
; resolution of the underlying continuum spectrum 
          ised1 = ised[igal]
          continuum = im_log_rebin(im_double(ised1.wave),ised1.flux,$
            vsc=velpixsize,minwave=minwave,maxwave=maxwave,outwave=obswave)

          obswave = exp(obswave)                          ; observed frame [A]
          logwave = alog10(obswave/(1+ispec[these[igal]].z)) ; rest-frame [A]
          npix = n_elements(obswave)

; get the *rest-frame* continuum flux around [OII]
          cindx = where(obswave gt 3727.0*(1+ispec[these[igal]].z)*0.995 and $
            obswave lt 3727.0*(1+ispec[these[igal]].z)*1.005)
          cflux = djs_median(continuum[cindx]*(1+ispec[these[igal]].z))
;         djs_plot, wave[these], cflux[these]

; mean intrinsic velocity width [km/s]       
          linesigma = djs_mean([ispec[these[igal]].oii_3727_1_sigma[0],$
            ispec[these[igal]].oii_3727_2_sigma[0]])
;         linesigma = sqrt(linesigma^2+ispec[these[igal]].line_inst_vdisp^2)

; rebuild the measured emission-line spectrum at the new velocity
; resolution, all relative to [OII] 3729; this step assumes that our 
; emission-line fitting is spot-on
          reflineflux = ispec[these[igal]].oii_3727_2_ew[0]*cflux ; [erg/s/cm2]
          emspectrum = build_emline(continuum*0,logwave=logwave,zshift=zshift,$
            lineflux=reflineflux,linesigma=linesigma,$
            linewave=ispec[these[igal]].oii_3727_2_wave)
       
          lineflux = reflineflux*ispec[these[igal]].oii_3727_1[0]/$
            ispec[these[igal]].oii_3727_2[0]
          emspectrum = build_emline(emspectrum,logwave=logwave,zshift=zshift,$
            lineflux=lineflux,linesigma=linesigma,$
            linewave=ispec[these[igal]].oii_3727_1_wave)

; now do the other nebular lines
          for ll = 0, nline-1 do emspectrum = build_emline(emspectrum,$
            logwave=logwave,zshift=zshift,lineflux=lineratio[ll]*reflineflux,$
            linesigma=linesigma,linewave=linewave[ll])

; shift to the observed frame and build the final spectrum
          emspectrum = emspectrum/(1+ispec[these[igal]].z)
          flux = continuum + emspectrum ; [erg/s/cm2/A]

          if keyword_set(debug) then begin
             mag = maggies2mag(zcat[these[igal]].maggies,magerr=magerr,$
               ivarmaggies=zcat[these[igal]].ivarmaggies)
             abflux = -2.5*alog10((flux*obswave^2.0/light_ang)>1D-50)-48.6
             yrange = [max(abflux),min(abflux)]
             djs_plot, obswave/1D4, abflux, psym=10, xsty=1, $
               ysty=3, xrange=xrange, yrange=yrange, $
               xtitle='Wavelength (\mu'+'m)', ytitle='Flux'
             oploterror, weff/1D4, mag, magerr, psym=6, $
               color=cgcolor('green'), errcolor=cgcolor('green')
             cc = get_kbrd(1)
          endif

          if igal eq 0 then outflux = fltarr(nthese,npix) ; output spectra

          outflux[igal,*] = flux
;         outinfo = ispec[these]

          outinfo = struct_trimtags(ispec[these],select=['GALAXY',$
            'ZCATINDX','MINWAVE','MAXWAVE','OBJNO','RA','DEC',$
            'MASK','SLIT','Z','CONTINUUM_SNR','OII_3727_*'])
       endfor ; close galaxy loop

       outfile = desipath+'desi_deep2_spectra/desi_deep2_chunk'+$
         string(ichunk+1,format='(I2.2)')+'.fits'
       
       mkhdr, hdr, outflux, /extend
       sxdelpar, hdr, 'DATE'
       sxdelpar, hdr, 'COMMENT'
       sxaddpar, hdr, 'OBJTYPE', 'ELG-DEEP2', 'object type'
       sxaddpar, hdr, 'FLUXUNIT', 'erg/s/cm^2/A', ' spectrum flux units'
       sxaddpar, hdr, 'CRVAL1', coeff0, ' reference log10(Angstrom)'
       sxaddpar, hdr, 'CDELT1', coeff1, ' delta log10(Angstrom)'
       sxaddpar, hdr, 'LOGLAM', 1, ' log10 spaced wavelengths'
       sxaddpar, hdr, 'VELSCALE', velpixsize, ' pixel size in km/s'

       im_mwrfits, transpose(outflux), outfile, hdr, /clobber, /nogzip

;      mkhdr, hdr, 0
;      sxdelpar, hdr, 'DATE'
;      sxdelpar, hdr, 'COMMENT'
;      sxaddpar, hdr, 'OBJTYPE', 'ELG-DEEP2', 'object type'
       
       im_mwrfits, outinfo, outfile, /append, /gzip

       
stop       
    endfor                      ; close chunk loop

return
end
