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

pro build_desi_deep2_templates, debug=debug
; jm13dec18siena - take the output of BUILD_DESI_DEEP2_TEMPLATE_SAMPLE
;   and generate the full-resolution spectra

    version = desi_deep2_template_version()
    templatepath = getenv('IM_PROJECTS_DIR')+'/desi/templates/'+version+'/'
    
    light = im_light(/kms)
    light_ang = im_light(/Ang)
    weff = k_lambda_eff(filterlist=deep2_filterlist())
    airtovac, weff
    synth_filters = (decam_filterlist())[[1,2,4]]

    oiifluxcut = 8D-17
    faintmag = 23.5
    pc10 = 3.085678D19          ; fiducial distance [10 pc in cm]

; wavelength of the [OII] 3729 reference line in vacuum
    oiirefwave = 3728.8149D
    airtovac, oiirefwave
    
; build the output *rest-wavelength* array; sample the rest-frame UV
; and optical wavelengths at high resolution, and the near-infrared at
; lower resolution (since it will only be used for filter
; convolutions) 
    
; simulation parameter defaults (rest-frame *vacuum* wavelengths!)
    velpixsize_hires = 20D ; [km/s]
;   velpixsize_lowres = 20D ; [km/s]
    velpixsize_lowres = 100D ; [km/s]

    pixsize_hires = velpixsize_hires/light/alog(10) ; [pixel size in log-10 A]
    pixsize_lowres = velpixsize_lowres/light/alog(10) ; [pixel size in log-10 A]

    minwave_hires = alog10(1500D)
    maxwave_hires = alog10(10500D)
    minwave_lowres = maxwave_hires+pixsize_lowres
    maxwave_lowres = alog10(3.5D4)

    npix_hires = round((maxwave_hires-minwave_hires)/pixsize_hires+1L)
    npix_lowres = round((maxwave_lowres-minwave_lowres)/pixsize_lowres+1L)

    restwave_hires = minwave_hires+dindgen(npix_hires)*pixsize_hires     ; log-10 spacing
    restwave_lowres = minwave_lowres+dindgen(npix_lowres)*pixsize_lowres ; log-10 spacing
    restwave = [restwave_hires,restwave_lowres]
    npix = n_elements(restwave)

; temporarily write out just the high-resolution templates
    restwave = restwave_hires
    npix = n_elements(restwave)
    
; read the fitted sample and select a fiducial subset 
    info = mrdfits(templatepath+'desi_deep2elg_template_sample_'+$
      version+'.fits.gz',1)
    keep = where(info.ugriz[2] lt faintmag,ngal)
    splog, 'Total number of templates', ngal

; note the minimum floor!    
    info[keep].sigma_kms = info[keep].sigma_kms > pixsize_hires*2.5 
    
; initialize the output data structure and header; the metadata header
; is written out below
    outinfo = {$
      objno:        0L,$
      ra:           0D,$
      dec:          0D,$
      z:           0.0,$
      sigma_kms:   0.0,$
      oii_3726:    0.0,$
      oii_3729:    0.0,$
      oii_3727:    0.0,$
      oii_3727_ew: 0.0,$ 
      mag_g:       0.0,$
      mag_r:       0.0,$
      mag_z:       0.0,$
      logmstar:    0.0,$
      logsfr:      0.0}
    outinfo = replicate(outinfo,ngal)
    outinfo.objno = info[keep].objno
    outinfo.ra = info[keep].ra
    outinfo.dec = info[keep].dec
    outinfo.z = info[keep].z
    outinfo.sigma_kms = info[keep].sigma_kms
    outinfo.oii_3726 = info[keep].oii_3727_1[0]
    outinfo.oii_3729 = info[keep].oii_3727_2[0]
    outinfo.oii_3727 = info[keep].oii_3727[0]
    outinfo.oii_3727_ew = info[keep].oii_3727_ew[0]

    outflux = fltarr(ngal,npix)
;   outflux = {flux: fltarr(ngal,npix), loglam: restwave}

    mkhdr, fluxhdr, outflux, /extend
    sxdelpar, fluxhdr, 'DATE'
    sxdelpar, fluxhdr, 'COMMENT'
    sxaddpar, fluxhdr, 'OBJTYPE', 'ELG-DEEP2', 'object type'
    sxaddpar, fluxhdr, 'DISPAXIS', 1, ' dispersion axis'
    sxaddpar, fluxhdr, 'CTYPE1', 'WAVE-WAV', ' meaning of wavelength array'
    sxaddpar, fluxhdr, 'CUNIT1', 'Angstrom', ' units of wavelength array'
    sxaddpar, fluxhdr, 'CRPIX1', 1, ' reference pixel number'
    sxaddpar, fluxhdr, 'CRVAL1', min(restwave), ' reference log10(Angstrom)'
    sxaddpar, fluxhdr, 'CDELT1', pixsize_hires, ' delta log10(Angstrom)'
    sxaddpar, fluxhdr, 'LOGLAM', 1, ' log10 spaced wavelengths'
    sxaddpar, fluxhdr, 'WAVEMIN', min(10D^restwave), ' minimum wavelength'
    sxaddpar, fluxhdr, 'WAVEMAX', max(10D^restwave), ' maximum wavelength'
    sxaddpar, fluxhdr, 'WAVEUNIT', 'Angstrom', ' wavelength units'
    sxaddpar, fluxhdr, 'AIRORVAC', 'vac', ' wavelengths in vacuum (vac) or air'
    sxaddpar, fluxhdr, 'VELSCALE', velpixsize_hires, ' pixel size in km/s'
    sxaddpar, fluxhdr, 'BUNIT', 'erg/s/cm2/A', ' spectrum flux units'
    sxaddpar, fluxhdr, 'FLUXUNIT', 'erg/s/cm2/A', ' spectrum flux units'

; divide the sample into chunks
;   splog, 'HACK!'
;   chunksize = 5
;   chunksize = 500L
    chunksize = ngal
    nchunk = ceil(ngal/float(chunksize))
    
    for ichunk = 0, nchunk-1 do begin
       i1 = ichunk*chunksize
       i2 = ((ichunk*chunksize+chunksize)<ngal)-1L
       nthese = i2-i1+1L
       these = lindgen(nthese)+i1

; read the iSEDfit continuum models into the rest-frame [erg/s/cm2/A],
; ignoring IGM attenuation
       ised = read_isedfit(templatepath+'desi_deep2_paramfile.par',/getmodels,$
         isedfit_dir=templatepath,montegrids_dir=templatepath+'montegrids/',$
         index=info[these].isedfit_id,/flam,/noigm,/restframe,thissfhgrid=1)
       outinfo[these].logmstar = ised.mstar_avg
       outinfo[these].logsfr = ised.sfr_avg
;      outinfo[these].logmstar = ised.mstar ; ised.mstar_50
;      outinfo[these].logsfr = ised.sfr     ; ised.sfr_50

       diff = fltarr(nthese)
       for igal = 0, nthese-1 do begin
          print, format='("Building template ",I0,"/",I0,A10,$)', igal, ngal, string(13b)
; resample to constant log10-lambda; note that we're ignoring the
; change in spectral resolution of the underlying continuum spectrum;
; the difference between IM_LOG_REBIN() and INTERPOL() is at the
; <+/-1% (peak-to-peak) level; therefore, since IM_LOG_REBIN() has
; trouble when the pixel spacing changes, just use INTERPOL()
          ised1 = ised[igal]
;         continuum_lowres = im_log_rebin(im_double(ised1.wave),ised1.flux,/log10,$  
;           minwave=10D^minwave_lowres,maxwave=10D^maxwave_lowres,vsc=velpixsize_lowres)
;         continuum_hires = im_log_rebin(im_double(ised1.wave),ised1.flux,/log10,$  
;           minwave=10D^minwave_hires,maxwave=10D^maxwave_hires,vsc=velpixsize_hires)
;         continuum = [continuum_lowres,continuum_hires]
          isedwave = ised1.wave
          airtovac = isedwave
          continuum = interpol(ised1.flux,alog10(isedwave),restwave)
          if n_elements(continuum) ne npix then message, 'Problem here!'
          
; rebuild the measured emission-line spectrum at the new velocity
; resolution, all relative to [OII] 3729; this step assumes that our 
; emission-line fitting is spot-on; for now use a fixed doublet ratio 
          linesigma = info[these[igal]].sigma_kms
          reflineflux = info[these[igal]].oii_3727_2[0] ; [erg/s/cm2]
          oiidoubletratio = info[these[igal]].oii_3727_1[0]/info[these[igal]].oii_3727_2[0]
          emspectrum = build_emline(continuum*0,logwave=restwave,zshift=zshift,$
            lineflux=reflineflux,linesigma=linesigma,linewave=oiirefwave)

; generate the theoretical emission-line and nebular continuum
; spectrum for this galaxy using the number of Lyman-continuum photons
; inferred from SED-fitting; 
          nebflux = isedfit_nebular(10D^(ised1.mstar+ised1.nlyc),inst_vsigma=0D,$
            vsigma=linesigma,oiiihb=oiiihb,wave=10D^restwave,line=line,$
            oiidoubletratio=oiidoubletratio,flam_line=flam_line,flam_cont=flam_cont,$
            /vacuum)

stop          
          
; NEBFLUX is defined at a distance of 10 pc; scale by the luminosity
; distance but leave the factor of (1+z) off so that the flux vector
; is in the rest frame
          dlumfactor = 1D/(10D^(lf_distmod(info[these[igal]].z,$
            omega0=0.3D,omegal0=0.7D)/5D)/0.7D)^2
;         dlumfactor = (pc10/dluminosity(info[these[igal]].z,/cm))^2D

          nebflux = nebflux*dlumfactor
          flam_line = flam_line*dlumfactor
          line.amp = line.amp*dlumfactor
          line.flux = line.flux*dlumfactor

; rescale the emission-line data structure by the desired [OII] flux
          isoii3729 = where(strtrim(line.name,2) eq '[OII]_3729')
          oiifactor = reflineflux/line[isoii3729].flux
          for ll = 0, n_elements(line)-1 do begin
             line[ll].flux = line[ll].flux*oiifactor
             line[ll].amp = line[ll].amp*oiifactor
          endfor

;; test to be sure we got the right total [OII] flux
;          testemspectrum = build_emline(emspectrum,logwave=restwave,lineflux=line[11].flux,$
;            zshift=zshift,linesigma=linesigma,linewave=line[11].wave)
;          diff[igal] = im_integral(10^restwave,testemspectrum)/info[these[igal]].oii_3727[0]
;          print, im_integral(10^restwave,testemspectrum), info[these[igal]].oii_3727[0], $
;            diff[igal], linesigma

; build the full emission-line spectrum and the final complete
; spectrum (all in the *rest-frame*!)
          for ll = 0, n_elements(line)-1 do emspectrum = build_emline(emspectrum,$
            logwave=restwave,lineflux=line[ll].flux,zshift=zshift,linesigma=linesigma,$
            linewave=line[ll].wave)

          flux = continuum + emspectrum ; [erg/s/cm2/A]
          outflux[these[igal],*] = flux
;         outflux.flux[igal,*] = flux

;; test code
;          pp = read_deep2(/fix,/ppxf)
;          match, pp.objno, info[these[igal]].objno, m1, m2
;          srt = sort(m2) & m1 = m1[srt] & m2 = m2[srt]
;          niceprint, pp[m1].oii_3727_ew[0], info[these[igal]].oii_3727_ew[0]
;          niceprint, pp[m1].oii_3727_continuum[0], info[these[igal]].cflux_3727_rest
;          factor = info[these[igal]].cflux_3727_rest/pp[m1].oii_3727_continuum[0]
;
;          jj = read_deep2_gandalf_specfit_dr4(info[these[igal]],/fix)
;
;          djs_plot, 10D^restwave/1D4, flux, ysty=3, xr=[0.37,0.38]
;          djs_oplot, 10D^restwave/1D4, nebflux, color='orange'          
;          ww = where(jj.wave ne 0)
;          djs_oplot, exp(jj.wave[ww])/1D4, factor[0]*(jj.continuum[ww]+jj.linefit[ww]), color='red'
;
;          djs_plot, 10D^restwave/1D4, emspectrum, ysty=3, xr=[0.37,0.375], psym=10
;          djs_oplot, 10D^restwave/1D4, nebflux, color='orange', psym=10
;          ww = where(jj.wave ne 0)
;          djs_oplot, exp(jj.wave[ww])/1D4, jj.linefit[ww], color='red', psym=10

; synthesize photometry and pack in the output structure
          grz = k_project_filters(k_lambda_to_edges(10D^restwave*(1+info[these[igal]].z)),$
            flux/(1+info[these[igal]].z),filterlist=synth_filters)
;         ugriz = -2.5*alog10(reform(k_project_filters(obswave_edges,flux,$
;           filterlist=[sdss_filterlist()])))
;         outinfo[these[igal]].lineratio = lineratio
;         outinfo[these[igal]].ugriz = ugriz
;         outinfo[these[igal]].mag_u = reform(ugriz[0,*])
          outinfo[these[igal]].mag_g = -2.5*alog10(reform(grz[0,0,0]))
          outinfo[these[igal]].mag_r = -2.5*alog10(reform(grz[0,0,1]))
          outinfo[these[igal]].mag_z = -2.5*alog10(reform(grz[0,0,2]))
          
          if keyword_set(debug) then begin
             mag = maggies2mag(info[these[igal]].maggies,magerr=magerr,$
               ivarmaggies=info[these[igal]].ivarmaggies)
             obswave = 10D^restwave*(1+info[these[igal]].z)
             abflux = -2.5*alog10((flux/(1+info[these[igal]].z)*$
               obswave^2.0/light_ang)>1D-50)-48.6
             yrange = [max(abflux),min(abflux)]
             djs_plot, obswave/1D4, abflux, psym=10, xsty=1, ysty=3, xrange=[0.2,1.5], $
               yrange=yrange, xtitle='Wavelength (\mu'+'m)', ytitle='Flux'
             oploterror, weff/1D4, mag, magerr, psym=6, symsize=2.0, $
               color=cgcolor('green'), errcolor=cgcolor('green')
             cc = get_kbrd(1)
          endif

          if igal eq 50 then stop
       endfor                   ; close galaxy loop

       if nchunk eq 1 then begin
          outfile = templatepath+'desi_deep2elg_templates_'+version+'.fits'
       endif else begin
          outfile = templatepath+'desi_deep2elg_templates_'+version+'_chunk'+$
            string(ichunk+1,format='(I2.2)')+'.fits'
       endelse

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
         'OII_3727,erg/s/cm2,total [OII]\lambda\lambda3726,29 emission-line flux',$
         'OII_3727_EW,Angstrom,rest-frame [OII]\lambda\lambda3726,29 emission-line equivalent width',$
;        'MAG_U,mag,synthesized DECam u-band AB mag',$
         'MAG_G,mag,synthesized DECam g-band AB mag',$
         'MAG_R,mag,synthesized DECam r-band AB mag',$
;        'MAG_I,mag,synthesized DECam i-band AB mag',$
         'MAG_Z,mag,synthesized DECam z-band AB mag',$
         'LOGMSTAR,Msun,log10(stellar mass) (Chabrier, h=0.7)',$
         'LOGSFR,Msun/yr,log10(star formation rate) (Chabrier, h=0.7)']
       metahdr =  im_update_header(outfile,units,nheader)

; finally write everything out       
       im_mwrfits, transpose(outflux), outfile, fluxhdr, /clobber, /nogzip
       im_mwrfits, outinfo, outfile, metahdr, /append, /gzip
stop       
    endfor                      ; close chunk loop

stop    
    
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
