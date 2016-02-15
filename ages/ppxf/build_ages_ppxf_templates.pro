function ages_convolve, oldflux, old_velscale, new_velscale
; convolve to the AGES instrumental resolution
    ntemp = (size(oldflux,/dim))[1]
    vdisp = sqrt(new_velscale^2 - old_velscale^2)
    smoothing = vdisp/old_velscale  ; [pixel]
;   kernel = psf_gaussian(npix=10.0*smoothing,$
;   fwhm=2.35*smoothing,/norm,ndim=1)
    nkpix = long(4.0*ceil(smoothing))*2L+3
    klam = findgen(nkpix)-float(nkpix-1.0)/2.0
    kernel = exp(-0.5*(klam/smoothing)^2)/sqrt(2.0*!dpi)/smoothing
    kernel = kernel/total(kernel)
    newflux = oldflux
    for ii = 0, ntemp-1 do newflux[*,ii] = $
      convol(oldflux[*,ii],kernel,/edge_truncate)
return, newflux
end

pro build_ages_ppxf_templates, bc03=bc03, debug=debug
; jm09nov12ucsd - write out the basic set of BC03 templates for the
;   PPXF template fitting
; jm16feb10siena - update to CKC14z templates (probably no longer backwards
;   compatible with BC03)

    outpath = ages_path(/ppxf)
    version = ages_version(/ppxf_templates)
    outfile = outpath+'ages_ppxf_ckc14z_templates_'+version+'.fits'

    light = 299792.458D
    lsun = 3.826D33         ; [erg/s]
    dist = 10.0*3.085678D18 ; 10 pc [cm]

;   minwave = 500.0D
;   maxwave = 5D4
;   dwave = 2.0D
;   npix = round((maxwave-minwave)/dwave+1L)
;   wave = minwave+dindgen(npix)*dwave
    
; resample to be constant in log10-lambda at the AGES pixel scale
    velscale = ages_ppxf_velscale()
    pixsize = velscale/light/alog(10) ; [pixel size in log-10 A]

;   minwave = 2000D
;   maxwave = 8000D
    minwave = 500D
    maxwave = 5D4
;   npix = round((maxwave-minwave)/pixsize+1L)
;   wave = minwave+dindgen(npix)*pixsize
    
    agegrid = [5.0,25.0,50.0,125.0,300.0,650.0,$
      1500.0,3000.0,7000.0,13000.0]

    if keyword_set(bc03) then begin
       bc03_Z004 = im_read_bc03(age=agegrid*1E6/1E9,$ ; chabrier IMF
         bc03_extras=extras,/silent,metallicity=2)
       bc03_Z02 = im_read_bc03(age=agegrid*1E6/1E9,$
         bc03_extras=extras,/silent,metallicity=4)
       bc03_Z05 = im_read_bc03(age=agegrid*1E6/1E9,$
         bc03_extras=extras,/silent,metallicity=5)
       bc03 = [bc03_Z004,bc03_Z02,bc03_Z05]

       outfile = outpath+'bc03_'+string(metal,format='(F5.3)')+$
         '_'+version+'.fits'
       for jj = 0, nmetal-1 do begin
          flux = fltarr(npix,nmodel)
          for ii = 0, nmodel-1 do begin
             linterp, bc03[jj].wave, bc03[jj].flux[*,ii], wave, newflux
             flux[*,ii] = lsun*newflux/extras[ii].m_/(4.0*!dpi*dist^2)
             if keyword_set(debug) then begin
                djs_plot, wave, newflux, xr=[1000,4E4], /xlog, /ylog
                djs_oplot, bc03[jj].wave, bc03[jj].flux[*,ii], color='orange'
                djs_oplot, wave, newflux1, color='blue'
                cc = get_kbrd(1)
             endif
          endfor
          
          info = {metallicity: metal[jj], age: fltarr(nmodel)}
          info.age = bc03[jj].age/1E9 ; [Gyr]
          
; instrumental velocity dispersion of the models    
          inst_vdisp = 60.0     ; [km/s]
          
; write out; generate an appropriate FITS header
          mkhdr, hdr, float(flux)
          sxdelpar, hdr, 'COMMENT'
          sxaddpar, hdr, 'CRVAL1', min(wave), ' central wavelength of first pixel [Angstrom]'
          sxaddpar, hdr, 'CRPIX1', 1, ' starting pixel (1-indexed)'
          sxaddpar, hdr, 'CD1_1', dwave, ' pixel size [Angstrom/pixel]'
          sxaddpar, hdr, 'VDISP', velscale, ' instrumental velocity dispersion [km/s]'
;         sxaddpar, hdr, 'VDISP', inst_vdisp, ' instrumental velocity dispersion [km/s]'
          sxaddpar, hdr, 'DATE', hogg_iso_date(), ' file creation date'
          
; write out
          splog, 'Writing '+outfile[jj]
          mwrfits, float(flux), outfile[jj], hdr, /create
          mwrfits, info, outfile[jj]
          spawn, 'gzip -f '+outfile[jj], /sh
       endfor
    endif else begin

       Zstr = ['Z0.0049','Z0.0190','Z0.0300']
       metal = [0.0049,0.019,0.03] ; metallicity grid
       nmodel = n_elements(agegrid)
       nmetal = n_elements(metal)

       ckc_velscale = 15D ; [km/s]
;      ckc_velscale = djs_median(alog10(ckc.wave)-shift(alog10(ckc.wave),1))*alog(10)*light
       ckc_pixsize = ckc_velscale/light/alog(10)
       
       outfile = outpath+'ckc14z_'+string(metal,format='(F6.4)')+$
         '_'+version+'.fits'
       for jj = 0, nmetal-1 do begin
          ckc = im_read_fsps(metallicity=Zstr[jj],/ckc,/flam)

; convolve to lower spectral resolution          
          these = where(ckc.wave gt minwave*0.98 and ckc.wave lt maxwave*1.02)
          ckcwave = ckc.wave[these]
          ckcflux = ckc.flux[these,*]*1D
;         ckcwave = ckc.wave
;         ckcflux = ckc.flux

          ageindx = findex(ckc.age/1D6,agegrid)
          ageflux = interpolate(ckcflux,ageindx)

          cflux = ages_convolve(ageflux,ckc_velscale,velscale)
          for ii = 0, nmodel-1 do begin
;            flux1 = im_log_rebin(ckc.wave+randomn(seed,n_elements(ckc.wave))*0.001,reform(cflux[*,ii]),$
             flux1 = im_log_rebin(ckcwave,reform(cflux[*,ii]),$
               vsc=velscale,outwave=lnwave,minwave=minwave,maxwave=maxwave)
             if ii eq 0 then begin
                npix = n_elements(flux1)
                flux = fltarr(npix,nmodel)
             endif
             flux[*,ii] = flux1
          endfor
;         flux = interpolate(flux,findex(alog10(ckc.wave),wave),lindgen(nmodel),/grid)

          if keyword_set(debug) then begin
             for ii = 0, nmodel-1 do begin
                djs_plot, ckcwave, ageflux[*,ii], xsty=1, ysty=1, xr=[minwave,maxwave], /xlog, /ylog ;, xr=[3000,9000];, xr=[3700,4000]
                djs_oplot, ckcwave, cflux[*,ii], color='blue'
                djs_oplot, exp(lnwave), flux[*,ii], color='orange'
                cc = get_kbrd(1)
             endfor
          endif
          
; divide by the stellar mass          
;         flux = interpolate(ckc.flux,findex(ckc.wave,wave),ageindx,/grid)
          flux = flux/rebin(reform(interpolate(ckc.mstar,ageindx),1,nmodel),npix,nmodel)
          
          info = {metallicity: metal[jj], age: fltarr(nmodel)}
          info.age = agegrid/1D3 ; [Gyr]
          
; instrumental velocity dispersion of the models    
;         inst_vdisp = 15.0     ; [km/s]
          
; write out; generate an appropriate FITS header
          mkhdr, hdr, float(flux)
          sxdelpar, hdr, 'COMMENT'
          sxaddpar, hdr, 'DISPAXIS', 1, ' dispersion axis'
          sxaddpar, hdr, 'CTYPE1', 'WAVE-WAV', ' meaning of wavelength array'
          sxaddpar, hdr, 'CUNIT1', 'Angstrom', ' units of wavelength array'
          sxaddpar, hdr, 'CRPIX1', 1, ' reference pixel number'
          sxaddpar, hdr, 'CRVAL1', min(lnwave), ' reference ln(Angstrom)'
          sxaddpar, hdr, 'CDELT1', pixsize*alog(10), ' delta ln(Angstrom)'
;         sxaddpar, hdr, 'CDELT1', pixsize, ' delta log10(Angstrom)'
          sxaddpar, hdr, 'LOGLAM', 1, ' log10 spaced wavelengths?'
          sxaddpar, hdr, 'WAVEMIN', min(exp(lnwave)), ' minimum wavelength (Angstrom)'
          sxaddpar, hdr, 'WAVEMAX', max(exp(lnwave)), ' maximum wavelength (Angstrom)'
          sxaddpar, hdr, 'WAVEUNIT', 'Angstrom', ' wavelength units'
          sxaddpar, hdr, 'AIRORVAC', 'AIR', ' wavelengths in vacuum (vac) or air'
          sxaddpar, hdr, 'DATE', hogg_iso_date(), ' file creation date'

; write out
          splog, 'Writing '+outfile[jj]
          mwrfits, float(flux), outfile[jj], hdr, /create
          mwrfits, info, outfile[jj]
          spawn, 'gzip -f '+outfile[jj], /sh
       endfor
    endelse
       
return
end
    
