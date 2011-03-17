function parse_ascii_header, ascii_header

    ascii_header = repstr(ascii_header[3:40],'# ','')
    keep = where(strcompress(ascii_header,/remove) ne '',nkeep)
    ascii_header = ascii_header[keep]

    hinfo = {keyword: '', value: '', comment: ''}
    hinfo = replicate(hinfo,nkeep)

    for i = 0L, nkeep-1L do begin
       hinfo[i].keyword = strmid(ascii_header[i],0,9)
       hinfo[i].value = repstr(strmid(ascii_header[i],11,30),"'",'')
       hinfo[i].comment = strmid(ascii_header[i],42,strlen(ascii_header[i])-1)
    endfor

return, hinfo
end    
    
pro unpack_nfgs, doplot=doplot, wfits=wfits
; jm04jan6uofa
; parse the raw NFGS ASCII spectra and write out gzipped FITS files
; and headers

    minwave = 3620.0
    maxwave = 7000.0
    
    datapath = nfgs_path(/analysis)
    outpath = nfgs_path(/spec1d)
    rawpath = nfgs_path(/original_spec1d)

    flist = file_basename(file_search(rawpath+'*.spec',count=ngalaxy))

; read the Jansen et al. 2000 database; we assume that JANSEN and
; FLIST are sorted properly

    jansen = read_00jansen()
;   niceprint, jansen.galaxy, jansen.nedgalaxy, jansen.nfgs_name, flist

; read the spectral normalizations
    
    readcol, datapath+'tab_contF5500.dat', contspec, type, junk, f5500_int, $
      int_prob, int_max, f5500_nuc, nuc_int, nuc_max, comment='#', $
      /silent, format='A,L,F,F,F,F,F,F,F'
    
;   for i = 30, 30 do begin
    for i = 0L, ngalaxy-1L do begin

       readfast, rawpath+flist[i], data, ascii_header, skipline=51
       hinfo = parse_ascii_header(ascii_header)
       wave = reform(data[0,*])
       
       if (total(data[3,*]) ne 0.0) then begin ; integrated spectrum

          spec = reform(data[3,*])*f5500_int[i]*1D-14
          sigspec = reform(data[4,*])*f5500_int[i]*1D-14

          good = where((spec gt 0.0) and (sigspec gt 0.0),ngood)
          if (ngood ne 0L) then begin
             spec = spec[good]
             sigspec = sigspec[good]
             oldwave = wave[good]
          endif
          npix = n_elements(spec)

; add systematic errors

;         sigspec = sqrt(sigspec^2.0 + (0.06*spec)^2.0) ; spectrophotometric error

          red = where(wave gt 6775.0,nred)
          if (nred ne 0L) then sigspec[red] = sqrt(sigspec[red]^2.0 + (0.05*spec[red])^2.0) ; blue contamination

; parse the header info          
          
          mkhdr, header, spec, /extend
          sxdelpar, header, 'COMMENT'

          for j = 0L, n_elements(hinfo)-1L do sxaddpar, header, $
            hinfo[j].keyword, strtrim(hinfo[j].value,2), ' '+$
            strtrim(hinfo[j].comment,2)
          z = float(sxpar(header,'Z-OBS'))
          z_err = float(sxpar(header,'EZ-OBS'))
          
          oldwave = oldwave*(1+z)
          
          dwave = fix(100.0*(djs_median(oldwave-shift(oldwave,1))))/100.0 ; [Angstrom/pixel]                   
          newwave = findgen((maxwave-minwave)/dwave+1)*dwave+minwave

          linterp, oldwave, spec, newwave, newspec, missing=missing
          linterp, oldwave, sigspec, newwave, newsigspec, missing=missing
          npix = n_elements(newspec)

; de-redden for foreground Galactic reddening

          glactc, 15.0D*im_hms2dec(sxpar(header,'RA')), im_hms2dec(sxpar(header,'DEC')), $
            1950.0, gl, gb, 1, /degree, /fk4
          ebv_mw = dust_getval(gl,gb,/interp)
          kl = k_lambda(newwave,/odonnell,R_V=3.1)
          newspec = newspec*10^(0.4*ebv_mw*kl)

; update the header          
          
          sxaddpar, header, 'RA', jansen[i].ra, ' right ascension [HMS]'
          sxaddpar, header, 'DEC', jansen[i].dec, ' declination [DMS]', after='RA'
          sxaddpar, header, 'EPOCH', 2000.0, ' coordinate epoch', format='(F6.1)', after='DEC'
          sxaddpar, header, 'EQUINOX', 2000.0, ' coordinate epoch', format='(F6.1)', after='EPOCH'

          galaxy = strtrim(jansen[i].ned_galaxy,2) ; <-- NOTE! Use the NED Galaxy Name!
;         galaxy = strtrim(jansen[i].galaxy,2)
          nedgal = strtrim(jansen[i].ned_galaxy,2)
          nfgsgal = strtrim(jansen[i].nfgs_galaxy,2)

; Rolf does not give a scan length and extraction aperture for this
; object because of a superposed star; however, an integrated spectrum
; exists, so assume the major- and minor-axis diameters          
          
          if strmatch(jansen[i].galaxy,'*UGC02296*') then begin
             jansen[i].nfgs_aperwid = 0.6
             jansen[i].nfgs_scanlen = 0.6
          endif
          
          sxaddpar, header, 'GALAXY', galaxy
;         sxaddpar, header, 'NEDGAL', nedgal
;         sxaddpar, header, 'NFGSGAL', nfgsgal
;         sxaddpar, header, 'OBJECT', object
;         sxaddpar, header, 'NFGSGAL', nfgsgal
          sxaddpar, header, 'NFGS_ID', jansen[i].nfgs_id
          sxaddpar, header, 'CRVAL1', minwave, ' central wavelength of first pixel'
          sxaddpar, header, 'CRPIX1', 1, ' starting pixel (1-indexed)'
          sxaddpar, header, 'CD1_1', dwave, ' dispersion [Angstrom/pixel]'
          sxaddpar, header, 'CDELT1', dwave, ' dispersion [Angstrom/pixel]'
          sxaddpar, header, 'CTYPE1', 'LINEAR', ' projection type'
;         sxaddpar, header, 'DC-FLAG', 0, ' log-linear flag'
          sxaddpar, header, 'Z', z, ' redshift'
          sxaddpar, header, 'Z_ERR', z_err, ' redshift error'
          sxaddpar, header, 'SCANLEN', jansen[i].nfgs_scanlen*60.0, $
                  ' scan length [arcsec]', format='(F12.2)'
          sxaddpar, header, 'POSANGLE', float(jansen[i].nfgs_pas), $
                  ' slit position angle [degrees]', format='(F12.1)'
          sxaddpar, header, 'APERWID', jansen[i].nfgs_aperwid*60.0, $
                  ' extraction aperture diameter [arcsec]', format='(F12.2)'
          sxaddpar, header1d, 'MEDSNR', median(newspec/newsigspec), $
            ' median signal-to-noise per pixel'

          fitsfile = strmid(flist[i],0,3)+'.'+strlowcase(galaxy)+'_int.fits'
          if (not keyword_set(wfits)) then splog, i+1, galaxy, nfgsgal, fitsfile, format='(I3,A15,A30)'

          if keyword_set(doplot) then begin
             ploterror, newwave, 1E15*newspec, 1E15*newsigspec, ps=10, xsty=3, ysty=3, $
;            plot, newwave, 1E15*newspec, ps=10, xsty=3, ysty=3, $; xrange=[6500,6620]*(1+z), $
               charsize=2.0, charthick=2.0, xthick=2.0, ythick=2.0, $
               xtitle='Wavelength', ytitle='Normalized Flux'
             print, median(newspec/newsigspec)
             cc = get_kbrd(1)
          endif
          
          if keyword_set(wfits) then begin

             splog, 'Writing '+outpath+fitsfile+'.'
             mwrfits, float(newspec), outpath+fitsfile, header, /create
             mwrfits, float(newsigspec), outpath+fitsfile
;            spawn, ['gzip -f '+outpath+fitsfile], /sh
             
          endif
          
       endif else begin

          galaxy = ''
          fitsfile = flist[i]
          splog, i+1, galaxy, fitsfile, ' -- no integrated data.', format='(I3,A15,A30,A25)'

       endelse

    endfor

return
end    
