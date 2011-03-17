pro alpha_write_spec1d, logwave, flux, ivar, skyflux, $
  info, cdelt=cdelt, outfile=outfile, std=std
; jm09jan07nyu - write out a 1D spectrum

    mkhdr, hdr, float(flux), /extend
    sxdelpar, hdr, 'COMMENT'
    sxdelpar, hdr, 'DATE'
    sxaddpar, hdr, 'OBJECT', info.obj
    sxaddpar, hdr, 'EXPTIME', float(info.exp)
    sxaddpar, hdr, 'AIRMASS', float(info.am)
    sxaddpar, hdr, 'RA', info.ra, format='(F11.7)'
    sxaddpar, hdr, 'DEC', info.dec, format='(F13.9)'
    if (not keyword_set(std)) then begin
       sxaddpar, hdr, 'GALAXY', info.galaxy
       sxaddpar, hdr, 'Z', float(info.z)
    endif
    sxaddpar, hdr, 'RAWFITS', info.img_root
    sxaddpar, hdr, 'CTYPE1', 'LOG'
    sxaddpar, hdr, 'CRPIX1', 1.0d
    sxaddpar, hdr, 'CRVAL1', min(logwave)
    sxaddpar, hdr, 'CDELT1', cdelt
    sxaddpar, hdr, 'DC-FLAG', 1
    
    splog, 'Writing '+outfile
    mwrfits, float(flux), outfile, hdr, /create
    mwrfits, float(ivar), outfile, hdr
    mwrfits, float(skyflux), outfile, hdr
    
return
end
    
