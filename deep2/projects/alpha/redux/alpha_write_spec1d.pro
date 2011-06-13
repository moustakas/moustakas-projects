pro alpha_write_spec1d, lnwave, flux, ivar, skyflux, $
  info, velpix=velpix, outfile=outfile, std=std
; jm09jan07nyu - write out a 1D spectrum

    mkhdr, hdr, float(flux), /extend
    sxdelpar, hdr, 'COMMENT'
    sxdelpar, hdr, 'DATE'
    sxaddpar, hdr, 'OBJECT', info.obj
    sxaddpar, hdr, 'EXPTIME', float(info.exp)
    sxaddpar, hdr, 'AIRMASS', float(info.am)
    sxaddpar, hdr, 'RA', info.ra, format='(F11.7)'
    sxaddpar, hdr, 'DEC', info.dec, format='(F13.9)'
    if (keyword_set(std) eq 0) then begin
       sxaddpar, hdr, 'GALAXY', info.galaxy
       sxaddpar, hdr, 'Z', float(info.z)
    endif
    sxaddpar, hdr, 'RAWFITS', strtrim(info.img_root,2)
    sxaddpar, hdr, 'CTYPE1', 'LN'
    sxaddpar, hdr, 'CRPIX1', 1D
    sxaddpar, hdr, 'CRVAL1', min(lnwave)
    sxaddpar, hdr, 'CD1_1', velpix
    sxaddpar, hdr, 'DC-FLAG', 1
    
    splog, 'Writing '+outfile
    mwrfits, float(flux), outfile, hdr, /create
    mwrfits, float(ivar), outfile, hdr
    mwrfits, float(skyflux), outfile, hdr
    
return
end
    
