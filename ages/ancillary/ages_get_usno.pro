pro ages_get_usno, out
; jm10feb02ucsd - grab all the USNO stars in the field and write out 

    usno = usno_read(218.0,34.5,3.0,catname='usno-a2.0') ; 'usno-b1.0')
    keep = where((usno.ra gt 216.0) and (usno.ra lt 220.0) and $
      (usno.dec gt 32.0) and (usno.dec lt 36.0),nkeep)
    splog, 'N = ', nkeep

    out = usno[keep]
    outfile = ages_path(/mycatalogs)+'ages_usno.fits'
    im_mwrfits, out, outfile, /clobber

return
end
    
