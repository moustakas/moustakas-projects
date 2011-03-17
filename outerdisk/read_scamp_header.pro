function read_scamp_header, fitsfile, ext=ext

    if (n_elements(ext) eq 0L) then ext = 1L
    allhdr = djs_readlines(fitsfile)
    len = (where(strmatch(allhdr,'END*',/fold) eq 1B))[0] + 1L
    iext = ext-1L
    hdr = allhdr[iext*len:(iext+1)*len-1L]

return, hdr
end
    
