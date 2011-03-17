pro ediscs_merge_apercor, out
; jm09jun17nyu 

    path = ediscs_path(/catalogs)+'apercor/'
    
    all = file_search(path+'*.dat',count=nall)
    out1 = {cluster: '', galaxy: '', frac: 0.0, apercor: 0.0}
    for ii = 0, nall-1 do begin
       readcol, all[ii], galaxy, frac, apercor, format='A,F,F', $
         comment='#', /silent
       ngal = n_elements(galaxy)
       thisout = replicate(out1,ngal)
       thisout.cluster = strmid(all[ii],0,11)
       thisout.galaxy = galaxy
       thisout.frac = frac
       thisout.apercor = apercor
       if (ii eq 0) then out = thisout else out = [out,thisout]
    endfor

    im_mwrfits, out, path+'kronapercorr.v2.3.fits'

return
end
    
