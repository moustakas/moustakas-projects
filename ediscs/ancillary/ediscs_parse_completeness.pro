pro ediscs_parse_completeness
; jm08may19nyu - parse Bianca's completeness files; see the
;                emails from Bianca in PATH for details

    path = ediscs_path(/catalogs)+'completeness/'
    mycatpath = ediscs_path(/mycatalogs)
    alltab = file_search(path+'*.tab',count=ntab)

    result_template = {galaxy: '', ra: 0.0D, dec: 0.0D, $
      z: '', spec_weight: 0.0, geom_weight: 0.0}
    
    for itab = 0L, ntab-1L do begin

       splog, 'Reading '+alltab[itab]
       readcol, alltab[itab], idnew, wmag, wgeom, /silent, $
         ra, dec, z, format='A,F,F,F,F,A', comment='#'
       nobj = n_elements(idnew)
       result1 = replicate(result_template,nobj)

       result1.galaxy = idnew
       result1.ra = 15.0D*ra
       result1.dec = dec
       result1.z = z
       result1.spec_weight = wmag
       result1.geom_weight = wgeom
       
       if (n_elements(result) eq 0L) then result = result1 else $
         result = [result1,result]
       
    endfor

; write out    
    outfile = mycatpath+'ediscs_spec_completeness.fits'
    im_mwrfits, result, outfile, clobber=clobber
    
return
end
    
