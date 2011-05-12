pro clash_merge_cat, out
; jm11apr11ucsd - merge the BCG and arc catalogs

    catpath = clash_path(/cat)
    filterlist = clash_filterlist(short_filter=filt) ; sort

    arc = rsex(catpath+'arc_a383.cat')
;   bcg = rsex(catpath+'bcg_a383.cat')
;   all = [bcg,arc]
    all = arc
    nobj = n_elements(all)

    out = struct_trimtags(all,select=['id','ra','dec',filt+'_mag*'])
    out = struct_addtags(out,replicate({galaxy: '', z: 0.0},nobj))
    out.galaxy = 'A383 Big Arc'
    out.z = 1.01
;   out.z = [0.189,1.01]
    
    outfile = catpath+'a383.fits'
    im_mwrfits, out, outfile, /clobber

;    bcg = im_empty_structure(arc,empty_value=-999)
;    out = [bcg,arc]
;
;    out[0].zb = 0.189 ; BCG redshift
;    out[0].ra = out[1].ra
;    out[0].dec = out[1].dec
;
;; pack in the BCG photometry    
;    allcat = file_search(catpath+'a383_*_bcg.cat',count=ncat)
;    root = file_basename(allcat)
;    rootfilt = repstr(repstr(root,'a383_',''),'_bcg.cat')
;    for ii = 0, ncat-1 do begin
;       cat1 = struct_trimtags(rsex(allcat[ii]),except=['fwhm_world'])
;       if (ii eq 0) then cat = cat1 else cat = [cat,cat1]
;    endfor
;    
;    filterlist = clash_filterlist(short_filter=filt) ; sort
;    match, rootfilt, filt, m1, m2
;    srt = sort(m2) & m1 = m1[srt] & m2 = m2[srt]
;    cat = cat[m1]
;    nfilt = n_elements(filt)
;
;    for ii = 0, nfilt-1 do begin
;       mag = tag_indx(out,filt[ii]+'_mag')
;
;stop       
;    endfor
;
;    out = struct_trimtags(cat,select=['mag*_auto'])
;
;    outfile = catpath+'a383_bcg.cat'
;    im_mwrfits, out, outfile, /clobber
    
return
end
    
