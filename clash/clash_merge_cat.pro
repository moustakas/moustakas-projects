pro clash_merge_cat
; jm11apr11ucsd - merge the catalogs

    catpath = clash_path(/cat)
    allcat = file_search(catpath+'*.cat',count=ncat)

    root = file_basename(allcat)
    rootfilt = repstr(repstr(root,'a383_',''),'_bcg.cat')

    for ii = 0, ncat-1 do begin
       cat1 = struct_trimtags(rsex(allcat[ii]),except=['fwhm_world'])
       if (ii eq 0) then cat = cat1 else cat = [cat,cat1]
    endfor

; sort
    filterlist = clash_filterlist(short_filter=filt)
    match, rootfilt, filt, m1, m2
    srt = sort(m2) & m1 = m1[srt] & m2 = m2[srt]
    
    cat = cat[m1]

    outfile = catpath+'a383_bcg.cat'
    out = 
    
stop    
    
return
end
    
