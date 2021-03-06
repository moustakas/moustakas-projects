function casjobs_remove_double, cat
; most of the fields in the casjobs catalog are unnecessarily double
; precision; change them to float
    tags = strlowcase(tag_names(cat))
    ntags = n_elements(tags)
    for ii = 0, ntags-1 do begin
       type = size(cat[0].(ii),/type)
       if (type eq 5) then begin
          if (tags[ii] eq 'alpha_j2000') or (tags[ii] eq 'delta_j2000') or $
            (tags[ii] eq 'ra') or (tags[ii] eq 'dec') then type = 5 else type = 4
       endif
       if (type eq 3) then begin
          if (tags[ii] eq 'ages_id') or (tags[ii] eq 'tile') or $
            (tags[ii] eq 'objid') then type = 3 else type = 2
       endif
       arr = (make_array(1,type=type))[0]
;      print, tags[ii], type, arr
       if (ii eq 0) then outcat = create_struct(tags[ii],arr) else $
         outcat = create_struct(outcat,tags[ii],arr)
    endfor
    outcat = replicate(outcat,n_elements(cat))
    for ii = 0, ntags-1 do outcat.(ii) = cat.(ii)
return, outcat
end

pro build_ages_galex, out_galex, query=query, gr=gr, clobber=clobber
; jm10apr30ucsd - build a line-matched GALEX catalog for the AGES
;   sample using the CasJobs output (see the README in the
;   mycatalogs/galex directory) 
; jm10jul23ucsd - updated to GR6
; jm13aug28siena - updated to GR6/7 + get more tags; what used to be
;   WRITE_AGES_GALEX_CASJOBS_INPUT is also part of this routine now;
;   see the README in the REDMAPPER directory for details

    if (n_elements(gr) eq 0) then gr = 'gr67'

    catpath = ages_path(/catalogs)
    path = ages_path(/mycatalogs)
    ages = mrdfits(catpath+'catalog.cat.noguidestars.fits.gz',1)
    ages.ra = ages.ra*15.0D
    ngal = n_elements(ages)
    ages_id = lindgen(ngal)

; build the input catalog
    if keyword_set(query) then begin
       out = struct_addtags(replicate({ages_id: 0L},ngal),$
         struct_trimtags(ages,select=['ra','dec']))
       out.ages_id = ages_id
       
       outfile = '~/tmp/ages_galex_'+gr+'_casjobs.dat'
       openw, lun, outfile, /get_lun
       printf, lun, '# ages_id ra dec'
       struct_print, out, lun=lun, ddigit=12, /no_head
       free_lun, lun
       return
    endif 
    
; output filename    
    outfile = path+'ages_galex_'+gr+'.fits'
    if file_test(outfile+'.gz') and (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+outfile+' exists; use /CLOBBER'
       return
    endif
    
; read the casjobs output; for some reason some input/output objects
; are repeated (not sure why), so remove them here
    infile = '~/tmp/ages_galex_'+gr+'_casjobs.fits'
    if file_test(infile) eq 0 then begin
       splog, 'Input file '+infile+' not found!'
       return
    endif
    incat = mrdfits(infile,1)
    incat = casjobs_remove_double(incat)
    incat = struct_addtags(replicate({galex_object_position: -999L},$
      n_elements(incat)),temporary(incat))
    incat.galex_object_position = lindgen(n_elements(incat))

    allid = strtrim(incat.ages_id,2)+strtrim(incat.objid,2)
    keep = uniq(allid,sort(allid))
    splog, 'Removing '+strtrim(n_elements(incat)-n_elements(keep),2)+$
      '/'+strtrim(n_elements(incat),2)+' repeated casjobs entries!'
    incat = incat[keep]

    outcat1 = galex_resolve_casjobs_duplicates(incat,idtag='ages_id')
    best = where(outcat1.isbest,nbest)
    outcat = outcat1[best]
    splog, 'Found '+strtrim(nbest,2)+'/'+strtrim(ngal,2)+$
      ' galaxies with good GALEX detections!'

    match, ages_id, outcat.ages_id, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    out_galex = im_empty_structure(outcat[0],$
      empty_value=-999.0,ncopies=ngal)
    out_galex[m1] = outcat[m2]

; add AGES_ID, RA, DEC from the original AGES structure to the output
; catalog 
    out_galex.ages_id = ages_id
    out_galex.ra = ages.ra
    out_galex.dec = ages.dec
    im_mwrfits, out_galex, outfile, /clobber

return
end
    
