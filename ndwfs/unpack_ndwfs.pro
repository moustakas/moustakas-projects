pro unpack_ndwfs
; jm09aug25ucsd - unpack the field-by-field NDWFS catalogs make a
;   merged, trimmed-down set of catalogs that are easier to manage

    indir = getenv('DATA_DIR')+'/data/ndwfs/'
;   indir = getenv('RESEARCHPATH')+'/data/ndwfs/'

    bands = ['Bw','R','I']
;   bands = ['Bw','R','I','K']
    field1 = ['32','33','34','35']
    field2 = ['33','34','35','36']

    tags = ['ndwfs_field','ndwfs_name','alpha_j2000','delta_j2000','flux_radius',$
      'class_star','mag_auto','magerr_auto','mag_aper*','magerr_aper*',$
      'flags','imaflags_iso','flag_duplicate','flag_splitmatch']
    
    for ii = 0, n_elements(field1)-1 do begin
       for jj = 0, n_elements(bands)-1 do begin
          file = indir+'NDWFS_'+bands[jj]+'_'+$
            field1[ii]+'_'+field2[ii]+'_cat_m.fits.gz'
          splog, 'Reading '+file
          cat1 = mrdfits(file,1);,range=[0,50])
          cat1 = struct_trimtags(cat1,select=tags)
          cat1 = im_struct_trimtags(cat1,select=tag_names(cat1),$
            newtags=bands[jj]+'_'+tag_names(cat1[0]))
          if (jj eq 0) then cat = cat1 else cat = $
            struct_addtags(temporary(cat),temporary(cat1))
       endfor
; only keep areas with photometry in Bw, R, *and* I
;      keep = where($
;        (cat.bw_mag_auto gt 0.0) and (cat.bw_mag_auto lt 90.0) and $
;        (cat.r_mag_auto gt 0.0) and (cat.r_mag_auto lt 90.0) and $
;        (cat.i_mag_auto gt 0.0) and (cat.i_mag_auto lt 90.0),nkeep)
       nkeep = n_elements(cat)
       keep = lindgen(nkeep)
       splog, 'NKEEP ', nkeep
       outfile = indir+'NDWFS_'+field1[ii]+'_'+field2[ii]+'.fits'
       im_mwrfits, cat[keep], outfile, /clobber

; also make an ubercatalog
       if (n_elements(bigcat) eq 0L) then $
         bigcat = temporary(cat) else $
         bigcat = [temporary(bigcat),temporary(cat)]
    endfor

; derive mean coordinates for each object in descending order of
; preference: I, R, Bw
    moretags = replicate({ra: -999.0D, dec: -999.0D},n_elements(bigcat))
    igood = where((moretags.ra lt 0.0) and (bigcat.i_alpha_j2000 lt 900.0))
    moretags[igood].ra = bigcat[igood].i_alpha_j2000
    moretags[igood].dec = bigcat[igood].i_delta_j2000

    rgood = where((moretags.ra lt 0.0) and (bigcat.r_alpha_j2000 lt 900.0))
    moretags[rgood].ra = bigcat[rgood].r_alpha_j2000
    moretags[rgood].dec = bigcat[rgood].r_delta_j2000
    
    bgood = where((moretags.ra lt 0.0) and (bigcat.bw_alpha_j2000 lt 900.0))
    moretags[bgood].ra = bigcat[bgood].bw_alpha_j2000
    moretags[bgood].dec = bigcat[bgood].bw_delta_j2000

    keep = where(moretags.ra gt 0.0,nkeep)
    splog, 'NKEEP ', nkeep
    bigcat = struct_addtags(moretags[keep],temporary(bigcat[keep]))
    
    im_mwrfits, bigcat, indir+'NDWFS_DR3.fits', /clobber

; remove duplicates and write out a new file
    ing = spheregroup(bigcat.ra,bigcat.dec,1.0/3600.0,$
      mult=mult,first=first,next=next)
    ngroups = max(ing)+1L
    cat = cat1[first[0:ngroups-1L]]
    im_mwrfits, cat, indir+'NDWFS_DR3.noduplicates.fits', /clobber
    
stop    

return    
end
