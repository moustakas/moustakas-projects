function sc1120_forage, speclist, datapath=datapath
; jm05jan17uofa - read the first header to initialize the data
;                 structure, and then match each tag in each
;                 subsequent object; i wrote it this way because not
;                 all objects have the same number of header entries,
;                 so we are missing *some* information
    
    if (n_elements(speclist) eq 0L) then begin
       print, 'Syntax - '
       return, -1
    endif

    if n_elements(datapath) eq 0L then datapath = cwd()

    pushd, datapath
    inspeclist = speclist
    speclist = file_search(speclist,count=nspec)

    forage = {$
      specfile:  '', $
      naxis:     0L, $
      naxis1:    0L, $
      object:    '', $
      type:      '', $
      galaxy:    '', $
      date:      '', $
      mask:      '', $ 
      quadrant:  '', $ 
      slit:      '', $ 
      ra:        '', $
      dec:       '', $
      jd:      0.0D, $
      gain:     0.0, $
      rdnoise:  0.0, $
      airmass:  0.0, $
      exptime:  0.0, $
      crval1:   0.0, $
      crpix1:   0.0, $
      cd1_1:    0.0, $
      crval2:   0.0, $
      crpix2:   0.0, $
      cd2_2:    0.0, $
      z:        0.0, $
      z_err:    0.0, $
      q:         0L, $
      yobj:      0L, $
      ysky:      0L}
      
    forage = replicate(forage,nspec)

    tags = tag_names(forage[0])
    ntags = n_elements(tags)
    
    for i = 0L, nspec-1L do begin
       
       h = headfits(speclist[i])
       for j = 0L, ntags-1L do begin
          val = sxpar(h,tags[j],count=count)
          if count ne 0L then forage[i].(j) = val
       endfor

       forage[i].specfile = speclist[i]
       forage[i].date = strcompress(sxpar(h,'DATE-OBS'),/remove)
       forage[i].type = strcompress(sxpar(h,'IMAGETYP'),/remove)
       
    endfor

; ##############    
; OLD CODE BELOW    
; ##############    
    
;   h = headfits(speclist[0])
;   data = vlt_header_forage(h)
;   ntags = n_tags(data)
;   tags = tag_names(data)
;   
;   forage = struct_append(data,mrd_struct(tags,replicate("''",ntags),nspec-1L))
;
; now loop on every remaining object    
;   
;   for j = 1L, nspec-1L do begin
;
;      h = headfits(speclist[j])
;      data = vlt_header_forage(h)
;
;      for itag = 0L, ntags-1L do begin
;
;         match = where(tags[itag] eq tag_names(data),nmatch)
;         if (nmatch eq 1L) then forage[j].(itag) = data.(match)
;         
;      endfor
;      
;   endfor
;   forage = reform(forage)

    popd    
    
return, forage
end    
