function ediscs_headgrab, speclist, datapath=datapath
; jm04feb16uofa
    
    nspec = n_elements(speclist)
    if nspec eq 0L then begin
       print, 'Syntax - '
       return, -1
    endif

    if n_elements(datapath) eq 0L then datapath = cwd()
    
    keep = ['NAXIS1','GALAXY','DATE_OBS','CRVAL1','CD1_1','CRPIX1','CTYPE1',$
      'CLUSTER','MASK','SLIT','SUBSLIT','Z','IMAG','TYPE']

    speclist = file_search(datapath+speclist,count=nspec)
    for j = 0L, nspec-1L do begin

       h = headfits(speclist[j])
       s = im_hdr2struct(h)
       forage1 = reform(struct_trimtags(s,select=keep))

       if (j eq 0L) then forage = forage1 else forage = [ [forage], [forage1] ]

    endfor
    forage = reform(forage)
    
return, forage
end    
