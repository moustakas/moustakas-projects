function ediscs_forage, speclist, datapath=datapath
; jm04jan15uofa

    if (n_elements(speclist) eq 0L) then begin
       print, 'Syntax - '
       return, -1
    endif

    if n_elements(datapath) eq 0L then datapath = cwd()

    pushd, datapath
    inspeclist = speclist
    speclist = file_search(speclist,count=nspec)

    keep = ['NAXIS','NAXIS1','OBNAME','TARGNAME','MJD_OBS','DATE_OBS','EXPTIME',$
      'RA','DEC','EQUINOX''AIRMSTAR','FWHMSTAR','PARANGEN','POSANG',$
      'PIXSCALE','MASKNAME','GRIS1NAM','GRIS1DIS','GRISWLE','FILT1ID','FILT1NAM','DETNAME',$
      'CD1_1','CDELT1','CRVAL1','CRPIX1','CTYPE1','OUT1GAIN','OUT1RON']

    for j = 0L, nspec-1L do begin

       h = headfits(speclist[j])
       s = im_hdr2struct(h)
       forage1 = struct_trimtags(s,select=keep)

; add tags here

       add = {file: speclist[j],id: 'EDCSJ'+repstr(strcompress(im_dec2hms(forage1.ra/15.0),/remove),'.','')+$
         repstr(strcompress(im_dec2hms(forage1.dec),/remove),'.','')}
       forage1 = struct_addtags(add,forage1)
       
       if (j eq 0L) then forage = forage1 else forage = [ [forage], [forage1] ]

    endfor
    forage = reform(forage)

    popd    
    
return, forage
end    
