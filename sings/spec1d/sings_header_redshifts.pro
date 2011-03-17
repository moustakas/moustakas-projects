pro sings_header_redshifts, update=update
; jm05jul25uofa - add redshifts and velocity dispersions to every
;                 object in SINGS
; jm06jan28uofa - use SINGS_INFO
    
    light = 2.99792458D5 ; speed of light [km/s]

; retrieve the sings data table

    datapath = sings_path(/spec1d)
    sings = sings_read_info()

    nsings = n_elements(sings)

    speclist = file_search(datapath+'*.ms.fits',count=nspec)
    speclist = file_basename(speclist)

    len = string(max(strlen(strtrim(sings.galaxy,2))),format='(I0)')
    
    for i = 0L, nspec-1L do begin

       h = headfits(datapath+speclist[i])
       galaxy = strtrim(sxpar(h,'GALAXY'),2)

       match = where(galaxy eq strtrim(sings.galaxy,2),nmatch)
       if (nmatch eq 0L) then continue
       zobj = sings[match].z
       
       print, speclist[i], galaxy, sings[match].galaxy, $
         zobj, format='(A35,A'+len+',2x,A'+len+',F10.5)'
       
       sxaddpar, h, 'Z', float(zobj), ' NED redshift', before='HISTORY'
;      sxaddhist, "'NED redshift added "+im_today()+"'", h

       if keyword_set(update) then modfits, datapath+speclist[i], 0, h
       
    endfor

return
end
