pro atlas1d_header_redshifts, update=update
; jm05jul22uofa - add redshifts and velocity dispersions to every
;                 object in ATLAS1D
; jm08apr14nyu - routine was re-written to use the INFO structure to
;                update the redshifts of the 1D spectra
    
    light = 2.99792458D5 ; speed of light [km/s]

; retrieve the atlas data table

    datapath = atlas_path(/atlas1d)
;   version = atlas_version(/ancillary)
;   atlas = mrdfits(atlas_path(/analysis)+'atlas_ancillary_data_'+version+'.fits.gz',1,/silent)
    atlas = atlas_read_info()
    natlas = n_elements(atlas)

    int = where(atlas.drift,nint)
    nuc = where(atlas.nuclear,nnuc)
    
    splog, 'Updating the redshifts for the INTEGRATED spectra.'
    for ii = 0L, nint-1L do begin

       specfile = strtrim(atlas[int[ii]].drift_file,2)
       h = headfits(datapath+specfile)
       struct_print, struct_trimtags(atlas[int[ii]],sel=['GALAXY',$
         'NED_GALAXY','DRIFT_FILE','Z']), /no_head
       sxaddpar, h, 'Z', float(atlas[int[ii]].z), ' NED redshift', before='HISTORY'
;      sxaddhist, "'NED redshift added "+im_today()+"'", h
       if keyword_set(update) then modfits, datapath+specfile, 0, h
       
    endfor

    splog, 'Updating the redshifts for the NUCLEAR spectra.'
    for ii = 0L, nnuc-1L do begin

       specfile = strtrim(atlas[nuc[ii]].nuclear_file,2)
       h = headfits(datapath+specfile)
       struct_print, struct_trimtags(atlas[nuc[ii]],sel=['GALAXY',$
         'NED_GALAXY','NUCLEAR_FILE','Z']), /no_head
       sxaddpar, h, 'Z', float(atlas[nuc[ii]].z), ' NED redshift', before='HISTORY'
;      sxaddhist, "'NED redshift added "+im_today()+"'", h
       if keyword_set(update) then modfits, datapath+specfile, 0, h
       
    endfor

;;    speclist = file_search(datapath+'*.ms.fits',count=nspec)
;;    speclist = file_basename(speclist)
;;
;;    len = string(max(strlen(strtrim(atlas.galaxy,2))),format='(I0)')
;;    for i = 0L, nspec-1L do begin
;;
;;       h = headfits(datapath+speclist[i])
;;       galaxy = strtrim(sxpar(h,'GALAXY'),2)
;;
;;       match = where(galaxy eq strtrim(atlas.galaxy,2),nmatch)
;;       if (nmatch eq 0L) then begin
;;          splog, 'WARNING! No match found for '+galaxy
;;          continue
;;       endif
;;       zobj = atlas[match].z
;;       
;;       print, speclist[i], galaxy, atlas[match].galaxy, $
;;         zobj, format='(A35,A'+len+',2x,A'+len+',F12.8)'
;;       
;;       sxaddpar, h, 'Z', float(zobj), ' NED redshift', before='HISTORY'
;;;      sxaddhist, "'NED redshift added "+im_today()+"'", h
;;
;;       if keyword_set(update) then modfits, datapath+speclist[i], 0, h
;;       
;;    endfor

return
end
