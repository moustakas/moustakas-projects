pro sc1120_write_catalogs, sc1120, write=write
; jm05feb01uofa
; generate an ASCII catalog for Vy, for the Spitzer proposal

    datapath = sc1120_path(/analysis)

    if (n_elements(sc1120) eq 0L) then sc1120 = read_sc1120()
    ngalaxy = n_elements(sc1120)
    
    tags = ['GALAXY','Z_OBJ','OII_3727','H_BETA','H_ALPHA',$
      'OII_3727_EW','H_BETA_EW','H_ALPHA_EW','LICK_HD_A_MODEL']
    format = ['A0','F7.5','2F10.5','2F8.3','2F10.5','2F8.3','2F10.5','2F8.3','2F8.3']    
    cat = struct_trimtags(sc1120,select=tags)

    scale = 1E17
    for i = 2L, 4L do begin
       oldflux = cat.(i)
       newflux = oldflux
       good = where(oldflux[1,*] gt 0.0,ngood)
       if (ngood ne 0L) then begin
          newflux[*,good] = scale*oldflux[*,good]
          cat.(i) = newflux
       endif
    endfor

    outcat = struct_trimtags(cat,select=tag_names(cat),format=format)

    if keyword_set(write) then begin
       openw, lun, datapath+'sc1120_flux_05feb01.dat', /get_lun
       struct_print, outcat, lun=lun
       free_lun, lun
    endif

stop    

return
end
    
