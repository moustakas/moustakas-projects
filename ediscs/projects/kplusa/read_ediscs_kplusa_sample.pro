function read_ediscs_kplusa_sample, cluster1=cluster1, field1=field1
; jm08may18nyu - 

    common ediscs_kplusa, kplusa_cluster1, kplusa_field1

    datapath = ediscs_path(/projects)+'kplusa/'

; cluster sample 1
    if keyword_set(cluster1) then begin
       thisfile = datapath+'ediscs_kplusa_cluster1.fits'
       if (size(kplusa_cluster1,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          kplusa_cluster1 = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, kplusa_cluster1
    endif

; field sample 1
    if keyword_set(field1) then begin
       thisfile = datapath+'ediscs_kplusa_field1.fits'
       if (size(kplusa_field1,/type) ne 8L) then begin
          if (not keyword_set(silent)) then splog, 'Reading '+thisfile+'.'
          kplusa_field1 = mrdfits(thisfile,1,silent=0)
       endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)+'.'
       return, kplusa_field1
    endif

end
