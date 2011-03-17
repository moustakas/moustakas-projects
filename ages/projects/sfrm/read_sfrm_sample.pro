function read_sfrm_sample, sdss=sdss, silent=silent
; jm10feb04ucsd - read the output from BUILD_SFRM_SAMPLE

    common sfrm_parent, ages_parent, sdss_parent

    datapath = ages_path(/projects)+'sfrm/'
    if keyword_set(sdss) then begin ; SDSS
       thisfile = datapath+'sdss_sfrm.fits.gz'
       if (size(sdss_parent,/type) ne 8L) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          sdss_parent = mrdfits(thisfile,1,silent=0)
; quality cuts: QUALITY=0 and NUV magnitude and at least *one* of JHK
;         good = where((sdss_parent.quality eq 0))
;         good = where((sdss_parent.quality eq 0) and $
;           (sdss_parent.maggies[1] gt 0) and $
;           (total(sdss_parent.maggies[7:9] gt 0.0,1) ge 1))
;         sdss_parent = sdss_parent[good]
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, sdss_parent
    endif else begin ; AGES
       thisfile = datapath+'ages_sfrm.fits.gz'
       if (size(ages_parent,/type) ne 8L) then begin
          if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
          ages_parent = mrdfits(thisfile,1,silent=0)
       endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
       return, ages_parent
    endelse

 end
