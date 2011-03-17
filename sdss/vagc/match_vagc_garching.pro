;+
; NAME:
;   MATCH_VAGC_GARCHING
;
; PURPOSE:
;   Match the Garching sample to the VAGC.
;
; INPUTS: 
;   None required.
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; COMMENTS:
;   This code reads the VAGC/object_sdss_imaging.fits file and matches
;   it to the *primary* (i.e., highest S/N; see BUILD_GARCHING_CATALOG)
;   objects in the Garching catalog.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Mar 06, UCSD
;
; Copyright (C) 2010, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro match_vagc_garching

    common com_match_vagc_garching, mpa, vagc
    
; read the garching catalog
    if (n_elements(mpa) eq 0) then begin
       mpafile = getenv('VAGC_REDUX')+'/garching_dr7/garching_dr7_catalog.fits.gz'
       splog, 'Reading '+mpafile
       mpa = hogg_mrdfits(mpafile,1,nrowchunk=20000)
    endif

; read the full SDSS imaging structure    
    if (n_elements(vagc) eq 0) then begin
       vagcfile = getenv('VAGC_REDUX')+'/object_sdss_imaging.fits.gz'
       splog, 'Reading '+vagcfile
       vagc = hogg_mrdfits(vagcfile,1,nrowchunk=20000)
    endif

    outfile = getenv('VAGC_REDUX')+'/garching_dr7/vagc_garching_dr7_catalog.fits'

; output structure    
    out = struct_addtags(replicate({object_position: -999L, $
      mpa_object_position: -999L},n_elements(mpa)),mpa)
    out.mpa_object_position = lindgen(n_elements(mpa))

; spherematch the VAGC *just* to the primary objects (i.e.,
; we will only keep track of how the *primary* objects match to the
; VAGC) 
    primary = where(mpa.garching_dr7_tag eq mpa.garching_dr7_tag_primary)
    spherematch, vagc.ra, vagc.dec, mpa[primary].ra, $
      mpa[primary].dec, 1.0/3600.0, m1, m2
    out[primary[m2]].object_position = m1

stop    
    
; write out the full garching catalog; to get the subset of *unique*
; matches to the VAGC you have to cut on OBJECT_POSITION NE -999
    im_mwrfits, out, outfile, /clobber

return
end
    
