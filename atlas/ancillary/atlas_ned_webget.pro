;+
; NAME:
;       ATLAS_NED_WEBGET
;
; PURPOSE:
;       Collect NED basic data, photometry, size and position angle
;       information from NED for the ATLAS sample.
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2008 Jan 30, NYU
;
; Copyright (C) 2008, John Moustakas
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

pro atlas_ned_webget, basic=basic, diam=diam, photo=photo

    path = atlas_path(/analysis)

    atlas = read_atlas_extra_info()
    galaxy = strtrim(atlas.nedgalaxy,2)
    ngalaxy = n_elements(galaxy)

;   prefix = 'test'
    prefix = 'atlas'
    
; basic data    
    
    if keyword_set(basic) then begin
       ned_webget_basic, galaxy, basic, $
         outfile=path+prefix+'_ned.fits', /write
;      struct_print, basic
    endif

; photometry
    
    if keyword_set(photo) then begin
       ned_webget_photometry, galaxy, photo, $
         outfile=path+prefix+'_ned_photo.fits', /write
    endif

; diameters
    
    if keyword_set(diam) then begin
       ned_webget_diameters, galaxy, diam, $
         outfile=path+prefix+'_diameters.fits', /write
    endif

return
end
    
