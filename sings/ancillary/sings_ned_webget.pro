;+
; NAME:
;       SINGS_NED_WEBGET
;
; PURPOSE:
;       Collect NED basic data, photometry, size and position angle
;       information from NED for the SINGS sample.
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007 Dec 24, NYU
;
; Copyright (C) 2007, John Moustakas
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

pro sings_ned_webget, basic=basic, diam=diam, photo=photo

    path = sings_path(/analysis)

    sings = read_sings_extra_info()
    galaxy = strtrim(sings.galaxy,2)
    ngalaxy = n_elements(galaxy)

    prefix = 'sings' ; 'test'
    
; basic data    
    
    if keyword_set(basic) then begin
       ned_webget_basic, galaxy, basic, $
         outfile=path+prefix+'_ned.fits', /write
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
    
