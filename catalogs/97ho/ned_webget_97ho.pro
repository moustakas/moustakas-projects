;+
; NAME:
;       NED_WEBGET_97HO
;
; PURPOSE:
;       Collect NED basic data, photometry, size and position angle
;       information from NED.
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2008 Feb 01, NYU
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

pro ned_webget_97ho

    path = getenv('CATALOGS_DIR')+'/97ho/'

    readcol, path+'galaxy_list.txt', galaxy, format='A', /silent
    galaxy = strcompress(galaxy,/remove)
    ngalaxy = n_elements(galaxy)

    prefix = '97ho' ; 'test'
    
; basic data    
    
    ned_webget_basic, galaxy, basic, outfile=path+prefix+'_ned.fits', /write
;   struct_print, basic

return
end
