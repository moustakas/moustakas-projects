;+
; NAME:
;   READ_NSA()
; PURPOSE:
;   Read the NSA catalog.
; MODIFICATION HISTORY:
;   J. Moustakas, 2012 Mar 06, UCSD 
;
; Copyright (C) 2012, John Moustakas
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

function read_nsa, wise=wise
    file = getenv('IM_DATA_DIR')+'/nsa/nsa_v0_1_2.fits.gz'
    nsa = mrdfits(file,1)
    if arg_present(wise) then begin
       file = getenv('IM_DATA_DIR')+'/nsa/nsa_v0_1_2_wise.fits.gz'
       wise = mrdfits(file,1)
    endif
return, nsa    
end
