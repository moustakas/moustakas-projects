;+
; NAME:
;   AGES_TARGET_WEIGHT()
;
; PURPOSE:
;   Return the AGES sparse-sampling (target) weight. 
;
; INPUTS: 
;   codes - structure containing the AGES targeting codes 
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   weight - sparse-sampling weight (TARGET_WEIGHT in Daniel's
;     catalog.spectweight) 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   From the 2005 documentation:
;
;  One can think of the sample as I<20 in which objects with more
;  interesting colors have been sparse-sampled:
;
;    if qshort > 0 then weight = 1
;    if qshort = 0 and gshort = 2048, weight = 0
;        [Notes: gshort = 2048 consists of all galaxies with I<=20
;         that were not included in the MAIN sample but were used as
;         filler. So they shouldn't be used for a properly weighted
;         sample.]
;    if qshort = 0 and gshort != 2048 and gbright > 0 then weight = 1
;        [Notes: the main sample targets all bright galaxies and then
;         a random set of fainter targets, so all gbright targets
;         should be set.]
;    if qshort = 0 and gshort != 2048 and gbright=0 and (grand&127) then
;                the weight is 1/0.3
;    else weight = 1/0.2
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 01, UCSD - based on code by R. Cool
;
; Copyright (C) 2009, John Moustakas
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

function ages_target_weight, codes

    weight = fltarr(n_elements(codes))

; set the default value
    weight[*] = 1./0.2
    
; now check for quasars
    kqshort = where(codes.qshort gt 0, ct)
    if (ct gt 0) then weight[kqshort] = 1.0

; check for fillers
    kfiller = where(codes.qshort eq 0 and $
      codes.gshort eq 2048,ct)
    if (ct gt 0) then weight[kfiller] = 0.0
    
; check for the bright samples
    kbright = where(codes.qshort eq 0 and $
      codes.gshort ne 2048 and $
      codes.gbright gt 0,ct)
    if (ct gt 0) then weight[kbright] = 1.0

; finally give weight to the random sample
    krand = where(codes.qshort eq 0 and $
      codes.gshort ne 2048 and $
      codes.gbright eq 0 and $
      ((codes.grand and 127) gt 0),ct)
    if (ct gt 0) then weight[krand] = 1.0/0.3

return, weight
end
