;+
; NAME:
;   SINGS_VERSION()
;
; PURPOSE:
;   Keep track of various versions of the SINGS data.
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;   ancillary - K-corrections, etc.
;   ispec     - emission- and absorption-line measurements
;   templates - version of the templates used in ispec 
;
; OUTPUTS: 
;   The version number.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   ancillary
;     v1.0 - earlier effort
;     v2.0 - gas masses added, among many other changes
;   specfit
;     v1.0 - earlier effort
;     v2.0 - bug! (overestimated errors due to use of SYSERR); better
;            treatment of broad-line AGN; new 2D reductions
;            (distortion maps, better sky subtraction, slightly
;            modified error maps); more BC03 templates 
;     v3.0 - rerun with new version of ISPECLINEFIT() and BC03
;            templates with the new stellar mass scaling; even
;            better broad-line deblending
;   templates
;     v1.0 - 7 templates, Salpeter, 3 Z
;     v1.1 - 10 templates, Salpeter, 3 Z
;     v1.2 - minor tweaks to the resolution and info structure to
;            match the new ISPECLINEFIT code
;     v2.0 - new mass scaling
;
;   ppxf_templates 
;     v1.0 - see SINGS_GANDALF_SPECFIT
;     v2.0 - expanded wavelength range, for M. Brown's photoz
;       template project
;   ppxf_specfit
;     v1.0 - see SINGS_GANDALF_SPECFIT
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Sep 08, NYU - written, based on earlier
;     code 

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

function sings_version, ancillary=ancillary, specfit=specfit, $
  templates=templates, ppxf_templates=ppxf_templates, $
  ppxf_specfit=ppxf_specfit

    version = -1.0

    if keyword_set(ancillary) then version = 'v2.0'
    if keyword_set(specfit) then version = 'v3.0'
    if keyword_set(templates) then version = 'v2.0'

    if keyword_set(ppxf_templates) then version = 'v2.0' ; 'v1.0'
    if keyword_set(ppxf_specfit) then version = 'v1.0'

return, version
end
