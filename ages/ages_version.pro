;+
; NAME:
;   AGES_VERSION()
;
; PURPOSE:
;   Keep track of various versions of the higher-order AGES data
;   products. 
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;   ancillary        - K-corrections, etc.
;   ispec_specfit    - emission- and absorption-line measurements
;                      derived using ispec 
;   unfluxed_specfit - emission-line measurements obtained from
;                      the unfluxed spectra
;   templates        - version of the templates used in ispec 
;
; OUTPUTS: 
;   The version number.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   ancillary ### OBSOLETE ###
;     v1.0 - earlier effort
;     v2.0 - mainly updated k-corrections
;     v2.1 - bug fix in k-correct masses
;     v2.2 - new k-correct code; removed synthesized magnitudes 
;     v2.3 - added original targeting codes
;     v2.4 - added R. Assef's and Hickox's AGN classification
;            results
;     v2.5 - added z-band data from Cool+07 and Vmax values from
;            Eisenstein+09
;     v2.6 - 
;     v3.0 - new band-merged catalogs; lots of other small changes 
;     v3.1 - updated to v2.1 photometry
;
;   kcorrect
;     v2.2 - earlier effort (branched from ANCILLARY)
;     v2.3 - now using BwRIzK and AGES_TO_MAGGIES 
;     v2.4 - additional tests with K-band and FLAMEX JKs
;     v2.5 - used the v1.0 photometry catalog; only fit GSHORT>0,
;            0.001<z<1.0 galaxies with good BwRI photometry
;     v3.0 - use the v2.0 photometry
;     v3.1 - use the v2.1 photometry
;     v3.2 - use the v2.2 photometry
;     v4.0 - use the v3.0 photometry
;     v5.0 - use the v4.0 photometry
; 
;   photometry 
;     v1.0 - new NDWFS zeropoints
;     v2.0 - BOOTES aperture-matched BwRIJHKsch[1-4] 
;     v2.1 - updated to GALEX/GR45 photometry and zeropoints
;     v2.2 - updated to GALEX/GR6
;     v3.0 - 2010b version of M. Brown's aperture-matched
;       catalogs; now includes U- and z-band photometry, improved PSF
;       magnitudes, background subtraction, and zeropoints
;     v3.1 - removed IMAFLAGS and SEGFLAGS tags and some of the
;       unnecessary aperture magnitudes
;     v4.0 - updated to 2011a Bootes photometry
;
;   ispec_specfit
;     v1.0 - earlier effort
;     v1.1 - the BC03 templates have now been normalized to the stellar
;            mass at the given age; hence the coefficients yield the
;            stellar mass [see AGES_RESTORE_BC03()]
;     v2.0 - better spectrophotometric tweaking and masking of the
;            red leak; do not fit Mg II; fit [OII] as a singlet
;
;   unfluxed_specfit
;     v1.0 - first effort
;     v2.0 - between smooth continuum modeling
;
;   templates
;     v1.0 - 7 templates, Salpeter, 3 Z
;     v1.1 - 10 templates, Salpeter, 3 Z
;     v1.2 - minor tweaks to the resolution and info structure to
;            match the new ISPECLINEFIT code
;
;   ppxf_templates 
;     v1.0 - see AGES_GANDALF_SPECFIT
;   ppxf_specfit
;     v1.0 - see AGES_GANDALF_SPECFIT
;     v2.0 - 
;     v2.1 - fix LOG_REBIN bug
;   ppxf_ancillary 
;     v1.0 - see BUILD_AGES_GANDALF_ANCILLARY
;     v2.0 - updated to use the v2.0 photometry and v3.0 K-corrections
;     v2.1 - updated to v2.1 photometry
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Mar 24, NYU - written
;   jm08apr02nyu - added TEMPLATES keyword; documented 
;   jm08jun18nyu - ancillary version pushed to v2.1
;   jm08sep02nyu - KCORRECT and ANCILLARY versions branched  
;   jm09nov09ucsd - new PHOTOMETRY version number

; Copyright (C) 2008-2009, John Moustakas
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

function ages_version, ancillary=ancillary, photometry=photometry, $
  ispec_specfit=ispec_specfit, unfluxed_specfit=unfluxed_specfit, $
  kcorrect=kcorrect, templates=templates, ppxf_templates=ppxf_templates, $
  ppxf_specfit=ppxf_specfit, unfluxed_ppxf_specfit=unfluxed_ppxf_specfit, $
  ppxf_ancillary=ppxf_ancillary

    version = -1.0

    if keyword_set(photometry) then version = 'v4.0'
    if keyword_set(kcorrect) then version = 'v5.0'

    if keyword_set(ispec_specfit) then version = 'v2.1'
    if keyword_set(templates) then version = 'v1.2'
    if keyword_set(unfluxed_specfit) then version = 'v2.1'

;   if keyword_set(ppxf_templates) then version = 'v1.0'
;   if keyword_set(ppxf_specfit) then version = 'v1.0'
    if keyword_set(ppxf_templates) then version = 'v2.0'
    if keyword_set(ppxf_specfit) then version = 'v2.1'

return, version
end
