;+
; NAME:
;   EDISCS_SFH_SMOOTH
;
; PURPOSE:
;   Smooth an EDisCS spectrum. 
;
; INPUTS: 
;   wave - input wavelength vector; assumed to have constant pixels in
;     km/s unless /LINEAR [NPIX]
;   flux - corresponding flux vector [NPIX]
;
; OPTIONAL INPUTS: 
;   ferr - error spectrum corresponding to FLUX [NPIX]
;   vdisp - galaxy velocity dispersion (default 0.0 km/s)
;   inst_vdisp - instrumental velocity dispersion (default 0.0 km/s)
;   final_vdisp - final velocity dispersion (default 350 km/s)
;
; KEYWORD PARAMETERS: 
;   linear - if this keyword is set then WAVE is assumed to have 
;     constant dispersion in *Angstroms* 
;
; OUTPUTS: 
;   smoothed_flux - smoothed version of FLUX
;
; OPTIONAL OUTPUTS:
;   smoothed_wave - if /LINEAR then this array contains the
;     interpolated wavelength vector corresponding to SMOOTHED_FLUX 
;   smoothed_ferr - smoothed version of FERR
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 May 07, UCSD
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

pro ediscs_sfh_smooth, wave, flux, smoothed_flux, ferr=ferr, $
  smoothed_wave=lnwave, smoothed_ferr=smoothed_ferr, linear=linear, $
  velscale=velscale, vdisp=vdisp, inst_vdisp=inst_vdisp, $
  final_vdisp=final_vdisp, debug=debug

    npix = n_elements(wave)

; if necessary, rebin the spectrum to have constant pixels in km/s
; (i.e., constant dispersion in ln-wavelength)
    if keyword_set(linear) then begin
       log_rebin, minmax(wave), flux, lnflux, $
         lnwave, velscale=velscale
       if (n_elements(ferr) ne 0) then begin
          log_rebin, minmax(wave), ferr^2, lnvar, $
            lnwave, velscale=velscale
          lnferr = sqrt(lnvar>0)
       endif
    endif else begin
       lnflux = flux
       lnwave = wave
       if (n_elements(ferr) ne 0) then lnferr = ferr
    endelse

    if (n_elements(vdisp) eq 0) then vdisp = 0.0 ; [km/s]
    if (n_elements(inst_vdisp) eq 0) then inst_vdisp = 0.0 ; [km/s]
    if (n_elements(final_vdisp) eq 0) then final_vdisp = 350.0 ; [km/s]
    vkernel = sqrt((final_vdisp^2-vdisp^2-inst_vdisp^2)>0) ; [km/s]

    if (vkernel gt 0.0) then begin
       smoothing = vkernel/velscale ; pixels
       npix = long(4.0*ceil(smoothing))*2L+3
       klam = findgen(npix)-float(npix-1.0)/2.0
       kernel = exp(-0.5*(klam/smoothing)^2)/sqrt(2.0*!dpi)/smoothing
       kernel = kernel/total(kernel)
       smoothed_flux = convol(lnflux,kernel,/edge_truncate)
       if (n_elements(ferr) ne 0) then begin
          smoothed_var = convol(lnferr^2,kernel,/edge_truncate)
          smoothed_ferr = sqrt(smoothed_var>0)
       endif
    endif else begin
       smoothed_flux = lnflux
       if (n_elements(ferr) ne 0) then smoothed_ferr = lnferr
    endelse
    
return
end
