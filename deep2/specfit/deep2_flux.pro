;+
; NAME:
;       DEEP2_FLUX
;
; PURPOSE:
;       Perform a rough flux calibration of the DEEP2 spectra.
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       Notes from N. Konidaris:
;
;         Converts e- counts to e- counts
;         The inverse of this function returns the throughput of DEIMOS
;       
;       Following factor accounts for (but deimos_correction does not account for):
;               Primary mirror size             pi*998cm^2
;               counts/hour -> counts/second    3600 sec/hour
;               Angstrom/Pixel                  0.33 angstrom/pix
;               Energy per e-                   2.00e-8 erg/ lambda(angstrom)
;               And slit corrections
;       
;       To convert to E, multiply by 5.3 (+-1.3) e-17 / lambda(angstrom)
;               erg/s/Ang/cm^2
;
;       W/out slit corrections, multiply by  2.16e-17 / lambda(angstrom)
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2007 Aug 04, NYU
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

pro deep2_flux, zcat1, debug=debug, write=write

    if (n_elements(zcat1) eq 0L) then zcat1 = read_deep2_zcat(/good)
    keep = where((zcat1.z gt 0.0) and (zcat1.zquality ge 3L))
    zcat = zcat1[keep]
    ngalaxy = n_elements(zcat)

    datapath = deep2_path(/dr2)

    if keyword_set(debug) then begin
       im_window, 0, xratio=0.6, yratio=0.6
       postthick1 = 2.0
       scale = 1D17
    endif

    for igal = 0L, ngalaxy-1L do begin

       print, format='("Fitting galaxy ",I0,"/",I0,".",A10,$)', igal+1, ngalaxy, string(13b)
       
       spec1 = mrdfits(datapath+strtrim(zcat[igal].file,2),1,/silent)
       spec2 = mrdfits(datapath+strtrim(zcat[igal].file,2),2,/silent)

       flux = [spec1.spec,spec2.spec]
       ivar = [spec1.ivar,spec2.ivar]
       wave = [spec1.lambda,spec2.lambda]

       good = where(ivar gt 0.0,ngood)
       wave = wave[good]
       flux1 = flux[good]
       ivar1 = ivar[good]
       ferr1 = 1.0/sqrt(ivar1)
       npix = n_elements(wave)

       flux = 5.3D-17 * flux1 / (-77.9026 + 0.0395916*wave - 7.49911e-06*wave^2 + $
         6.29692e-10*wave^3 - 1.97967e-14*wave^4) / wave

       
       
       
       z = zcat[igal].z

; debugging plot       

       if keyword_set(debug) then begin

          cstats = im_stats(scale*flux,sigrej=5.0)
          yrange = cstats.median_rej + [-5.0,8.0]*cstats.sigma_rej

          title = 'ID = '+strtrim(zcat[igal].zcatindx,2)+', '+$
            'z = '+strtrim(string(z,format='(F12.4)'),2)+', '+$
            'Q = '+string(zcat[igal].zquality,format='(I0)')
          
          djs_plot, wave, smooth(scale*flux,3), xsty=3, ysty=3, xthick=postthick1, ythick=postthick1, $
            charthick=postthick1, charsize=2.0, ps=10, yrange=yrange, title=title, $;, color='grey', $
            xtitle='Wavelength (\AA)', ytitle='Flux (10^{-17} '+flam_units()+')'
          cc = get_kbrd(1)

       endif 

       if keyword_set(write) then begin

          
          

       endif
       
    endfor 
    
return
end
    
    
