;+
; NAME:
;   READ_09RIEKE()
;
; PURPOSE:
;   Read the Rieke et al. (2009) infrared models and put them into a 
;   standard format.
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;   local - read Table 3 of the local galaxy SEDs used to construct
;     the models and return
;
; OUTPUTS: 
;   data - structure containing the following fields:
;     lir - total infrared luminosity (8-1000 micron) [NMODEL, L_sun] 
;     l24 - Spitzer/24-micron luminosity [NMODEL, L_sun]
;     wave - wavelength vector [NPIX, A]
;     flux - flux density [NPIX,NMODEL, erg/s/A]
;
; COMMENTS:
;   The models are read from the electronic Table 4 in the Rieke et
;   al. paper.
;
;   The models are oversampled by a factor of 5 to ensure an accurate
;   convolution with the 24-micron filter.
; 
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 17, UCSD
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

function qpint1d_func, x, wave=wave, flux=flux
return, interpol(flux,wave,x)
end
 
function read_09rieke, local=local
    
    path = getenv('CATALOGS_DIR')+'/09rieke/'
    if (file_test(path,/dir) eq 0) then begin
       splog, 'Data path '+path+' does not exist!'
       return, -1
    endif

; optionally return the local SEDs and return    
    if keyword_set(local) then begin
       file = 'table3_localseds.txt' 
       data = im_read_fmr(path+file)
       return, data
    endif

; read the models and pack them into a more convenient format; store
; the nominal L(IR) and invert eq. (A6) to get L(24)
    temp1 = im_read_fmr(path+'table4_avgtemplates.txt')
    npix = n_elements(temp1.wavelength)
    nmodel = n_tags(temp1)-1
    temp = {lir: fltarr(nmodel), l24: fltarr(nmodel), $
      wave: temp1.wavelength, flux: dblarr(npix,nmodel)}
    for ii = 0, nmodel-1 do temp.flux[*,ii] = temp1.(ii+1)
    tags = (tag_names(temp1))[1:nmodel]
    temp.lir = strmid(tags,6)/100.0
    temp.l24 = (temp.lir-1.445)/0.945
    
; now put the models into a standard format; the models are in Jy at
; an assumed distance of 10 Mpc
    filt24 = 'spitzer_mips_24.par'
    weff = (k_lambda_eff(filterlist=filt24))[0]
    lsun = 3.826D33
    light = 2.9979246D18
    dist = 3.085678D18*10.0*1D6 ; [10 Mpc]

; oversample the models, otherwise the filter convolution is fubar     
    nsamp = 5
    owave = temp.wave*1D4 ; [A]
    wave = range(min(owave),max(owave),npix*nsamp,/log) ; [A]
    
    data = {$
      lir:  dblarr(nmodel),$
      l24:  dblarr(nmodel),$
      wave:           wave,$
      flux: dblarr(npix*nsamp,nmodel)}
    
; convert the models from [Jy] to [erg/s/A]
    wave_edges = k_lambda_to_edges(wave)
    for ii = 0, nmodel-1 do begin
       fnu = temp.flux[*,ii]*1D-23*4.0*!dpi*dist^2 ; [erg/s/Hz]
       flam = fnu*(light/owave^2)                  ; [erg/s/A]
       data.flux[*,ii] = interpol(flam,owave,wave)
; integrate to get L(IR)=L(8-1000)
       data.lir[ii] = qpint1d('qpint1d_func',8D4,1000D4,$
         functargs={wave:wave,flux:data.flux[*,ii]})/lsun
; convolve the 24-micron filter curve to get L(24); we need the factor
; of 10^40 to trick k_project_filters() from running into numerical
; problems
       f24 = k_project_filters(wave_edges,$
         data.flux[*,ii]/1D40,filterlist=filt24)
       data.l24[ii] = f24*10.0^(-0.4*48.6)*1D40*(light/weff^2)*weff/lsun
    endfor

return, data
end
