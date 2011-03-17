;+
; NAME:
;   READ_02DALE()
;
; PURPOSE:
;   Read the Dale & Helou (2002) infrared models and put them into a
;   standard format.
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   data - structure containing the following fields:
;     lir - total infrared luminosity (8-1000 micron) [NMODEL, L_sun] 
;     l24 - Spitzer/24-micron luminosity [NMODEL, L_sun]
;     wave - wavelength vector [NPIX, A]
;     flux - flux density [NPIX,NMODEL, erg/s/A]
;
; COMMENTS:
;   The version of the Dale & Helou (2002) models used here were
;   downloaded from http://www.its.caltech.edu/~rchary and were
;   converted into a luminosity-dependent form by Ranga-Ram:
;
;	from	Ranga-Ram Chary <rchary@caltech.edu>
;               Fri, Feb 19, 2010 at 11:17 AM
;
;     "DH02 give SEDs as a function of alpha which is related to the
;     60/100 micron colors. However, the 60/100 micron color is
;     related to LIR among the IRAS BGS galaxies. So I went from alpha
;     to 60/100 micron to LIR to do the conversion since I need
;     luminosity dependent templates."
; 
;   The models are oversampled by a factor of 10 to ensure an accurate
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
 
function read_02dale
    
    path = getenv('CATALOGS_DIR')+'/02dale/'
    if (file_test(path,/dir) eq 0) then begin
       splog, 'Data path '+path+' does not exist!'
       return, -1
    endif

; restore Ranga-Ram's IDL save set; my values of L(IR) match his
; reasonably well, although it's not clear what wavelength
; interval he used (I use 8-1000 micron), so the discrepancies are
; luminosity-dependent 
    restore, path+'Edited_Dale_Spec.save'

    filt24 = 'spitzer_mips_24.par'
    weff = (k_lambda_eff(filterlist=filt24))[0]
    lsun = 3.826D33
    light = 2.9979246D18

    sz = size(nulnuinlsun,/dim)
    nmodel = sz[0]
    npix = sz[1]

; oversample the models, otherwise the filter convolution is fubar     
    nsamp = 10
    owave = lambda*1D4 ; [A]
    wave = range(min(owave),max(owave),npix*nsamp,/log) ; [A]
    
    data = {$
      lir:  dblarr(nmodel),$
      l24:  dblarr(nmodel),$
      wave:           wave,$
      flux: dblarr(npix*nsamp,nmodel)}
    
; convert the models from [nu*L_nu/L_sun] to [erg/s/A]
    wave_edges = k_lambda_to_edges(wave)
    for ii = 0, nmodel-1 do begin
       fnu = lsun*nulnuinlsun[ii,*]*(owave/light) ; [erg/s/Hz]
       flam = fnu*(light/owave^2)                 ; [erg/s/A]
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

;; old code
;    root = '02dale'
;    path = getenv('CATALOGS_DIR')+'/'+root+'/'
;    file = root+'.fits'
;    dale = mrdfits(path+file,1,/silent)
;    out = {$
;      alpha:   dale.alpha,$
;      r60_100: dale.r60_100,$
;      loglum:  dale.loglum,$
;      lir:     fltarr(nmodel),$
;      l24:     fltarr(nmodel),$
;      wave:    dblarr(npix*nsamp),$
;      flux:    dblarr(npix*nsamp,nmodel)}

return, data
end
