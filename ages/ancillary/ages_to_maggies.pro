;+
; NAME:
;   AGES_TO_MAGGIES
;
; PURPOSE:
;   Convert the AGES photometry to extinction-corrected AB maggies.
;
; INPUTS: 
;   cat - input photometric catalog (see AGES_MERGE_CATALOGS) [NGAL] 
;
; KEYWORD PARAMETERS:
;   sdss - use SDSS photometry of AGES sources
;   twomass - use TWOMASS photometry of AGES sources
;
; OUTPUTS: 
;   maggies - output maggies [6,NGAL]
;   ivarmaggies - corresponding inverse variance array [6,NGAL] 
;
; COMMENTS:
;   A minimum error of 0.05 mag is applied to every bandpass. 
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Mar 03, NYU - written
;   jm09may15nyu - added TWOMASS keyword
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

pro ages_to_maggies, cat, maggies, ivarmaggies, sdss=sdss, $
  filterlist=filterlist, psf=psf, use_aper=use_aper, $
  totalmag=totalmag, itot=itot, bands=bands

    nobj = n_elements(cat)
    if (nobj eq 0) then begin
       doc_library, 'ages_to_maggies'
       return
    endif

; convert the SDSS/2MASS photometry to maggies, ignoring the other
; multiwavelength photometry 
    if keyword_set(sdss) then begin
       nband = 8
       maggies = dblarr(nband,nobj)
       ivarmaggies = dblarr(nband,nobj)
       filterlist = ['sdss_'+['u0','g0','r0','i0','z0'],$
         'twomass_'+['J','H','Ks']]+'.par'
; SDSS
       good = where(cat.sdss_ra gt -900)
       sdss_to_maggies, maggies1, ivarmaggies1, $
         calibobj=cat[good], flux='model'
       maggies[0:4,good] = maggies1
       ivarmaggies[0:4,good] = ivarmaggies1
; 2MASS
       good = where(cat.twomass_match)
       ages_twomass_to_maggies, cat[good], maggies1, ivarmaggies1
       maggies[5:7,good] = maggies1
       ivarmaggies[5:7,good] = ivarmaggies1
       return
    endif

; GALEX
    im_galex_to_maggies, cat, galex_maggies, galex_ivarmaggies
    bands = ['FUV','NUV']

; BOOTES
    bootes_to_maggies, cat, bootes_maggies, bootes_ivarmaggies, $
      filterlist=bootes_filterlist, psf=psf, use_aper=use_aper, $
      totalmag=totalmag, itot=itot, bands=bootes_bands
    bands = [bands,bootes_bands]

; the z-band photometry of objects near the edges of the individual
; mosaics are suspect, so toss them out here
    zfile = ages_path(/mycat)+'zbootes/zbootes_field_centers.fits.gz'
;   splog, 'Reading '+zfile
    field = mrdfits(zfile,1)
    nfield = n_elements(field)
    zmask = intarr(nobj) ; all start out bad
    pad = 15.0/3600.0 ; 15" padding
    
    for ii = 0, nfield-1 do zmask = zmask or $
      ((cat.ra gt field[ii].ra-field[ii].dra/2.0+pad) and $
      (cat.ra lt field[ii].ra+field[ii].dra/2.0-pad) and $
      (cat.dec gt field[ii].dec-field[ii].ddec/2.0+pad) and $      
      (cat.dec lt field[ii].dec+field[ii].ddec/2.0-pad))
;   good = where(zmask);,comp=bad)
    zfilt = (where(strmatch(bootes_filterlist,'*bok*')))[0]
    bad = where((zmask eq 0) and (bootes_maggies[zfilt,*] gt 0.0))

;   djs_plot, cat.ra, cat.dec, ps=3, ysty=3               
;   djs_oplot, cat[bad].ra, cat[bad].dec, ps=3, color='red'
;   bootes_maggies[zfilt,bad] = 0.0
    bootes_ivarmaggies[zfilt,bad] = 0.0

; fix some other crummy photometry
    
    
    
; final photometry     
    maggies = [galex_maggies,bootes_maggies]
    ivarmaggies = [galex_ivarmaggies,bootes_ivarmaggies]
    filterlist = strtrim([galex_filterlist(),bootes_filterlist],2)

return
end
