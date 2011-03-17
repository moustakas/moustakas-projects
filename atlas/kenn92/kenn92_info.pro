;+
; NAME:
;       KENN92_INFO
;
; PURPOSE:
;       Merge all the ancillary and spectral data on each unique
;       galaxy with an integrated spectrum in the Kennicutt (1992)
;       spectral atlas.  The output data structure is the master
;       binary FITS table for the sample.   
;
; CALLING SEQUENCE:
;       kenn92_info, kenn92, /write
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       write - write out
;
; OUTPUTS:
;       kenn92 - data structure
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       SPLOG, IFORAGE(), MATCH, STRUCT_ADDTAGS(), STRUCT_TRIMTAGS(),
;       IM_HMS2DEC(), MWRFITS
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Aug 02, U of A - written
;-

pro kenn92_info, kenn92, write=write

    analysis_path = kenn92_path(/analysis)
    datapath = kenn92_path(/spec1d)
    
; forage the 1D spectra for additional spectral information

    splog, 'Foraging 1D spectral headers.'
    info = iforage(file_search('*.fits*'),datapath=datapath)

; construct a unique set of galaxy names

    galaxy = strtrim(info.galaxy,2)
    uindx = uniq(galaxy,sort(galaxy))
    ugalaxy = galaxy[uindx]
    ngalaxy = n_elements(ugalaxy)

    splog, 'Found '+string(ngalaxy,format='(I0)')+' unique galaxies with integrated spectra.'

; initialize the output information structure    
    
    kenn92 = {$
      galaxy:            ' ', $
      ned_galaxy:        ' ', $
      alt_galaxy:        ' ', $
      nice_galaxy:       ' ', $
      ra:                ' ', $
      dec:               ' ', $

      drift:              0L, $ ; drift scanned spectrum
      drift_file:        ' ', $ 
      drift_ra:          ' ', $ 
      drift_dec:         ' ', $ 
      drift_ap:          0.0, $ ; extraction aperture [arcsec]
      drift_scan:        0.0, $ ; drift scan length [arcsec]
      drift_strap:       ' ', $ ; aperture x scanlength
      drift_posangle:    0.0, $ ; position angle [degree]
      drift_exptime:     0.0, $ ; effective exposure time [s]
      drift_agnflag:      0L, $ ; 1 = broad-line AGN 
      drift_lowresflag:   0L}   ; 1 = low resolution spectrum
    kenn92 = replicate(kenn92,ngalaxy)

    kenn92.galaxy = ugalaxy

; read additional spectral information here

    extra = read_kenn92_extra_info()
    extra_galaxy = strtrim(extra.galaxy,2)

; read and append the ancillary data

    splog, 'Reading the ancillary data.'
    ancillary = mrdfits(analysis_path+'kenn92_ancillary_data.fits.gz',1,/silent)
    ancillary_galaxy = strtrim(ancillary.galaxy,2)

    match, ugalaxy, ancillary_galaxy, indx1, match
    if (monotonic(indx1) ne 1L) then message, 'This index should be monotonic.'

    kenn92 = struct_addtags(kenn92,struct_trimtags(ancillary[match],$
      except_tags=['GALAXY','NED_GALAXY','RA','DEC']))

    kenn92.ned_galaxy  = ancillary[match].ned_galaxy
    kenn92.ra          = ancillary[match].ra
    kenn92.dec         = ancillary[match].dec

; loop on each object and parse all the header data

    for j = 0L, ngalaxy-1L do begin

       match = where(ugalaxy[j] eq galaxy,nmatch)
       if (nmatch eq 0L) then message, 'Open the pod bay doors Hal.'

; integrated spectra       
       
       w = where(info[match].scanlen gt 0.0,nw)
       if (nw ne 0L) then begin
          
          kenn92[j].drift           = 1L
          kenn92[j].drift_file      = strtrim(info[match[w]].file,2)
          kenn92[j].drift_ra        = strtrim(info[match[w]].ra,2)
          kenn92[j].drift_dec       = strtrim(info[match[w]].dec,2)
          kenn92[j].drift_ap        = info[match[w]].aperwid
          kenn92[j].drift_scan      = info[match[w]].scanlen
          kenn92[j].drift_posangle  = info[match[w]].posangle
          kenn92[j].drift_exptime   = info[match[w]].exptime
          
          kenn92[j].drift_strap     = string(kenn92[j].drift_ap,format='(I0)')+' x '+$
            string(kenn92[j].drift_scan,format='(I0)')

       endif

; extra information       
       
       indx = where(ugalaxy[j] eq strtrim(extra_galaxy,2),nindx)
       if (nindx ne 0L) then begin

          kenn92[j].nice_galaxy      = strtrim(extra[indx].nicegalaxy,2)
          kenn92[j].alt_galaxy       = strtrim(extra[indx].altgalaxy[indx],2)
          kenn92[j].lit_type         = strtrim(extra[indx].morph_type,2)
          kenn92[j].drift_agnflag    = extra[indx].agnflag
          kenn92[j].drift_lowresflag = extra[indx].lowresflag

       endif else splog, 'No entry in the EXTRAS file for '+ugalaxy[j]+'!'

    endfor

    srtra = sort(im_hms2dec(kenn92.ra))
    kenn92 = kenn92[srtra]

    if keyword_set(write) then begin

       outname = 'kenn92_info.fits'
           
       splog, 'Writing '+analysis_path+outname+'.'
       mwrfits, kenn92, analysis_path+outname, /create
       spawn, ['gzip -f '+analysis_path+outname], /sh

    endif

return
end
