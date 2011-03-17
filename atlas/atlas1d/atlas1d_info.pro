;+
; NAME:
;       ATLAS1D_INFO
;
; PURPOSE:
;       Merge all the ancillary and spectral data on each unique
;       galaxy with an integrated spectrum in the spectral atlas.
;       The output data structure is the master binary FITS table on
;       the sample.
;
; CALLING SEQUENCE:
;       atlas1d_info, atlas, /write
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       write - write out
;
; OUTPUTS:
;       atlas - data structure
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
;       J. Moustakas, Spring 2003, U of A
;       jm05jul22uofa - updated
;-

pro atlas1d_info, atlas, write=write

    version = atlas_version(/ancillary)

    analysis_path = atlas_path(/analysis)
    datapath = atlas_path(/atlas1d)
    
; forage the 1D spectra for additional spectral information

    splog, 'Foraging 1D spectral headers.'
    info = iforage(file_search(datapath+'*.ms.fits*'))

; read additional spectral information and the ancillary data 

    extra = read_atlas_extra_info()
    extra_galaxy = strtrim(extra.galaxy)

    splog, 'Reading the ancillary data.'
    ancillary = mrdfits(analysis_path+'atlas_ancillary_data_'+version+'.fits.gz',1,/silent)

; construct a unique set of galaxy names using only the integrated
; spectra; the corresponding nuclear spectra will be identified below;
; this has to happen because we have many (~40) nuclear spectra
; without the corresponding integrated spectra

    galaxy = strtrim(info.galaxy,2)
    integrated = where(info.scanlen gt 0.0,nintegrated)
    intgalaxy = strtrim(info[integrated].galaxy,2)
    uindx = uniq(intgalaxy,sort(intgalaxy))
    ugalaxy = intgalaxy[uindx]
    if (n_elements(ugalaxy) ne n_elements(intgalaxy)) then $
      message, 'UGALAXY and GALAXY must have the same number of elements!'
    ngalaxy = n_elements(ugalaxy)

    splog, 'Found '+string(ngalaxy,format='(I0)')+' unique galaxies with integrated spectra.'
    
; initialize the output information structure    
    
    atlas = {$
      atlas_id:            0L, $
      galaxy:             ' ', $
      ned_galaxy:         ' ', $
      leda_galaxy:        ' ', $
      alt_galaxy:         ' ', $
      nice_galaxy:        ' ', $
      ra:                 ' ', $
      dec:                ' ', $
                          
      drift:               0L, $ ; drift scanned spectrum
      drift_file:         ' ', $ 
      drift_2dfile:       ' ', $ 
      drift_asciifile:    ' ', $ 
      drift_ra:           ' ', $ 
      drift_dec:          ' ', $ 
      drift_date:         ' ', $ 
      drift_ap:           0.0, $ ; extraction aperture [arcsec]
      drift_scan:         0.0, $ ; drift scan length [arcsec]
      drift_strap:        ' ', $ ; aperture x scanlength
      drift_posangle:     0.0, $ ; position angle [degree]
      drift_exptime:      0.0, $ ; effective exposure time [s]
      drift_oexptime:     0.0, $ ; original exposure time [s]
      drift_minwave:      0.0, $ ; minimum wavelength [Angstroms]
      drift_maxwave:      0.0, $ ; maximum wavelength [Angstroms]
      drift_zptshift:     0.0, $
      drift_note:         ' ', $
      drift_photflag:     ' ', $
      drift_absphoterr:   0.0, $ ; absolute spectrophotometric error
      drift_agnflag:       0L, $ ; 1 = broad-line AGN 
      drift_code:         -1L, $ ; 1 = integrated; 2 = missing LSB features; 3 = not truly integrated
      drift_comments:     ' ', $
                          
      nuclear:             0L, $ ; nuclear spectrum
      nuclear_file:       ' ', $ 
      nuclear_2dfile:     ' ', $ 
      nuclear_asciifile:  ' ', $ 
      nuclear_ra:         ' ', $ 
      nuclear_dec:        ' ', $ 
      nuclear_date:       ' ', $ 
      nuclear_ap:         0.0, $ ; extraction aperture [arcsec]
      nuclear_scan:       2.5, $ ; slit width [arcsec]
      nuclear_strap:      ' ', $ ; aperture x scanlength
      nuclear_posangle:   0.0, $ ; position angle [degree]
      nuclear_exptime:    0.0, $ ; exposure time [s]
      nuclear_minwave:    0.0, $ ; minimum wavelength [Angstroms]
      nuclear_maxwave:    0.0, $ ; maximum wavelength [Angstroms]
      nuclear_zptshift:   0.0, $
      nuclear_note:       ' ', $
      nuclear_photflag:   ' ', $
      nuclear_absphoterr: 0.0, $ ; absolute spectrophotometric error
      nuclear_agnflag:     0L, $ ; 1 = broad-line AGN 
      nuclear_comments:   ' '}
    atlas = replicate(atlas,ngalaxy)

    atlas.galaxy = strtrim(ugalaxy,2)

; cross-match UGALAXY against EXTRA.GALAXY to get the NED_GALAXY names
; which were used to retrieve the ANCILLARY data

    match, strtrim(ugalaxy,2), strtrim(extra.galaxy,2), match_ugalaxy, match_extra
    if (n_elements(match_ugalaxy) ne n_elements(match_extra)) then $
      message, 'These should match!'
    if (monotonic(match_ugalaxy) eq 0L) then message, 'This should be monotonic!'
    extra1 = extra[match_extra]
;   niceprint, atlas.galaxy, extra1.galaxy

    match, strtrim(extra1.nedgalaxy,2), strtrim(ancillary.galaxy,2), $
      match_extra_ned_galaxy, match_ancillary_ned_galaxy
    if (n_elements(match_ancillary_ned_galaxy) ne n_elements(match_extra_ned_galaxy)) then $
      message, 'These should match!'
;   if (monotonic(match_extra_ned_galaxy) eq 0L) then message, 'This should be monotonic!'

    srt = sort(match_extra_ned_galaxy)
    ancillary1 = ancillary[match_ancillary_ned_galaxy[srt]]
;   niceprint, atlas.galaxy, ancillary1.galaxy, ancillary1.ned_galaxy

    atlas = struct_addtags(atlas,struct_trimtags(ancillary1,$
      except_tags=['GALAXY','NED_GALAXY','LEDA_GALAXY','RA','DEC']))

    atlas.ned_galaxy  = ancillary1.ned_galaxy
    atlas.leda_galaxy = ancillary1.leda_galaxy
    atlas.ra          = ancillary1.ra
    atlas.dec         = ancillary1.dec

; loop on each object and parse all the header data

    for j = 0L, ngalaxy-1L do begin

       match = where(ugalaxy[j] eq galaxy,nmatch)
       if (nmatch eq 0L) then message, 'Open the pod bay doors Hal.'

; integrated spectra       
       
       w = where(info[match].scanlen gt 0.0,nw)
       if (nw ne 0L) then begin
          
          atlas[j].drift            = 1L
          atlas[j].drift_file       = file_basename(strtrim(info[match[w]].file,2))
          atlas[j].drift_2dfile     = strtrim(info[match[w]].spec2d,2)
          atlas[j].drift_asciifile  = repstr(repstr(repstr(strtrim(atlas[j].drift_file,2),$
            '.fits','.txt'),'.ms',''),'.gz','')
          atlas[j].drift_ra         = strtrim(info[match[w]].ra,2)
          atlas[j].drift_dec        = strtrim(info[match[w]].dec,2)
          atlas[j].drift_date       = strtrim(strmid(info[match[w]].date,0,10),2)
          atlas[j].drift_ap         = info[match[w]].aperwid
          atlas[j].drift_scan       = info[match[w]].scanlen
          atlas[j].drift_posangle   = info[match[w]].posangle
          atlas[j].drift_exptime    = info[match[w]].exptime
          atlas[j].drift_oexptime   = info[match[w]].oexptime
          atlas[j].drift_minwave    = info[match[w]].crval1
          atlas[j].drift_maxwave    = info[match[w]].crval1+(info[match[w]].naxis1-1)*info[match[w]].cd1_1
          atlas[j].drift_zptshift   = info[match[w]].zptshift
          atlas[j].drift_note       = strtrim(info[match[w]].object,2)
          atlas[j].drift_photflag   = strtrim(info[match[w]].photflag,2)
          atlas[j].drift_absphoterr = info[match[w]].abserror
          
          atlas[j].drift_strap     = string(atlas[j].drift_ap,format='(I0)')+' x '+$
            string(atlas[j].drift_scan,format='(I0)')

       endif

; nuclear spectra       
       
       w = where(info[match].scanlen eq 0.0,nw)
       if (nw ne 0L) then begin
          
          atlas[j].nuclear            = 1L
          atlas[j].nuclear_file       = file_basename(strtrim(info[match[w]].file,2))
          atlas[j].nuclear_2dfile     = strtrim(info[match[w]].spec2d,2)
          atlas[j].nuclear_asciifile  = repstr(repstr(repstr(strtrim(atlas[j].nuclear_file,2),$
            '.fits','.txt'),'.ms',''),'.gz','')
          atlas[j].nuclear_ra         = strtrim(info[match[w]].ra,2)
          atlas[j].nuclear_dec        = strtrim(info[match[w]].dec,2)
          atlas[j].nuclear_date       = strtrim(strmid(info[match[w]].date,0,10),2)
          atlas[j].nuclear_ap         = info[match[w]].aperwid
          atlas[j].nuclear_scan       = info[match[w]].aperture
          atlas[j].nuclear_posangle   = info[match[w]].posangle
          atlas[j].nuclear_exptime    = info[match[w]].exptime
          atlas[j].nuclear_minwave    = info[match[w]].crval1
          atlas[j].nuclear_maxwave    = info[match[w]].crval1+(info[match[w]].naxis1-1)*info[match[w]].cd1_1
          atlas[j].nuclear_zptshift   = info[match[w]].zptshift
          atlas[j].nuclear_note       = strtrim(info[match[w]].object,2)
          atlas[j].nuclear_photflag   = strtrim(info[match[w]].photflag,2)
          atlas[j].nuclear_absphoterr = info[match[w]].abserror
          
          atlas[j].nuclear_strap     = string(atlas[j].nuclear_ap,format='(F3.1)')+' x '+$
            string(atlas[j].nuclear_scan,format='(F3.1)')

       endif

; extra information       
       
       indx = where(ugalaxy[j] eq strtrim(extra_galaxy,2),nindx)
;      if (nindx eq 0L) then message, 'Open the pod bay doors Hal.'
       if (nindx ne 0L) then begin

          atlas[j].nice_galaxy     = strtrim(extra[indx].nicegalaxy,2)
          atlas[j].alt_galaxy      = strtrim(extra[indx].altgalaxy[indx],2)
          atlas[j].drift_code      = strtrim(extra[indx].drift_code,2)
          atlas[j].drift_comments  = strtrim(extra[indx].comments,2)
          atlas[j].lit_type        = strtrim(extra[indx].morph_type,2)
          atlas[j].lit_t           = strtrim(extra[indx].morph_t,2)
          atlas[j].drift_agnflag   = extra[indx].int_agnflag
          atlas[j].nuclear_agnflag = extra[indx].nuc_agnflag

       endif else splog, 'No entry in the EXTRAS file for '+ugalaxy[j]+'!'

    endfor

; compute the gas fraction using relations in McGaugh & de Block
; (1997); use the literature morphological types stored in
; ATLAS1D_INFO

    good = where(atlas.lit_t gt -900.0,ngood)
    if (ngood ne 0L) then begin
       type = atlas[good].lit_t
       gas_ratio = ((3.7-0.8*type+0.043*type^2)>0.0)<4.0 ; gas_ratio --> 0 as type --> 10
       atlas[good].eta_gas = 1.4*(1+gas_ratio)
    endif

; only use the LEDA HI masses    
    
    good = where((atlas.eta_gas gt -900.0) and (atlas.mass_HI_leda gt -900.0),ngood) 
    if (ngood ne 0L) then begin
       atlas[good].mass_gas     = atlas[good].mass_HI_leda + alog10(atlas[good].eta_gas)
       atlas[good].mass_gas_err = atlas[good].mass_HI_leda_err
    endif
    
; sort by RA and store the running ID

    srtra = sort(im_hms2dec(atlas.ra))
    atlas = atlas[srtra]
    atlas.atlas_id = lindgen(ngalaxy)+1L
    
    if keyword_set(write) then begin

       outname = 'atlas1d_info_'+version+'.fits'
           
       splog, 'Writing '+analysis_path+outname+'.'
       mwrfits, atlas, analysis_path+outname, /create
       spawn, ['gzip -f '+analysis_path+outname], /sh

; generate a galaxy list of the integrated spectra

;      w = where(atlas.drift)
;      outfile = 'integrated_galaxy_list.txt'
;      splog, 'Writing '+analysis_path+outfile+'.'
;      openw, lun, analysis_path+outfile, /get_lun
;      struct_print, struct_trimtags(atlas[w],select='galaxy'), /no_head, lun=lun
;      free_lun, lun    

    endif

return
end
