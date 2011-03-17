;+
; NAME:
;       SINGS_INFO
;
; PURPOSE:
;       Merge all the ancillary and spectral data on each unique
;       galaxy with an integrated spectrum in the spectral sings.
;       The output data structure is the master binary FITS table on
;       the sample.
;
; CALLING SEQUENCE:
;       sings_info, sings, /write
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       write - write out
;
; OUTPUTS:
;       sings - data structure
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
;       J. Moustakas, 2005 Jul 26, U of A - updated
;-

pro sings_info, sings, write=write
    
    analysis_path = sings_path(/analysis)
    datapath = sings_path(/spec1d)

    version = sings_version(/ancillary)
    
; forage the 1D spectra for additional spectral information

    splog, 'Foraging 1D spectral headers.'
    info = iforage(file_search(datapath+'*.ms.fits*'))
    info_galaxy = strtrim(info.galaxy,2)
    
; define the sample of 75 galaxies using the ancillary data, *not* by
; the existence (or not) of a spectrum (contrast: ATLAS1D_INFO)

    splog, 'Reading the ancillary data.'
    ancillary = mrdfits(analysis_path+'sings_ancillary_data_'+version+'.fits.gz',1,/silent)
    ugalaxy = strtrim(ancillary.galaxy,2)
    ngalaxy = n_elements(ugalaxy)

    splog, 'SINGS: '+string(ngalaxy,format='(I0)')+' galaxies.'

; initialize the output information structure    
    
    sings = {$
      sings_id:                   0L, $
      galaxy:                    ' ', $
      ned_galaxy:                ' ', $
      nice_galaxy:               ' ', $
      alt_galaxy:                ' ', $
      ra:                        ' ', $
      dec:                       ' ', $
                              
      drift56:                    0L, $ ; is there a drift56 spectrum?
      drift56_file:              ' ', $ ; 1D spectrum file name
      drift56_2dfile:            ' ', $ ; 2D spectrum file name
      drift56_asciifile:         ' ', $ ; ascii file name - see SINGS_WRITE_ASCII
      drift56_ra:                ' ', $ ; RA
      drift56_dec:               ' ', $ ; DEC
      drift56_date:              ' ', $ ; date observed
      drift56_observat:          ' ', $ ; observatory
      drift56_ap:                0.0, $ ; extraction aperture [arcsec]
      drift56_scan:              0.0, $ ; drift scan length [arcsec]
      drift56_strap:             ' ', $ ; aperture x scanlength
      drift56_latexap:           ' ', $ ; aperture x scanlength
      drift56_posangle:       -999.0, $ ; position angle
      drift56_exptime:        -999.0, $ ; effective exposure time [s]
      drift56_oexptime:       -999.0, $ ; original exposure time [s]
      drift56_minwave:           0.0, $ ; minimum wavelength [Angstroms]
      drift56_maxwave:           0.0, $ ; maximum wavelength [Angstroms]
      drift56_photflag:          ' ', $ ; photometric?
      drift56_absphoterr:        0.0, $ ; absolute spectrophotometric error
      drift56_agnflag:            0L, $ ; 1 = broad-line AGN 
      drift56_comments:          ' ', $ ; other comments
      drift56_fraction_b:     -999.0, $ ; B-band light fraction
      drift56_fraction_b_err: -999.0, $ ; error
      
      drift20:                    0L, $ ; is there a drift20 spectrum?
      drift20_file:              ' ', $ ; 1D spectrum file name
      drift20_2dfile:            ' ', $ ; 2D spectrum file name
      drift20_asciifile:         ' ', $ ; ascii file name - see SINGS_WRITE_ASCII
      drift20_ra:                ' ', $ ; RA
      drift20_dec:               ' ', $ ; DEC
      drift20_date:              ' ', $ ; date observed
      drift20_observat:          ' ', $ ; observatory
      drift20_ap:                0.0, $ ; extraction aperture [arcsec]
      drift20_scan:              0.0, $ ; drift scan length [arcsec]
      drift20_strap:             ' ', $ ; aperture x scanlength
      drift20_latexap:           ' ', $ ; aperture x scanlength
      drift20_posangle:       -999.0, $ ; position angle
      drift20_exptime:        -999.0, $ ; effective exposure time [s]
      drift20_oexptime:       -999.0, $ ; original exposure time [s]
      drift20_minwave:           0.0, $ ; minimum wavelength [Angstroms]
      drift20_maxwave:           0.0, $ ; maximum wavelength [Angstroms]
      drift20_photflag:          ' ', $ ; photometric?
      drift20_absphoterr:        0.0, $ ; absolute spectrophotometric error
      drift20_agnflag:            0L, $ ; 1 = broad-line AGN 
      drift20_comments:          ' ', $ ; other comments
      drift20_fraction_b:     -999.0, $ ; B-band light fraction
      drift20_fraction_b_err: -999.0, $ ; error
                              
      nuclear:                    0L, $ ; is there a nuclear spectrum?
      nuclear_file:              ' ', $ ; 1D spectrum file name
      nuclear_2dfile:            ' ', $ ; 2D spectrum file name
      nuclear_asciifile:         ' ', $ ; ascii file name - see SINGS_WRITE_ASCII
      nuclear_ra:                ' ', $ ; RA
      nuclear_dec:               ' ', $ ; DEC
      nuclear_date:              ' ', $ ; date observed
      nuclear_observat:          ' ', $ ; observatory
      nuclear_ap:                0.0, $ ; extraction aperture [arcsec]
      nuclear_scan:              2.5, $ ; slit width [arcsec]
      nuclear_strap:             ' ', $ ; aperture x scanlength
      nuclear_latexap:           ' ', $ ; aperture x scanlength
      nuclear_posangle:       -999.0, $ ; position angle
      nuclear_exptime:        -999.0, $ ; effective exposure time [s]
      nuclear_minwave:        -999.0, $ ; minimum wavelength [Angstroms]
      nuclear_maxwave:           0.0, $ ; maximum wavelength [Angstroms]
      nuclear_photflag:          ' ', $ ; photometric?
      nuclear_absphoterr:        0.0, $ ; absolute spectrophotometric error
      nuclear_agnflag:            0L, $ ; 1 = broad-line AGN 
      nuclear_comments:          ' ', $ ; other comments
      nuclear_fraction_b:     -999.0, $ ; B-band light fraction
      nuclear_fraction_b_err: -999.0}   ; error
    sings = replicate(sings,ngalaxy)
    sings.galaxy = strtrim(ugalaxy,2)
    
; read any additional information here

    extra = read_sings_extra_info()
    extra_galaxy = strtrim(extra.galaxy,2)

; read and parse D. Dale's light fractions

; total    
    
    readcol, analysis_path+'flux.all.dat', dale_galaxy1, flux_total, $
      skip=2L, format='A,X,X,F', /silent
    readcol, analysis_path+'uncertainty.all.dat', eflux_total, $
      skip=2L, format='X,X,X,F', /silent

; drift56    
    
    readcol, analysis_path+'flux.drift56.dat', dale_galaxy, flux_drift56, $
      skip=1L, format='A,X,X,F', /silent
    readcol, analysis_path+'eflux.drift56.dat', eflux_drift56, $
      skip=1L, format='X,X,X,F', /silent
    good = where((flux_drift56 gt 0.0) and (eflux_drift56 gt 0.0),ngood)
    sings[good].drift56_fraction_b     = flux_drift56[good]/flux_total[good]
    sings[good].drift56_fraction_b_err = im_compute_error(flux_drift56[good],$
      eflux_drift56[good],flux_total[good],eflux_total[good],/quotient)

; drift20
    
    readcol, analysis_path+'flux.drift20.dat', dale_galaxy, flux_drift20, $
      skip=1L, format='A,X,X,F', /silent
    readcol, analysis_path+'eflux.drift20.dat', eflux_drift20, $
      skip=1L, format='X,X,X,F', /silent
    good = where((flux_drift20 gt 0.0) and (eflux_drift20 gt 0.0),ngood)
    sings[good].drift20_fraction_b     = flux_drift20[good]/flux_total[good]
    sings[good].drift20_fraction_b_err = im_compute_error(flux_drift20[good],$
      eflux_drift20[good],flux_total[good],eflux_total[good],/quotient)

; nuclear    
    
    readcol, analysis_path+'flux.nuc.dat', dale_galaxy, flux_nuclear, $
      skip=1L, format='A,X,X,F', /silent
    readcol, analysis_path+'eflux.nuc.dat', eflux_nuclear, $
      skip=1L, format='X,X,X,F', /silent
    good = where((flux_nuclear gt 0.0) and (eflux_nuclear gt 0.0),ngood)
    sings[good].nuclear_fraction_b     = flux_nuclear[good]/flux_total[good]
    sings[good].nuclear_fraction_b_err = im_compute_error(flux_nuclear[good],$
      eflux_nuclear[good],flux_total[good],eflux_total[good],/quotient)

; append the ancillary data
    
    sings = struct_addtags(sings,struct_trimtags(ancillary,$
      except_tags=['GALAXY','NED_GALAXY','LEDA_GALAXY','RA','DEC']))

    sings.ned_galaxy  = ancillary.ned_galaxy
    sings.ra          = ancillary.ra
    sings.dec         = ancillary.dec

; loop on each object and parse all the header data

    for j = 0L, ngalaxy-1L do begin

       match = where(ugalaxy[j] eq info_galaxy,nmatch)
       if (nmatch eq 0L) then begin

          splog, 'No spectral data for '+ugalaxy[j]+'.'

       endif else begin

          w56 = where(info[match].scanlen ge 55.0,nw56)
          if (nw56 ne 0L) then begin
             
             sings[j].drift56            = 1L
             sings[j].drift56_file       = file_basename(strtrim(info[match[w56]].file,2))
             sings[j].drift56_2dfile     = strtrim(info[match[w56]].spec2d,2)
             sings[j].drift56_asciifile  = repstr(repstr(repstr(strtrim(sings[j].drift56_file,2),$
               '.fits','.txt'),'.ms',''),'.gz','')
             sings[j].drift56_ra         = strtrim(info[match[w56]].ra,2)
             sings[j].drift56_dec        = strtrim(info[match[w56]].dec,2)
             sings[j].drift56_date       = strtrim(strmid(info[match[w56]].date,0,10),2)
             sings[j].drift56_observat   = strtrim(info[match[w56]].observat,2)
             sings[j].drift56_ap         = info[match[w56]].aperwid
             sings[j].drift56_scan       = info[match[w56]].scanlen
             sings[j].drift56_posangle   = info[match[w56]].posangle
             sings[j].drift56_exptime    = info[match[w56]].exptime
             sings[j].drift56_oexptime   = info[match[w56]].oexptime
             sings[j].drift56_minwave    = info[match[w56]].crval1
             sings[j].drift56_maxwave    = info[match[w56]].crval1+(info[match[w56]].naxis1-1)*info[match[w56]].cd1_1
             sings[j].drift56_photflag   = strtrim(info[match[w56]].photflag,2)
             sings[j].drift56_absphoterr = info[match[w56]].abserror
             
             sings[j].drift56_strap = string(sings[j].drift56_ap,format='(I0)')+' x '+$
               string(sings[j].drift56_scan,format='(I0)')
             sings[j].drift56_latexap = '$'+string(sings[j].drift56_ap,format='(I0)')+'\arcsec\times'+$
               string(sings[j].drift56_scan,format='(I0)')+'\arcsec$'

          endif

          w20 = where(info[match].scanlen eq 20.0,nw20)
          if (nw20 ne 0L) then begin

             sings[j].drift20           = 1L
             sings[j].drift20_file      = file_basename(strtrim(info[match[w20]].file,2))
             sings[j].drift20_2dfile    = strtrim(info[match[w20]].spec2d,2)
             sings[j].drift20_asciifile = repstr(repstr(repstr(strtrim(sings[j].drift20_file,2),$
               '.fits','.txt'),'.ms',''),'.gz','')
             sings[j].drift20_ra        = strtrim(info[match[w20]].ra,2)
             sings[j].drift20_dec       = strtrim(info[match[w20]].dec,2)
             sings[j].drift20_date      = strtrim(strmid(info[match[w20]].date,0,10),2)
             sings[j].drift20_observat  = strtrim(info[match[w20]].observat,2)
             sings[j].drift20_ap        = info[match[w20]].aperwid
             sings[j].drift20_scan      = info[match[w20]].scanlen
             sings[j].drift20_posangle  = info[match[w20]].posangle
             sings[j].drift20_exptime   = info[match[w20]].exptime
             sings[j].drift20_oexptime  = info[match[w20]].oexptime
             sings[j].drift20_minwave   = info[match[w20]].crval1
             sings[j].drift20_maxwave   = info[match[w20]].crval1+(info[match[w20]].naxis1-1)*info[match[w20]].cd1_1
             sings[j].drift20_photflag  = strtrim(info[match[w20]].photflag,2)
             sings[j].drift20_absphoterr = info[match[w20]].abserror
             
             sings[j].drift20_strap     = string(sings[j].drift20_ap,format='(I0)')+' x '+$
               string(sings[j].drift20_scan,format='(I0)')
             sings[j].drift20_latexap = '$'+string(sings[j].drift20_ap,format='(I0)')+'\arcsec\times'+$
               string(sings[j].drift20_scan,format='(I0)')+'\arcsec$'

          endif

          wnuc = where(info[match].scanlen eq 0.0,nwnuc)
          if (nwnuc ne 0L) then begin

             sings[j].nuclear           = 1L
             sings[j].nuclear_file      = file_basename(strtrim(info[match[wnuc]].file,2))
             sings[j].nuclear_2dfile    = strtrim(info[match[wnuc]].spec2d,2)
             sings[j].nuclear_asciifile = repstr(repstr(repstr(strtrim(sings[j].nuclear_file,2),$
               '.fits','.txt'),'.ms',''),'.gz','')
             sings[j].nuclear_ra        = strtrim(info[match[wnuc]].ra,2)
             sings[j].nuclear_dec       = strtrim(info[match[wnuc]].dec,2)
             sings[j].nuclear_date      = strtrim(strmid(info[match[wnuc]].date,0,10),2)
             sings[j].nuclear_observat  = strtrim(info[match[wnuc]].observat,2)
             sings[j].nuclear_ap        = info[match[wnuc]].aperwid
             sings[j].nuclear_scan      = info[match[wnuc]].aperture
             sings[j].nuclear_posangle  = info[match[wnuc]].posangle
             sings[j].nuclear_exptime   = info[match[wnuc]].exptime
             sings[j].nuclear_minwave   = info[match[wnuc]].crval1
             sings[j].nuclear_maxwave   = info[match[wnuc]].crval1+(info[match[wnuc]].naxis1-1)*info[match[wnuc]].cd1_1
             sings[j].nuclear_photflag  = strtrim(info[match[wnuc]].photflag,2)
             sings[j].nuclear_absphoterr = info[match[wnuc]].abserror
             
             sings[j].nuclear_strap     = string(sings[j].nuclear_ap,format='(F3.1)')+' x '+$
               string(sings[j].nuclear_scan,format='(F3.1)')
             sings[j].nuclear_latexap = '$'+repstr(string(sings[j].nuclear_ap,format='(F3.1)'),'.','\farcs')+'\times'+$
               repstr(string(sings[j].nuclear_scan,format='(F3.1)'),'.','\farcs')+'$'

          endif

       endelse
          
; extra information       

       indx = where(ugalaxy[j] eq strtrim(extra_galaxy,2),nindx)
       if (nindx ne 0L) then begin

          sings[j].nice_galaxy      = strtrim(extra[indx].nicegalaxy,2)
          sings[j].alt_galaxy       = strtrim(extra[indx].altgalaxy,2)
          sings[j].drift56_comments = strtrim(extra[indx].drift56_comments,2)
          sings[j].drift20_comments = strtrim(extra[indx].drift20_comments,2)
          sings[j].nuclear_comments = strtrim(extra[indx].nuclear_comments,2)
          sings[j].lit_type         = strtrim(extra[indx].type,2)   ; from the SINGS team!
          sings[j].lit_t            = morph2type(repstr(repstr(sings[j].lit_type,'~',''),'pec','')) ; jm06oct10nyu
          sings[j].drift56_agnflag  = extra[indx].drift56_agnflag
          sings[j].drift20_agnflag  = extra[indx].drift20_agnflag
          sings[j].nuclear_agnflag  = extra[indx].nuclear_agnflag

       endif else splog, 'No entry in the EXTRAS file for '+ugalaxy[j]+'!'
       
    endfor
    
; sort by RA

    srtra = sort(im_hms2dec(sings.ra))
    sings = sings[srtra]
    sings.sings_id = lindgen(ngalaxy)+1L

    if keyword_set(write) then begin
       outname = 'sings_info_'+version+'.fits'
       splog, 'Writing '+analysis_path+outname+'.'
       mwrfits, sings, analysis_path+outname, /create
       spawn, ['gzip -f '+analysis_path+outname], /sh
    endif

return
end
