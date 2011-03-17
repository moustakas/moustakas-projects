;+
; NAME:
;       NFGS_INFO
;
; PURPOSE:
;       Merge all the ancillary and spectral data on each unique
;       galaxy with an integrated spectrum in the NFGS.  The output
;       data structure is the master binary FITS table for the sample.  
;
; CALLING SEQUENCE:
;       nfgs_info, nfgs, /write
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       write - write out
;
; OUTPUTS:
;       nfgs - data structure
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
;       J. Moustakas, 2005 July 24, U of A - written
;-

pro nfgs_info, nfgs, write=write

    analysis_path = nfgs_path(/analysis)
    datapath = nfgs_path(/spec1d)
    
; forage the 1D spectra for additional spectral information

    splog, 'Foraging 1D spectral headers.'
    info = iforage(file_search('*int*.fits*'),datapath=datapath)
    ngalaxy = n_elements(info)
    splog, 'Found '+string(ngalaxy,format='(I0)')+' galaxies with integrated spectra.'
    
; initialize the output information structure    
    
    nfgs = {$
      nfgs_id:           -1L, $
      galaxy:            ' ', $
      ned_galaxy:        ' ', $
      leda_galaxy:       ' ', $
      nfgs_galaxy:       ' ', $
;     alt_galaxy:        ' ', $
      nice_galaxy:       ' ', $
      ra:                ' ', $
      dec:               ' ', $

      drift:              0L, $ ; drift scanned spectrum
      drift_file:        ' ', $ 
      drift_asciifile:   ' ', $ 
      drift_ra:          ' ', $ 
      drift_dec:         ' ', $ 
      drift_date:        ' ', $ 
      drift_ap:          0.0, $ ; extraction aperture [arcsec]
      drift_scan:        0.0, $ ; drift scan length [arcsec]
      drift_strap:       ' ', $ ; aperture x scanlength
      drift_posangle:    0.0, $ ; position angle [degree]
      drift_exptime:     0.0, $ ; effective exposure time [s]
      drift_oexptime:    0.0, $ ; original exposure time [s]
;     drift_zptshift:    0.0, $
;     drift_note:        ' ', $
;     drift_photflag:    ' ', $
      drift_agnflag:      0L}   ; 1 = broad-line AGN 
;     drift_code:        -1L, $ ; 1 = integrated; 2 = missing LSB features; 3 = not truly integrated
;     drift_comments:    ' '}

;     nuclear:            0L, $ ; nuclear spectrum
;     nuclear_file:      ' ', $ 
;     nuclear_asciifile: ' ', $ 
;     nuclear_ra:        ' ', $ 
;     nuclear_dec:       ' ', $ 
;     nuclear_date:      ' ', $ 
;     nuclear_ap:        0.0, $ ; extraction aperture [arcsec]
;     nuclear_scan:      2.5, $ ; slit width [arcsec]
;     nuclear_strap:     ' ', $ ; aperture x scanlength
;     nuclear_posangle:  0.0, $ ; position angle [degree]
;     nuclear_exptime:   0.0, $ ; exposure time [s]
;     nuclear_zptshift:  0.0, $
;     nuclear_note:      ' ', $
;     nuclear_photflag:  ' ', $
;     nuclear_agnflag:    0L, $ ; 1 = broad-line AGN 
;     nuclear_comments:  ' '}

    nfgs = replicate(nfgs,ngalaxy)
    nfgs.galaxy = info.galaxy

; read and append the ancillary data
    
    splog, 'Reading the ancillary data.'
    ancillary = mrdfits(analysis_path+'nfgs_ancillary_data.fits.gz',1,/silent)
    ancillary_galaxy = strtrim(ancillary.ned_galaxy,2) ; <-- NOTE!

    match, strtrim(nfgs.galaxy,2), ancillary_galaxy, indx1, match
    
    nfgs = struct_addtags(nfgs,struct_trimtags(ancillary[match],$
      except_tags=['GALAXY','NED_GALAXY','LEDA_GALAXY','RA','DEC']))

    nfgs[indx1].ned_galaxy  = ancillary[match].ned_galaxy
    nfgs[indx1].leda_galaxy = ancillary[match].leda_galaxy
    nfgs[indx1].ra          = ancillary[match].ra
    nfgs[indx1].dec         = ancillary[match].dec

; read and store additional spectral information here

    splog, 'Reading the extra information.'
    extra = read_nfgs_extra_info()
    match, strtrim(nfgs.galaxy,2), strtrim(extra.galaxy,2), indx1, match

    nfgs[indx1].nice_galaxy     = strtrim(extra[match].nicegalaxy,2)
    nfgs[indx1].lit_type        = strtrim(extra[match].morph_type,2)
    nfgs[indx1].lit_t           = strtrim(extra[match].morph_t,2)
    nfgs[indx1].drift_agnflag   = extra[match].int_agnflag

; concatenate the Jansen et al. (2000) broadband photometry

    jansen = read_00jansen()
    match, strtrim(jansen.ned_galaxy,2), strtrim(nfgs.ned_galaxy,2), indx1, match

    nfgs[match].nfgs_galaxy   = jansen[indx1].nfgs_galaxy
    nfgs[match].nfgs_id       = jansen[indx1].nfgs_id
;   niceprint, nfgs.nfgs_id, nfgs.nfgs_galaxy, nfgs.galaxy, nfgs.ned_galaxy

    nfgs[match].u_obs         = jansen[indx1].u
    nfgs[match].u_obs_err     = jansen[indx1].u_err
    nfgs[match].u_lum_obs     = jansen[indx1].u_lum
    nfgs[match].u_lum_obs_err = jansen[indx1].u_lum_err
    nfgs[match].m_u_obs       = jansen[indx1].m_u
    nfgs[match].m_u_obs_err   = jansen[indx1].m_u_err
    
    nfgs[match].b_obs         = jansen[indx1].b
    nfgs[match].b_obs_err     = jansen[indx1].b_err
    nfgs[match].b_lum_obs     = jansen[indx1].b_lum
    nfgs[match].b_lum_obs_err = jansen[indx1].b_lum_err
    nfgs[match].m_b_obs       = jansen[indx1].m_b
    nfgs[match].m_b_obs_err   = jansen[indx1].m_b_err

    nfgs[match].v_obs         = jansen[indx1].v
    nfgs[match].v_obs_err     = jansen[indx1].v_err
    nfgs[match].v_lum_obs     = jansen[indx1].v_lum
    nfgs[match].v_lum_obs_err = jansen[indx1].v_lum_err
    nfgs[match].m_v_obs       = jansen[indx1].m_v
    nfgs[match].m_v_obs_err   = jansen[indx1].m_v_err
    
    nfgs[match].r_obs         = jansen[indx1].r
    nfgs[match].r_obs_err     = jansen[indx1].r_err
    nfgs[match].r_lum_obs     = jansen[indx1].r_lum
    nfgs[match].r_lum_obs_err = jansen[indx1].r_lum_err
    nfgs[match].m_r_obs       = jansen[indx1].m_r
    nfgs[match].m_r_obs_err   = jansen[indx1].m_r_err
    
; loop on each object and parse all the header data; INFO and NFGS are
; sorted in the same way

;   niceprint, nfgs.galaxy, info.galaxy
    
    for j = 0L, ngalaxy-1L do begin

       nfgs[j].drift           = 1L
       nfgs[j].drift_file      = strtrim(info[j].file,2)
       nfgs[j].drift_asciifile = repstr(repstr(repstr(strtrim(nfgs[j].drift_file,2),$
         '.fits','.txt'),'.ms',''),'.gz','')
       nfgs[j].drift_ra        = strtrim(info[j].ra,2)
       nfgs[j].drift_dec       = strtrim(info[j].dec,2)
       nfgs[j].drift_date      = strtrim(strmid(info[j].date,0,10),2)
       nfgs[j].drift_ap        = info[j].aperwid
       nfgs[j].drift_scan      = info[j].scanlen
       nfgs[j].drift_posangle  = info[j].posangle
       nfgs[j].drift_exptime   = info[j].exptime
       nfgs[j].drift_oexptime  = info[j].oexptime
;      nfgs[j].drift_zptshift  = info[j].zptshift
;      nfgs[j].drift_note      = strtrim(info[j].object,2)
;      nfgs[j].drift_photflag  = strtrim(info[j].photflag,2)
       
       nfgs[j].drift_strap     = string(nfgs[j].drift_ap,format='(I0)')+' x '+$
         string(nfgs[j].drift_scan,format='(I0)')

    endfor 

; compute the gas fraction using relations in McGaugh & de Block
; (1997); use the literature morphological types stored in
; NFGS_INFO

    good = where(nfgs.lit_t gt -900.0,ngood)
    if (ngood ne 0L) then begin
       type = nfgs[good].lit_t
       gas_ratio = ((3.7-0.8*type+0.043*type^2)>0.0)<4.0 ; gas_ratio --> 0 as type --> 10
       nfgs[good].eta_gas = 1.4*(1+gas_ratio)
    endif

; only use the LEDA HI masses    
    
    good = where((nfgs.eta_gas gt -900.0) and (nfgs.mass_HI_leda gt -900.0),ngood) 
    if (ngood ne 0L) then begin
       nfgs[good].mass_gas     = nfgs[good].mass_HI_leda + alog10(nfgs[good].eta_gas)
       nfgs[good].mass_gas_err = nfgs[good].mass_HI_leda_err
    endif
    
; sort by RA

    srtra = sort(im_hms2dec(nfgs.ra))
    nfgs = nfgs[srtra]

    if keyword_set(write) then begin

       outname = 'nfgs_info.fits'
           
       splog, 'Writing '+analysis_path+outname+'.'
       mwrfits, nfgs, analysis_path+outname, /create
       spawn, ['gzip -f '+analysis_path+outname], /sh

    endif

return
end

;    galaxy = strtrim(info.galaxy,2)
;    integrated = where(info.scanlen gt 0.0,nintegrated)
;    
;    intgalaxy = strtrim(info[integrated].galaxy,2)
;    uindx = uniq(intgalaxy,sort(intgalaxy))
;    ugalaxy = intgalaxy[uindx]
;    if (n_elements(ugalaxy) ne n_elements(intgalaxy)) then $
;      message, 'UGALAXY and GALAXY must have the same number of elements!'
;    ngalaxy = n_elements(ugalaxy)

