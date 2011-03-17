;+
; NAME:
;   PARSE_EDISCS_GANDALF_SPECFIT
;
; PURPOSE:
;   Parse the RUN12 & RUN34 output from EDISCS_GANDALF_SPECFIT and
;   write out the final spectral fits and line-measurements.
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;   debug - render some plots on the screen for debugging
;
; OUTPUTS: 
;   Various.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Nov 13, UCSD
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

pro do_parse_ediscs_gandalf_specfit, run34=run34, debug=debug, $
  solar=solar, specdata=specdata, specfit=specfit

    delvarx, specdata, specfit  ; delete memory!!

    if (n_elements(fluxscale) eq 0) then fluxscale = 1D-17
    snrcut_line = 3.0  ; see GANDALF_CLEAN_EMISSION
    light = 2.99792458D5 ; speed of light [km/s]
    velscale = ediscs_ppxf_velscale(run34=run34)

; path names and emission-line file names
    version = ediscs_version(/ppxf_specfit)
    specfitpath = ediscs_path(/ppxf)
    indexfile = ediscs_path(/ppxf)+'lick_indexlist_'+version+'.dat'
    linefile = ediscs_path(/ppxf)+'gandalf_elinelist_'+version+'.dat'
    linefile_broad = repstr(linefile,'elinelist_','elinelist_broad_')

; parse the EDISCS_GANDALF_SPECFIT output
    if keyword_set(run34) then suffix = 'run34' else suffix = 'run12'
    if keyword_set(solar) then suffix = 'solar_'+suffix
    rawfile = specfitpath+'ediscs_specdata_raw_'+suffix+'_'+version+'.fits.gz'
    rawfile_broad = repstr(rawfile,'raw_','raw_broad_')

; read the model templates so that we can compute the
; luminosity-weighted age for each object
    junk = read_ediscs_ppxf_templates(tempinfo=tempinfo,solar=solar)
    
; initialize the output data structures
    nfinalpix = 2550
    specfit_template = {$
      ediscs_id:                    0L,$
      galaxy:                       '',$
      specfile:                     '',$
      cluster:                      '',$
      run:                          '',$
      ra:                         0.0D,$
      dec:                        0.0D,$
      z:                           0.0,$
      zabs:                        0.0,$
      vdisp:                       0.0,$
      wave:          dblarr(nfinalpix),$
      flux:          fltarr(nfinalpix),$
      ferr:          fltarr(nfinalpix),$
      linefit:       fltarr(nfinalpix),$
      continuum:     fltarr(nfinalpix),$
      smooth_continuum: fltarr(nfinalpix)}

    t1 = systime(1)
    if (file_test(rawfile,/reg) eq 0) then begin
       splog, 'File '+rawfile+' not found'
       return
    endif

; read the data-files    
    splog, 'Reading '+rawfile
    raw = mrdfits(rawfile,1,/silent)
    nobj = n_elements(raw)

; is there a data structure of broad-line AGN?
    if file_test(rawfile_broad,/reg) then begin
       splog, 'Reading '+rawfile_broad
       raw_broad = mrdfits(rawfile_broad,1,/silent)
       broad = 1
    endif else broad = 0

; loop on each object
    t0 = systime(1)
;   for iobj = 49, 49 do begin
    for iobj = 0, nobj-1 do begin
       print, "Object ", iobj, " of ", nobj, string(13b), $
         format='(A0,I0,A0,I0,A1,$)'
       if broad then begin
          match = where(strtrim(raw[iobj].galaxy,2) eq $
            strtrim(raw_broad.galaxy,2),nmatch)
          if (nmatch ne 0) then begin
             thisraw = raw_broad[match]
             thislinefile = linefile_broad
             thisbroad = 1
          endif else begin
             thisraw = raw[iobj]
             thislinefile = linefile 
             thisbroad = 0
          endelse
       endif else begin
          thisraw = raw[iobj]
          thislinefile = linefile 
          thisbroad = 0
       endelse
       
       notzero = where((thisraw.wave gt 0.0),nnotzero)
       if (nnotzero eq 0) then message, 'Problem here!'
       wave = thisraw.wave[notzero]         ; rest
       flux = thisraw.flux[notzero]*1.0D    ; note!
       ferr = thisraw.ferr[notzero]*1.0D 
       continuum = thisraw.continuum[notzero]
       smooth_continuum = thisraw.smooth_continuum[notzero]
       npix = n_elements(wave)
       
       zabs = thisraw.zabs
       vsys = alog(1.0+zabs)*light ; systemic velocity

; rebuild LINEPARS and the emission-line templates
       junk = read_gandalf_elinelist(thislinefile,$
         fitlines=linepars,/actionfit)
       junk = gandalf_mask_emission_lines(npix,vsys,$
         linepars,velscale,wave[0]+alog(1.0+zabs),$
         velscale/light,sigma=300.0)
       these = where(linepars.kind eq 'l',nlines)
       fitlines = where((linepars[these].action ne 'i'),nfitlines)

       sol = thisraw.sol[*,0:nfitlines-1]
       esol = thisraw.esol[*,0:nfitlines-1]
       
       etemplates = thisraw.etemplates[notzero,0:nfitlines-1]
       if (size(etemplates,/n_dim) eq 1) then $
         linefit = etemplates else $
           linefit = total(etemplates,2)

; remove undetected emission lines from both the matrix of
; best-fitting parameters *and* the emission-line templates; demand
; S/N>3 on the amplitude of each line
       bestfit = continuum + linefit + smooth_continuum
       mask = ferr lt 1E5
       new_sol = gandalf_clean_emission(wave,flux/fluxscale,$
         bestfit/fluxscale,mask,etemplates/fluxscale,linepars,sol,$
         esol,snrcut=snrcut_line,velscale=velscale,debug=debug,$
         new_etemplates=new_etemplates,new_linepars=new_linepars)
       new_etemplates = fluxscale*new_etemplates
       
       if (size(new_etemplates,/n_dim) eq 1) then $
         bestlinefit = new_etemplates else $
           bestlinefit = total(new_etemplates,2)
; parse
       elinefit = gandalf_parse_elinefit(new_sol,esol,$
         new_linepars,vsys=vsys,fluxscale=fluxscale)

; measure raw and emission-line corrected spectral indices; do not do
; the inverse-variance weighting!
       ivar = 1.0/ferr^2
       badpixels = where((ferr gt 1E5) or (abs(flux) gt 100.0),nbad)
       if (nbad ne 0) then ivar[badpixels] = 0.0
       indices_raw = spectral_indices(exp(wave),flux,ivar=ivar,$
         debug=0,indexfile=indexfile,_extra=extra,/silent)

       cleanflux = flux - bestlinefit; - smooth_continuum
       indices_cor = spectral_indices(exp(wave),cleanflux,ivar=ivar,$
         debug=0,indexfile=indexfile,_extra=extra,/silent)

; if the continuum coefficients are all zero then just return the
; empty structure
       if (total(thisraw.continuum_coeff,/double) eq 0) then $
         empty_structure = 1 else empty_structure = 0
       modelivar = ferr*0.0+1.0/median(ferr)^2 ; the model errors are illustrative
       indices_model = spectral_indices(exp(wave),continuum,ivar=modelivar,$
         debug=0,indexfile=indexfile,empty_structure=empty_structure,$ ; no weighting
         _extra=extra,/silent)

       newindices = [strtrim(indices_raw.indices,2)+'_raw',$ ; index names
         strtrim(indices_cor.indices,2)+'_cor',$
         strtrim(indices_model.indices,2)+'_model']
       indices_raw = struct_trimtags(temporary(indices_raw),except='INDICES')
       indices_raw = im_struct_trimtags(indices_raw,$
         select=tag_names(indices_raw),newtags=tag_names(indices_raw)+'_raw')

       indices_cor = struct_trimtags(temporary(indices_cor),except='INDICES')
       indices_cor = im_struct_trimtags(indices_cor,$
         select=tag_names(indices_cor),newtags=tag_names(indices_cor)+'_cor')

       indices_model = struct_trimtags(temporary(indices_model),except='INDICES')
       indices_model = im_struct_trimtags(indices_model,$
         select=tag_names(indices_model),newtags=tag_names(indices_model)+'_model')

       indices = struct_addtags(struct_addtags(struct_addtags({indices: newindices},$
         temporary(indices_raw)),temporary(indices_cor)),temporary(indices_model))

; pack into a structure and write out; for now, toss *out* the
; information on the broad Balmer lines; this should be easy enough to
; change if necessary
       specdata1 = struct_trimtags(thisraw,except=['wave','flux',$
         'ferr','continuum','smooth_continuum','etemplates','sol','esol'])
       specdata1 = create_struct(specdata1,'isbroad',0,'lumage',0.0)

; luminosity-weighted age          
       total_weight = total(specdata1.continuum_coeff,/double)
       if (total_weight gt 0) then $
         specdata1.lumage = total(specdata1.continuum_coeff*$
         tempinfo.age,/double)/total_weight ; [Gyr]
       
       if thisbroad then begin
          specdata1.isbroad = 1
; toss out the broad Balmer lines, for now             
          keep = where(strmatch(elinefit.linename,'*broad*',/fold) eq 0)
          linename = strtrim(elinefit.linename[keep],2)
          elinefit = struct_trimtags(elinefit,except=['h_*_broad*','*linename*'])
          elinefit = create_struct('linename',linename,elinefit)
       endif
       specdata1 = struct_addtags(temporary(specdata1),elinefit)
       specdata1 = struct_addtags(temporary(specdata1),indices)
       if (n_elements(specdata) eq 0) then specdata = specdata1 else $
         specdata = [temporary(specdata),specdata1]

       specfit1 = specfit_template
       struct_assign, specdata1, specfit1, /nozero
       specfit1.wave[0:npix-1] = wave ; rest frame
       specfit1.flux[0:npix-1] = flux
       specfit1.ferr[0:npix-1] = ferr
       specfit1.linefit[0:npix-1] = bestlinefit
       specfit1.continuum[0:npix-1] = continuum
       specfit1.smooth_continuum[0:npix-1] = smooth_continuum
       if (n_elements(specfit) eq 0) then specfit = specfit1 else $
         specfit = [temporary(specfit),specfit1]

;      djs_plot, exp(wave), flux, psym=10, xsty=3, ysty=3
;      djs_oplot, exp(wave), continuum, psym=10, color='red'
;      djs_oplot, exp(wave), continuum+bestlinefit, psym=10, color='blue'

    endfor
    splog, 'Total time for all files = ', (systime(1)-t0)/60.0, ' minutes'

return
end    

pro parse_ediscs_gandalf_specfit, debug=debug, solar=solar

; parse RUN12 and R23 separately; then merge and line-match to the
; info structure
    do_parse_ediscs_gandalf_specfit, run34=0, debug=debug, $
      solar=solar, specdata=specdata12, specfit=specfit12
    do_parse_ediscs_gandalf_specfit, run34=1, debug=debug, $
      solar=solar, specdata=specdata34, specfit=specfit34

; line-match
    ediscs = read_ediscs(/spec1d)
    specdata = [specdata12,specdata34]
    specfit = [specfit12,specfit34]
    match, ediscs.ediscs_id, specdata.ediscs_id, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    specdata = specdata[m2]
    specfit = specfit[m2]
    niceprint, ediscs[m1].galaxy, specdata.galaxy, $
      specdata.cluster, ediscs[m1].cluster

; write out!
    version = ediscs_version(/ppxf_specfit)
    specfitpath = ediscs_path(/ppxf)

    if keyword_set(solar) then suffix = '_solar' else suffix = ''
    specdatafile = specfitpath+'ediscs_specdata'+suffix+'_'+version+'.fits'
    specfitfile = repstr(specdatafile,'specdata','specfit')
    
    im_mwrfits, specdata, specdatafile, /clobber
    im_mwrfits, specfit, specfitfile, /clobber

return
end