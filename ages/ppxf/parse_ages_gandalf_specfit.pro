;+
; NAME:
;   PARSE_AGES_GANDALF_SPECFIT
;
; PURPOSE:
;   Parse the output from AGES_GANDALF_SPECFIT and write out the final
;   spectral fits and line-measurements.
;
; INPUTS: 
;   pass - pass number to parse
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;   test - test parsing
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

pro parse_ages_gandalf_specfit, pass1, solar=solar, debug=debug

    snrcut_line = 3.0  ; see GANDALF_CLEAN_EMISSION
    light = 2.99792458D5 ; speed of light [km/s]
    velscale = ages_ppxf_velscale()
    fluxscale = ages_ppxf_fluxscale()

; path names and emission-line file name
    version = ages_version(/ppxf_specfit)
    specfitpath = ages_path(/spec1d)+'fluxed/tweak/'+version+'/' 
    indexfile = ages_path(/ppxf)+'lick_indexlist_'+version+'.dat'
    linefile = ages_path(/ppxf)+'gandalf_elinelist_'+version+'.dat'
    linefile_broad = repstr(linefile,'elinelist_','elinelist_broad_')

; parse the full list of files written by AGES_GANDALF_SPECFIT, or a
; subset
    if keyword_set(solar) then suffix = '_solar' else suffix = ''
    if (n_elements(pass1) eq 0) then begin
       rawfiles = file_search(specfitpath+'specdata_raw_???'+$
         suffix+'.fits.gz')
    endif else begin
       rawfiles = specfitpath+'specdata_raw_'+$
         string(pass1,format='(I3.3)')+suffix+'.fits.gz'
    endelse
    nraw = n_elements(rawfiles)

    rawfiles_broad = repstr(rawfiles,'raw_','raw_broad_')
    
; output file names
    specdatafile = repstr(repstr(rawfiles,'_raw',''),'.gz','')
    specfitfile = repstr(specdatafile,'specdata','specfit')

; read the model templates so that we can compute the
; luminosity-weighted age for each object
    junk = read_ages_ppxf_templates(tempinfo=tempinfo,solar=solar)
    
; initialize the output data structures
    nfinalpix = 5000
    specfit_template = {$
      ages_id:   0L, $
      pass:       0, $
      aper:       0, $
      z:        0.0, $
      zabs:     0.0, $
      vdisp:    0.0, $
      wave:          dblarr(nfinalpix), $
      flux:          fltarr(nfinalpix), $
      ferr:          fltarr(nfinalpix), $
      linefit:       fltarr(nfinalpix), $
      continuum:     fltarr(nfinalpix), $
      smooth_continuum: fltarr(nfinalpix)}

; loop on each plate
    t1 = systime(1)
    for ii = 0, nraw-1 do begin
       if (file_test(rawfiles[ii],/reg) eq 0) then begin
          splog, 'File '+rawfiles[ii]+' not found'
          continue
       endif
       splog, 'Reading '+rawfiles[ii]
       raw = mrdfits(rawfiles[ii],1,/silent)
       nobj = n_elements(raw)

; is there a data structure of broad-line AGN?
       if file_test(rawfiles_broad[ii],/reg) then begin
          splog, 'Reading '+rawfiles_broad[ii]
          raw_broad = mrdfits(rawfiles_broad[ii],1,/silent)
          broad = 1
       endif else broad = 0

; loop on each object
       t0 = systime(1)
;      for iobj = 178, 178 do begin
;      for iobj = 3, 3 do begin
       for iobj = 0, nobj-1 do begin
;         print, iobj, nobj
          if broad then begin
             match = where(raw[iobj].aper eq raw_broad.aper,nmatch)
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
          wave = thisraw.wave[notzero] ; rest
          flux = thisraw.flux[notzero]*1.0D ; note!
          ferr = thisraw.ferr[notzero]*1.0D 
          continuum = thisraw.continuum[notzero]
          smooth_continuum = thisraw.smooth_continuum[notzero]
          npix = n_elements(wave)
          
          zabs = thisraw.zabs
          vsys = alog(1.0+zabs)*light ; systemic velocity

; rebuild LINEPARS and the emission-line templates
          junk = read_gandalf_elinelist(thislinefile,$
            fitlines=linepars,/actionfit)
          nfitlines = fix(total(thisraw.fitlines ne -1))
          fitlines = thisraw.fitlines[0:nfitlines-1]

          linepars.action = 'i'
          linepars[fitlines].action = 'f'
          tied = where(strmatch(linepars.fit_iter2,'t*') and strmatch(linepars.kind,'d*'),ntied)
          if (ntied ne 0) then linepars[tied].action = $
            linepars[strmid(linepars[tied].fit_iter2,1)].action
          
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
            new_etemplates=new_etemplates,new_linepars=new_linepars,$
            fitlines=fitlines)
          new_etemplates = fluxscale*new_etemplates
          
          if (size(new_etemplates,/n_dim) eq 1) then $
            bestlinefit = new_etemplates else $
              bestlinefit = total(new_etemplates,2)
; parse
          elinefit = gandalf_parse_elinefit(new_sol,esol,$
            new_linepars,vsys=vsys,fluxscale=fluxscale)
;         help, new_linepars.name, elinefit.linename & print, iobj; & stop

; measure spectral indices; if the continuum coefficients are all zero
; then just return the empty structure
          cleanflux = flux - bestlinefit - smooth_continuum
          mask = ferr lt 1E5
          indices = spectral_indices(exp(wave),cleanflux,$
            ivar=mask/ferr^2,debug=0,indexfile=indexfile,$
            _extra=extra,/silent,/nolick)

          if (total(thisraw.continuum_coeff,/double) eq 0) then $
            empty_structure = 1 else empty_structure = 0
          modelivar = ferr*0.0+1.0/median(ferr)^2 ; the model errors are illustrative
          model_indices = spectral_indices(exp(wave),continuum,$
            ivar=modelivar,debug=0,indexfile=indexfile,$
            empty_structure=empty_structure,_extra=extra,/silent,/nolick)

          newindices = [indices.indices,strtrim(model_indices.indices,2)+'_model']
          indices = struct_addtags({indices: newindices},$
            struct_trimtags(indices,except='INDICES'))
          indices = struct_addtags({indices: newindices},$
            struct_addtags(struct_trimtags(indices,$
            except='INDICES'),struct_trimtags(im_struct_trimtags($
            model_indices,select=tag_names(model_indices),$
            newtags=tag_names(model_indices)+'_model'),except='INDICES_MODEL')))

; pack into a structure and write out; for now, toss *out* the
; information on the broad Balmer lines; this should be easy enough to
; change if necessary
          specdata1 = struct_trimtags(thisraw,except=['fitlines','wave','flux',$
            'ferr','continuum','smooth_continuum','etemplates','sol','esol'])
          specdata1 = create_struct(specdata1,'isbroad',0,'lumage',0.0)

; luminosity-weighted age          
          total_weight = total(specdata1.continuum_coeff,/double)
          if (total_weight gt 0) then $
            specdata1.lumage = total(specdata1.continuum_coeff*$
            tempinfo.age,/double)/total_weight ; [Gyr]
             
          if thisbroad then begin
             specdata1.isbroad = 1
; rename the Mg II tags
             tags = tag_names(elinefit)
             mgii = where(strmatch(tags,'*mgii*broad*',/fold),nmgii)
             if (nmgii ne 0) then begin
                newtags = tags
                newtags[mgii] = repstr(newtags[mgii],'_BROAD','') ; case-sensitive
                elinefit = im_struct_trimtags(elinefit,select=tags,newtags=newtags)
                mgname = where(strmatch(elinefit.linename,'*mgii*broad*',/fold))
                elinefit.linename[mgname] = repstr(elinefit.linename[mgname],'_broad','')
             endif
; toss out the Balmer lines, for now             
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
          
;         djs_plot, exp(wave), flux, psym=10, xsty=3, ysty=3
;         djs_oplot, exp(wave), continuum, psym=10, color='red'
;         djs_oplot, exp(wave), continuum+bestlinefit, psym=10, color='blue'
       endfor 

       im_mwrfits, specdata, specdatafile[ii], /clobber
       im_mwrfits, specfit, specfitfile[ii], /clobber
       splog, 'Total time = ', (systime(1)-t0)/60.0, ' minutes'

       delvarx, specdata, specfit ; delete memory!!

    endfor 
    splog, 'Total time for all files = ', (systime(1)-t1)/60.0, ' minutes'
       
return
end
