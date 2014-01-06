;+
; NAME:
;   PARSE_DEEP2_GANDALF_SPECFIT_DR4
;
; PURPOSE:
;   Parse the output from DEEP2_GANDALF_SPECFIT_DR4 and write out
;   the final spectral fits and line-measurements.
;
; INPUTS: 
;   mask - mask number to parse
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
;   J. Moustakas, 2013 Dec 20, Siena 
;
; Copyright (C) 2013, John Moustakas
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

pro parse_deep2_gandalf_specfit_dr4, thismask=thismask, $
  firstmask=firstmask, lastmask=lastmask, debug=debug 

    snrcut_line = 3.0  ; see GANDALF_CLEAN_EMISSION
    light = 2.99792458D5 ; speed of light [km/s]
    velscale = deep2_ppxf_velscale()

; path names and emission-line file name
    version = deep2_version(/ppxf)
    specfitpath = deep2_path(/ppxf,/dr4)
    linefile = specfitpath+'gandalf_elinelist_'+version+'.dat'
    linefile_broad = repstr(linefile,'elinelist_','elinelist_broad_')

; parse the full list of files written by DEEP2_GANDALF_SPECFIT, or a
; subset
    if (n_elements(thismask) eq 0) then begin
       rawfiles = file_search(specfitpath+'specdata_raw_????.fits.gz')
    endif else begin
       rawfiles = specfitpath+'specdata_raw_'+$
         string(thismask,format='(I0)')+'.fits.gz'
    endelse
    nraw = n_elements(rawfiles)
    rawfiles_broad = repstr(rawfiles,'raw_','raw_broad_')

    if (n_elements(firstmask) eq 0L) then firstmask = 0L
    if (n_elements(lastmask) eq 0L) then lastmask = nraw-1L

; output file names
    specdatafile = repstr(repstr(rawfiles,'_raw',''),'.gz','')
    specfitfile = repstr(specdatafile,'specdata','specfit')

; loop on each plate
    t1 = systime(1)
    for ii = firstmask, lastmask do begin
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
;      for iobj = 18, 18 do begin
;      for iobj = 2, 2 do begin
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
          npix = n_elements(wave)
          
          zobj = thisraw.z
          vsys = alog(1.0+zobj)*light ; systemic velocity

; rebuild LINEPARS and the emission-line templates
          junk = read_gandalf_elinelist(thislinefile,$
            fitlines=linepars,/actionfit)

          nfitlines = fix(total(thisraw.fitlines ne -1))
          fitlines = thisraw.fitlines[0:nfitlines-1]
          linepars.action = 'i'
          linepars[fitlines].action = 'f'

          sol = thisraw.sol[*,0:nfitlines-1]
          esol = thisraw.esol[*,0:nfitlines-1]
          etemplates = thisraw.etemplates[notzero,0:nfitlines-1]
          if (size(etemplates,/n_dim) eq 1) then $
            linefit = etemplates else $
              linefit = total(etemplates,2)

; remove undetected emission lines from both the matrix of
; best-fitting parameters *and* the emission-line templates; demand
; S/N>3 on the amplitude of each line
          bestfit = continuum + linefit
          mask = ferr lt 1E5
          new_sol = gandalf_clean_emission(wave,flux,$
            bestfit,mask,etemplates,linepars,sol,$
            esol,snrcut=snrcut_line,velscale=velscale,debug=0,$
            new_etemplates=new_etemplates,new_linepars=new_linepars,$
            fitlines=fitlines)
          new_etemplates = new_etemplates

          if (size(new_etemplates,/n_dim) eq 1) then $
            bestlinefit = new_etemplates else $
              bestlinefit = total(new_etemplates,2)
; parse
          elinefit = gandalf_parse_elinefit(new_sol,esol,$
            new_linepars,vsys=vsys,fluxscale=1.0)

; pack into a structure and write out; for now, toss *out* the
; information on the broad Balmer lines; this should be easy enough to
; change if necessary
          specdata1 = struct_trimtags(thisraw,except=['fitlines','wave','flux',$
            'ferr','continuum','etemplates','sol','esol'])
          specdata1 = create_struct(specdata1,'isbroad',0)

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

          if (n_elements(specdata) eq 0) then specdata = specdata1 else $
            specdata = [temporary(specdata),specdata1]

          nfinalpix = n_elements(thisraw.wave)
          specfit1 = struct_addtags(struct_trimtags(specdata1,$
            select=['galaxy','file','objno','mask','z']),$
            {wave: dblarr(nfinalpix), flux: dblarr(nfinalpix), $
            ferr: dblarr(nfinalpix), linefit: dblarr(nfinalpix), $
            continuum: dblarr(nfinalpix)})

          specfit1.wave[0:npix-1] = wave ; rest frame
          specfit1.flux[0:npix-1] = flux
          specfit1.ferr[0:npix-1] = ferr
          specfit1.linefit[0:npix-1] = bestlinefit
          specfit1.continuum[0:npix-1] = continuum
          if (n_elements(specfit) eq 0) then specfit = specfit1 else $
            specfit = [temporary(specfit),specfit1]
          
;         djs_plot, exp(wave), flux, psym=10, xsty=3, ysty=3
;         djs_oplot, exp(wave), continuum, psym=10, color='red'
;         djs_oplot, exp(wave), continuum+bestlinefit, psym=10, color='blue'
       endfor 

       im_mwrfits, specdata, specdatafile[ii], /clobber
       im_mwrfits, specfit, specfitfile[ii], /clobber
       splog, 'Total time = ', systime(1)-t0, ' sec'

       delvarx, specdata, specfit ; delete memory!!

    endfor 
    splog, 'Total time for all files = ', (systime(1)-t1)/60.0, ' minutes'
       
return
end