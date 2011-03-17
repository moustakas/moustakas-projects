;+
; NAME:
;   PARSE_MPA_LINEFIT()
;
; PURPOSE:
;   Parse a MPA formatted data structure into an iSPEC1d
;   formatted data structure.
;
; INPUTS:
;   line - emission-line data structure outputted from
;          platefit, Christy's SDSS spectral fitting code  
;
; OPTIONAL INPUTS:
;   galindx - absorption-line index data structure outputted from
;             platefit  
;   scale   - flux normalization factor (default 1E-17)
;
; KEYWORD PARAMETERS:
;   silent - suppress messages to STDOUT
;
; OUTPUTS:
;   newline - data structure that is compatible with iSPEC1d 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The [O III] 4959 emission-line flux and equivalent width is
;   set to one-third of the value in the [O III] 5007 line.  
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2005 Jan 26, U of A - written, based on earlier
;                                       code in PARSE_SDSS()
;   jm05jun03uofa - check that TARGETTYPE is a valid structure tag 
;   jm06apr20uofa - significantly streamlined; when computing
;                   EW's, do not propagate the continuum
;                   uncertainty, or the extra flux error scale
;                   factor; the continuum flux errors, especially
;                   for the [OII] doublet may need to be reduced
;                   by sqrt(N) 
;   jm08feb27nyu  - GALINDX was made an optional input
;   jm08mar21nyu  - bug! do not assume that 4959 is 0.33*5007
;                   the *observed* ratio depends weakly on the
;                   reddening, so do this in IM_ABUNDANCE()
;   jm08sep24nyu  - renamed PARSE_MPA_LINEFIT() from
;     PARSE_TREMONTI_LINEFIT() 
;   jm10mar07ucsd - add the _sub tags (and rename them _cor), which
;     are the indices measured after subtracting off emission lines
;   jm10dec14ucsd - properly deal with upper limits
;
; Copyright (C) 2005-2006, 2008, 2010, John Moustakas
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

function parse_mpa_linefit, line, galindx=galindx, scale=scale, $
  sdss_id=sdss_id, silent=silent

    ngal = n_elements(line)
    if (ngal eq 0L) then begin
       doc_library, 'parse_mpa_linefit'
       return, -1L
    endif

    if arg_present(galindx) then begin
       if (ngal ne n_elements(galindx)) then begin
          splog, 'LINE and GALINDX have incompatible dimensions'
          return, -1L
       endif
    endif

    light = im_light() ; [km/s]
    oratio = 2.984 ; intrinsic [O III] doublet ratio

    if (n_elements(scale) eq 0L) then scale = 1D-17
    
    linetags = [$
      'oii_3726','oii_3729','neiii_3869','h_delta','h_gamma','oiii_4363','h_beta',$
      'oiii_4959','oiii_5007','oi_6300','nii_6548','h_alpha','nii_6584',$
      'sii_6717','sii_6731'$
      ]
    nline = n_elements(linetags)

    fluxtags = [linetags+'_flux',linetags+'_flux_err']
    conttags = [linetags+'_cont',linetags+'_cont_err']
    instrestags = linetags+'_inst_res'

    newfluxtags = [$
      'oii_3727','neiii_3869','h_delta','h_gamma','oiii_4363','h_beta',$
      'oiii_4959','oiii_5007','oi_6300','nii_6548','h_alpha','nii_6584',$
      'sii_6716','sii_6731']
    newconttags = newfluxtags+'_continuum'
    newewtags = newfluxtags+'_ew'
    newlimittags = newfluxtags+'_limit'
    newewlimittags = newfluxtags+'_ew_limit'
    
    linewave = [$
      3727.4235,3869.06,4101.734,4340.464,4363.21,4861.325,$
      4958.911,5006.843,6300.304,6548.04,6562.8,6583.46,$
      6716.14,6730.81]
;   niceprint, newfluxtags, linewave

    fluxcut = struct_trimtags(line,select=fluxtags)
    contcut = struct_trimtags(line,select=conttags)
    instrescut = struct_trimtags(line,select=instrestags)

    ntags = n_elements(newfluxtags)

; --------------------------------------------------
; RELEGATED!    
;   newfluxtags2 = [newfluxtags,'oiii_4959'] ; used below
;   nnewfluxtags2 = n_elements(newfluxtags2)
; --------------------------------------------------
    linenamestruct = replicate({linename: newfluxtags},ngal)

; combine the [O II] doublet    
    
    indx = 0L
    newflux = mrd_struct(newfluxtags,replicate('fltarr(2)-1.0',ntags),ngal)
    newcont = mrd_struct(newconttags,replicate('fltarr(2)-1.0',ntags),ngal)
    newlimit = mrd_struct(newlimittags,replicate('-1.0',ntags),ngal)
    newewlimit = mrd_struct(newewlimittags,replicate('-1.0',ntags),ngal)
    newew = mrd_struct(newewtags,replicate('fltarr(2)-1.0',ntags),ngal)

    if (not keyword_set(silent)) then splog, 'Parsing the SDSS fluxes and equivalent widths'

; assumed parsing rules: if the flux *error* is <0 then the line was
; not measured (error=-2); if the line-amplitude is <3*rms, where rms
; is the error in the continuum around the line, then drop the line
; (error=-1); for the latter case compute the limiting 1-sigma flux
; assuming a Gaussian line-profile
    sigma_avg = (line.sigma_forbidden+line.sigma_balmer)/2.0
;   zero = where(sigma_avg le 0.0,nzero)
;   if (nzero ne 0L) then sigma_avg[zero] = 100.0 ; [km/s]
    
    snrcut = 3.0
    for i = 1L, nline-1L do begin ; offset from the first [O II] line

; special case: [O II]
       if (i eq 1L) then begin 
          flux1 = fluxcut.(i)    & ferr1 = fluxcut.(i+nline)
          flux2 = fluxcut.(i-1L) & ferr2 = fluxcut.((i-1L)+nline)

          cont1 = contcut.(i)    & conterr1 = contcut.(i+nline)      
          cont2 = contcut.(i-1L) & conterr2 = contcut.((i-1L)+nline) 

          instres = (instrescut.(i)+instrescut.(i-1))/2.0 ; instrumental resolution [km/s]
          
          miss = where((ferr1 le 0.0) or (ferr2 le 0.0) or $ ; not measured
            (conterr1 le 0.0) or (conterr2 le 0.0),nmiss)            
          if (nmiss ne 0L) then newflux[miss].(indx) = [0.0,-2.0]          

          good = where((ferr1 gt 0.0) and (ferr2 gt 0.0) and $
            (conterr1 gt 0.0) and (conterr2 gt 0.0),ngood)
          if (ngood ne 0L) then begin
;            check = where((conterr1[good] le 0) or (conterr2[good] le 0))
;            if (check[0] ne -1L) then message, 'Bad, bad, bad'

             newflux[good].(indx) = scale*transpose([[(flux1[good]+flux2[good])],$ ; total
               [(sqrt(ferr1[good]^2+ferr2[good]^2))]])
             newcont[good].(indx) = scale*transpose([[(cont1[good]+cont2[good])/2.0],$ ; average
               [(sqrt(conterr1[good]^2+conterr2[good]^2))/sqrt(2.0)]]) 

             cnorm = (newcont[good].(indx))[0,*] ; no extra uncertainty
             newew[good].(indx) = newflux[good].(indx)/rebin(cnorm,2,ngood)
          endif

       endif else begin ; all other lines

          flux1 = fluxcut.(i)    
          ferr1 = fluxcut.(i+nline)

          cont1 = contcut.(i)    
          conterr1 = contcut.(i+nline)      

          instres = instrescut.(i) ; instrumental resolution [km/s]
          
          miss = where((ferr1 le 0.0) or (conterr1 le 0.0),nmiss) ; not measured
          if (nmiss ne 0L) then newflux[miss].(indx) = [0.0,-2.0]

          good = where((ferr1 gt 0.0) and (conterr1 gt 0.0),ngood)
          if (ngood ne 0L) then begin
;            check = where((conterr1[good] le 0))
;            if (check[0] ne -1L) then message, 'Bad, bad, bad'

             newflux[good].(indx) = scale*transpose([[flux1[good]],[ferr1[good]]])
             newcont[good].(indx) = scale*transpose([[cont1[good]],[conterr1[good]]])

             cnorm = (newcont[good].(indx))[0,*] ; no extra uncertainty
             newew[good].(indx) = newflux[good].(indx)/rebin(cnorm,2,ngood)
          endif
       endelse   

; now check for significant detections based on the line-amplitude;
; should we keep the original flux measurements!?!
       good = where(((newflux.(indx))[1,*] ne -2.0),ngood)
       if (ngood ne 0L) then begin
          lflux = reform((newflux[good].(indx))[0,*])
          conterr = reform((newcont[good].(indx))[1,*])
          
          sigma = sqrt(sigma_avg[good]^2.0+instres[good]^2)*linewave[i-1]/light ; [A]
          amp = lflux/(sqrt(2.0*!dpi)*sigma) ; [erg/s/cm2/A]

          det = where(amp gt snrcut*conterr,ndet,comp=lim,ncomp=nlim)
          if (nlim ne 0L) then newflux[good[lim]].(indx) = [0.0,-1.0] 

          newlimit[good].(indx) = sqrt(2.0*!dpi)*conterr*sigma ; [erg/s/cm2]
          newewlimit[good].(indx) = newlimit[good].(indx)/reform(cnorm) ; [A]
       endif
       indx++ ; increment the line
    endfor

; --------------------------------------------------
; RELEGATED!    
;; append the [O III] 4959 emission-line flux
;
;   newflux = struct_addtags(temporary(newflux),replicate({oiii_4959: [0.0,-2.0]},ngal))
;   newflux.oiii_4959 = newflux.oiii_5007
;   good = where(newflux.oiii_5007[1] gt 0.0,ngood)
;   if (ngood ne 0L) then newflux[good].oiii_4959 = newflux[good].oiii_4959/oratio
;
;; append the [O III] 4959 emission-line continuum; assume that the
;; continuum between [O III] 5007 and [O III] 4959 is identical
;
;    newcont = struct_addtags(temporary(newcont),replicate({oiii_4959_continuum: [0.0,-2.0]},ngal))
;    newcont.oiii_4959_continuum = newcont.oiii_5007_continuum
;
;; append the [O III] 4959 emission-line EW; assume that the continuum
;; between [O III] 5007 and [O III] 4959 is identical
;
;    newew = struct_addtags(temporary(newew),replicate({oiii_4959_ew: [0.0,-2.0]},ngal))
;    newew.oiii_4959_ew = newew.oiii_5007_ew
;
;    good = where(newew.oiii_5007_ew[1] gt 0.0,ngood)
;    if (ngood ne 0L) then newew[good].oiii_4959_ew = newew[good].oiii_4959_ew/oratio
; --------------------------------------------------

; increase the line flux errors according to
; http://www.mpa-garching.mpg.de/SDSS/DR4/raw_data.html

    errorfactor = [$
      2.199,$ ; [O II] 3727
      1.731,$ ; [Ne III] 3869
      1.400,$ ; H-delta (private communication)
      1.600,$ ; H-gamma (private communication)
      1.306,$ ; [O III] 4363
      1.882,$ ; H-beta
      1.573,$ ; [O III] 4959
      1.566,$ ; [O III] 5007
      1.378,$ ; [O I] 6300
      2.039,$ ; [N II] 6548 (assumed)
      2.473,$ ; H-alpha
      2.039,$ ; [N II] 6584
      1.621,$ ; [S II] 6716
      1.621]  ; [S II] 6731
    errorfactor = errorfactor*0.0 + 1.0 ; actually, don't do it; jm06nov16nyu

    tags = tag_names(newflux)

    if (not keyword_set(silent)) then splog, 'Scaling the flux uncertainties'
    
    for i = 0L, n_elements(newfluxtags)-1L do begin
;   for i = 0L, n_elements(newfluxtags2)-1L do begin

       factor = rebin([1.0,errorfactor[i]],2,ngal)
       match = tag_exist(newflux[0],newfluxtags[i],index=tagindx)
;      match = tag_exist(newflux[0],newfluxtags2[i],index=tagindx)

       theflux = newflux.(tagindx)
       good = where(theflux[1,*] gt 0.0,ngood)
       if (ngood ne 0L) then begin
          theflux[*,good] = theflux[*,good]*factor[*,good]
          newflux.(tagindx) = theflux
       endif
    endfor

; parse the Balmer absorption-line measurements, if they exist
    if (not keyword_set(silent)) then splog, 'Parsing Balmer absorption-line equivalent widths'

    linetags = tag_names(line[0])
    babstags = ['H_DELTA','H_GAMMA','H_BETA','H_ALPHA']
    babswave = [4101.73,4340.46,4861.33,6562.80]
    nbabs = n_elements(babstags)

    babs = {babs_linename: babstags}
    for iabs = 0L, nbabs-1L do babs = create_struct(babs,$
      'BABS_'+babstags[iabs]+'_WAVE',      babswave[iabs], $
      'BABS_'+babstags[iabs],              [0.0,-2.0], $
      'BABS_'+babstags[iabs]+'_EW',        [0.0,-2.0], $
      'BABS_'+babstags[iabs]+'_CONTINUUM', [0.0,-2.0])
    babs = replicate(babs,ngal)

    for iabs = 0L, nbabs-1L do begin

       if tag_exist(line,babstags[iabs]+'_CONT') and $
         tag_exist(line,babstags[iabs]+'_CONT_ERR') then begin

          cont = scale*line.(where(babstags[iabs]+'_CONT' eq linetags))
          cont_err = scale*line.(where(babstags[iabs]+'_CONT_ERR' eq linetags))

          zero = where((cont eq 0.0) or (cont_err eq 0.0),nzero)
          if (nzero ne 0L) then begin
             cont[zero] = 0.0
             cont_err[zero] = -1.0
          endif

          reqw = line.(where(babstags[iabs]+'_REQW' eq linetags))
          reqw_err = line.(where(babstags[iabs]+'_REQW_ERR' eq linetags))

          eqw = line.(where(babstags[iabs]+'_EQW' eq linetags))
          eqw_err = line.(where(babstags[iabs]+'_EQW_ERR' eq linetags))

          ew = reqw - eqw
          ew_err = sqrt(reqw_err^2 + eqw_err^2)

          ewtag = where('BABS_'+babstags[iabs]+'_EW' eq tag_names(babs[0]))
          conttag = where('BABS_'+babstags[iabs]+'_CONTINUUM' eq tag_names(babs[0]))
          fluxtag = where('BABS_'+babstags[iabs] eq tag_names(babs[0]))
          
          babs.(ewtag) = transpose([ [ew], [ew_err] ])
          babs.(conttag) = transpose([ [cont], [cont_err] ])
          babs.(fluxtag) = babs.(ewtag)*babs.(conttag)
       endif
    endfor

; wavelength structure    
    wavestr = mrd_struct(newfluxtags+'_wave',replicate('0.0',ntags),ngal)
    for i = 0L, ntags-1L do wavestr.(i) = linewave[i]
;   wavestr = mrd_struct(newfluxtags2+'_wave',replicate('0.0',nnewfluxtags2),ngal)
;   for i = 0L, nnewfluxtags2-1L do wavestr.(i) = linewave[i]

; finally put everything together
;   splog, 'Setting LINE to zero!'
    line = 0                    ; NOTE!
;   newline = struct_addtags(struct_addtags(struct_addtags(temporary(linenamestruct),$
;     struct_addtags(struct_addtags(temporary(newflux),temporary(newcont)),$
;     temporary(newew))),temporary(wavestr)),temporary(babs))

; leave off the continuum and Balmer-absorption tags
    newline = struct_addtags(temporary(linenamestruct),temporary(newflux))
    newline = struct_addtags(temporary(newline),temporary(newcont))
    newline = struct_addtags(temporary(newline),temporary(newew))
    newline = struct_addtags(temporary(newline),temporary(newlimit))
    newline = struct_addtags(temporary(newline),temporary(newewlimit))
    newline = struct_addtags(temporary(newline),temporary(wavestr))

; ---------------------------------------------------------------------------    
; parse the lick index measurements
    if arg_present(galindx) then begin
       
       if (not keyword_set(silent)) then splog, 'Parsing Lick indices'

       lick = {$
         D4000:              [0.0,0.0], $
         D4000_cor:          [0.0,0.0], $
         D4000_model:        [0.0,0.0], $

         D4000_narrow:       [0.0,0.0], $
         D4000_narrow_cor:   [0.0,0.0], $
         D4000_narrow_model: [0.0,0.0], $

         lick_hg_a:          [0.0,0.0], $
         lick_hg_a_cor:      [0.0,0.0], $
         lick_hg_a_model:    [0.0,0.0], $

         lick_hd_a:          [0.0,0.0], $
         lick_hd_a_cor:      [0.0,0.0], $
         lick_hd_a_model:    [0.0,0.0]}
       lick = replicate(lick,ngal)
       
       if tag_exist(galindx,'D4000') and tag_exist(galindx,'D4000_ERR') then $
         lick.D4000 = transpose([ [galindx.D4000],[galindx.D4000_err] ])
       if tag_exist(galindx,'D4000_sub') and tag_exist(galindx,'D4000_sub_ERR') then $
         lick.D4000_cor = transpose([ [galindx.D4000_sub],[galindx.D4000_sub_err] ])
       if tag_exist(galindx,'D4000_MODEL') and tag_exist(galindx,'D4000_ERR') then $
         lick.D4000_model = transpose([ [galindx.D4000_model],[galindx.D4000_err] ])

       if tag_exist(galindx,'D4000_N') and tag_exist(galindx,'D4000_N_ERR') then $
         lick.D4000_narrow = transpose([ [galindx.D4000_N],[galindx.D4000_N_err] ])
       if tag_exist(galindx,'D4000_N_sub') and tag_exist(galindx,'D4000_N_sub_ERR') then $
         lick.D4000_narrow_cor = transpose([ [galindx.D4000_N_sub],[galindx.D4000_N_sub_err] ])
       if tag_exist(galindx,'D4000_N_MODEL') and tag_exist(galindx,'D4000_N_ERR') then $
         lick.D4000_narrow_model = transpose([ [galindx.D4000_N_model],[galindx.D4000_N_err] ])

       if tag_exist(galindx,'LICK_HG_A') and tag_exist(galindx,'LICK_HG_A_ERR') then $
         lick.lick_hg_a = transpose([ [galindx.lick_hg_a],[galindx.lick_hg_a_err] ])
       if tag_exist(galindx,'LICK_HG_A_sub') and tag_exist(galindx,'LICK_HG_A_sub_ERR') then $
         lick.lick_hg_a_cor = transpose([ [galindx.lick_hg_a_sub],[galindx.lick_hg_a_sub_err] ])
       if tag_exist(galindx,'LICK_HG_A_MODEL') and tag_exist(galindx,'LICK_HG_A_ERR') then $
         lick.lick_hg_a_model = transpose([ [galindx.lick_hg_a_model],[galindx.lick_hg_a_err] ])

       if tag_exist(galindx,'LICK_HD_A') and tag_exist(galindx,'LICK_HD_A_ERR') then $
         lick.lick_hd_a = transpose([ [galindx.lick_hd_a],[galindx.lick_hd_a_err] ])
       if tag_exist(galindx,'LICK_HD_A_sub') and tag_exist(galindx,'LICK_HD_A_sub_ERR') then $
         lick.lick_hd_a_cor = transpose([ [galindx.lick_hd_a_sub],[galindx.lick_hd_a_sub_err] ])
       if tag_exist(galindx,'LICK_HD_A_MODEL') and tag_exist(galindx,'LICK_HD_A_ERR') then $
         lick.lick_hd_a_model = transpose([ [galindx.lick_hd_a_model],[galindx.lick_hd_a_err] ])
       
       newline = struct_addtags(temporary(newline),temporary(lick))
    endif
       
; prepend any other tags of interest

;   prepend = {$
;     galaxy:   '', $
;     sdss_id: -1L, $
;     linename: newfluxtags2 $
;     }
;   prepend = replicate(prepend,ngal)

;   if (n_elements(sdss_id) ne 0L) then prepend.sdss_id = sdss_id
    
;   newline = struct_addtags(prepend,newline)
;   if tag_exist(line[0],'TARGETTYPE') then newline.galaxy = strcompress(line.targettype,/remove)
    
return, newline
end    
