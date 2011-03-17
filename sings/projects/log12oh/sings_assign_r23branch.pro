;+
; NAME:
;   SINGS_ASSIGN_R23BRANCH()
;
; PURPOSE:
;   Given the output from IM_ABUNDANCE() assign the appropriate
;   R23 branch. 
;
; INPUTS:
;   abund - output structure from IM_ABUNDANCE() or SINGS_ABUNDANCE() 
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;   kk04 - consider the Kobulnicky & Kewley (2004) abundances
;          (default) 
;   pt05 - consider the Pilyugin & Thuan (2005) abundances
;   m91  - consider the McGaugh (1991) abundances
;
; OUTPUTS:
;   result - 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2006 Apr 18, U of A - written
;   jm07sep18nyu - added JUSTEW and JUSTFLUX keywords
;   jm08feb11nyu - added FRACCUT optional input and SILENT keyword
;   jm08oct22nyu - FRACCUT is now obsolete because the formal
;     1-sigma calculation is done (see SINGS paper)
;
; Copyright (C) 2006-2008, John Moustakas
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

function sings_assign_r23branch, abund, r23branch=r23branch, justew=justew, $
  justflux=justflux, kk04=kk04, pt05=pt05, m91=m91, silent=silent, $
  debug=debug, _extra=extra

    nobj = n_elements(abund)
    if (nobj eq 0L) then begin
       doc_library, 'sings_assign_r23branch'
       return, -1L
    endif
    
    if (n_elements(kk04) eq 0L) and (n_elements(pt05) eq 0L) and $
      (n_elements(m91) eq 0L) then kk04 = 1L
    if keyword_set(kk04) then begin
       suffix = 'KK04' & newsuffix = suffix
    endif
    if keyword_set(pt05) then begin
       suffix = 'PT05' & newsuffix = suffix
    endif
    if keyword_set(m91) then begin
       suffix = 'M91' & newsuffix = suffix
    endif

; default: compute abundances based on both EWS and FLUXES 
    if (keyword_set(justew) eq 0) and (keyword_set(justflux) eq 0) then begin
       justew = 1
       justflux = 1
    endif
    
    result = {$
      ewalpha:                1.0,$
      r23branch:              '?',$
      r23branch_ew:           '?',$
      zstrong_niiha:       -999.0,$
      zstrong_niiha_err:   -999.0,$
      zstrong_niioii:      -999.0,$
      zstrong_niioii_err:  -999.0,$
      zstrong_o32:         -999.0,$
      zstrong_o32_err:     -999.0,$
      zstrong_ew_o32:      -999.0,$
      zstrong_ew_o32_err:  -999.0,$
      zstrong_r23:         -999.0,$
      zstrong_r23_err:     -999.0,$
      zstrong_ew_r23:      -999.0,$
      zstrong_ew_r23_err:  -999.0,$
      zstrong_p:           -999.0,$
      zstrong_p_err:       -999.0,$
      zstrong_ew_p:        -999.0,$
      zstrong_ew_p_err:    -999.0,$
      zstrong_logu:        -999.0,$
      zstrong_logu_err:    -999.0,$
      zstrong_ew_logu:     -999.0,$
      zstrong_ew_logu_err: -999.0,$
      zstrong_12oh:        -999.0,$
      zstrong_12oh_err:    -999.0,$
      zstrong_ew_12oh:     -999.0,$
      zstrong_ew_12oh_err: -999.0}
    result = replicate(result,nobj)

    struct_assign, abund, result, /nozero
    ngalaxy = n_elements(result)

; assign the appropriate suffix to the output tags    
    finaltags = tag_names(result[0])
    alter = where(strmatch(finaltags,'*12oh*',/fold) or $
      strmatch(finaltags,'*logu*',/fold) or $
      strmatch(finaltags,'*r23branch*',/fold))
    finaltags[alter] = repstr(finaltags[alter]+'_'+$
      newsuffix,'_ERR_'+newsuffix,'_'+newsuffix+'_ERR')

    data = struct_trimtags(abund,select=['ZSTRONG_12OH_'+$
      suffix+'_*','ZSTRONG_EW_12OH_'+suffix+'_*',$
      'ZSTRONG_LOGU_'+suffix+'_*','ZSTRONG_EW_LOGU_'+suffix+'_*'])
    data = im_struct_trimtags(data,select=tag_names(data),$
      newtags=repstr(tag_names(data),suffix+'_',''))

; #########################
; fluxes    
    if keyword_set(justflux) then begin
       if (keyword_set(silent) eq 0) then splog, 'Assigning R23 branches from fluxes:'
       up = where((data.zstrong_12oh_upper gt -900.0) and (r23branch eq 'U'),nup)
       if (nup ne 0L) then begin
          result[up].r23branch        = 'U'
          result[up].zstrong_12oh     = data[up].zstrong_12oh_upper
          result[up].zstrong_12oh_err = data[up].zstrong_12oh_upper_err
          if keyword_set(kk04) then begin
             result[up].zstrong_logu     = data[up].zstrong_logu_upper
             result[up].zstrong_logu_err = data[up].zstrong_logu_upper_err
          endif
       endif
       
       lo = where((data.zstrong_12oh_lower gt -900.0) and (r23branch eq 'L'),nlo)
       if (nlo ne 0L) then begin
          result[lo].r23branch        = 'L'
          result[lo].zstrong_12oh     = data[lo].zstrong_12oh_lower
          result[lo].zstrong_12oh_err = data[lo].zstrong_12oh_lower_err
          if keyword_set(kk04) then begin
             result[lo].zstrong_logu     = data[lo].zstrong_logu_lower
             result[lo].zstrong_logu_err = data[lo].zstrong_logu_lower_err
          endif
       endif
       
; adopt the "average" lower/upper abundance of objects with ambiguous
; abundances, near the turn-around region
       ambig = where((data.zstrong_12oh_upper gt -900) and (data.zstrong_12oh_lower gt -900) and $
         ((data.zstrong_12oh_upper lt data.zstrong_12oh_lower)),nambig) ; off the R23 calibration
       if (keyword_set(silent) eq 0) then splog, '   Ambiguous: (O/H)_lower>(O/H)_upper: '+$
         string(nambig,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' ('+$
         strtrim(string(100.0*nambig/float(ngalaxy),format='(F12.1)'),2)+'%)'

; reject objects whose upper- and lower-branch abundances differ by
; more than 1-sigma, where 1-sigma is given by the width of the
; overlapping histogram (12oh_avg)
       if (nambig ne 0L) then begin
          good = where($
            (data[ambig].zstrong_12oh_upper+data[ambig].zstrong_12oh_upper_err) gt $
            (data[ambig].zstrong_12oh_lower-data[ambig].zstrong_12oh_lower_err),ngood,$
            comp=reject,ncomp=nreject)
          if (keyword_set(silent) eq 0) then splog, '     Retain: '+$
            string(ngood,format='(I0)')+'/'+string(nambig,format='(I0)')+' ('+$
            strtrim(string(100.0*ngood/float(nambig),format='(F12.1)'),2)+'%)'
          if (ngood ne 0L) then begin
             result[ambig[good]].r23branch        = 'A'
             result[ambig[good]].zstrong_12oh     = data[ambig[good]].zstrong_12oh_avg
             result[ambig[good]].zstrong_12oh_err = data[ambig[good]].zstrong_12oh_avg_err
             if keyword_set(kk04) then begin
                result[ambig[good]].zstrong_logu     = data[ambig[good]].zstrong_logu_avg
                result[ambig[good]].zstrong_logu_err = data[ambig[good]].zstrong_logu_avg_err
             endif
          endif
          if (keyword_set(silent) eq 0) then splog, '     Reject: '+$
            string(nreject,format='(I0)')+'/'+string(nambig,format='(I0)')+' ('+$
            strtrim(string(100.0*nreject/float(nambig),format='(I0)'),2)+'%)'
          if (nreject ne 0L) then begin
             result[ambig[reject]].r23branch        = 'Rejected'
             result[ambig[reject]].zstrong_12oh     = -999.0
             result[ambig[reject]].zstrong_12oh_err = -999.0
             result[ambig[reject]].zstrong_logu     = -999.0
             result[ambig[reject]].zstrong_logu_err = -999.0
          endif
       endif

; build the QAplot, if requested
       if keyword_set(debug) then begin
          djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
            yrange=[7,9.5], xrange=[-0.5,1.3], _extra=extra
          linestyle = [0,2,3,4,5]
          if keyword_set(kk04) then begin ; KK04 models
             model_logq = alog10([1.2D8,8D7,4D7,2D7,1D7]) ; alog10(4E7)
             model_logr23 = range(-0.5,1.1,1500)
             for iq = 0L, n_elements(model_logq)-1L do begin
                model_logoh_upper = 9.72D - 0.777*model_logr23 - $
                  0.951*model_logr23^2 - 0.072*model_logr23^3 - $
                  0.811*model_logr23^4 - model_logq[iq]*(0.0737 - $
                  0.0713*model_logr23 - 0.141*model_logr23^2 + $
                  0.0373*model_logr23^3 - 0.058*model_logr23^4)
                model_logoh_lower = 9.40D + 4.65D*model_logr23 - $
                  3.17D*model_logr23^2 - model_logq[iq]*$
                  (0.272D + 0.547D*model_logr23 - 0.513D*model_logr23^2)
                model_good1 = where((model_logoh_upper gt model_logoh_lower))
                model_good2 = where((model_logoh_lower[model_good1] gt 7.5))
                djs_oplot, model_logr23[model_good1], model_logoh_upper[model_good1], $
                  linestyle=linestyle[iq]
                djs_oplot, model_logr23[model_good1], model_logoh_lower[model_good1], $
                  linestyle=linestyle[iq]
             endfor
             legend, 'log(U)='+string(model_logq-alog10(im_light(/cm)),format='(F5.2)'), $
               /right, /bottom, box=0, charsize=1.3, line=linestyle, pspacing=1.4
          endif
          if keyword_set(pt05) then begin ; PT05 models
             model_r23 = range(0.1,10.0,1500)
             model_p = [0.1,0.3,0.5,0.8,1.0]
;            model_p = range(0.1,0.85,5)
             for ip = 0L, n_elements(model_p)-1L do begin
                model_logoh_lower = (model_r23 + 106.4 + 106.8*model_p[ip] - 3.40*model_p[ip]^2) / $
                  (17.72 + 6.60*model_p[ip] + 6.95*model_p[ip]^2 - 0.302*model_r23)
                model_logoh_upper = (model_r23 + 726.1 + 842.2*model_p[ip] + 337.5*model_p[ip]^2) / $
                  (85.96 + 82.76*model_p[ip] + 43.98*model_p[ip]^2 + 1.793*model_r23)
                model_good1 = where((model_logoh_upper gt model_logoh_lower))
                model_good2 = where((model_logoh_lower[model_good1] gt 7.0))
                djs_oplot, alog10(model_r23[model_good1]), model_logoh_upper[model_good1], $
                  linestyle=linestyle[ip]
                djs_oplot, alog10(model_r23[model_good1]), model_logoh_lower[model_good1], $
                  linestyle=linestyle[ip]
             endfor
             legend, 'P='+string(model_p,format='(F4.2)'), /left, /top, $
               box=0, charsize=1.3, line=linestyle, pspacing=1.4
          endif
; plot the data
          if (nup ne 0L) then oploterror, alog10(result[up].zstrong_r23), $
            result[up].zstrong_12oh, result[up].zstrong_r23_err/result[up].zstrong_r23/alog(10), $
            result[up].zstrong_12oh_err, psym=6, color=djs_icolor('orange'), $
            errcolor=djs_icolor('orange')
          if (nlo ne 0L) then oploterror, alog10(result[lo].zstrong_r23), $
            result[lo].zstrong_12oh, result[lo].zstrong_r23_err/result[lo].zstrong_r23/alog(10), $
            result[lo].zstrong_12oh_err, psym=6, color=djs_icolor('blue'), $
            errcolor=djs_icolor('blue')
          if (nambig ne 0L) then begin
             if (ngood ne 0L) then oploterror, alog10(result[ambig[good]].zstrong_r23), $
               result[ambig[good]].zstrong_12oh, result[ambig[good]].zstrong_r23_err/$
               result[ambig[good]].zstrong_r23/alog(10), $
               result[ambig[good]].zstrong_12oh_err, psym=symcat(16), $
               color=djs_icolor('red'), errcolor=djs_icolor('red')
             if (nreject ne 0L) then begin
                for jj = 0, nreject-1 do begin
                   oplot, alog10(result[ambig[reject[jj]]].zstrong_r23)*[1,1], $
                     [data[ambig[reject[jj]]].zstrong_12oh_upper,$
                     data[ambig[reject[jj]]].zstrong_12oh_lower], thick=3, $
                     psym=-6, color=djs_icolor('dark green');, symsize=1.4
                endfor
;               stop
             endif 
          endif
          legend, ['Upper','Lower','Ambig','Reject'], /right, /top, $
            box=0, psym=[6,6,symcat(16),6], $
            color=djs_icolor(['orange','blue','red','dark green']), $
            textcolor=djs_icolor(['orange','blue','red','dark green']), $
            charsize=1.4, thick=5
       endif 
    endif                       ; close JUSTFLUX condition
       
; #########################
; equivalent widths
    if keyword_set(justew) and tag_exist(data,'ZSTRONG_EW_12OH_LOWER') then begin
       if (keyword_set(silent) eq 0) then splog, 'Assigning R23 branches from EWs:'
       up = where((data.zstrong_ew_12oh_upper gt -900.0) and (r23branch eq 'U'),nup)
       if (nup ne 0L) then begin
          result[up].r23branch_ew        = 'U'
          result[up].zstrong_ew_12oh     = data[up].zstrong_ew_12oh_upper
          result[up].zstrong_ew_12oh_err = data[up].zstrong_ew_12oh_upper_err
          if keyword_set(kk04) then begin
             result[up].zstrong_ew_logu     = data[up].zstrong_ew_logu_upper
             result[up].zstrong_ew_logu_err = data[up].zstrong_ew_logu_upper_err
          endif
       endif
       
       lo = where((data.zstrong_ew_12oh_lower gt -900.0) and (r23branch eq 'L'),nlo)
       if (nlo ne 0L) then begin
          result[lo].r23branch_ew        = 'L'
          result[lo].zstrong_ew_12oh     = data[lo].zstrong_ew_12oh_lower
          result[lo].zstrong_ew_12oh_err = data[lo].zstrong_ew_12oh_lower_err
          if keyword_set(kk04) then begin
             result[lo].zstrong_ew_logu     = data[lo].zstrong_ew_logu_lower
             result[lo].zstrong_ew_logu_err = data[lo].zstrong_ew_logu_lower_err
          endif
       endif

; adopt the "average" lower/upper abundance of objects with ambiguous
; abundances, near the turn-around region; reject objects with a
; "frac" value less than fraccut (see IM_ABUNDANCE for details)
       ambig = where((data.zstrong_ew_12oh_upper gt -900) and (data.zstrong_ew_12oh_lower gt -900) and $
         ((data.zstrong_ew_12oh_upper lt data.zstrong_ew_12oh_lower)),nambig) ; off the R23 calibration
       if (keyword_set(silent) eq 0) then splog, '   Ambiguous: (O/H)_lower>(O/H)_upper: '+$
         string(nambig,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' ('+$
         strtrim(string(100.0*nambig/float(ngalaxy),format='(F12.1)'),2)+'%)'

; reject objects whose upper- and lower-branch abundances differ by
; more than 1-sigma, where 1-sigma is given by the width of the
; overlapping histogram (12oh_avg)
       if (nambig ne 0L) then begin
          good = where($
            (data[ambig].zstrong_ew_12oh_upper+data[ambig].zstrong_ew_12oh_upper_err) gt $
            (data[ambig].zstrong_ew_12oh_lower-data[ambig].zstrong_ew_12oh_lower_err),ngood,$
            comp=reject,ncomp=nreject)
          if (keyword_set(silent) eq 0) then splog, '     Retain: '+$
            string(ngood,format='(I0)')+'/'+string(nambig,format='(I0)')+' ('+$
            strtrim(string(100.0*ngood/float(nambig),format='(F12.1)'),2)+'%)'
          if (ngood ne 0L) then begin
             result[ambig[good]].r23branch_ew        = 'A'
             result[ambig[good]].zstrong_ew_12oh     = data[ambig[good]].zstrong_ew_12oh_avg
             result[ambig[good]].zstrong_ew_12oh_err = data[ambig[good]].zstrong_ew_12oh_avg_err
             if keyword_set(kk04) then begin
                result[ambig[good]].zstrong_ew_logu     = data[ambig[good]].zstrong_ew_logu_avg
                result[ambig[good]].zstrong_ew_logu_err = data[ambig[good]].zstrong_ew_logu_avg_err
             endif
          endif
          if (keyword_set(silent) eq 0) then splog, '     Reject: '+$
            string(nreject,format='(I0)')+'/'+string(nambig,format='(I0)')+' ('+$
            strtrim(string(100.0*nreject/float(nambig),format='(F12.1)'),2)+'%)'
          if (nreject ne 0L) then begin
             result[ambig[reject]].r23branch_ew        = 'Rejected'
             result[ambig[reject]].zstrong_ew_12oh     = -999.0
             result[ambig[reject]].zstrong_ew_12oh_err = -999.0
             result[ambig[reject]].zstrong_ew_logu     = -999.0
             result[ambig[reject]].zstrong_ew_logu_err = -999.0
          endif
       endif
    endif                       ; close JUSTEW condition
       
; rename the tags and return    
    result = im_struct_trimtags(result,$
      select=tag_names(result),newtags=finaltags)

return, result
end
