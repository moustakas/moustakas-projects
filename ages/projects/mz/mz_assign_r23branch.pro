;+
; NAME:
;   MZ_ASSIGN_R23BRANCH()
;
; PURPOSE:
;   Given the output from IM_ABUNDANCE() assign the appropriate
;   R23 branch. 
;
; INPUTS:
;   abund - output structure from IM_ABUNDANCE() or MZ_ABUNDANCE() 
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;   kk04 - consider the Kobulnicky & Kewley (2004) abundances
;          (default) 
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

function mz_assign_r23branch, abund, r23branch=r23branch, justew=justew, $
  justflux=justflux, kk04=kk04, m91=m91, silent=silent, debug=debug, _extra=extra

    nobj = n_elements(abund)
    if (nobj eq 0L) then begin
       doc_library, 'mz_assign_r23branch'
       return, -1
    endif
    
    if (keyword_set(kk04) eq 0) and (keyword_set(m91) eq 0) then begin
       splog, 'Must specify one of KK04 or M91'
       return, -1
    endif
    if keyword_set(kk04) then suffix = 'KK04'
    if keyword_set(m91) then suffix = 'M91'

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
      zstrong_niiha_limit: -999.0,$
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
    ngal = n_elements(result)

; assign the appropriate suffix to the output tags    
    finaltags = tag_names(result[0])
    alter = where(strmatch(finaltags,'*12oh*',/fold) or $
      strmatch(finaltags,'*logu*',/fold))
;     strmatch(finaltags,'*r23branch*',/fold))
    finaltags[alter] = repstr(finaltags[alter]+'_'+$
      suffix,'_ERR_'+suffix,'_'+suffix+'_ERR')

    data = struct_trimtags(abund,select=['ZSTRONG_12OH_'+$
      suffix+'_*','ZSTRONG_EW_12OH_'+suffix+'_*',$
      'ZSTRONG_LOGU_'+suffix+'_*','ZSTRONG_EW_LOGU_'+suffix+'_*'])
    data = im_struct_trimtags(data,select=tag_names(data),$
      newtags=repstr(tag_names(data),suffix+'_',''))

; special case
    if keyword_set(kk04) then begin
       converge_upper = abund.zstrong_converge_kk04_upper
       converge_lower = abund.zstrong_converge_kk04_lower
       converge_ew_upper = abund.zstrong_ew_converge_kk04_upper
       converge_ew_lower = abund.zstrong_ew_converge_kk04_lower
    endif else begin
       converge_upper = intarr(ngal)+1
       converge_lower = intarr(ngal)+1
       converge_ew_upper = intarr(ngal)+1
       converge_ew_lower = intarr(ngal)+1
    endelse
    
;; get R23-branches based on [NII]/Ha, assuming a default upper branch
;    r23branch = replicate('U',nobj) ; default
;    lodetect = (result.zstrong_niiha gt -900) and (result.zstrong_niiha lt -1.1)
;    lolimit = (result.zstrong_niiha_limit gt -900) and (result.zstrong_niiha_limit lt -1.1)
;    lo = where(lodetect or lolimit,nlo)
;    if (nlo ne 0L) then r23branch[lo] = 'L'

; #########################
; fluxes    
    if keyword_set(justflux) then begin
       if (keyword_set(silent) eq 0) then splog, 'Assigning R23 branches from fluxes:'

       up = where((data.zstrong_12oh_upper gt -900.0) and (r23branch eq 'U') and $
         (converge_upper eq 1),nup)
       if (nup ne 0L) then begin
          result[up].r23branch        = 'U'
          result[up].zstrong_12oh     = data[up].zstrong_12oh_upper
          result[up].zstrong_12oh_err = data[up].zstrong_12oh_upper_err
          if keyword_set(kk04) then begin
             result[up].zstrong_logu     = data[up].zstrong_logu_upper
             result[up].zstrong_logu_err = data[up].zstrong_logu_upper_err
          endif
       endif
       
       lo = where((data.zstrong_12oh_lower gt -900.0) and (r23branch eq 'L') and $
         (converge_lower eq 1),nlo)
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
         string(nambig,format='(I0)')+'/'+string(ngal,format='(I0)')+' ('+$
         strtrim(string(100.0*nambig/float(ngal),format='(F12.1)'),2)+'%)'

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
            yrange=[7,9.5], xrange=[-0.3,1.3], _extra=extra
          if keyword_set(kk04) then mzoplot_kk04_models
          if keyword_set(m91) then mzoplot_m91_models
          if (nobj gt 1000) then psym = 3 else psym = 6
; plot the data
          if (nup ne 0L) then oploterror, alog10(result[up].zstrong_r23), $
            result[up].zstrong_12oh, result[up].zstrong_r23_err/result[up].zstrong_r23/alog(10), $
            result[up].zstrong_12oh_err, psym=psym, color=djs_icolor('orange'), $
            errcolor=djs_icolor('orange')
          if (nlo ne 0L) then oploterror, alog10(result[lo].zstrong_r23), $
            result[lo].zstrong_12oh, result[lo].zstrong_r23_err/result[lo].zstrong_r23/alog(10), $
            result[lo].zstrong_12oh_err, psym=psym, color=djs_icolor('blue'), $
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
            box=0, psym=[psym,psym,symcat(16),6], $
            color=djs_icolor(['orange','blue','red','dark green']), $
            textcolor=djs_icolor(['orange','blue','red','dark green']), $
            charsize=1.4, thick=5
          cc = get_kbrd(1)
       endif 
    endif                       ; close JUSTFLUX condition
       
; #########################
; equivalent widths
    if keyword_set(justew) and tag_exist(data,'ZSTRONG_EW_12OH_LOWER') then begin
       if (keyword_set(silent) eq 0) then splog, 'Assigning R23 branches from EWs:'
       up = where((data.zstrong_ew_12oh_upper gt -900.0) and (r23branch eq 'U') and $
         (converge_ew_upper eq 1),nup)
       if (nup ne 0L) then begin
          result[up].r23branch_ew        = 'U'
          result[up].zstrong_ew_12oh     = data[up].zstrong_ew_12oh_upper
          result[up].zstrong_ew_12oh_err = data[up].zstrong_ew_12oh_upper_err
          if keyword_set(kk04) then begin
             result[up].zstrong_ew_logu     = data[up].zstrong_ew_logu_upper
             result[up].zstrong_ew_logu_err = data[up].zstrong_ew_logu_upper_err
          endif
       endif
       
       lo = where((data.zstrong_ew_12oh_lower gt -900.0) and (r23branch eq 'L') and $
         (converge_ew_lower eq 1),nlo)
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
         string(nambig,format='(I0)')+'/'+string(ngal,format='(I0)')+' ('+$
         strtrim(string(100.0*nambig/float(ngal),format='(F12.1)'),2)+'%)'

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

; build the QAplot, if requested
       if keyword_set(debug) then begin
          djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
            yrange=[7,9.5], xrange=[-0.5,1.3], _extra=extra
          if keyword_set(kk04) then mzoplot_kk04_models
          if keyword_set(m91) then mzoplot_m91_models
          if (nobj gt 1000) then psym = 3 else psym = 6
; plot the data
          if (nup ne 0L) then oploterror, alog10(result[up].zstrong_ew_r23), $
            result[up].zstrong_ew_12oh, result[up].zstrong_ew_r23_err/result[up].zstrong_ew_r23/alog(10), $
            result[up].zstrong_ew_12oh_err, psym=psym, color=djs_icolor('orange'), $
            errcolor=djs_icolor('orange')
          if (nlo ne 0L) then oploterror, alog10(result[lo].zstrong_ew_r23), $
            result[lo].zstrong_ew_12oh, result[lo].zstrong_ew_r23_err/result[lo].zstrong_ew_r23/alog(10), $
            result[lo].zstrong_ew_12oh_err, psym=psym, color=djs_icolor('blue'), $
            errcolor=djs_icolor('blue')
          if (nambig ne 0L) then begin
             if (ngood ne 0L) then oploterror, alog10(result[ambig[good]].zstrong_ew_r23), $
               result[ambig[good]].zstrong_ew_12oh, result[ambig[good]].zstrong_ew_r23_err/$
               result[ambig[good]].zstrong_ew_r23/alog(10), $
               result[ambig[good]].zstrong_ew_12oh_err, psym=symcat(16), $
               color=djs_icolor('red'), errcolor=djs_icolor('red')
             if (nreject ne 0L) then begin
                for jj = 0, nreject-1 do begin
                   oplot, alog10(result[ambig[reject[jj]]].zstrong_ew_r23)*[1,1], $
                     [data[ambig[reject[jj]]].zstrong_ew_12oh_upper,$
                     data[ambig[reject[jj]]].zstrong_ew_12oh_lower], thick=3, $
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
          cc = get_kbrd(1)
       endif 
    endif                       ; close JUSTEW condition
       
; rename the tags and return    
    result = im_struct_trimtags(result,$
      select=tag_names(result),newtags=finaltags)

return, result
end
