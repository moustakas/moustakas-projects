pro sings_log12oh_class, class, snrcut_class=snrcut_class, silent=silent
; jm06oct14uofa - assign spectral classifications; this should only
;   need to be run only after the spectra have been re-fitted  
; jm08feb11nyu - use a S/N>5 cut on [NII], otherwise you end up with
;   spurious classifications, especially with the nuclear spectra 
; jm10mar10ucsd - updated to final form

    version = sings_log12oh_version()
    outpath = sings_path(/projects)+'log12oh/'

    allnuclear = read_sings_gandalf(/nuclear)
    alldrift20 = read_sings_gandalf(/drift20)
    alldrift56 = read_sings_gandalf(/drift56)

    sings = sings_read_info()
    ngal = n_elements(sings)

    if (n_elements(snrcut_class) eq 0) then snrcut_class = 2.0
    
    class = {$
      galaxy: '', $
      ned_galaxy: '', $
      nice_galaxy: '', $
      nuclear_broad:         0, drift20_broad:         0, drift56_broad:         0, $
      nuclear_bpt_d:       0.0, drift20_bpt_d:       0.0, drift56_bpt_d:       0.0, $
      nuclear_bpt_phi:     0.0, drift20_bpt_phi:     0.0, drift56_bpt_phi:     0.0, $
      nuclear_n2flag:        0, drift20_n2flag:        0, drift56_n2flag:        0, $
      nuclear_class: '\nodata', drift20_class: '\nodata', drift56_class: '\nodata', $
      ho_class:      '\nodata', $
      class:         '\nodata', $
      class_remarks: '\nodata'}
    class = replicate(class,ngal)

    dale = rsex(outpath+'sings_dale_class.txt')

    class.galaxy = strtrim(sings.galaxy,2)
    class.ned_galaxy = strtrim(sings.ned_galaxy,2)
    class.nice_galaxy = strtrim(sings.nice_galaxy,2)
    class.class = dale.final_class
    class.class_remarks = repstr(dale.remarks,'__',' ')

; redo the classifications using a higher S/N cut than the default (3
; vs 1)

    if (not keyword_set(silent)) then splog, 'Nuclear:'
    nucflags = iclassification(allnuclear,ratios=nuc,$
      snrcut_class=snrcut_class,niihacut=-0.25)
    if (not keyword_set(silent)) then print

    if (not keyword_set(silent)) then splog, 'Drift20:'
    d20flags = iclassification(alldrift20,ratios=d20,$
      snrcut_class=snrcut_class,niihacut=-0.25)
    if (not keyword_set(silent)) then print

    if (not keyword_set(silent)) then splog, 'Drift56:'
    d56flags = iclassification(alldrift56,ratios=d56,$
      snrcut_class=snrcut_class,niihacut=-0.25)
    if (not keyword_set(silent)) then print

    for igal = 0, ngal-1 do begin

       match = where(strtrim(sings[igal].ned_galaxy,2) eq strtrim(allnuclear.ned_galaxy,2),nmatch)
       if (nmatch eq 0L) then begin
          class[igal].nuclear_class = '\nodata'
       endif else begin
          class[igal].nuclear_bpt_d   = nuc[match].bpt_d
          class[igal].nuclear_bpt_phi = nuc[match].bpt_phi
          class[igal].nuclear_broad = allnuclear[match].isbroad
          class[igal].nuclear_class = repstr(repstr(nuc[match].final_class,'UNKNOWN','?'),'NII_','')
          if strmatch(nuc[match].final_class,'NII_*',/fold) then class[igal].nuclear_n2flag = 1
;         if strmatch(class[igal].galaxy,'*ngc1512*',/fold) then stop
; BROAD-LINE AGN
          if strmatch(class[igal].galaxy,'*ngc1566*',/fold) then begin
             class[igal].nuclear_class = 'AGN'
             class[igal].class_remarks = 'Broad-line__AGN.'
          endif
          if strmatch(class[igal].galaxy,'*ngc3031*',/fold) then begin
             class[igal].nuclear_class = 'AGN'
             class[igal].class_remarks = 'Broad-line__AGN.'
          endif
          if strmatch(class[igal].galaxy,'*ngc4579*',/fold) then begin
             class[igal].nuclear_class = 'AGN'
             class[igal].class_remarks = 'Broad-line__AGN.'
          endif
          if strmatch(class[igal].galaxy,'*ngc4594*',/fold) then begin
             class[igal].nuclear_class = 'AGN'
             class[igal].class_remarks = 'Broad-line__AGN.'
          endif
          if strmatch(class[igal].galaxy,'*ngc5033*',/fold) then begin
             class[igal].nuclear_class = 'AGN'
             class[igal].class_remarks = 'Broad-line__AGN.'
          endif
; FIX SOME SPURIOUS CLASSIFICATIONS (jm07oct05nyu)
          if strmatch(class[igal].galaxy,'*ngc3034*',/fold) then begin
             class[igal].nuclear_class = 'SF'
             class[igal].class_remarks = 'Spurious__nuclear__classification.'
          endif
          if strmatch(class[igal].galaxy,'*ngc6946*',/fold) then begin
             class[igal].nuclear_class = 'SF'
             class[igal].class_remarks = 'Spurious__nuclear__classification.'
          endif
       endelse 

       match = where(strtrim(sings[igal].ned_galaxy,2) eq strtrim(alldrift20.ned_galaxy,2),nmatch)
       if (nmatch eq 0L) then begin
          class[igal].drift20_class = '\nodata'
       endif else begin
          class[igal].drift20_bpt_d   = d20[match].bpt_d
          class[igal].drift20_bpt_phi = d20[match].bpt_phi
          class[igal].drift20_broad = alldrift20[match].isbroad
          class[igal].drift20_class = repstr(repstr(d20[match].final_class,'UNKNOWN','?'),'NII_','')
          if strmatch(d20[match].final_class,'NII_*',/fold) then class[igal].drift20_n2flag = 1
; FIX SOME SPURIOUS CLASSIFICATIONS (jm07oct05nyu)
          if strmatch(class[igal].galaxy,'*ngc1377*',/fold) then begin
             class[igal].drift20_class = 'SF' ; AGN --> SF
             class[igal].class_remarks = 'Spurious__circumnuclear__classification.'
          endif
          if strmatch(class[igal].galaxy,'*ngc2403*',/fold) then begin
             class[igal].drift20_class = 'SF'
             class[igal].class_remarks = 'Spurious__circumnuclear__classification.'
          endif
          if strmatch(class[igal].galaxy,'*ngc3034*',/fold) then begin
             class[igal].drift20_class = 'SF'
             class[igal].class_remarks = 'Spurious__circumnuclear__classification.'
          endif
          if strmatch(class[igal].galaxy,'*ngc3198*',/fold) then begin
             class[igal].drift20_class = 'SF'
             class[igal].class_remarks = 'Spurious__circumnuclear__classification.'
          endif
          if strmatch(class[igal].galaxy,'*ngc6946*',/fold) then begin
             class[igal].drift20_class = 'SF'
             class[igal].class_remarks = 'Spurious__circumnuclear__classification.'
          endif
       endelse

; #drift56
       match = where(strtrim(sings[igal].ned_galaxy,2) eq strtrim(alldrift56.ned_galaxy,2),nmatch)
       if (nmatch eq 0L) then begin
          class[igal].drift56_class = '\nodata'
       endif else begin
          class[igal].drift56_bpt_d   = d56[match].bpt_d
          class[igal].drift56_bpt_phi = d56[match].bpt_phi
          class[igal].drift56_broad = alldrift56[match].isbroad
          class[igal].drift56_class = repstr(repstr(d56[match].final_class,'UNKNOWN','?'),'NII_','')
          if strmatch(d56[match].final_class,'NII_*',/fold) then class[igal].drift56_n2flag = 1
; BROAD-LINE AGN
          if strmatch(class[igal].galaxy,'*ngc1566*',/fold) then begin
             class[igal].drift56_class = 'AGN'
             class[igal].class_remarks = 'Broad-line__AGN.'
          endif
; FIX SOME SPURIOUS CLASSIFICATIONS (jm07aug10nyu)
          if strmatch(class[igal].galaxy,'*ngc1097*',/fold) then begin
             class[igal].drift56_class = 'SF' ; SF/AGN --> SF
             class[igal].class_remarks = 'Spurious__radial__strip__classification.'
          endif
          if strmatch(class[igal].galaxy,'*ngc1482*',/fold) then begin
             class[igal].drift56_class = 'SF' ; SF/AGN --> SF
             class[igal].class_remarks = 'Spurious__radial__strip__classification.'
          endif
          if strmatch(class[igal].galaxy,'*ngc1512*',/fold) then begin
             class[igal].drift56_class = 'SF' ; SF/AGN --> SF
             class[igal].class_remarks = 'Spurious__radial__strip__classification.'
          endif
       endelse

    endfor

; Ho et al. classifications; do not apply any S/N cut on the emission
; lines 
    ho = read_97ho()
    hoflags = iclassification(ho,snrcut_class=0.0,/silent,ratios=ho_class)

    match, strtrim(class.ned_galaxy,2), strtrim(ho.ned_galaxy,2), class_match, ho_match
    nho_match = n_elements(ho_match)
;   splog, 'Identified '+string(nho_match,format='(I0)')+'/'+string(ngal,format='(I0)')+$
;     ' SINGS galaxies in the Ho et al. 1997 atlas.'
;   niceprint, class[class_match].galaxy, ho[ho_match].galaxy, ho[ho_match].class
    
;   class[class_match].ho_class = repstr(ho_class[ho_match].bpt_nii_mixture_class,'Unknown','?')
    class[class_match].ho_class = repstr(repstr(ho_class[ho_match].final_class,'UNKNOWN','?'),'NII_')

    hoclass = where((strtrim(class[class_match].ho_class,2) ne '\nodata') and $
      (strtrim(class[class_match].ho_class,2) ne '?'),nhoclass,comp=ho_noclass,$
      ncomp=nho_noclass)
    
    splog, '------------------------------'
    splog, 'Ho et al. 1997:'
    splog, '   Matching  : '+string(nho_match,format='(I0)')+'/'+string(ngal,format='(I0)')+$
      ' = '+string(100.0*nho_match/ngal,format='(F4.1)')+'%'
    splog, '   Classified: '+string(nhoclass,format='(I0)')+'/'+string(nho_match,format='(I0)')+$
      ' = '+string(100.0*nhoclass/nho_match,format='(F5.1)')+'%'

; print out statistics on classified galaxies

;   nuclear = where((strtrim(class.nuclear_class,2) ne '\nodata') and $
;     (strtrim(class.nuclear_class,2) ne '?'),nnuclear,comp=nuclear_noclass,$
;     ncomp=nnuclear_noclass)
    nuclear = where((strtrim(class.nuclear_class,2) ne '\nodata'),nnuclear)
    nuclear_noclass = where((strtrim(class[nuclear].nuclear_class,2) eq '?'),nnuclear_noclass)
    nuclear_sf = where((strtrim(class[nuclear].nuclear_class,2) eq 'SF'),nnuclear_sf)
    nuclear_sfagn = where((strtrim(class[nuclear].nuclear_class,2) eq 'SF/AGN'),nnuclear_sfagn)
    nuclear_agn = where((strtrim(class[nuclear].nuclear_class,2) eq 'AGN'),nnuclear_agn)

;   drift20 = where((strtrim(class.drift20_class,2) ne '\nodata') and $
;     (strtrim(class.drift20_class,2) ne '?'),ndrift20,comp=drift20_noclass,$
;     ncomp=ndrift20_noclass)
    drift20 = where((strtrim(class.drift20_class,2) ne '\nodata'),ndrift20)
    drift20_noclass = where((strtrim(class[drift20].drift20_class,2) eq '?'),ndrift20_noclass)
    drift20_sf = where((strtrim(class[drift20].drift20_class,2) eq 'SF'),ndrift20_sf)
    drift20_sfagn = where((strtrim(class[drift20].drift20_class,2) eq 'SF/AGN'),ndrift20_sfagn)
    drift20_agn = where((strtrim(class[drift20].drift20_class,2) eq 'AGN'),ndrift20_agn)

;   drift56 = where((strtrim(class.drift56_class,2) ne '\nodata') and $
;     (strtrim(class.drift56_class,2) ne '?'),ndrift56,comp=drift56_noclass,$
;     ncomp=ndrift56_noclass)
    drift56 = where((strtrim(class.drift56_class,2) ne '\nodata'),ndrift56)
    drift56_noclass = where((strtrim(class[drift56].drift56_class,2) eq '?'),ndrift56_noclass)
    drift56_sf = where((strtrim(class[drift56].drift56_class,2) eq 'SF'),ndrift56_sf)
    drift56_sfagn = where((strtrim(class[drift56].drift56_class,2) eq 'SF/AGN'),ndrift56_sfagn)
    drift56_agn = where((strtrim(class[drift56].drift56_class,2) eq 'AGN'),ndrift56_agn)

    splog, '------------------------------'
    splog, 'Nuclear:'
    splog, '   SF    : '+string(nnuclear_sf,format='(I0)')+'/'+string(nnuclear,format='(I0)')+$
      ' = '+string(100.0*nnuclear_sf/nnuclear,format='(F4.1)')+'%'
    splog, '   SF/AGN:  '+string(nnuclear_sfagn,format='(I0)')+'/'+string(nnuclear,format='(I0)')+$
      ' = '+string(100.0*nnuclear_sfagn/nnuclear,format='(F4.1)')+'%'
    splog, '   AGN   : '+string(nnuclear_agn,format='(I0)')+'/'+string(nnuclear,format='(I0)')+$
      ' = '+string(100.0*nnuclear_agn/nnuclear,format='(F4.1)')+'%'
    splog, '   ?     :  '+string(nnuclear_noclass,format='(I0)')+'/'+string(nnuclear,format='(I0)')+$
      ' = '+string(100.0*nnuclear_noclass/nnuclear,format='(F4.1)')+'%'

    splog, '------------------------------'
    splog, 'Circumnuclear:'
    splog, '   SF    : '+string(ndrift20_sf,format='(I0)')+'/'+string(ndrift20,format='(I0)')+$
      ' = '+string(100.0*ndrift20_sf/ndrift20,format='(F4.1)')+'%'
    splog, '   SF/AGN:  '+string(ndrift20_sfagn,format='(I0)')+'/'+string(ndrift20,format='(I0)')+$
      ' = '+string(100.0*ndrift20_sfagn/ndrift20,format='(F4.1)')+'%'
    splog, '   AGN   : '+string(ndrift20_agn,format='(I0)')+'/'+string(ndrift20,format='(I0)')+$
      ' = '+string(100.0*ndrift20_agn/ndrift20,format='(F4.1)')+'%'
    splog, '   ?     :  '+string(ndrift20_noclass,format='(I0)')+'/'+string(ndrift20,format='(I0)')+$
      ' = '+string(100.0*ndrift20_noclass/ndrift20,format='(F4.1)')+'%'

    splog, '------------------------------'
    splog, 'Radial Strip:'
    splog, '   SF    : '+string(ndrift56_sf,format='(I0)')+'/'+string(ndrift56,format='(I0)')+$
      ' = '+string(100.0*ndrift56_sf/ndrift56,format='(F4.1)')+'%'
    splog, '   SF/AGN: '+string(ndrift56_sfagn,format='(I0)')+'/'+string(ndrift56,format='(I0)')+$
      ' = '+string(100.0*ndrift56_sfagn/ndrift56,format='(F4.1)')+'%'
    splog, '   AGN   :  '+string(ndrift56_agn,format='(I0)')+'/'+string(ndrift56,format='(I0)')+$
      ' = '+string(100.0*ndrift56_agn/ndrift56,format='(F4.1)')+'%'
    splog, '   ?     : '+string(ndrift56_noclass,format='(I0)')+'/'+string(ndrift56,format='(I0)')+$
      ' = '+string(100.0*ndrift56_noclass/ndrift56,format='(F4.1)')+'%'

; write out
    splog, 'Writing '+outpath+'sings_class_'+version+'.fits'
    struct_print, class, file=outpath+'sings_class_'+version+'.dat'
    mwrfits, class, outpath+'sings_class_'+version+'.fits', /create

return
end
    
