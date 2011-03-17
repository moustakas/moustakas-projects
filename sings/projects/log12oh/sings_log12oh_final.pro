pro sings_log12oh_final, result, clobber=clobber
; jm10mar11ucsd - build the final list of oxygen abundances for the
; SINGS galaxies

; read the data    
    version = sings_log12oh_version()
    outpath = sings_path(/projects)+'log12oh/'

    outfile = outpath+'sings_log12oh_final_'+version+'.fits'
    if file_test(outfile+'.gz') and (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+outfile+'.gz exists; use /CLOBBER'
       return
    endif

    galinfo = mrdfits(outpath+'sings_log12oh_'+version+'.fits.gz',1)
    hiiinfo = mrdfits(outpath+'sings_log12oh_hiiregions_'+version+'.fits.gz',1)
    ngal = n_elements(galinfo)

    rr25_char = 0.4 ; definition of "characteristic" metallicity
    rr25_avg = 0.5

; there was some discussion with Christy whether I should use the
; central metallicity, which is defined at R/R25=0.1, or the
; metallicity at R=0; the latter is ~0.04 dex higher than the former;
; examining the O/H vs R/R25 "final" plot (Fig 8 in the paper) it
; looks like the metallicity at R=0 might be more appropriate given
; the average R/R25 values of our nuclear and circumnuclear spectra 
    usenuclear = 1
    if usenuclear then begin
       rr25_cent = 0.0          ; definition of "central" metallicity 
       rr25_cut = 0.1
    endif else begin
; old values!
       rr25_cent = 0.1          ; definition of "central" metallicity 
       rr25_cut = 0.2
    endelse
    
    result = {$
      sings_id:                                 0, $
      galaxy:                                  '', $
      ned_galaxy:                              '', $
      nice_galaxy:                             '', $
      ra:                                      '', $
      dec:                                     '', $
      type:                                    '', $
      t:                                      0.0, $
      mb:                                  -999.0, $
      mb_err:                              -999.0, $
      r23_branch:                             '?', $
      class:                            '\nodata', $ ; final spectral class
      class_remarks:                    '\nodata', $
      
      log12oh_pt05_char:          [-999.0,-999.0], $
      log12oh_pt05_char_comment:            '...', $
      log12oh_pt05_central:       [-999.0,-999.0], $
      log12oh_pt05_central_comment:         '...', $
      log12oh_pt05_lz:            [-999.0,-999.0], $
      
      log12oh_kk04_char:          [-999.0,-999.0], $
      log12oh_kk04_char_comment:            '...', $
      log12oh_kk04_central:       [-999.0,-999.0], $
      log12oh_kk04_central_comment:         '...', $
      log12oh_kk04_lz:            [-999.0,-999.0]}
    result = replicate(result,ngal)
    struct_assign, galinfo, result, /nozero
   
; KK04; collect all six possible metallicities for each object
    splog, '##################################################'
    splog, 'KK04'

    final1 = {galaxy: '', oh: 0.0, oh_err: 0.0, $
      rr25: 0.0, comment: '', class: ''}

    nuc = where((strtrim(galinfo.nuclear_class,2) ne 'AGN') and $
      (galinfo.nuclear_log12oh_kk04[0] gt -900.0),nnuc)
    nucfinal = replicate(final1,nnuc)
    nucfinal.galaxy = galinfo[nuc].galaxy
    nucfinal.oh = galinfo[nuc].nuclear_log12oh_kk04[0]
    nucfinal.oh_err = galinfo[nuc].nuclear_log12oh_kk04[1]
    nucfinal.rr25 = galinfo[nuc].nuclear_rr25
    nucfinal.class = galinfo[nuc].nuclear_class
    nucfinal.comment = 'Nuclear'

    circum = where((strtrim(galinfo.drift20_class,2) ne 'AGN') and $
      (galinfo.drift20_log12oh_kk04[0] gt -900.0),ncircum)
    circumfinal = replicate(final1,ncircum)
    circumfinal.galaxy = galinfo[circum].galaxy
    circumfinal.oh = galinfo[circum].drift20_log12oh_kk04[0]
    circumfinal.oh_err = galinfo[circum].drift20_log12oh_kk04[1]
    circumfinal.rr25 = galinfo[circum].drift20_rr25
    circumfinal.class = galinfo[circum].drift20_class
    circumfinal.comment = 'Circum'

    strip = where((strtrim(galinfo.drift56_class,2) ne 'AGN') and $
      (galinfo.drift56_log12oh_kk04[0] gt -900.0),nstrip)
    stripfinal = replicate(final1,nstrip)
    stripfinal.galaxy = galinfo[strip].galaxy
    stripfinal.oh = galinfo[strip].drift56_log12oh_kk04[0]
    stripfinal.oh_err = galinfo[strip].drift56_log12oh_kk04[1]
    stripfinal.rr25 = galinfo[strip].drift56_rr25
    stripfinal.class = galinfo[strip].drift56_class
    stripfinal.comment = 'Strip'

    char = where((hiiinfo.hii_kk04_log12oh_char[0] gt -900.0),nchar)
    charfinal = replicate(final1,nchar)
    charfinal.galaxy = galinfo[char].galaxy
    charfinal.oh = hiiinfo[char].hii_kk04_log12oh_char[0]
    charfinal.oh_err = hiiinfo[char].hii_kk04_log12oh_char[1]
    charfinal.rr25 = replicate(rr25_char,nchar)
    charfinal.comment = 'Char' ; string(rr25_char,format='(G0)')+'*R25'

    if usenuclear then begin
; central (R=0) metallicities
       cent = where((hiiinfo.hii_kk04_log12oh_nuclear[0] gt -900.0),ncent)
       centfinal = replicate(final1,ncent)
       centfinal.galaxy = galinfo[cent].galaxy
       centfinal.oh = hiiinfo[cent].hii_kk04_log12oh_nuclear[0]
       centfinal.oh_err = hiiinfo[cent].hii_kk04_log12oh_nuclear[1]
       centfinal.rr25 = replicate(rr25_cent,ncent)
       centfinal.comment = 'Cent' ; string(rr25_cent,format='(G0)')+'*R25'
    endif else begin
; central (R/R25=0.1) metallicities!
       cent = where((hiiinfo.hii_kk04_log12oh_central[0] gt -900.0),ncent)
       centfinal = replicate(final1,ncent)
       centfinal.galaxy = galinfo[cent].galaxy
       centfinal.oh = hiiinfo[cent].hii_kk04_log12oh_central[0]
       centfinal.oh_err = hiiinfo[cent].hii_kk04_log12oh_central[1]
       centfinal.rr25 = replicate(rr25_cent,ncent)
       centfinal.comment = 'Cent' ; string(rr25_cent,format='(G0)')+'*R25'
    endelse

; only use the average metallicity if the characteristic metallicity
; doesn't exist (so we're not double-counting HII regions)
    avg = where(((hiiinfo.hii_kk04_log12oh_avg[0] gt -900.0) and $
      (hiiinfo.hii_kk04_log12oh_char[0] lt -900.0)),navg)
    avgfinal = replicate(final1,navg)
    avgfinal.galaxy = galinfo[avg].galaxy
    avgfinal.oh = hiiinfo[avg].hii_kk04_log12oh_avg[0]
    avgfinal.oh_err = hiiinfo[avg].hii_kk04_log12oh_avg[1]
    avgfinal.rr25 = replicate(rr25_avg,navg)
;   avgfinal.rr25 = hiiinfo[avg].hii_kk04_rr25_avg[0]
    avgfinal.comment = 'HIIAvg'

;; electron temperature abundance: just keep the dwarf galaxies
;    ohte = where((hiiinfo.hii_te_log12oh_avg[0] gt -900.0) and $
;      (strtrim(hiiinfo.galaxy,2) ne 'NGC0628') and $
;      (strtrim(hiiinfo.galaxy,2) ne 'NGC0925') and $
;      (strtrim(hiiinfo.galaxy,2) ne 'NGC5194') and $
;      (strtrim(hiiinfo.galaxy,2) ne 'NGC6946') and $
;      (strtrim(hiiinfo.galaxy,2) ne 'NGC7793'),nohte)
;    ohtefinal = replicate(final1,nohte)
;    ohtefinal.galaxy = galinfo[ohte].galaxy
;    ohtefinal.oh = hiiinfo[ohte].hii_te_log12oh_avg[0]
;    ohtefinal.oh_err = hiiinfo[ohte].hii_te_log12oh_avg[1]
;    ohtefinal.comment = 'Te'

;   gg = [nucfinal,circumfinal,stripfinal]
;   hh = [avgfinal,charfinal,centfinal]
;   splog, mean(hh.oh_err), mean(gg.oh_err)
    
; now put everything together and average in two bins of R/R25
    final_kk04 = [nucfinal,circumfinal,stripfinal,charfinal,centfinal,avgfinal]
    final_kk04.galaxy = strtrim(final_kk04.galaxy,2)
    final_kk04 = final_kk04[sort(final_kk04.galaxy)]
;   struct_print, final_kk04

    openw, lun, repstr(outfile,'.fits','_kk04.txt'), /get_lun
    for igal = 0, ngal-1 do begin
       printf, lun, '--------------------------------------------------'
       printf, lun, '### '+result[igal].galaxy
       these = where(strtrim(result[igal].galaxy,2) eq $
         strtrim(final_kk04.galaxy,2),nthese)
       if (nthese ne 0) then begin
;         struct_print, final_kk04[these]
; identify metallicities corresponding to the inner and outer regions,
; which will correspond to our "central" and "characteristic"
; metallicities; if only either "inner" or "outer" metallicities
; exist, then, in the absence of any other information, make the 
; central and characteristic metallicities equal
          inner = where(final_kk04[these].rr25 lt rr25_cut,ninner,$
            comp=outer,ncomp=nouter)
          if (ninner ne 0) then begin
             printf, lun, ' '
             for ii = 0, ninner-1 do printf, lun, 'Inner: ', final_kk04[these[inner[ii]]].oh, $
               final_kk04[these[inner[ii]]].oh_err, final_kk04[these[inner[ii]]].class, $
               final_kk04[these[inner[ii]]].comment, format='(A7,F7.5,2x,F7.5,2x,A9,A0)'
             oh = final_kk04[these[inner]].oh
             oh_err = final_kk04[these[inner]].oh_err
             wavg = sings_weighted_mean(oh,oh_err,wmean_err=werr)
             result[igal].log12oh_kk04_central = [wavg,werr>0.01]
             result[igal].log12oh_kk04_central_comment = strjoin(final_kk04[these[inner]].comment,'/')
             printf, lun, '  Avg: ', result[igal].log12oh_kk04_central[0], result[igal].log12oh_kk04_central[1], $
               result[igal].log12oh_kk04_central_comment, format='(A7,F7.5,2x,F7.5,2x,A0)'
             if (nouter eq 0) then begin
                result[igal].log12oh_kk04_char = result[igal].log12oh_kk04_central
                result[igal].log12oh_kk04_char_comment = 'Central'
             endif
;            if strmatch(result[igal].galaxy,'*0925*') then stop
          endif
          if (nouter ne 0) then begin
             printf, lun, ' '
             for ii = 0, nouter-1 do printf, lun, 'Outer: ', final_kk04[these[outer[ii]]].oh, $
               final_kk04[these[outer[ii]]].oh_err, final_kk04[these[outer[ii]]].class, $
               final_kk04[these[outer[ii]]].comment, format='(A7,F7.5,2x,F7.5,2x,A9,A0)'
             oh = final_kk04[these[outer]].oh
             oh_err = final_kk04[these[outer]].oh_err
             wavg = sings_weighted_mean(oh,oh_err,wmean_err=werr)
             result[igal].log12oh_kk04_char = [wavg,werr>0.01]
             result[igal].log12oh_kk04_char_comment = strjoin(final_kk04[these[outer]].comment,'/')
             printf, lun, '  Avg: ', result[igal].log12oh_kk04_char[0], result[igal].log12oh_kk04_char[1], $
               result[igal].log12oh_kk04_char_comment, format='(A7,F7.5,2x,F7.5,2x,A0)'
             if (ninner eq 0) then begin
                result[igal].log12oh_kk04_central = result[igal].log12oh_kk04_char
                result[igal].log12oh_kk04_central_comment = 'Char'
             endif
          endif 
; if (O/H)_central<(O/H)_char then take the weighted average of the
; two values
          if (result[igal].log12oh_kk04_char[0] gt result[igal].log12oh_kk04_central[0]) then begin
             printf, lun, ' '
             printf, lun, ' ** (O/H)_central<(O/H)_char **'
             oh = [result[igal].log12oh_kk04_char[0],result[igal].log12oh_kk04_central[0]]
             oh_err = [result[igal].log12oh_kk04_char[1],result[igal].log12oh_kk04_central[1]]
             wavg = sings_weighted_mean(oh,oh_err,wmean_err=werr)

             result[igal].log12oh_kk04_char = [wavg,werr>0.01]
             result[igal].log12oh_kk04_central = result[igal].log12oh_kk04_char
             comment = 'cent>char:'+result[igal].log12oh_kk04_central_comment+$
               '/'+result[igal].log12oh_kk04_char_comment
             result[igal].log12oh_kk04_char_comment = comment
             result[igal].log12oh_kk04_central_comment = comment
             printf, lun, 'Final: ', result[igal].log12oh_kk04_char[0], result[igal].log12oh_kk04_char[1], $
               result[igal].log12oh_kk04_char_comment, format='(A7,F7.5,2x,F7.5,2x,A0)'
          endif
; average the abundances of certain special cases
          if $
            strmatch(result[igal].galaxy,'*ngc0337*',/fold) or $
            strmatch(result[igal].galaxy,'*ngc1512*',/fold) or $
            strmatch(result[igal].galaxy,'*ic2574*',/fold) or $
            strmatch(result[igal].galaxy,'*mrk0033*',/fold) or $
            strmatch(result[igal].galaxy,'*ngc4536*',/fold) or $
            strmatch(result[igal].galaxy,'*ngc4631*',/fold) or $
            strmatch(result[igal].galaxy,'*ngc0024*',/fold) then begin
             printf, lun, ' '
             printf, lun, ' ** Special case **'
             oh = [result[igal].log12oh_kk04_char[0],result[igal].log12oh_kk04_central[0]]
             oh_err = [result[igal].log12oh_kk04_char[1],result[igal].log12oh_kk04_central[1]]
             wavg = sings_weighted_mean(oh,oh_err,wmean_err=werr)
             result[igal].log12oh_kk04_char = [wavg,werr>0.01]
             result[igal].log12oh_kk04_central = result[igal].log12oh_kk04_char
             print, result[igal].galaxy
;            if strmatch(result[igal].galaxy,'*ngc2915*',/fold) then stop
          
             comment = 'specialcase:'+result[igal].log12oh_kk04_central_comment+$
               '/'+result[igal].log12oh_kk04_char_comment
             result[igal].log12oh_kk04_char_comment = comment
             result[igal].log12oh_kk04_central_comment = comment
             printf, lun, 'Final: ', result[igal].log12oh_kk04_char[0], result[igal].log12oh_kk04_char[1], $
               result[igal].log12oh_kk04_char_comment, format='(A7,F7.5,2x,F7.5,2x,A0)'
          endif
       endif else printf, lun, '  No abundances!'
;      if strmatch(sings[igal].galaxy,'*2915*') then stop
       printf, lun, ' '
    endfor
    free_lun, lun

; PT05; collect all six possible metallicities for each object
    splog, '##################################################'
    splog, 'PT05'

    final1 = {galaxy: '', oh: 0.0, oh_err: 0.0, $
      rr25: 0.0, comment: '', class: ''}

    nuc = where((strtrim(galinfo.nuclear_class,2) ne 'AGN') and $
      (galinfo.nuclear_log12oh_pt05[0] gt -900.0),nnuc)
    nucfinal = replicate(final1,nnuc)
    nucfinal.galaxy = galinfo[nuc].galaxy
    nucfinal.oh = galinfo[nuc].nuclear_log12oh_pt05[0]
    nucfinal.oh_err = galinfo[nuc].nuclear_log12oh_pt05[1]
    nucfinal.rr25 = galinfo[nuc].nuclear_rr25
    nucfinal.class = galinfo[nuc].nuclear_class
    nucfinal.comment = 'Nuclear'

    circum = where((strtrim(galinfo.drift20_class,2) ne 'AGN') and $
      (galinfo.drift20_log12oh_pt05[0] gt -900.0),ncircum)
    circumfinal = replicate(final1,ncircum)
    circumfinal.galaxy = galinfo[circum].galaxy
    circumfinal.oh = galinfo[circum].drift20_log12oh_pt05[0]
    circumfinal.oh_err = galinfo[circum].drift20_log12oh_pt05[1]
    circumfinal.rr25 = galinfo[circum].drift20_rr25
    circumfinal.class = galinfo[circum].drift20_class
    circumfinal.comment = 'Circum'

    strip = where((strtrim(galinfo.drift56_class,2) ne 'AGN') and $
      (galinfo.drift56_log12oh_pt05[0] gt -900.0),nstrip)
    stripfinal = replicate(final1,nstrip)
    stripfinal.galaxy = galinfo[strip].galaxy
    stripfinal.oh = galinfo[strip].drift56_log12oh_pt05[0]
    stripfinal.oh_err = galinfo[strip].drift56_log12oh_pt05[1]
    stripfinal.rr25 = galinfo[strip].drift56_rr25
    stripfinal.class = galinfo[strip].drift56_class
    stripfinal.comment = 'Strip'

    char = where((hiiinfo.hii_pt05_log12oh_char[0] gt -900.0),nchar)
    charfinal = replicate(final1,nchar)
    charfinal.galaxy = galinfo[char].galaxy
    charfinal.oh = hiiinfo[char].hii_pt05_log12oh_char[0]
    charfinal.oh_err = hiiinfo[char].hii_pt05_log12oh_char[1]
    charfinal.rr25 = replicate(rr25_char,nchar)
    charfinal.comment = 'Char' ; string(rr25_char,format='(G0)')+'*R25'

    if usenuclear then begin
; central (R=0) metallicities
       cent = where((hiiinfo.hii_pt05_log12oh_nuclear[0] gt -900.0),ncent)
       centfinal = replicate(final1,ncent)
       centfinal.galaxy = galinfo[cent].galaxy
       centfinal.oh = hiiinfo[cent].hii_pt05_log12oh_nuclear[0]
       centfinal.oh_err = hiiinfo[cent].hii_pt05_log12oh_nuclear[1]
       centfinal.rr25 = replicate(rr25_cent,ncent)
       centfinal.comment = 'Cent' ; string(rr25_cent,format='(G0)')+'*R25'
    endif else begin
; central (R/R25=0.1) metallicities!
       cent = where((hiiinfo.hii_pt05_log12oh_central[0] gt -900.0),ncent)
       centfinal = replicate(final1,ncent)
       centfinal.galaxy = galinfo[cent].galaxy
       centfinal.oh = hiiinfo[cent].hii_pt05_log12oh_central[0]
       centfinal.oh_err = hiiinfo[cent].hii_pt05_log12oh_central[1]
       centfinal.rr25 = replicate(rr25_cent,ncent)
       centfinal.comment = 'Cent' ; string(rr25_cent,format='(G0)')+'*R25'
    endelse

; only use the average metallicity if the characteristic metallicity
; doesn't exist (so we're not double-counting HII regions)
    avg = where(((hiiinfo.hii_pt05_log12oh_avg[0] gt -900.0) and $
      (hiiinfo.hii_pt05_log12oh_char[0] lt -900.0)),navg)
    avgfinal = replicate(final1,navg)
    avgfinal.galaxy = galinfo[avg].galaxy
    avgfinal.oh = hiiinfo[avg].hii_pt05_log12oh_avg[0]
    avgfinal.oh_err = hiiinfo[avg].hii_pt05_log12oh_avg[1]
    avgfinal.rr25 = replicate(rr25_avg,navg)
;   avgfinal.rr25 = hiiinfo[avg].hii_pt05_rr25_avg[0]
    avgfinal.comment = 'HIIAvg'

;; electron temperature abundance: just keep the dwarf galaxies
;    ohte = where((hiiinfo.hii_te_log12oh_avg[0] gt -900.0) and $
;      (strtrim(hiiinfo.galaxy,2) ne 'NGC0628') and $
;      (strtrim(hiiinfo.galaxy,2) ne 'NGC0925') and $
;      (strtrim(hiiinfo.galaxy,2) ne 'NGC5194') and $
;      (strtrim(hiiinfo.galaxy,2) ne 'NGC6946') and $
;      (strtrim(hiiinfo.galaxy,2) ne 'NGC7793'),nohte)
;    ohtefinal = replicate(final1,nohte)
;    ohtefinal.galaxy = galinfo[ohte].galaxy
;    ohtefinal.oh = hiiinfo[ohte].hii_te_log12oh_avg[0]
;    ohtefinal.oh_err = hiiinfo[ohte].hii_te_log12oh_avg[1]
;    ohtefinal.comment = 'Te'

;   gg = [nucfinal,circumfinal,stripfinal]
;   hh = [avgfinal,charfinal,centfinal]
;   splog, mean(hh.oh_err), mean(gg.oh_err)
    
; now put everything together and average in two bins of R/R25
    final_pt05 = [nucfinal,circumfinal,stripfinal,charfinal,centfinal,avgfinal]
    final_pt05.galaxy = strtrim(final_pt05.galaxy,2)
    final_pt05 = final_pt05[sort(final_pt05.galaxy)]
;   struct_print, final_pt05

    openw, lun, repstr(outfile,'.fits','_pt05.txt'), /get_lun
    for igal = 0, ngal-1 do begin
       printf, lun, '--------------------------------------------------'
       printf, lun, '### '+result[igal].galaxy
       these = where(strtrim(result[igal].galaxy,2) eq $
         strtrim(final_pt05.galaxy,2),nthese)
       if (nthese ne 0) then begin
;         struct_print, final_pt05[these]
; identify metallicities corresponding to the inner and outer regions,
; which will correspond to our "central" and "characteristic"
; metallicities; if only either "inner" or "outer" metallicities
; exist, then, in the absence of any other information, make the 
; central and characteristic metallicities equal
          inner = where(final_pt05[these].rr25 lt rr25_cut,ninner,$
            comp=outer,ncomp=nouter)
          if (ninner ne 0) then begin
             printf, lun, ' '
             for ii = 0, ninner-1 do printf, lun, 'Inner: ', final_pt05[these[inner[ii]]].oh, $
               final_pt05[these[inner[ii]]].oh_err, final_pt05[these[inner[ii]]].class, $
               final_pt05[these[inner[ii]]].comment, format='(A7,F7.5,2x,F7.5,2x,A9,A0)'
             oh = final_pt05[these[inner]].oh
             oh_err = final_pt05[these[inner]].oh_err
             wavg = sings_weighted_mean(oh,oh_err,wmean_err=werr)
             result[igal].log12oh_pt05_central = [wavg,werr>0.01]
             result[igal].log12oh_pt05_central_comment = strjoin(final_pt05[these[inner]].comment,'/')
             printf, lun, '  Avg: ', result[igal].log12oh_pt05_central[0], result[igal].log12oh_pt05_central[1], $
               result[igal].log12oh_pt05_central_comment, format='(A7,F7.5,2x,F7.5,2x,A0)'
             if (nouter eq 0) then begin
                result[igal].log12oh_pt05_char = result[igal].log12oh_pt05_central
                result[igal].log12oh_pt05_char_comment = 'Central'
             endif
;            if strmatch(result[igal].galaxy,'*0925*') then stop
          endif
          if (nouter ne 0) then begin
             printf, lun, ' '
             for ii = 0, nouter-1 do printf, lun, 'Outer: ', final_pt05[these[outer[ii]]].oh, $
               final_pt05[these[outer[ii]]].oh_err, final_pt05[these[outer[ii]]].class, $
               final_pt05[these[outer[ii]]].comment, format='(A7,F7.5,2x,F7.5,2x,A9,A0)'
             oh = final_pt05[these[outer]].oh
             oh_err = final_pt05[these[outer]].oh_err
             wavg = sings_weighted_mean(oh,oh_err,wmean_err=werr)
             result[igal].log12oh_pt05_char = [wavg,werr>0.01]
             result[igal].log12oh_pt05_char_comment = strjoin(final_pt05[these[outer]].comment,'/')
             printf, lun, '  Avg: ', result[igal].log12oh_pt05_char[0], result[igal].log12oh_pt05_char[1], $
               result[igal].log12oh_pt05_char_comment, format='(A7,F7.5,2x,F7.5,2x,A0)'
             if (ninner eq 0) then begin
                result[igal].log12oh_pt05_central = result[igal].log12oh_pt05_char
                result[igal].log12oh_pt05_central_comment = 'Char'
             endif
          endif 
; if (O/H)_central<(O/H)_char then take the weighted average of the
; two values
          if (result[igal].log12oh_pt05_char[0] gt result[igal].log12oh_pt05_central[0]) then begin
             printf, lun, ' '
             printf, lun, ' ** (O/H)_central<(O/H)_char **'
             oh = [result[igal].log12oh_pt05_char[0],result[igal].log12oh_pt05_central[0]]
             oh_err = [result[igal].log12oh_pt05_char[1],result[igal].log12oh_pt05_central[1]]
             wavg = sings_weighted_mean(oh,oh_err,wmean_err=werr)

             result[igal].log12oh_pt05_char = [wavg,werr>0.01]
             result[igal].log12oh_pt05_central = result[igal].log12oh_pt05_char
             comment = 'cent>char:'+result[igal].log12oh_pt05_central_comment+$
               '/'+result[igal].log12oh_pt05_char_comment
             result[igal].log12oh_pt05_char_comment = comment
             result[igal].log12oh_pt05_central_comment = comment
             printf, lun, 'Final: ', result[igal].log12oh_pt05_char[0], result[igal].log12oh_pt05_char[1], $
               result[igal].log12oh_pt05_char_comment, format='(A7,F7.5,2x,F7.5,2x,A0)'
          endif
; average the abundances of certain special cases
          if $
            strmatch(result[igal].galaxy,'*ngc3773*',/fold) or $
            strmatch(result[igal].galaxy,'*ic2574*',/fold) or $
            strmatch(result[igal].galaxy,'*mrk0033*',/fold) or $
            strmatch(result[igal].galaxy,'*ngc2976*',/fold) or $
            strmatch(result[igal].galaxy,'*ngc4536*',/fold) or $
            strmatch(result[igal].galaxy,'*ngc1512*',/fold) then begin
             printf, lun, ' '
             printf, lun, ' ** Special case **'
             oh = [result[igal].log12oh_pt05_char[0],result[igal].log12oh_pt05_central[0]]
             oh_err = [result[igal].log12oh_pt05_char[1],result[igal].log12oh_pt05_central[1]]
             wavg = sings_weighted_mean(oh,oh_err,wmean_err=werr)
             result[igal].log12oh_pt05_char = [wavg,werr>0.01]
             result[igal].log12oh_pt05_central = result[igal].log12oh_pt05_char
             comment = 'specialcase:'+result[igal].log12oh_pt05_central_comment+$
               '/'+result[igal].log12oh_pt05_char_comment
             result[igal].log12oh_pt05_char_comment = comment
             result[igal].log12oh_pt05_central_comment = comment
             printf, lun, 'Final: ', result[igal].log12oh_pt05_char[0], result[igal].log12oh_pt05_char[1], $
               result[igal].log12oh_pt05_char_comment, format='(A7,F7.5,2x,F7.5,2x,A0)'
          endif
       endif else printf, lun, '  No abundances!'
;      if strmatch(sings[igal].galaxy,'*3627*') then stop
       printf, lun, ' '
    endfor
    free_lun, lun

; use the characteristic abundances to fit the B-band LZ relation
;   mbpivot = -18.0
    mbpivot = 0.0
    kk04_lzfit = sings_lzfit(result.mb,result.log12oh_kk04_char[0],pivot=mbpivot)
    pt05_lzfit = sings_lzfit(result.mb,result.log12oh_pt05_char[0],pivot=mbpivot)
    lzfit = [kk04_lzfit,pt05_lzfit]
    struct_print, lzfit
    
; now assign LZ abundances to the full sample
    result.log12oh_kk04_lz[0] = poly(result.mb-mbpivot,kk04_lzfit.coeff)
    result.log12oh_pt05_lz[0] = poly(result.mb-mbpivot,pt05_lzfit.coeff)
    result.log12oh_kk04_lz[1] = 0.2 ; kk04_lzfit.scatter
    result.log12oh_pt05_lz[1] = 0.2 ; pt05_lzfit.scatter

; print out some info    
    out = im_struct_trimtags(result,select=['GALAXY','TYPE','mb',$
      'LOG12OH_PT05_CENTRAL','LOG12OH_KK04_CENTRAL','LOG12OH_KK04_CENTRAL_COMMENT',$
      'LOG12OH_PT05_CHAR','LOG12OH_KK04_CHAR','LOG12OH_KK04_CHAR_COMMENT'],$
      newtags=['GALAXY','type','mb','PT05_12OH_CENT','KK04_12OH_CENT','CENT_COMMENT',$
      'PT05_12OH_CHAR','KK04_12OH_CHAR','CHAR_COMMENT'])
    need = where((result.log12oh_kk04_char[0] lt -900.0),nneed)
    struct_print, out & print
;   struct_print, out[need]

; write out
    splog, 'Writing '+outfile
    mwrfits, result, outfile, /create
    mwrfits, final_kk04, outfile
    mwrfits, final_pt05, outfile
    spawn, 'gzip -f '+outfile, /sh

    outfile = outpath+'lzfit_'+version+'.fits'
    im_mwrfits, lzfit, outfile, clobber=clobber

return
end
