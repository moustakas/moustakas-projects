pro build_sings_table_galaxies, emulateapj=emulateapj
; build Table 10 in the paper on all the observed HII regions

    if keyword_set(emulateapj) then emapj = '.emapj' else emapj = ''
    
    scale = 1D15
    sigcut = 1.0
    
    metpath = sings_path(/projects)+'log12oh/'
    texpath = sings_path(/papers)+'log12oh/'

    version = sings_log12oh_version()
    classinfo = mrdfits(metpath+'sings_class_'+version+'.fits',1,/silent)

    sings = sings_read_info()
    ngal = n_elements(sings)
    sings = struct_addtags(sings,replicate({empty1: '', empty2: '', empty3: ''},ngal))

    result = mrdfits(metpath+'sings_log12oh_'+version+'.fits.gz',1,/silent)
    result = struct_addtags(result,replicate({empty1: '', empty2: '', empty3: ''},ngal))

    final = mrdfits(metpath+'sings_log12oh_final_'+version+'.fits.gz',1,/silent)
    final = struct_addtags(final,replicate({empty1: '', empty2: '', empty3: ''},ngal))

    allnuclear = read_sings_gandalf(/nuclear)
    alldrift20 = read_sings_gandalf(/drift20)
    alldrift56 = read_sings_gandalf(/drift56)
    srt = sort(allnuclear.ra) & allnuclear = allnuclear[srt]
    srt = sort(alldrift20.ra) & alldrift20 = alldrift20[srt]
    srt = sort(alldrift56.ra) & alldrift56 = alldrift56[srt]

; ---------------------------------------------------------------------------    
; final abundances
    select = ['nice_galaxy','log12oh_kk04_central','log12oh_kk04_char','log12oh_kk04_lz',$
      'empty1','log12oh_pt05_central','log12oh_pt05_char','log12oh_pt05_lz']
    format = ['A0','F12.2','F12.2','F12.2','A0','F12.2','F12.2','F12.2']
    texcenter = ['l','c','c','c','c','c','c','c']

    colhead1 = '\colhead{} & \multicolumn{3}{c}{12+log(O/H)$_{\rm KK04}$} & \colhead{} & '+$
      '\multicolumn{3}{c}{12+log(O/H)$_{\rm PT05}$} \\'
    colhead2 = '\cline{2-4} \cline{6-8}'
    colhead3 = '\colhead{'+[ ['Galaxy','Central','Characteristic','$L-Z$','',$
      'Central','Characteristic']+'} & ',['$L-Z$']+'} ']
    colunits = '\colhead{'+[ ['(1)','(2)','(3)','(4)','','(5)','(6)']+'} & ',['(7)']+'} ']

    table = html_structure_parse(final,tags=select,tagformats=format,$
      tagindices=match,/keep_error,blank='\nodata')
    ntags = n_tags(table)
    tags = tag_names(table)

    tablenotetext = [$
      '{a}{The \pagel{} branch of NGC~5474 based on its circumnuclear spectrum is '+$
      'ambiguous, according to the criteria defined in \S\ref{sec:intnucoh}.}']

; notes
;   splog, 'Add the HoIX note and electron temperatures here!!!  Also note that Croxall+09 get 8.65 for HoIX!!!!'
;   hoix = where(strmatch(sings.galaxy,'*HolmbergIX*',/fold),nhoix)
;   if (nhoix ne 0L) then table[hoix].nice_galaxy = strtrim(table[hoix].nice_galaxy,2)+'\tablenotemark{a}'
;   tablenotetext = [$
;     '{a}{\citet{miller95a} derive $\logoh\simeq8.2\pm0.2$ from modeling a supernova '+$
;     'remnant in Ho~IX.}']

; comments
    tablecomments = ['Final central and characteristic nebular oxygen abundances for the SINGS galaxies based on '+$
      'the KK04 and PT05 abundance calibrations.  The central abundance '+$
      'characterizes the inner, or nuclear metallicity of each galaxy, while the characteristic abundance is a '+$
      'globally-averaged metallicity.  Note that the uncertainties on the oxygen abundances include statistical measurement '+$
      'errors only.  We also provide the oxygen abundances inferred from the $B$-band '+$
      'luminosity-metallicity ($L-Z$) relation for every object; the uncertainty in the $L-Z$ abundance is '+$
      'at least $0.2$~dex (see \S\ref{sec:final} and Fig.~\ref{fig:lz}).']
    tablecomments = strjoin(tablecomments,' ')
    
; replace the comments with numbers    
;   table.log12oh_kk04_char_comment = repstr(repstr(repstr(repstr(repstr(table.log12oh_kk04_char_comment,'RadialStrip','RS'),$
;     'LZRelation','LZ'),'R=0.4*R25','Char'),'Circumnuclear','Circum'),'HII-RegionAverage','Avg')

; write out

    texfile = 'finaloh'+emapj+'.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
       printf, lun, '\pagestyle{empty}'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;   if keyword_set(emulateapj) then printf, lun, '\tablefontsize{\footnotesize}'
;   printf, lun, '\tabletypesize{\small}'
    printf, lun, '\tablecaption{Final Oxygen Abundances\label{table:final}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    niceprintf, lun, colhead3
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0L, ngal-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else if (i lt ngal-1L) then suffix = ' \\ ' else suffix = ''
          if strmatch(tags[j],'*galaxy*',/fold) then line[j] = table[i].(j)+suffix
          if strmatch(tags[j],'*empty*',/fold) then line[j] = ' '+suffix
          if (strmatch(tags[j],'*log12oh*',/fold) eq 1B) and (strmatch(tags[j],'*log12oh*lz*',/fold) eq 0B) then begin
             if (strtrim(table[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else $
               line[j] = '$'+strtrim(table[i].(j)[0],2)+'\pm'+strtrim(table[i].(j)[1],2)+'$'+suffix
          endif
          if (strmatch(tags[j],'*log12oh*lz*',/fold) eq 1B) then begin
             if (strtrim(table[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else $
               line[j] = '$'+strtrim(table[i].(j)[0],2)+'$'+suffix
          endif
       endfor
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
;   niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; ---------------------------------------------------------------------------    
; abundances

    nucselect = ['nice_galaxy','r23_branch','nuclear_r23','nuclear_o32','nuclear_p','nuclear_logu_kk04','nuclear_log12oh_kk04','nuclear_log12oh_pt05']
    d20select = ['nice_galaxy','r23_branch','drift20_r23','drift20_o32','drift20_p','drift20_logu_kk04','drift20_log12oh_kk04','drift20_log12oh_pt05']
    d56select = ['nice_galaxy','r23_branch','drift56_r23','drift56_o32','drift56_p','drift56_logu_kk04','drift56_log12oh_kk04','drift56_log12oh_pt05']
      
    format = ['A0','A0','F12.2','F12.2','F12.2','F12.2','F12.2','F12.2']
    texcenter = ['l','c','c','c','c','c','c','c']

    colhead1 = '\colhead{Galaxy} & \colhead{$R_{23}$ Branch} & \colhead{$R_{23}$} & \colhead{$O_{32}$} & '+$
      '\colhead{$P$} & \colhead{$\log\,(U)_{\rm KK04}$} & \colhead{$12+\log\,({\rm O/H})_{\rm KK04}$} & \colhead{$12+\log\,({\rm O/H})_{\rm PT05}$} \\'
    colunits = '\colhead{'+[ ['(1)','(2)','(3)','(4)','(5)','(6)','(7)']+'} & ',['(8)']+'} ']

    tablecomments = [$
      '(1) Galaxy name; ',$
      '(2-5) Adopted $R_{23}$ branch (U=upper, L=lower, A=ambiguous), \pagel{} parameter, and excitation parameters \ioniz{} and $P$ (see \S\ref{sec:calib} and \S\ref{sec:intnucoh}); ',$
      '(6-8) Ionization parameter and oxygen abundances based on the KK04 and PT05 strong-line calibrations.  Note that the '+$
      'uncertainties on the ionization parameter and oxygen abundances include statistical measurement errors only.']
    tablecomments = strjoin(tablecomments,' ')

; parse the measurements, but only keep objects with at least one
; abundance measurement    
    nuctable = html_structure_parse(result,tags=nucselect,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error)
    keep = where((strmatch(nuctable.nuclear_r23[0],'*nodata*',/fold) eq 0),nnuc)
    nuctable = nuctable[keep]
    nuctags = tag_names(nuctable)

    d20table = html_structure_parse(result,tags=d20select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error)
    keep = where((strmatch(d20table.drift20_r23[0],'*nodata*',/fold) eq 0),nd20)
    d20table = d20table[keep]
    d20tags = tag_names(d20table)

    d56table = html_structure_parse(result,tags=d56select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error)
    keep = where((strmatch(d56table.drift56_r23[0],'*nodata*',/fold) eq 0),nd56)
    d56table = d56table[keep]
    d56tags = tag_names(d56table)

    ntags = n_tags(nuctable)

; write out the integrated abundances table

    texfile = 'intnucoh'+emapj+'.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       if keyword_set(emulateapj) then begin
          printf, lun, '\documentclass{emulateapj}' 
;         printf, lun, '\usepackage{lscape}' 
       endif else begin
          printf, lun, '\documentclass[10pt,preprint]{aastex}'
       endelse
       printf, lun, '\begin{document}'
       printf, lun, '\pagestyle{empty}'
    endif
    if keyword_set(emulateapj) then begin
       printf, lun, '\LongTables'
;      printf, lun, '\begin{landscape}'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;   printf, lun, '\tabletypesize{\scriptsize}'
    printf, lun, '\tabletypesize{\small}'
    if (keyword_set(emulateapj) eq 0) then printf, lun, '\rotate'
    printf, lun, '\tablecaption{Nuclear, Circumnuclear, \& Radial Strip Oxygen Abundances\label{table:intnucoh}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
; nuclear    
    printf, lun, '\cline{1-8}'
    printf, lun, '\multicolumn{8}{c}{Nuclear} \\'
    printf, lun, '\cline{1-8}'
;   printf, lun, '\cutinhead{Nuclear}'
    for i = 0, nnuc-1 do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else begin
             if (i le nnuc-1) then suffix = ' \\ '
          endelse
          if strmatch(nuctags[j],'*galaxy*',/fold) then line[j] = nuctable[i].(j)+suffix
          if strmatch(nuctags[j],'*branch*',/fold) then line[j] = nuctable[i].(j)+suffix
          if (strmatch(nuctags[j],'*galaxy*',/fold) eq 0B) and (strmatch(nuctags[j],'*branch*',/fold) eq 0B) and $
            (strmatch(nuctags[j],'*log12oh*',/fold) eq 0B) and (strmatch(nuctags[j],'*logu*',/fold) eq 0B) then begin
             if (strtrim(nuctable[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else $
               line[j] = '$'+strtrim(nuctable[i].(j)[0],2)+'\pm'+strtrim(nuctable[i].(j)[1],2)+'$'+suffix
          endif
          if strmatch(nuctags[j],'*log12oh*',/fold) then begin
             if (strtrim(nuctable[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else begin
                pre1 = '' & suf1 = ''
;               if strmatch(log12oh[i].class,'*AGN*',/fold) then begin
;                  pre1 = '(' & suf1 = ')' ; put parantheses around AGNs
;               endif
                note = ''
                line[j] = pre1+'$'+strtrim(nuctable[i].(j)[0],2)+'\pm'+strtrim(nuctable[i].(j)[1],2)+'$'+suf1+note+suffix
             endelse
          endif
          if strmatch(nuctags[j],'*logu*',/fold) then begin
             if (strtrim(nuctable[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else $
               line[j] = '$'+strtrim(nuctable[i].(j)[0],2)+'\pm'+strtrim(nuctable[i].(j)[1],2)+'$'+suffix
          endif 
       endfor
       printf, lun, line
    endfor
; drift20
    printf, lun, '\cline{1-8}'
    printf, lun, '\multicolumn{8}{c}{Circumnuclear} \\'
    printf, lun, '\cline{1-8}'
;   printf, lun, '\cutinhead{Circumnuclear}'
    for i = 0, nd20-1 do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else begin
             if (i le nd20-1) then suffix = ' \\ '
          endelse
          if strmatch(d20tags[j],'*galaxy*',/fold) then line[j] = d20table[i].(j)+suffix
          if strmatch(d20tags[j],'*branch*',/fold) then begin
             if strmatch(d20table[i].nice_galaxy,'*ngc 5474*',/fold) then begin ; ambiguous R23 branch
                thisbranch = 'A'
                note = '\tablenotemark{b}' 
             endif else begin
                thisbranch = d20table[i].(j)
                note = ''
             endelse
             line[j] = thisbranch+note+suffix
          endif else note = ''
          if (strmatch(d20tags[j],'*galaxy*',/fold) eq 0B) and (strmatch(d20tags[j],'*branch*',/fold) eq 0B) and $
            (strmatch(d20tags[j],'*log12oh*',/fold) eq 0B) and (strmatch(d20tags[j],'*logu*',/fold) eq 0B) then begin
             if (strtrim(d20table[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else $
               line[j] = '$'+strtrim(d20table[i].(j)[0],2)+'\pm'+strtrim(d20table[i].(j)[1],2)+'$'+suffix
          endif
          if strmatch(d20tags[j],'*log12oh*',/fold) then begin
             if (strtrim(d20table[i].(j)[1],2) lt 0.0) then begin
                note = ''
;               if strmatch(d20table[i].nice_galaxy,'*5474*',/fold) and $
;                 strmatch(d20tags[j],'*drift20_log12oh_pt05*',/fold) then note = '\tablenotemark{a}'
                line[j] = '\nodata'+note+suffix 
                note = ''
             endif else begin
                pre1 = '' & suf1 = ''
;               if strmatch(log12oh[i].class,'*AGN*',/fold) then begin
;                  pre1 = '(' & suf1 = ')' ; put parantheses around AGNs
;               endif
                line[j] = pre1+'$'+strtrim(d20table[i].(j)[0],2)+'\pm'+strtrim(d20table[i].(j)[1],2)+'$'+suf1+note+suffix
             endelse
          endif 
          if strmatch(d20tags[j],'*logu*',/fold) then begin
             if (strtrim(d20table[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else $
               line[j] = '$'+strtrim(d20table[i].(j)[0],2)+'\pm'+strtrim(d20table[i].(j)[1],2)+'$'+suffix
          endif 
       endfor
       printf, lun, line
    endfor
; drift56
    printf, lun, '\cline{1-8}'
    printf, lun, '\multicolumn{8}{c}{Radial Strip} \\'
    printf, lun, '\cline{1-8}'
;   printf, lun, '\cutinhead{Radial Strip}'
    for i = 0L, nd56-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else begin
             if (i lt nd56-1) then suffix = ' \\ ' else suffix = ''
          endelse
          if strmatch(d56tags[j],'*galaxy*',/fold) then line[j] = d56table[i].(j)+suffix
          if strmatch(d56tags[j],'*branch*',/fold) then line[j] = d56table[i].(j)+suffix
          if (strmatch(d56tags[j],'*galaxy*',/fold) eq 0B) and (strmatch(d56tags[j],'*branch*',/fold) eq 0B) and $
            (strmatch(d56tags[j],'*log12oh*',/fold) eq 0B) and (strmatch(d56tags[j],'*logu*',/fold) eq 0B) then begin
             if (strtrim(d56table[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else $
               line[j] = '$'+strtrim(d56table[i].(j)[0],2)+'\pm'+strtrim(d56table[i].(j)[1],2)+'$'+suffix
          endif
          if strmatch(d56tags[j],'*log12oh*',/fold) then begin
             if (strtrim(d56table[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else begin
                pre1 = '' & suf1 = ''
;               if strmatch(log12oh[i].class,'*AGN*',/fold) then begin
;                  pre1 = '(' & suf1 = ')' ; put parantheses around AGNs
;               endif
                note = ''
;               if strmatch(sings[i].galaxy,'*ngc4536*',/fold) and strmatch(d56tags[j],'*drift56_log12oh_kk04*',/fold) then note = '\tablenotemark{c}'
;               if strmatch(sings[i].galaxy,'*ngc4536*',/fold) and strmatch(d56tags[j],'*drift56_log12oh_pt05*',/fold) then note = '\tablenotemark{c}'
                line[j] = pre1+'$'+strtrim(d56table[i].(j)[0],2)+'\pm'+strtrim(d56table[i].(j)[1],2)+'$'+suf1+note+suffix
             endelse
          endif
          if strmatch(d56tags[j],'*logu*',/fold) then begin
             if (strtrim(d56table[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else $
               line[j] = '$'+strtrim(d56table[i].(j)[0],2)+'\pm'+strtrim(d56table[i].(j)[1],2)+'$'+suffix
          endif 
       endfor
       printf, lun, line
    endfor
   
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(emulateapj) then begin
;      printf, lun, '\clearpage'
;      printf, lun, '\end{landscape}'
    endif
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; ---------------------------------------------------------------------------    
; reddening
    select = ['nice_galaxy',$
      'nuclear_hahb','nuclear_ebv',$
      'empty1',$
      'drift20_hahb','drift20_ebv',$
      'empty2',$
      'drift56_hahb','drift56_ebv']
    format = ['A0',$
      'F12.2','F12.2',$
      'A0',$
      'F12.2','F12.2',$
      'A0',$
      'F12.2','F12.2']
    texcenter = ['l',$
      'c','c',$
      'c',$
      'c','c',$
      'c',$
      'c','c']

    colhead1 = '\colhead{} & \multicolumn{2}{c}{Nuclear} & \colhead{} & \multicolumn{2}{c}{Circumnuclear} & \colhead{} & \multicolumn{2}{c}{Radial Strip} \\'
    colhead2 = '\cline{2-3} \cline{5-6} \cline{8-9}'
    colhead3 = '\colhead{Galaxy} & \colhead{H$\alpha$/H$\beta$} & \colhead{$E(B-V)$} & \colhead{} & '+$
      '\colhead{H$\alpha$/H$\beta$} & \colhead{$E(B-V)$} & \colhead{} & \colhead{H$\alpha$/H$\beta$} & \colhead{$E(B-V)$}'
    colunits = '\colhead{'+[ ['(1)','(2)','(3)','','(4)','(5)','','(6)']+'} & ',['(7)']+'} ']

    table = html_structure_parse(result,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error)
    ntags = n_tags(table)
    tags = tag_names(table)

; only keep objects with at least one reddening measurement
    keep = where($
      (strmatch(table.nuclear_hahb[0],'*nodata*',/fold) eq 0) or $
      (strmatch(table.drift20_hahb[0],'*nodata*',/fold) eq 0) or $
      (strmatch(table.drift56_hahb[0],'*nodata*',/fold) eq 0),nobj)
    table = table[keep]
    
    tablecomments = [$
      'Balmer decrement, H$\alpha$/H$\beta$, and corresponding color excess, $E(B-V)$, in magnitudes, asuming '+$
      'the \citet{odonnell94a} Milky Way extinction curve and an intrinsic (unreddened) decrement of $2.86^{+0.18}_{-0.11}$.  '+$
      'Note that the tabulated uncertainties only include the statistical measurement errors of the emission lines and the '+$
      'uncertainty in the adopted intrinsic Balmer decrement; they do not include possible systematic errors arising from, '+$
      'for example, a breakdown of our assumption of a foreground dust screen geometry (e.g., \citealt{witt00a}, but see \citealt{kenn09a}).']
    tablecomments = strjoin(tablecomments,' ')

; write out the reddening table
    texfile = 'reddening'+emapj+'.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       if keyword_set(emulateapj) then begin
          printf, lun, '\documentclass{emulateapj}' 
;         printf, lun, '\usepackage{lscape}' 
       endif else begin
          printf, lun, '\documentclass[10pt,preprint]{aastex}'
       endelse
       printf, lun, '\begin{document}'
       printf, lun, '\pagestyle{empty}'
    endif
    if keyword_set(emulateapj) then begin
       printf, lun, '\LongTables'
;      printf, lun, '\begin{landscape}'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\small}'
;   if (not keyword_set(emulateapj)) then printf, lun, '\rotate'
    printf, lun, '\tablecaption{Nebular Dust Reddening\label{table:reddening}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    niceprintf, lun, colhead3
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0, nobj-1 do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else begin
             if (i lt nobj-1) then suffix = ' \\ ' else suffix = ''
          endelse
          if strmatch(tags[j],'*galaxy*',/fold) then line[j] = table[i].(j)+suffix
          if strmatch(tags[j],'*empty*',/fold) then line[j] = ' '+suffix
          if (strmatch(tags[j],'*galaxy*',/fold) eq 0B) and (strmatch(tags[j],'*empty*',/fold) eq 0B) then begin
             if (strtrim(table[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else begin
                line[j] = '$'+strtrim(table[i].(j)[0],2)+'\pm'+strtrim(table[i].(j)[1],2)+'$'+suffix
             endelse
          endif
       endfor
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
;   niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(emulateapj) then begin
;      printf, lun, '\clearpage'
;      printf, lun, '\end{landscape}'
    endif
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; ---------------------------------------------------------------------------    
; emission-line equivalent widths

    select = ['nice_galaxy','oii_3727_ew','h_gamma_ew',$
      'h_beta_ew','oiii_5007_ew','h_alpha_ew',$
      'nii_6584_ew','sii_6716_ew','sii_6731_ew']
    select_fluxes = ['nice_galaxy',['oii_3727','h_gamma',$
      'h_beta','oiii_5007','h_alpha',$
      'nii_6584','sii_6716','sii_6731']]
    format = [$
      'A0','F15.7','F15.7',$
      'F15.7','F15.7','F15.7',$
      'F15.7','F15.7','F15.7']
    texcenter = [$
      'l','c','c',$
      'c','c','c',$
      'c','c','c']
    colhead = '\colhead{'+[ [$
      'Galaxy','[O~{\sc ii}]$~\lambda3727$','H$\gamma~\lambda4340$',$
      'H$\beta~\lambda4861$','[O~{\sc iii}]$~\lambda5007$','H$\alpha~\lambda6563$',$
      '[N~{\sc ii}]$~\lambda6584$','[S~{\sc ii}]$~\lambda6716$']+'} & ',['[S~{\sc ii}]$~\lambda6731$']+'}']

    ntags = n_elements(select)
    
    tablecomments = ['Rest-frame emission-line equivalent widths in \AA.  Note that the '+$
      'uncertainties include only statistical measurement uncertainties; they do not include '+$
      'systematic errors due to, for example, imperfect continuum subtraction.  For most applications we recommend that a minimum ' +$
      '${\rm S/N}>2$ signal-to-noise ratio cut be applied to the EWs listed in this table.']
    tablecomments = strjoin(tablecomments,' ')

;   tablenotetext = [$
;     '{a}{This Balmer emission line is significantly affected by underlying stellar absorption; therefore '+$
;     'the measured equivalent width may be problematic.}']

    nuclear = struct_trimtags(allnuclear,select=select)
    for i = 0L, n_elements(allnuclear)-1L do begin 
       for iline = 1L, n_tags(nuclear)-1L do begin ; offset from GALAXY
          if (nuclear[i].(iline)[0]/nuclear[i].(iline)[1] le sigcut) then $
            nuclear[i].(iline)[1] = -1.0 ; kludge!
       endfor 
    endfor
    nuclear = html_structure_parse(nuclear,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error,/keep_upperlimits)

    drift20 = struct_trimtags(alldrift20,select=select)
    for i = 0L, n_elements(alldrift20)-1L do begin 
       for iline = 1L, n_tags(drift20)-1L do begin ; offset from GALAXY
          if (drift20[i].(iline)[0]/drift20[i].(iline)[1] le sigcut) then $
            drift20[i].(iline)[1] = -1.0
       endfor 
    endfor
    drift20 = html_structure_parse(drift20,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error,/keep_upperlimits)

    drift56 = struct_trimtags(alldrift56,select=select)
    for i = 0L, n_elements(alldrift56)-1L do begin 
       for iline = 1L, n_tags(drift56)-1L do begin ; offset from GALAXY
          if (drift56[i].(iline)[0]/drift56[i].(iline)[1] le sigcut) then $
            drift56[i].(iline)[1] = -1.0 ; kludge!
       endfor 
    endfor
    drift56 = html_structure_parse(drift56,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error,/keep_upperlimits)

; write out the EW table

    texfile = 'ews'+emapj+'.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       if keyword_set(emulateapj) then begin
          printf, lun, '\documentclass{emulateapj}' 
          printf, lun, '\usepackage{lscape}' 
       endif else begin
          printf, lun, '\documentclass[10pt,preprint]{aastex}'
       endelse
       printf, lun, '\begin{document}'
       printf, lun, '\pagestyle{empty}'
    endif
    if keyword_set(emulateapj) then begin
       printf, lun, '\LongTables'
;      printf, lun, '\begin{landscape}'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\small}'
    if (keyword_set(emulateapj) eq 0) then printf, lun, '\rotate'
    printf, lun, '\tablecaption{Optical Emission-Line Equivalent Widths\label{table:ews}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
    printf, lun, '\cline{1-9}'
    printf, lun, '\multicolumn{9}{c}{Nuclear} \\'
    printf, lun, '\cline{1-9}'
;   printf, lun, '\cutinhead{Nuclear}'
; nuclear    
    for i = 0L, n_elements(nuclear)-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j eq 0L) then line[j] = nuclear[i].(j)+' & ' else begin ; GALAXY tag
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '
             if (strtrim(nuclear[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix
             if (strtrim(nuclear[i].(j)[1],2) gt 0.0) then begin
                flux = nuclear[i].(j)[0]
                ferr = nuclear[i].(j)[1]
                note = ''
                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+note+suffix
             endif
          endelse
       endfor
       printf, lun, line
    endfor
; drift20    
    printf, lun, '\cline{1-9}'
    printf, lun, '\multicolumn{9}{c}{Circumnuclear} \\'
    printf, lun, '\cline{1-9}'
;   printf, lun, '\cutinhead{Circumnuclear}'
    for i = 0L, n_elements(drift20)-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j eq 0L) then line[j] = drift20[i].(j)+' & ' else begin ; GALAXY tag
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '
             if (strtrim(drift20[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix
             if (strtrim(drift20[i].(j)[1],2) gt 0.0) then begin
                flux = drift20[i].(j)[0]
                ferr = drift20[i].(j)[1]
                note = ''
                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+note+suffix
             endif
          endelse
       endfor
       printf, lun, line
    endfor
; drift56
    printf, lun, '\cline{1-9}'
    printf, lun, '\multicolumn{9}{c}{Radial Strip} \\'
    printf, lun, '\cline{1-9}'
;   printf, lun, '\cutinhead{Radial Strip}'
    for i = 0L, n_elements(drift56)-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j eq 0L) then line[j] = drift56[i].(j)+' & ' else begin ; GALAXY tag
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '
             if (strtrim(drift56[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix
             if (strtrim(drift56[i].(j)[1],2) gt 0.0) then begin
                flux = drift56[i].(j)[0]
                ferr = drift56[i].(j)[1]
                note = ''
                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+note+suffix
             endif
          endelse
       endfor
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
;   niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(emulateapj) then begin
;      printf, lun, '\clearpage'
;      printf, lun, '\end{landscape}'
    endif
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; ---------------------------------------------------------------------------    
; spectral classifications

    select = ['nice_galaxy','nuclear_class','drift20_class','drift56_class','ho_class','class']
    format = ['A0','A0','A0','A0','A0','A0']
    texcenter = ['l','c','c','c','c','c']

    colhead = '\colhead{'+[ ['Galaxy','Nuclear','Circumnuclear','Radial Strip','Ho97']+'} & ',['Adopted']+'} \\ ']
    colunits = '\colhead{'+[ ['(1)','(2)','(3)','(4)','(5)']+'} & ',['(6)']+'} ']

    tablecomments = ['Optical spectral classifications based on our nuclear, circumnuclear, and '+$
      'radial-strip spectra and the \nii/\ha{} vs.~\oiii/\hb{} emission-line diagnostic diagram (see Fig.~\ref{fig:bpt}).  '+$
      'For comparison we also tabulate the classifications based on the \citet[hereafter Ho97]{ho97a} nuclear spectra.  '+$
      'Ellipses indicate no data are available, while a question mark indicates that one or more emission line '+$
      'failed our ${\rm S/N}>2$ requirement (see \S\ref{sec:class}).  We list the final, adopted classification for each galaxy, usually based '+$
      'on the nuclear or circumnuclear spectrum, in the last column.  See \S\ref{sec:class} for additional details.']
    tablecomments = strjoin(tablecomments,' ')

    table = html_structure_parse(classinfo,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error)
    ntags = n_tags(table)
    tags = tag_names(table)

; these comments were figured out "by hand" and are described in more
; detail in the text
    ngc0584 = speclinefit_locate(sings,['ngc0584'])
    table[ngc0584].class = strtrim(table[ngc0584].class,2)+'\tablenotemark{a}'
    
;   n2flag = where(classinfo.nuclear_n2flag or classinfo.drift20_n2flag or classinfo.drift56_n2flag)
;   table[n2flag].class = strtrim(table[n2flag].class,2)+'\tablenotemark{q}'
    n2flag = where(classinfo.nuclear_n2flag); & help, n2flag
    table[n2flag].nuclear_class = strtrim(table[n2flag].nuclear_class,2)+'\tablenotemark{b}'
    n2flag = where(classinfo.drift20_n2flag); & help, n2flag
    table[n2flag].drift20_class = strtrim(table[n2flag].drift20_class,2)+'\tablenotemark{b}'
    n2flag = where(classinfo.drift56_n2flag); & help, n2flag
    table[n2flag].drift56_class = strtrim(table[n2flag].drift56_class,2)+'\tablenotemark{b}'
;   niceprint, table[n2flag].nice_galaxy, table[n2flag].class
    
    ngc1377 = speclinefit_locate(sings,['ngc1377'])
    table[ngc1377].drift20_class = strtrim(table[ngc1377].drift20_class,2)+'\tablenotemark{c}'
    
    ngc1404 = speclinefit_locate(sings,['ngc1404'])
    table[ngc1404].class = strtrim(table[ngc1404].class,2)+'\tablenotemark{d}'

    ngc = speclinefit_locate(sings,['ngc1097','ngc1482','ngc1512'])
    table[ngc].drift56_class = table[ngc].drift56_class+'\tablenotemark{e}'
    ngc = speclinefit_locate(sings,['ngc2403','ngc3198'])
    table[ngc].drift20_class = table[ngc].drift20_class+'\tablenotemark{e}'

    broad = where(classinfo.nuclear_broad or classinfo.drift20_broad or classinfo.drift56_broad)
    table[broad].class = strtrim(table[broad].class,2)+'\tablenotemark{f}'
    
    dwarfs = speclinefit_locate(sings,['holmbergii','m81dwa','holmbergi','holmbergix','m81dwb',$
      'ddo053','ddo154','ddo165','ngc5408','ic4710','ngc6822'],/exact)
    table[dwarfs].class = strtrim(table[dwarfs].class,2)+'\tablenotemark{g}'

    ngc3034 = speclinefit_locate(sings,['ngc3034'])
    table[ngc3034].nuclear_class = strtrim(table[ngc3034].nuclear_class,2)+'\tablenotemark{h}'
    table[ngc3034].drift20_class = strtrim(table[ngc3034].drift20_class,2)+'\tablenotemark{h}'

    ngc6946 = speclinefit_locate(sings,['ngc6946'])
    table[ngc6946].nuclear_class = strtrim(table[ngc6946].nuclear_class,2)+'\tablenotemark{i}'
    table[ngc6946].drift20_class = strtrim(table[ngc6946].drift20_class,2)+'\tablenotemark{i}'

    ngc7552 = speclinefit_locate(sings,['ngc7552'])
    table[ngc7552].class = strtrim(table[ngc7552].class,2)+'\tablenotemark{j}'

    tablenotetext = [$
      '{a}{The lower limit on the nuclear [N~{\sc ii}]/H$\alpha$ ratio for NGC~0584 suggests that it hosts a weak AGN.}',$
      '{b}{These spectra were classified as SF or AGN on the basis of the \nii/\ha{} ratio alone because '+$
      'of an unconstrained \oiii/\hb{} ratio, where we adopt $\log\,(\nii/\ha)=-0.25$ as the boundary between the two classes.}',$
      '{c}{NGC~1377 is a ``nascent starburst" with a highly dust-obscured nucleus and AGN-like emission-line ratios '+$
      '\citep{roussel06a}; therefore, although the formal classification based on our circumnuclear spectrum is an AGN, '+$
      'we have changed it to SF.}',$
      '{d}{\citet{tajer05a} report an X-ray luminosity for NGC~1404 of $\log\,(L_{X})=41.46$~erg~s$^{-1}$, '+$
      'suggesting the presence of a weak AGN in this elliptical galaxy.}',$
      '{e}{The radial strip spectra of NGC~1097, NGC~1482, and NGC~1512 formally yield a SF/AGN classification; however, '+$
      'the more physically realistic SF classification is consistent with the errors in the emission-line ratios, and with '+$
      'the uncertainty in the classification method.  For the same reasons we have changed the circumnuclear classifications of '+$
      'NGC~2403 and NGC~3198 from SF/AGN to SF.}',$
      '{f}{One or more spectra of these objects exhibits broad Balmer emission lines.}',$
      '{g}{Ho~II, M~81~Dw~A, DDO~053, Ho~I, Ho~IX, M~81~DwB, DDO~154, DDO~165, NGC~5408, IC~4710, '+$
      'and NGC~6822 are clearly star-forming galaxies; '+$
      'they are diffuse, low-luminosity ($M_{B}\gtrsim-18$~mag) galaxies with no well-defined bulge or nucleus.}',$
      '{h}{Our nuclear and circumnuclear spectra of the starburst galaxy NGC~3034=M~82 exhibit enhanced '+$
      '[N~{\sc ii}]/H$\alpha$ ratios, likely from shock-ionized gas, resulting in AGN-like spectral classifications; '+$
      'therefore, we have changed these classifications to SF, although we note that recent evidence suggests that '+$
      'M~82 may host a low-luminosity AGN \citep{matsumoto01a}.}',$
      '{i}{NGC~6946 is undergoing a powerful nuclear starburst, resulting in a shock-enhanced [N~{\sc ii}]/H$\alpha$ ratio '+$
      'in our nuclear and circumnuclear spectra; therefore, the formal SF/AGN classifications indicated by these two spectra '+$
      'have been changed to SF.}',$
      '{j}{We classify NGC~7552 as SF using the nuclear emission-line flux ratios published by \citet{kewley01a}, '+$
      'although \citet{durret88a} argue that this object may host a weak AGN.}']

; write out

    texfile = 'class'+emapj+'.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
       printf, lun, '\pagestyle{empty}'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\small}'
    printf, lun, '\tablecaption{Optical Spectral Classifications\label{table:class}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0L, ngal-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else if (i lt ngal-1L) then suffix = ' \\ ' else suffix = ''
          line[j] = repstr(table[i].(j),'HII','SF')+suffix 
;         line[j] = repstr(table[i].(j),'HII','H~{\sc ii}')+suffix
       endfor
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; ---------------------------------------------------------------------------    
; emission-line fluxes

    select = ['nice_galaxy',['oii_3727','h_gamma',$
      'h_beta','oiii_5007','h_alpha',$
      'nii_6584','sii_6716','sii_6731']]
    format = [$
      'A0','F15.7','F15.7',$
      'F15.7','F15.7','F15.7',$
      'F15.7','F15.7','F15.7']
    texcenter = [$
      'l','c','c',$
      'c','c','c',$
      'c','c','c']
    colhead = '\colhead{'+[ [$
      'Galaxy','[O~{\sc ii}]$~\lambda3727$','H$\gamma~\lambda4340$',$
      'H$\beta~\lambda4861$','[O~{\sc iii}]$~\lambda5007$','H$\alpha~\lambda6563$',$
      '[N~{\sc ii}]$~\lambda6584$','[S~{\sc ii}]$~\lambda6716$']+'} & ',['[S~{\sc ii}]$~\lambda6731$']+'}']

    ntags = n_elements(select)
    tags = strupcase(select)
    
    tablecomments = ['Optical emission-line fluxes in units of $10^{-15}$~erg~s$^{-1}$~cm$^{-2}$, '+$
      'corrected for foreground Galactic extinction \citep[$R_{V}=3.1$;][]{odonnell94a, schlegel98a}. '+$
      'Note that the uncertainties include only '+$
      'statistical measurement uncertainties, and do not include systematic '+$
      'errors due to, for example, imperfect continuum subtraction.  For most applications we recommend that a ' +$
      'minimum ${\rm S/N}>2$ signal-to-noise ratio cut be applied to the fluxes listed in this table.']
    tablecomments = strjoin(tablecomments,' ')

; multiply the line fluxes and errors by SCALE    
    nuclear = struct_trimtags(allnuclear,select=select)
    for i = 0L, n_elements(allnuclear)-1L do begin 
       for iline = 1L, n_tags(nuclear)-1L do begin ; offset from GALAXY
          flux = nuclear[i].(iline)[0]
          ferr = nuclear[i].(iline)[1]
          if (flux/ferr gt sigcut) then begin
             nuclear[i].(iline)[0] = scale*flux
             nuclear[i].(iline)[1] = scale*ferr
          endif else nuclear[i].(iline)[1] = -1.0 ; kludge!
       endfor
    endfor
    nuclear = html_structure_parse(nuclear,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error,/keep_upperlimits)

    drift20 = struct_trimtags(alldrift20,select=select)
    for i = 0L, n_elements(alldrift20)-1L do begin 
       for iline = 1L, n_tags(drift20)-1L do begin ; offset from ATLAS_ID
          flux = drift20[i].(iline)[0]
          ferr = drift20[i].(iline)[1]
          if (flux/ferr gt sigcut) then begin
             drift20[i].(iline)[0] = scale*flux
             drift20[i].(iline)[1] = scale*ferr
          endif else drift20[i].(iline)[1] = -1.0 ; kludge!
       endfor
    endfor
    drift20 = html_structure_parse(drift20,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error,/keep_upperlimits)

    drift56 = struct_trimtags(alldrift56,select=select)
    for i = 0L, n_elements(alldrift56)-1L do begin 
       for iline = 1L, n_tags(drift56)-1L do begin ; offset from ATLAS_ID
          flux = drift56[i].(iline)[0]
          ferr = drift56[i].(iline)[1]
          if (flux/ferr gt sigcut) then begin
             drift56[i].(iline)[0] = scale*flux
             drift56[i].(iline)[1] = scale*ferr
          endif else drift56[i].(iline)[1] = -1.0 ; kludge!
       endfor
    endfor
    drift56 = html_structure_parse(drift56,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error,/keep_upperlimits)

; write out the line flux table
    texfile = 'fluxes'+emapj+'.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       if keyword_set(emulateapj) then begin
          printf, lun, '\documentclass{emulateapj}' 
          printf, lun, '\usepackage{lscape}' 
       endif else begin
          printf, lun, '\documentclass[10pt,preprint]{aastex}'
       endelse
       printf, lun, '\begin{document}'
       printf, lun, '\pagestyle{empty}'
    endif
    if keyword_set(emulateapj) then begin
       printf, lun, '\LongTables'
;      printf, lun, '\begin{landscape}'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\small}'
    if (keyword_set(emulateapj) eq 0) then printf, lun, '\rotate'
    printf, lun, '\tablecaption{Optical Emission-Line Fluxes\label{table:fluxes}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
; nuclear    
    printf, lun, '\cline{1-9}'
    printf, lun, '\multicolumn{9}{c}{Nuclear} \\'
    printf, lun, '\cline{1-9}'
;   printf, lun, '\cutinhead{Nuclear}'
    for i = 0L, n_elements(nuclear)-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j eq 0L) then line[j] = nuclear[i].(j)+' & ' else begin ; GALAXY tag
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '
             if (strtrim(nuclear[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix
             if (strtrim(nuclear[i].(j)[1],2) gt 0.0) then begin
                note = ''
;               if strmatch(tags[j],'*h_gamma*',/fold) and $
;                 (allnuclear[i].h_gamma_ew[0] lt allnuclear[i].babs_h_gamma_ew[0]) then $
;                   note = '\tablenotemark{a}'
;               if strmatch(tags[j],'*h_beta*',/fold) and $
;                 (allnuclear[i].h_beta_ew[0] lt allnuclear[i].babs_h_beta_ew[0]) then $
;                   note = '\tablenotemark{a}'
;               if strmatch(tags[j],'*h_alpha*',/fold) and $
;                 (allnuclear[i].h_alpha_ew[0] lt allnuclear[i].babs_h_alpha_ew[0]) then $
;                   note = '\tablenotemark{a}'
                flux = nuclear[i].(j)[0]
                ferr = nuclear[i].(j)[1]
                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+note+suffix
             endif
          endelse
       endfor
       printf, lun, line
    endfor
; drift20    
    printf, lun, '\cline{1-9}'
    printf, lun, '\multicolumn{9}{c}{Circumnuclear} \\'
    printf, lun, '\cline{1-9}'
;   printf, lun, '\cutinhead{Circumnuclear}'
    for i = 0L, n_elements(drift20)-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j eq 0L) then line[j] = drift20[i].(j)+' & ' else begin ; GALAXY tag
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '
             if (strtrim(drift20[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix
             if (strtrim(drift20[i].(j)[1],2) gt 0.0) then begin
                flux = drift20[i].(j)[0]
                ferr = drift20[i].(j)[1]
                note = ''
                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+note+suffix
             endif
          endelse
       endfor
       printf, lun, line
    endfor
; drift56
    printf, lun, '\cline{1-9}'
    printf, lun, '\multicolumn{9}{c}{Radial Strip} \\'
    printf, lun, '\cline{1-9}'
;   printf, lun, '\cutinhead{Radial Strip}'
    for i = 0L, n_elements(drift56)-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j eq 0L) then line[j] = drift56[i].(j)+' & ' else begin ; GALAXY tag
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '
             if (strtrim(drift56[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix
             if (strtrim(drift56[i].(j)[1],2) gt 0.0) then begin
                flux = drift56[i].(j)[0]
                ferr = drift56[i].(j)[1]
                note = ''
;               if strmatch(tags[j],'*h_gamma*',/fold) and $
;                 (alldrift56[i].h_gamma_ew[0] lt alldrift56[i].babs_h_gamma_ew[0]) then $
;                   note = '\tablenotemark{a}'
;               if strmatch(tags[j],'*h_beta*',/fold) and $
;                 (alldrift56[i].h_beta_ew[0] lt alldrift56[i].babs_h_beta_ew[0]) then $
;                   note = '\tablenotemark{a}'
;               if strmatch(tags[j],'*h_alpha*',/fold) and $
;                 (alldrift56[i].h_alpha_ew[0] lt alldrift56[i].babs_h_alpha_ew[0]) then $
;                   note = '\tablenotemark{a}'
                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+note+suffix
             endif
          endelse
       endfor 
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
;   niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(emulateapj) then begin
;      printf, lun, '\clearpage'
;      printf, lun, '\end{landscape}'
    endif
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; ---------------------------------------------------------------------------    
; journal of observations

    select = ['nice_galaxy','nuclear_latexap','nuclear_posangle','nuclear_exptime','nuclear_photflag','nuclear_comments',$
      'empty1','drift20_latexap','drift20_posangle','drift20_exptime','drift20_photflag','drift20_comments',$
      'empty2','drift56_latexap','drift56_posangle','drift56_exptime','drift56_photflag','drift56_comments']
    format = ['A0','A0','I0','I0','A0','A0','A0','A0','I0','I0','A0','A0','A0','A0','I0','I0','A0','A0']
    texcenter = ['l','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c']

    colhead1 = ['\colhead{} & ','\multicolumn{5}{c}{Nuclear} & ','\colhead{} & ','\multicolumn{5}{c}{Circumnuclear} & ',$
      '\colhead{} & ','\multicolumn{5}{c}{Radial Strip} \\']
    colhead2 = '\cline{2-6}\cline{8-12}\cline{14-18}'
    colhead3 = ['\colhead{Galaxy} & ',$
      '\colhead{Aperture\tablenotemark{a}} & ','\colhead{$\theta_{\rm slit}$\tablenotemark{b}} & ',$
      '\colhead{$t_{\rm exp}$\tablenotemark{c}} & ','\colhead{Phot\tablenotemark{d}} & ','\colhead{Notes\tablenotemark{e}} & ',$
      '\colhead{} & ',$
      '\colhead{Aperture\tablenotemark{a}} & ','\colhead{$\theta_{\rm slit}$\tablenotemark{b}} & ',$
      '\colhead{$t_{\rm exp}$\tablenotemark{c}} & ','\colhead{Phot\tablenotemark{d}} & ','\colhead{Notes\tablenotemark{e}} & ',$
      '\colhead{} & ',$
      '\colhead{Aperture\tablenotemark{a}} & ','\colhead{$\theta_{\rm slit}$\tablenotemark{b}} & ',$
      '\colhead{$t_{\rm exp}$\tablenotemark{c}} & ','\colhead{Phot\tablenotemark{d}} & ','\colhead{Notes\tablenotemark{e}}']

    tablecomments = ['Summary of our nuclear, circumnuclear, and radial-strip observations of the '+$
      'SINGS galaxies.']

;    
;The
;\emph{effective} exposure time of the integrated spectra, $t_{\rm
;  eff}$, or the time spent on a spatially fixed location in each
;galaxy, is given by $t\times(\Delta_{\rm slit}/\Delta_{\rm scan})$,
;where $\Delta_{\rm slit}$ is the slit width in arcseconds and
;$\Delta_{\rm scan}$ is the drift-scan length perpendicular to the slit
;in the same units \citep{moustakas06a}.      
    
    tablenotetext = [$
      '{a}{Spectrophotometric extraction aperture along and perpendicular to the slit.}',$
      '{b}{Long-slit position angle, measured positive from North to East, in degrees.}',$
      '{c}{Effective exposure time, as defined in \citet{moustakas06a}, in seconds.}',$
      '{d}{Flag indicating whether the spectra were obtained during photometric (Y) or non-photometric conditions (N).}',$
      '{e}{Extraction notes: (1) Extended spatial profile; (2) One or more foreground stars subtracted; '+$
      '(3) Multiple pointings stitched together; '+$
;     '(4) No well-defined nucleus ; '+$
      '(4) Extraction centered on bright \ion{H}{2} region; '+$
;     '(6) Centered on IRS spectroscopic aperture at $\alpha_{\rm J2000} = $09:56:05.2, $\delta_{\rm J2000} = $+69:41:22.0 ; '+$
      '(5) Residual foreground stellar contamination.}']

    table = html_structure_parse(sings,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error)
    ntags = n_tags(table)
    tags = tag_names(table)

; only keep objects that were actually observed!
    keep = where((strtrim(table.nuclear_latexap,2) ne '') and $
      (strtrim(table.drift20_latexap,2) ne '') and $
      (strtrim(table.drift56_latexap,2) ne ''),nobj)
    table = table[keep]
    
; add comments
    table.nuclear_comments = ''
    table.drift20_comments = ''
    table.drift56_comments = ''

    w = where(strmatch(sings.nuclear_comments,'*extended*',/fold),nw)
    if (nw ne 0L) then table[w].nuclear_comments = '1'
;   w = where(strmatch(sings.nuclear_comments,'*no nucle*',/fold),nw) 
;   if (nw ne 0L) then table[w].nuclear_comments = '4'
    w = where(strmatch(sings.nuclear_comments,'*centered*',/fold),nw) 
    if (nw ne 0L) then table[w].nuclear_comments = '4' ; '5'
    w = where(strmatch(sings.drift20_comments,'*extended*',/fold),nw) 
    if (nw ne 0L) then table[w].drift20_comments = '1'
    w = where(strmatch(sings.drift20_comments,'*centered*',/fold),nw)
    if (nw ne 0L) then table[w].drift20_comments = '4' ; '5'
    
    w = where((strmatch(sings.drift56_comments,'*foregrou*',/fold) eq 1B) and $
      (strmatch(sings.drift56_comments,'*stitch*',/fold) eq 0B),nw)
    if (nw ne 0L) then table[w].drift56_comments = '2'
    w = where((strmatch(sings.drift56_comments,'*foregrou*',/fold) eq 1B) and $
      (strmatch(sings.drift56_comments,'*stitch*',/fold) eq 1B),nw) 
    if (nw ne 0L) then table[w].drift56_comments = '2,3'
    w = where((strmatch(sings.drift56_comments,'*contamin*',/fold) eq 1B),nw) 
    if (nw ne 0L) then table[w].drift56_comments = '5' ; '7' ; tol89
    
;   w = where(strmatch(sings.galaxy,'*ngc3034*',/fold)) & table[w].drift56_comments = '6'

; write out
    texfile = 'journal'+emapj+'.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
       printf, lun, '\pagestyle{empty}'
    endif
    if keyword_set(emulateapj) then begin
       printf, lun, '\LongTables'
       printf, lun, '\begin{landscape}'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\scriptsize}'
;   printf, lun, '\tabletypesize{\small}'
    if (keyword_set(emulateapj) eq 0) then printf, lun, '\rotate'
    printf, lun, '\tablecaption{Journal of Observations\label{table:journal}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    niceprintf, lun, colhead3
;   niceprintf, lun, colhead4
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0, nobj-1 do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else begin
             if (i lt nobj-1) then suffix = ' \\ ' else suffix = ''
          endelse
          if (strcompress(table[i].(j),/remove) eq '') and (strmatch(tags[j],'*empty*',/fold) eq 0) then $
            line[j] = '\nodata'+suffix else line[j] = table[i].(j)+suffix
       endfor
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(emulateapj) then begin
       printf, lun, '\clearpage'
       printf, lun, '\end{landscape}'
    endif
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; ---------------------------------------------------------------------------    
; general properties    

    select = ['nice_galaxy','type','ebv_mw','r25','incl','pa','mb','bv','distance','distance_method','distance_texref']
    format = ['A0','A0','F12.3','F12.2','I0','I0','F12.2','F12.3','F12.3','A0','A0']
    texcenter = ['l','l','c','c','c','c','c','c','c','c','c']

;   colhead1 = '\colhead{'+[ ['','','$E(B-V)$','$\rho_{25}$','$i$','$\theta$',$
;     '$D$','Distance']+'} & ',['Distance']+'} \\ ']
;   colhead2 = '\colhead{'+[ ['Galaxy','Type','(mag)','(arcmin)','(deg)','(deg)',$
;     '(Mpc)','Method']+'} & ',['Ref.']+'} \\']
    colhead1 = '\colhead{'+[ ['','','$E(B-V)$','$\rho_{25}$',$
      '$i$','$\theta$','$M_{B}$','$B-V$','$D$','Distance']+'} & ',['Distance']+'} \\ ']
    colhead2 = '\colhead{'+[ ['Galaxy','Type','(mag)','(arcmin)','(deg)','(deg)','(mag)','(mag)',$
      '(Mpc)','Method']+'} & ',['Ref.']+'} \\']
    colunits = '\colhead{'+[ ['(1)','(2)','(3)','(4)','(5)','(6)','(7)','(8)','(9)','(10)']+'} & ',['(11)']+'} ']

;   niceprint, select, format, texcenter, colhead1, colhead2, colunits
    table = html_structure_parse(result,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata')
    ntags = n_tags(table)
    tags = tag_names(table)

    table.nice_galaxy = strtrim(table.nice_galaxy,2)

    table.distance_method = repstr(table.distance_method,'FLOW','Flow')
    table.distance_method = repstr(table.distance_method,'STARS','BS')
    table.distance_method = repstr(table.distance_method,'CEPH','Ceph')
    table.distance_method = repstr(table.distance_method,'CMULT','Mult')
    table.distance_method = repstr(table.distance_method,'MEM','Mem')
    
; parse the distance references and assign unique reference numbers
    refs = strsplit(strjoin(strtrim(result.distance_texref,2),','),',',/extract)
    nrefs = n_elements(refs)

    refnum = lonarr(nrefs)-1L
    refcount = 1L
    for i = 0L, nrefs-1L do begin
       if (refnum[i] eq -1L) then begin
          match = where(refs[i] eq refs)
          refnum[match] = refcount
          refcount = refcount + 1L
       endif
    endfor

;   niceprint, refs, refnum
    urefnum = uniq(refnum,sort(refnum))
    tablerefs = strjoin('('+string(refnum[urefnum],format='(I0)')+') \citet{'+$
      string(refs[urefnum],format='(A0)')+'}','; ')+'.'

; footnotes
    
;   rep = where((sings.twomass_posangle lt -900.0) and (sings.optical_posangle gt -900),nrep)
;   if (nrep ne 0L) then table[rep].pa = reform(string(sings[rep].optical_posangle,format='(I0)'),1,nrep)+'\tablenotemark{a}'

    m51 = where(strmatch(sings.galaxy,'*5194*'),nm51)
    if (nm51 ne 0L) then table[m51].incl = strtrim(table[m51].incl,2)+'\tablenotemark{c}' ; '\tablenotemark{i}'

    n6822 = where(strmatch(sings.galaxy,'*6822*'),nn6822)
    if (nn6822 ne 0L) then begin
       note = '\tablenotemark{e}' ; '\tablenotemark{k}'
       table[n6822].incl = strtrim(table[n6822].incl,2)+note
       table[n6822].pa = strtrim(table[n6822].pa,2)+note
    endif

; generate the object tags

;   table.nice_galaxy = reform('\object['+strtrim(table.nice_galaxy,2)+']{'+strtrim(table.nice_galaxy,2)+'}',1,ngal)
;   table.nice_galaxy = reform('\object['+strtrim(sings.ned_galaxy,2)+']{'+strtrim(table.nice_galaxy,2)+'}',1,ngal)

    tablecomments = ['(1) SINGS galaxy name; ',$
      '(2) Morphological type from \citet{kenn03b}; ',$
      '(3) Color excess used to correct the optical spectroscopy and photometry for '+$
      'foreground Galactic reddening \citep{schlegel98a}; ',$
      '(4) Radius of the major axis at the $\mu_{B}=25$~mag~arcsec$^{-2}$ isophote \citep{devac91a}; ',$
      '(5-6) Galaxy inclination and position angles \citep{devac91a, jarrett03a}.  We '+$
      'compute the inclination angle, $i$, following \citet{tully98a}: $\sin\,(i) = 1.02\times\sqrt{1-(b/a)^2}$, '+$
      'where $(b/a)$ is the minor-to-major axis ratio at $\mu_{B} = '+$
      '25$~mag~arcsec$^{-2}$ \citep{devac91a}, or the isophotal axis ratio at '+$
      '$\mu_{K_{s}} = 20$~mag~arcsec$^{-2}$ \citep{jarrett03a};',$
      '(7-8) Absolute $B$-band magnitude and $B-V$ color, both relative to Vega \citep{dale07a, munoz09a} ',$
      '(see \S\ref{sec:sample} for more details); ',$
      '(9-11) Distance, distance method, and corresponding reference to the literature.  The distance '+$
      'methods have the following meaning: Flow = Hubble distance assuming $H_{0} = '+$
      '70$~\kms~Mpc$^{-1}$ corrected for peculiar motions; '+$
      'SBF = surface-brightness fluctuations; BS = bright stars (supergiants); Ceph = Cepheid variables; '+$
      'Mult = an average of multiple, cross-calibrated techniques; TRGB = tip of the red-giant branch; '+$
      'Mem = group/cluster membership; SN = supernova; PNLF = planetary nebula luminosity function.']
    tablecomments = strjoin(tablecomments,' ')

    tablenotetext = [$
      '{a}{For these objects we adopt the \emph{HST} Key Project Cepheid distance uncorrected for metallicity variations '+$
      '\citep[see][]{freedman01a}.}',$
      '{b}{NGC~4254, NGC~4450, NGC~4569, and NGC~4579 have been placed at the distance of the Virgo cluster '+$
      '\citep[$16.5\pm0.6$~Mpc;][]{mei07a}.}',$
      '{c}{For NGC~5194 we adopt the kinematic inclination angle determined by \citet{tully74a}.}',$
      '{d}{The interacting pair NGC~5194 (M~51~a) and NGC~5195 have been placed at a common distance.}',$
      '{e}{For NGC~6822 we adopt the kinematic '+$
      'inclination and position angle derived by \citet{brandenburg98a}.  We have also adjusted the '+$
      'apparent $B$- and $V$-band magnitudes for this object from \citet{dale07a} by $-0.37$~mag '+$
      'to better match the $B$-band magnitude from \citet{kara04a}.}']

;   tablenotetext = [$
;     '{a}{Foreground Galactic reddening value \citep{schlegel98}.}',$
;     '{b}{Diameter of the major axis at the $\mu_{B} = 25$~mag~arcsec$^{-2}$ isophote.}',$
;     '{c}{Derived inclination angle.  Note that our inclination angles are only '+$
;     'strictly appropriate for spiral galaxies '+$
;     'assuming they are oblate spheroids with an intrinsic axial ratio of $0.2$.}',$
;     '{d}{Projected position angle, measured positive from North to East.}',$
;     '{e}{Absolute $B$-band (Vega) magnitude ($H_{0} = 70~{\rm km~s}^{-1}~{\rm Mpc}^{-1}$) and '+$
;     '$B-V$ color, corrected for Galactic '+$
;     'extinction \citep[$R_{V}=3.1$;][]{odonnell94} but not for inclination or internal extinction}',$
;     '{f}{Distance, distance method, and corresponding reference to the literature: '+$
;     'Flow = distance based on an analysis of the local peculiar velocity field; '+$
;     'SBF = surface-brightness fluctuations; BS = bright stars (supergiants); Ceph = Cepheid variables; '+$
;     'Mult = an average of multiple, cross-calibrated techniques; TRGB = tip of the red-giant branch; '+$
;     'Mem = group/cluster membership; SN = supernova; PNLF = planetary nebula luminosity function.}',$
;     '{g}{We have adopted the \emph{HST} Key Project Cepheid distance uncorrected for metallicity variations '+$
;     '\citep[see][]{freedman01}.}',$
;     '{h}{NGC~4254, NGC~4450, NGC~4569, and NGC~4579 have been placed at the distance of the Virgo cluster '+$
;     '\citep[$16.5\pm0.6$~Mpc;][]{mei07}.}',$
;     '{i}{We have adopted the kinematic inclination angle determined by \citet{tully74} for NGC~5194.}',$
;     '{j}{NGC~5194 and NGC~5195 have been placed at a common distance.}',$
;     '{k}{Following \citet{lee06}, for NGC~6822 we have adopted the kinematic '+$
;     'inclination and position angle derived by \citet{brandenburg98}.}'] 

; write out the general properties table

    texfile = 'properties'+emapj+'.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       if keyword_set(emulateapj) then begin
          printf, lun, '\documentclass{emulateapj}' 
          printf, lun, '\usepackage{lscape}' 
       endif else begin
          printf, lun, '\documentclass[10pt,preprint]{aastex}'
       endelse
       printf, lun, '\begin{document}'
       printf, lun, '\pagestyle{empty}'
    endif
    if keyword_set(emulateapj) then begin
       printf, lun, '\LongTables'
;      printf, lun, '\begin{landscape}'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;   printf, lun, '\tabletypesize{\scriptsize}'
    printf, lun, '\tabletypesize{\small}'
    if (keyword_set(emulateapj) eq 0) then printf, lun, '\rotate'
    printf, lun, '\tablecaption{Properties of the SINGS Galaxies\label{table:properties}}'
    printf, lun, '\tablewidth{0pt}'
;   printf, lun, '\tablewidth{460pt}'
;   printf, lun, '\tablecolumns{11}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0L, ngal-1L do begin
       line = strarr(ntags)
       for j = 0, ntags-1 do begin
          if (j lt ntags-1) then suffix = ' & ' else begin
             if (i lt ngal-1) then suffix = ' \\ ' else suffix = ''
          endelse
          if (strcompress(table[i].(j),/remove) eq '') then $ ; this needs to be first
            line[j] = '\nodata'+suffix else line[j] = table[i].(j)+suffix
          if strmatch(tags[j],'*MB*',/fold) then begin
             flux = result[i].mb & ferr = result[i].mb_err
             note = ''
             if strmatch(sings[i].galaxy,'*6822*') then note = '\tablenotemark{e}'
             if (flux lt -900.0) then line[j] = '\nodata'+suffix else $
               line[j] = '$'+string(flux,format='(F7.2)')+'\pm'+string(ferr,format='(F5.2)')+'$'+note+suffix
          endif
          if strmatch(tags[j],'BV',/fold) then begin
             flux = result[i].bv & ferr = result[i].bv_err
             note = ''
             if strmatch(sings[i].galaxy,'*6822*') then note = '\tablenotemark{e}'
             if (flux lt -900.0) then line[j] = '\nodata'+suffix else $
               line[j] = '$'+string(flux,format='(F7.2)')+'\pm'+string(ferr,format='(F7.2)')+'$'+note+suffix
          endif
          if strmatch(tags[j],'*texref*',/fold) then begin
             theserefs = strsplit(strtrim(table[i].distance_texref,2),',',/extract)
             these = cmset_op(refs,'and',theserefs,/index)
             thesenum = refnum[these[uniq(refnum[these],sort(refnum[these]))]]
             line[j] = strjoin(string(thesenum,format='(I0)'),',')+suffix
          endif
          if strmatch(tags[j],'*distance',/fold) then begin
             flux = result[i].distance
             ferr = result[i].distance_err
             format_flux_error, flux, ferr, newflux, newferr
             note = ''
             if strmatch(sings[i].galaxy,'*0925*') or $ ; strmatch(sings[i].galaxy,'*2403*') or $
               strmatch(sings[i].galaxy,'*3031*') or strmatch(sings[i].galaxy,'*3198*') or $
               strmatch(sings[i].galaxy,'*3351*') or strmatch(sings[i].galaxy,'*3621*') or $
               strmatch(sings[i].galaxy,'*3627*') or strmatch(sings[i].galaxy,'*4321*') or $
               strmatch(sings[i].galaxy,'*4536*') or strmatch(sings[i].galaxy,'*4725*') or $
               strmatch(sings[i].galaxy,'*7331*') then note = '\tablenotemark{a}' ; '\tablenotemark{g}'
             if strmatch(sings[i].galaxy,'*4254*') or strmatch(sings[i].galaxy,'*4450*') or $
               strmatch(sings[i].galaxy,'*4569*') or strmatch(sings[i].galaxy,'*4579*') then $
                 note = '\tablenotemark{b}' ; '\tablenotemark{h}'
             if strmatch(sings[i].galaxy,'*5194*') or strmatch(sings[i].galaxy,'*5195*') then $
               note = '\tablenotemark{d}' ; '\tablenotemark{j}'
             line[j] = '$'+newflux+'\pm'+newferr+'$'+note+suffix
          endif
       endfor
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\tablerefs{'+tablerefs+'}'
    printf, lun, '\end{deluxetable}'
    if keyword_set(emulateapj) then begin
;      printf, lun, '\clearpage'
;      printf, lun, '\end{landscape}'
    endif
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

return
end    
