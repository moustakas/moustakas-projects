pro build_sings_table_hiiregions, emulateapj=emulateapj
; build Table 10 in the paper on all the observed HII regions
    
    if keyword_set(emulateapj) then emapj = '.emapj' else emapj = ''

    metpath = sings_path(/projects)+'log12oh/'
    texpath = sings_path(/papers)+'log12oh/'

    version = sings_log12oh_version()
    hiifile = metpath+'sings_log12oh_hiiregions_'+version+'.fits.gz'
    splog, 'Reading '+hiifile
    hiiinfo = mrdfits(hiifile,1)
    hii = mrdfits(hiifile,2)
    nhii = n_elements(hii)
    ngal = n_elements(hiiinfo)

; ---------------------------------------------------------------------------    
; HII-region oxygen abundances
    kk04select = ['nice_galaxy','r23_branch',$
      'hii_kk04_slope','hii_kk04_log12oh_central','hii_kk04_log12oh_char',$
      'hii_kk04_log12oh_avg','hii_kk04_nhii_used','hii_kk04_texrefs']
    pt05select = ['nice_galaxy','r23_branch',$
      'hii_pt05_slope','hii_pt05_log12oh_central','hii_pt05_log12oh_char',$
      'hii_pt05_log12oh_avg','hii_pt05_nhii_used','hii_pt05_texrefs']

    format = ['A0','A0',$
      'F12.2','F12.2','F12.2','F12.2','I0','A0']
      
    texcenter = ['l','c',$
      'c','c','c','c','c','c']

    colhead1 = '\colhead{} & \colhead{} & \colhead{Gradient} & \colhead{12+log(O/H)} & \colhead{12+log(O/H)} & \colhead{12+log(O/H)} & \colhead{} & \colhead{} \\'
    colhead2 = '\colhead{Galaxy} & \colhead{$R_{23}$ Branch} & \colhead{(dex~$\rho_{25}^{-1}$)} & \colhead{at $\rho=0$} & \colhead{at $\rho=0.4\,\rho_{25}$} & \colhead{Average} & '+$
      '\colhead{N(H~{\sc ii})} & \colhead{Refs.} \\'
    colunits = '\colhead{'+[ ['(1)','(2)','(3)','(4)','(5)','(6)','(7)']+'} & ',['(8)']+'} ']

; makarova02a, weisz08a
    tablenotetext = [$
      '{a}{Although Ho~IX is among the faintest galaxies in our sample, the \nii/\ha{} and \nii/\oii{} line-ratios '+$
      'clearly place it on the upper \pagel{} branch \citep{croxall09a}.}']

    tablecomments = [$
      '(1) Galaxy name; ',$
      '(2) Adopted $R_{23}$ branch (U=upper, L=lower); ',$
      '(3) Slope of the radial abundance gradient; ',$
      '(4) Central oxygen abundance based on the derived abundance gradient; ',$
      '(5) Characteristic oxygen abundance, defined as the metallicity at $\rho=0.4\rho_{25}$, based on the derived abundance gradient; ',$
      '(6) Unweighted average and standard deviation of all the H~{\sc ii}-region abundances (see \S\ref{sec:hiiavg}); ',$
      '(7) Number of H~{\sc ii} regions with measured abundances; ',$
      '(8) References to the literature from which the \hii-region line-ratios were taken (see Appendix~\ref{sec:hiiappendix}).']
    tablecomments = strjoin(tablecomments,' ')

;   splog, 'Computing N(HII)!'
;   rep = where((log12oh.hii_nhii gt 0L) and (log12oh.hii_kk04_nhii_used eq 0L),nrep)
;   if (nrep ne 0L) then log12oh[rep].hii_kk04_nhii_used = log12oh[rep].hii_nhii

; only keep objects with abundance measurements
    kk04table = html_structure_parse(hiiinfo,tags=kk04select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error)
    keep = where(kk04table.hii_kk04_nhii_used gt 0,nkk04)
    splog, 'KK04 ', nkk04
    kk04table = kk04table[keep]
    kk04tags = tag_names(kk04table)

    pt05table = html_structure_parse(hiiinfo,tags=pt05select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error)
    keep = where(pt05table.hii_pt05_nhii_used gt 0,npt05)
    splog, 'PT05 ', npt05
    pt05table = pt05table[keep]
    pt05tags = tag_names(pt05table)

    ntags = n_tags(kk04table)

; parse the HII-region references and assign unique reference numbers

    kk04_refs = strsplit(strjoin(strtrim(hiiinfo.hii_kk04_texrefs,2),','),',',/extract)
    nkk04_refs = n_elements(kk04_refs)

    kk04_refnum = lonarr(nkk04_refs)-1L
    kk04_refcount = 1L
    for i = 0L, nkk04_refs-1L do begin
       if (kk04_refnum[i] eq -1L) then begin
          kk04_match = where(kk04_refs[i] eq kk04_refs)
          kk04_refnum[kk04_match] = kk04_refcount
          kk04_refcount = kk04_refcount + 1L
       endif
    endfor

;   niceprint, refs, refnum
    kk04_urefnum = uniq(kk04_refnum,sort(kk04_refnum))
    kk04_tablerefs = strjoin('('+string(kk04_refnum[kk04_urefnum],format='(I0)')+') \citet{'+$
      string(kk04_refs[kk04_urefnum],format='(A0)')+'}','; ')+'.'

    finaltablerefs = strarr(ngal)
    for i = 0, nkk04-1 do begin
       if (kk04table[i].hii_kk04_nhii_used gt 0) then begin
          kk04_theserefs = strsplit(strtrim(kk04table[i].hii_kk04_texrefs,2),',',/extract)
          kk04_these = cmset_op(kk04_refs,'and',kk04_theserefs,/index)
          kk04_thesenum = kk04_refnum[kk04_these[uniq(kk04_refnum[kk04_these],sort(kk04_refnum[kk04_these]))]]
          finaltablerefs[i] = strjoin(string(kk04_thesenum,format='(I0)'),',')
       endif else finaltablerefs[i] = '\nodata'
    endfor

; write out
    texfile = 'hiioh'+emapj+'.tex'
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
       printf, lun, '\begin{landscape}'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\small}'
    if (keyword_set(emulateapj) eq 0) then printf, lun, '\rotate'
    printf, lun, '\tablecaption{H~{\sc ii}-Region Oxygen Abundances\label{table:hiioh}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    printf, lun, colhead1
    printf, lun, colhead2
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
; KK04
    printf, lun, '\cline{1-8}'
    printf, lun, '\multicolumn{8}{c}{KK04} \\'
    printf, lun, '\cline{1-8}'
;   printf, lun, '\cutinhead{KK04}'
    for i = 0, nkk04-1 do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else if (i le nkk04-1) then suffix = ' \\ '; else suffix = ''
          if strmatch(kk04table[i].nice_galaxy,'*ho ix*',/fold) then note = '\tablenotemark{a}' else note = ''
          line[j] = kk04table[i].(j)+suffix
          if strmatch(kk04tags[j],'*branch*',/fold) then line[j] = kk04table[i].(j)+note+suffix ; R23 branch tag
          if strmatch(kk04tags[j],'*nhii_used*',/fold) then $
            if (kk04table[i].(j) eq 0L) then line[j] = '\nodata'+suffix else line[j] = kk04table[i].(j)+suffix ; N(HII) tag
          if strmatch(kk04tags[j],'*texrefs*',/fold) then line[j] = finaltablerefs[i]+suffix
          if strmatch(kk04tags[j],'*log12oh*',/fold) or strmatch(kk04tags[j],'*slope*',/fold) then begin
             if (strtrim(kk04table[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else begin
                note = ''
;               if strmatch(sings[i].galaxy,'*ngc6822*',/fold) and strmatch(kk04tags[j],'*hii_kk04_slope*',/fold) then note = '\tablenotemark{a}'
                line[j] = '$'+strtrim(kk04table[i].(j)[0],2)+'\pm'+strtrim(kk04table[i].(j)[1],2)+'$'+note+suffix
             endelse
          endif
       endfor 
       printf, lun, line
    endfor
; PT05
    printf, lun, '\cline{1-8}'
    printf, lun, '\multicolumn{8}{c}{PT05} \\'
    printf, lun, '\cline{1-8}'
;   printf, lun, '\cutinhead{PT05}'
    for i = 0L, npt05-1 do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else if (i lt npt05-1) then suffix = ' \\ ' else suffix = ''
          note = ''
          if strmatch(pt05table[i].nice_galaxy,'*ho ix*',/fold) then note = '\tablenotemark{a}'
          thisbranch = pt05table[i].r23_branch
;         if strmatch(pt05table[i].nice_galaxy,'*ngc 2915*',/fold) or $ ; ambiguous R23 branch
;           strmatch(pt05table[i].nice_galaxy,'*mrk 33*',/fold) then begin
;            thisbranch = 'A'
;            note = '\tablenotemark{b}'
;         endif else thisbranch = pt05table[i].r23_branch
          line[j] = pt05table[i].(j)+suffix
          if strmatch(kk04tags[j],'*branch*',/fold) then line[j] = thisbranch+note+suffix ; R23 branch tag
;         if strmatch(kk04tags[j],'*branch*',/fold) then line[j] = pt05table[i].(j)+suffix ; R23 branch tag
          if strmatch(pt05tags[j],'*nhii_used*',/fold) then $
            if (pt05table[i].(j) eq 0L) then line[j] = '\nodata'+suffix else line[j] = pt05table[i].(j)+suffix ; N(HII) tag
          if strmatch(pt05tags[j],'*texrefs*',/fold) then line[j] = finaltablerefs[i]+suffix
          if strmatch(pt05tags[j],'*log12oh*',/fold) or strmatch(pt05tags[j],'*slope*',/fold) then begin
             if (strtrim(pt05table[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else begin
                note = ''
;               if strmatch(sings[i].galaxy,'*ngc6822*',/fold) and strmatch(pt05tags[j],'*hii_pt05_slope*',/fold) then note = '\tablenotemark{a}'
;               if strmatch(pt05table[i].nice_galaxy,'*ho ix*',/fold) and $
;                 strmatch(pt05tags[j],'*hii_pt05_log12oh_avg*',/fold) then note = '\tablenotemark{a}'
                line[j] = '$'+strtrim(pt05table[i].(j)[0],2)+'\pm'+strtrim(pt05table[i].(j)[1],2)+'$'+note+suffix
             endelse
          endif
       endfor 
       printf, lun, line
    endfor
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\tablerefs{'+kk04_tablerefs+'}'
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
; HII-Region Database

    select = ['number','nice_galaxy','region',$
      'rc3_rr25','r23','o32','p','log12oh_kk04','log12oh_pt05','texref']
    format = ['I3.3','A0','A0','F12.2','F12.2','F12.2',$
      'F12.2','F12.2','F12.2','A0']
    texcenter = ['c','l','l','c','c','c','c','c','c','c']

    colhead1 = '\colhead{'+[ ['No.','Galaxy','Region','$\rho/\rho_{25}$',$
      '$R_{23}$','$O_{32}$','$P$','$12+\log\,({\rm O/H})_{\rm KK04}$',$
      '$12+\log\,({\rm O/H})_{\rm PT05}$']+'} & ',['Ref.']+'} \\ ']
    colunits = '\colhead{'+[ ['(1)','(2)','(3)','(4)','(5)','(6)','(7)',$
      '(8)','(9)']+'} & ',['(10)']+'} ']

; apply a minimum error    
    for ii = 0L, n_elements(hii)-1 do begin
       if (hii[ii].log12oh_kk04[1] gt 0.0) then $
         hii[ii].log12oh_kk04[1] = hii[ii].log12oh_kk04[1]>0.01
       if (hii[ii].log12oh_pt05[1] gt 0.0) then $
         hii[ii].log12oh_pt05[1] = hii[ii].log12oh_pt05[1]>0.01
       if (hii[ii].r23[1] gt 0.0) then hii[ii].r23[1] = hii[ii].r23[1]>0.01
       if (hii[ii].o32[1] gt 0.0) then hii[ii].o32[1] = hii[ii].o32[1]>0.01
       if (hii[ii].p[1] gt 0.0) then hii[ii].p[1] = hii[ii].p[1]>0.01
    endfor

    table = html_structure_parse(hii,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error)
    ntags = n_tags(table)
    tags = tag_names(table)

    table.region = repstr(table.region,'_',' ')

; parse the references and assign unique reference numbers
    refs = reform(strtrim(table.texref,2))
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

    tablecomments = [$
      '(1-2) Unique \hii-region identification number and SINGS galaxy name;',$
      '(3) H~{\sc ii}-region name as listed in the originating reference, if available;',$
      '(4) Deprojected galactocentric radius normalized by the $\rho_{25}$ radius of the galaxy (see Table~\ref{table:properties}); ',$
      '(5-7) Abundance-sensitive \pagel{} parameter and excitation parameters \ioniz{} and $P$ (see \S\ref{sec:calib})',$
      '(8-9) Oxygen abundances derived using the KK04 and PT05 strong-line calibrations '+$
      'assuming the \pagel{} branch listed in Table~\ref{table:intnucoh}, unless otherwise noted; ',$
      '(10) References to the literature from which the \hii-region line-ratios were taken (see Appendix~\ref{sec:hiiappendix}).']
    tablecomments = strjoin(tablecomments,' ')

    tablenotetext = [$
      '{a}{The \pagel{} branch of this \hii{} region is ambiguous '+$
      'according to the criteria defined in \S\ref{sec:intnucoh}.}',$
      '{b}{The oxygen abundance of this \hii{} region is not defined '+$
      'according to the quantitative criteria established in \S\ref{sec:intnucoh}.}']
    
; write out the hii-region database table
    texfile = 'hiidata'+emapj+'.tex'
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
       printf, lun, '\begin{landscape}'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\small}'
    if (keyword_set(emulateapj) eq 0) then printf, lun, '\rotate'
    printf, lun, '\tablecaption{Oxygen Abundances of the H~{\sc ii} Regions in the SINGS Galaxies\label{table:hiidata}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0L, nhii-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else begin
             if (i lt nhii-1L) then suffix = ' \\ ' else suffix = ''
          endelse
          if (strcompress((table[i].(j))[0],/remove) eq '') then $ ; this needs to be first
            line[j] = '\nodata'+suffix else line[j] = table[i].(j)+suffix
          if (strmatch(tags[j],'*number*',/fold) eq 0) and (strmatch(tags[j],'*galaxy*',/fold) eq 0) and $
            (strmatch(tags[j],'*region*',/fold) eq 0) and (strmatch(tags[j],'*rc3_rr25*',/fold) eq 0) and $
            (strmatch(tags[j],'*texref*',/fold) eq 0) then begin
             if (strtrim(table[i].(j)[1],2) lt 0.0) then begin
                note = ''
                if strmatch(tags[j],'*pt05*',/fold) then $
                  branch = hii[i].r23_branch_pt05 else $
                  branch = hii[i].r23_branch_kk04
                if strmatch(branch,'*rej*',/fold) and strmatch(tags[j],'*log12oh*',/fold) then $
                  note = '\tablenotemark{b}'
                line[j] = '\nodata'+note+suffix
             endif else begin
                note = ''
                if strmatch(tags[j],'*pt05*',/fold) then $
                  branch = hii[i].r23_branch_pt05 else $
                  branch = hii[i].r23_branch_kk04
                if (strtrim(branch,2) eq 'A') and strmatch(tags[j],'*log12oh*',/fold) then $
                  note = '\tablenotemark{a}' ; NOTE MINIMUM ERROR!
                line[j] = '$'+strtrim(table[i].(j)[0],2)+'\pm'+$
                  strtrim(table[i].(j)[1],2)+'$'+note+suffix
             endelse
          endif
          if strmatch(tags[j],'*texref*',/fold) then begin
             theserefs = strtrim(table[i].texref,2)
             these = cmset_op(refs,'and',theserefs,/index)
             thesenum = refnum[these[uniq(refnum[these],sort(refnum[these]))]]
             line[j] = strjoin(string(thesenum,format='(I0)'),',')+suffix
          endif
       endfor
       if (strtrim(hii[(i-1L)>0L].galaxy,2) ne strtrim(hii[i].galaxy,2)) then printf, lun, '\cline{1-10}'
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\tablerefs{'+tablerefs+'}'
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

return
end
    

    

;; ---------------------------------------------------------------------------    
;; Table: HII-Region Database
;
;    select = ['number','galaxy','region','raoffset','deoffset',$
;      'rc3_rr25','r23','o32','p','log12oh_kk04','log12oh_pt05','texref']
;    format = ['I3.3','A0','A0','I0','I0','F12.2','F12.2','F12.2',$
;      'F12.2','F12.2','F12.2','A0']
;    texcenter = ['c','l','l','c','c','c','c','c','c','c','c','c']
;
;    colhead1 = '\colhead{'+[ ['No.','Galaxy','Region','$\Delta\alpha$',$
;      '$\Delta\delta$','$\rho/\rho_{25}$',$
;      '$R_{23}$','$O_{32}$','$P$','$12+\log\,({\rm O/H})_{\rm KK04}$',$
;      '$12+\log\,({\rm O/H})_{\rm PT05}$']+'} & ',['Ref.']+'} \\ ']
;    colunits = '\colhead{'+[ ['(1)','(2)','(3)','(4)','(5)','(6)','(7)',$
;      '(8)','(9)','(10)','(11)']+'} & ',['(12)']+'} ']
;
;    table = html_structure_parse(hii,tags=select,tagformats=format,$
;      tagindices=match,blank='\nodata',/keep_error)
;    ntags = n_tags(table)
;    tags = tag_names(table)
;
;    table.region = repstr(table.region,'_',' ')
;
;; parse the references and assign unique reference numbers
;    refs = reform(strtrim(table.texref,2))
;    nrefs = n_elements(refs)
;
;    refnum = lonarr(nrefs)-1L
;    refcount = 1L
;    for i = 0L, nrefs-1L do begin
;       if (refnum[i] eq -1L) then begin
;          match = where(refs[i] eq refs)
;          refnum[match] = refcount
;          refcount = refcount + 1L
;       endif
;    endfor
;
;;   niceprint, refs, refnum
;    urefnum = uniq(refnum,sort(refnum))
;    tablerefs = strjoin('('+string(refnum[urefnum],format='(I0)')+') \citet{'+$
;      string(refs[urefnum],format='(A0)')+'}','; ')+'.'
;
;    tablecomments = [$
;      'Col. (1) Unique running identification number;',$
;      'Col. (2) SINGS galaxy name;',$
;      'Col. (3) H~{\sc ii}-region name as listed in the originating reference, if available;',$
;      'Cols. (4-5) Relative position in right ascension and declination, respectively, relative to the galaxy nucleus in arcseconds;',$ 
;      'Col. (6) Deprojected galactocentric radius normalized by the $\rho_{25}$ radius of the galaxy (see Table~\ref{table:properties}); ',$
;      'Cols. (7-9) Abundance-sensitive \pagel{} parameter and excitation parameters \ioniz{} and $P$ (see \S\ref{sec:calib})',$
;      'Cols. (10-11) Oxygen abundance derived using the KK04 and PT05 strong-line calibration, '+$
;      'respectively, assuming the \pagel{} branch listed in Table~\ref{table:intnucoh}; ',$
;      'Cols. (12) Corresponding references to the literature.']
;    tablecomments = strjoin(tablecomments,' ')
;
;; write out the hii-region database table
;
;    texfile = 'hiidata.tex'
;    splog, 'Writing '+texpath+texfile+'.'
;    openw, lun, texpath+texfile, /get_lun
;    if keyword_set(textest) then begin
;       if keyword_set(emulateapj) then begin
;          printf, lun, '\documentclass{emulateapj}' 
;          printf, lun, '\usepackage{lscape}' 
;       endif else begin
;          printf, lun, '\documentclass[10pt,preprint]{aastex}'
;       endelse
;       printf, lun, '\begin{document}'
;       printf, lun, '\pagestyle{empty}'
;    endif
;    if keyword_set(emulateapj) then begin
;       printf, lun, '\LongTables'
;       printf, lun, '\begin{landscape}'
;    endif
;    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;    printf, lun, '\tabletypesize{\small}'
;    if (keyword_set(emulateapj) eq 0) then printf, lun, '\rotate'
;    printf, lun, '\tablecaption{Oxygen Abundances of the H~{\sc ii} Regions in the SINGS Galaxies\label{table:hiidata}}'
;    printf, lun, '\tablewidth{0pt}'
;    printf, lun, '\tablehead{'
;    niceprintf, lun, colhead1
;    niceprintf, lun, colunits
;    printf, lun, '}'
;    printf, lun, '\startdata'
;
;    for i = 0L, nhii-1L do begin
;       line = strarr(ntags)
;       for j = 0L, ntags-1L do begin
;          if (j lt ntags-1L) then suffix = ' & ' else begin
;             if (i lt nhii-1L) then suffix = ' \\ ' else suffix = ''
;          endelse
;          if (strcompress((table[i].(j))[0],/remove) eq '') then $ ; this needs to be first
;            line[j] = '\nodata'+suffix else line[j] = table[i].(j)+suffix
;          if (strmatch(tags[j],'*number*',/fold) eq 0B) and (strmatch(tags[j],'*galaxy*',/fold) eq 0B) and $
;            (strmatch(tags[j],'*region*',/fold) eq 0B) and (strmatch(tags[j],'*raoffset*',/fold) eq 0B) and $
;            (strmatch(tags[j],'*deoffset*',/fold) eq 0B) and (strmatch(tags[j],'*rc3_rr25*',/fold) eq 0B) and $
;            (strmatch(tags[j],'*texref*',/fold) eq 0B) then begin
;             if (strtrim(table[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else begin
;                note = ''
;                line[j] = '$'+strtrim(table[i].(j)[0],2)+'\pm'+strtrim(table[i].(j)[1],2)+'$'+note+suffix
;             endelse
;          endif
;          if strmatch(tags[j],'*texref*',/fold) then begin
;             theserefs = strtrim(table[i].texref,2)
;             these = cmset_op(refs,'and',theserefs,/index)
;             thesenum = refnum[these[uniq(refnum[these],sort(refnum[these]))]]
;             line[j] = strjoin(string(thesenum,format='(I0)'),',')+suffix
;          endif
;       endfor
;       if (strtrim(hii[(i-1L)>0L].galaxy,2) ne strtrim(hii[i].galaxy,2)) then printf, lun, '\cline{1-12}'
;       printf, lun, line
;    endfor
;    
;    printf, lun, '\enddata'
;    printf, lun, '\tablecomments{'+tablecomments+'}'
;;   niceprintf, lun, '\tablenotetext'+tablenotetext
;    printf, lun, '\tablerefs{'+tablerefs+'}'
;    printf, lun, '\end{deluxetable}'
;    if keyword_set(emulateapj) then begin
;       printf, lun, '\clearpage'
;       printf, lun, '\end{landscape}'
;    endif
;    if keyword_set(textest) then printf, lun, '\end{document}'
;    free_lun, lun
;    if keyword_set(textest) then begin
;       pushd, texpath
;       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
;       popd
;    endif
;
