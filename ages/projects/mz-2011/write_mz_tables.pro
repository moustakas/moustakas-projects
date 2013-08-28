pro write_mz_tables, agesdust, agesancillary, textest=textest, emulateapj=emulateapj
; jm07apr11nyu - based on WRITE_ATLAS_TABLES

    emulateapj = 1L
    nstub = 15L

    texpath = ages_path(/papers)+'mz/'
    mrtpath = ages_path(/papers)+'mz/MRT/'
    if keyword_set(textest) then texpath = texpath+'textest/'

; ---------------------------------------------------------------------------    
; Table: LZ coefficients
; ---------------------------------------------------------------------------    

; read this from a file!    

    table = {$
      sample:     '', $
      z:          '', $
      selection:  '', $
      bandpass:   '', $
      slope:     0.0, $
      intercept: 0.0, $
      scatter:   0.0, $
      npts:       0L}
    table = replicate(table,3L)
    table.sample    = ['AGES','AGES','SDSS']
    table.z         = ['$0.01-0.8$','$0.01-0.8$','$0.033-0.25$']
    table.selection = ['$I<19.95$','$I<19.95$','$r<17.77$']
    table.bandpass  = ['$M_{0.1g}$','$M_{0.1r}$','$M_{0.1g}$']

    ntable = n_elements(table)
    ntags = n_tags(table)
    
    colhead1 = '\colhead{'+[ ['Sample','$z$','Selection\tablenotemark{b}','Bandpass',$
      'Slope','Intercept','Scatter']+'} & ',['N']+'} \\ ']
    colunits = '\colhead{'+[ ['(1)','(2)','(3)','(4)','(5)','(6)','(7)']+'} & ',['(8)']+'} ']

    texcenter = ['c','c','c','c','c','c','c','c']

    tablenotetext = [$
      '{a}{Fits were performed using a ``robust'' ordinaryleast-squares bisector fit.}',$
      '{b}{$I$-band magnitudes are on the Vega system, while $g$- and $r$-band magnitudes are AB-relative.}']
    
; write it out    
    
    texfile = 'ages_lzcoeff.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       if keyword_set(emulateapj) then begin
          printf, lun, '\documentclass[onecolumn]{emulateapj}' 
;         printf, lun, '\usepackage{lscape}' 
       endif else begin
          printf, lun, '\documentclass[10pt,preprint]{aastex}'
       endelse
       printf, lun, '\begin{document}'
       printf, lun, '\pagestyle{empty}'
    endif
    if keyword_set(emulateapj) then begin
;      printf, lun, '\LongTables'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;   printf, lun, '\tabletypesize{\small}'
    printf, lun, '\tablecaption{Luminosity-Metallicity Coefficient Fitting Results\tablenotemark{a}\label{table:lzcoeff}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0L, ntable-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else if (i lt ntable-1L) then suffix = ' \\ ' else suffix = ''
          if size(table[i].(j),/type) eq 7L then line[j] = table[i].(j)+suffix else $
            line[j] = '$'+strtrim(string(table[i].(j),format='(G0.0)'),2)+'$'+suffix
              
       endfor
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then begin
;      printf, lun, '\bibliographystyle{apj}'
;      printf, lun, '\bibliography{/Users/ioannis/ay/bin/bibtex/im}'
       printf, lun, '\end{document}'
    endif
    free_lun, lun

    if keyword_set(textest) then begin
       pushd, texpath
       spawn, 'mkps '+repstr(texfile,'.tex',''), /sh
;      spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

stop    
    
; ---------------------------------------------------------------------------    
; Table: measured properties
; ---------------------------------------------------------------------------    

    if (n_elements(agesdust) eq 0L) then agesdust = read_ages_mz_sample(agesancillary=agesancillary)
    ngalaxy = n_elements(agesdust)
    srt = sort(agesancillary.ages_id)

    ages = replicate({id: '', ra: '', dec: '', z: '', weight: '', $
      m_b: '', m_v: '', m_r: '', mass: '', $
      oii: '', hbeta: '', oiii: ''},ngalaxy)
    ntags = n_tags(ages)

    ages.id     = string(agesancillary[srt].ages_id,format='(I0)')
    ages.ra     = im_dec2hms(agesancillary[srt].ra/15D,/colon)
    ages.dec    = im_dec2hms(agesancillary[srt].dec,/colon)
    pos = where(agesancillary[srt].dec gt 0.0) & ages[pos].dec = '+'+ages[pos].dec

    ages.z      = string(agesancillary[srt].z,format='(F6.4)')
    ages.weight = string(agesancillary[srt].spec_weight,format='(F8.3)')

    ages.m_b  = '$'+string(agesancillary[srt].m_b,format='(F6.2)')+'\pm'+string(agesancillary[srt].m_b_err,format='(F4.2)')+'$'
    ages.m_v  = '$'+string(agesancillary[srt].m_v,format='(F6.2)')+'\pm'+string(agesancillary[srt].m_v_err,format='(F4.2)')+'$'
    ages.m_r  = '$'+string(agesancillary[srt].m_r,format='(F6.2)')+'\pm'+string(agesancillary[srt].m_r_err,format='(F4.2)')+'$'
    ages.mass = '$'+string(agesancillary[srt].kcorr_mass,format='(F5.2)')+'\pm'+string(agesancillary[srt].kcorr_mass_err>0.01,format='(F4.2)')+'$'
    
    format_flux_error, agesdust[srt].oii_3727_ew[0],  agesdust[srt].oii_3727_ew[1],  f_oii,  e_oii  & ages.oii   = '$'+f_oii +'\pm'+e_oii+'$'
    format_flux_error, agesdust[srt].oiii_5007_ew[0], agesdust[srt].oiii_5007_ew[1], f_oiii, e_oiii & ages.oiii  = '$'+f_oiii+'\pm'+e_oiii+'$'
    format_flux_error, agesdust[srt].h_beta_ew[0],    agesdust[srt].h_beta_ew[1],    f_hb,   e_hb   & ages.hbeta = '$'+f_hb  +'\pm'+e_hb+'$'

; replace negative [O III] values with upper limits

    neg = where(agesdust[srt].oiii_5007_ew[0] lt 0.0,nneg)
    if (nneg ne 0L) then ages[neg].oiii = '$>'+string(agesdust[srt[neg]].oiii_5007_ew_limit,format='(F4.1)')+'$' ; make sure the limits are positive!
;   format_flux_error, abs(agesdust[srt[neg]].oiii_5007_ew[0]), agesdust[srt[neg]].oiii_5007_ew[1], f, e
    
    colhead = strjoin('\colhead{'+['ID\tablenotemark{a}','$\alpha_{\rm J2000}$\tablenotemark{a}','$\delta_{\rm J2000}$\tablenotemark{a}',$
      'Redshift\tablenotemark{a}','Weight\tablenotemark{a}',$
      '$M_{B}$\tablenotemark{b}','$M_{V}$\tablenotemark{b}','$M_{R}$\tablenotemark{b}','$\log\,(M/M_{\sun})$\tablenotemark{b}',$
      'EW([O~{\sc ii}~$\lambda3727$])\tablenotemark{c}','EW(H$\beta$)\tablenotemark{c}','EW([O~{\sc iii}~$\lambda5007$])\tablenotemark{c}']+'}',' & ')+' \\'
    colunits = strjoin('\colhead{('+string(lindgen(ntags)+1L,format='(I0)')+')}',' & ') ;+' \\'
    colcenter = replicate('c',ntags)

    tablecomments = ['(1)','(2)']
    tablecomments = strjoin(tablecomments,' ')
    tablenotetext = [$
      '{a}{AGES identification numbers, celestial coordinates, spectroscopic redshifts, and statistical weights are '+$
      'taken from C.~S.~Kochanek et~al. (2007, in prep.) and D.~J. Eisenstein et~al. (2007, in prep.).}',$
      '{b}{Absolute magnitudes (Vega; $h=0.7$) and stellar masses (\citealt{salpeter55} IMF; $0.1-100~M_{\sun}$) have '+$
      'been derived using {\sc k-correct} \citep[ver.~4.1.3;][]{blanton06b}.}',$
      '{c}{Absorption-corrected, rest-frame emission-line equivalent widths in \AA.}']
    
    texfile = 'ages_mzdata.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\usepackage{lscape}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\begin{landscape}'
    printf, lun, '\begin{deluxetable}{'+strjoin(colcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
;   printf, lun, '\rotate'
    printf, lun, '\tablecaption{AGES: Spectrophotometric Properties and Oxygen Abundances\label{table:ages_mzdata}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0L, 24L do begin
;   for i = 0L, ngalaxy-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else begin
             if (i lt ngalaxy-1L) then suffix = ' \\ ' else suffix = ''
          endelse
          if (strcompress(ages[i].(j),/remove) eq '') then $
            line[j] = '\nodata'+suffix else line[j] = ages[i].(j)+suffix
       endfor
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    printf, lun, '\clearpage'
    printf, lun, '\end{landscape}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

stop    

return
end
