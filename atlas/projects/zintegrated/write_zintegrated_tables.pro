pro write_zintegrated_tables, textest=textest
; jm05sep19uofa - write data tables for the ZINTEGRATED paper 
; jm06mar26uofa - re-written, based on WRITE_SINGS_TABLES

; read the results from ZINTEGRATED_GRADIENTS

    datapath = atlas_path(/projects)+'zintegrated/'
    texpath = atlas_path(/papers)+'zintegrated/'
    if keyword_set(textest) then texpath = texpath+'textest/'

    gradients1 = mrdfits(datapath+'zintegrated_gradients.fits.gz',1,/silent)
    stats = mrdfits(datapath+'zintegrated_gradients.fits.gz',2,/silent)

;   if (n_elements(gradients1) eq 0L) then gradients = mrdfits(datapath+$
;     'zintegrated_gradients.fits.gz',1,/silent) else gradients = gradients1
    ngalaxy = n_elements(gradients1)
    gradients = struct_addtags(gradients1,replicate({empty1: '', empty2: '', empty3: ''},ngalaxy))

    voffset = '0.5'
    voffset2 = '0.8'

; ---------------------------------------------------------------------------    
; Table: correlation coefficient statistics
; ---------------------------------------------------------------------------    

; replicate a three-element structure, corresponding to the slope of
; the abundance gradient, E(B-V), and the inclination, and to the
; observed, reddening-corrected, and EW abundances

    spearman = replicate({variables: '', m91_coeff: -1.0, m91_prob: -1.0, empty1: '', $
      pt05_coeff: -1.0, pt05_prob: -1.0},12)

    spearman.variables = [$
      '$\Delta\log({\rm O/H})_{\rm obs}$ vs. Slope',$
      '$\Delta\log({\rm O/H})_{\rm obs}$ vs. $E(B-V)$',$
      '$\Delta\log({\rm O/H})_{\rm obs}$ vs. $i$',$
      '$\Delta\log({\rm O/H})_{\rm obs}$ vs. $\Sigma({\rm H}\alpha)$',$

      '$\Delta\log({\rm O/H})_{\rm cor}$ vs. Slope',$
      '$\Delta\log({\rm O/H})_{\rm cor}$ vs. $E(B-V)$',$
      '$\Delta\log({\rm O/H})_{\rm cor}$ vs. $i$',$
      '$\Delta\log({\rm O/H})_{\rm cor}$ vs. $\Sigma({\rm H}\alpha)$',$

      '$\Delta\log({\rm O/H})_{\rm EW}$ vs. Slope',$
      '$\Delta\log({\rm O/H})_{\rm EW}$ vs. $E(B-V)$',$
      '$\Delta\log({\rm O/H})_{\rm EW}$ vs. $i$',$
      '$\Delta\log({\rm O/H})_{\rm EW}$ vs. $\Sigma({\rm H}\alpha)$']

; observed
    
    spearman[0+0*4L].m91_coeff   = stats.int_obs_char_m91_slope_coeff
    spearman[0+0*4L].m91_prob    = stats.int_obs_char_m91_slope_prob
    spearman[0+0*4L].pt05_coeff  = stats.int_obs_char_pt05_slope_coeff
    spearman[0+0*4L].pt05_prob   = stats.int_obs_char_pt05_slope_prob
;   print, stats.int_obs_char_m91_slope_coeff, stats.int_obs_char_m91_slope_prob, $
;     stats.int_obs_char_pt05_slope_coeff, stats.int_obs_char_pt05_slope_prob

    spearman[1+0*4L].m91_coeff   = stats.int_obs_char_m91_ebv_coeff
    spearman[1+0*4L].m91_prob    = stats.int_obs_char_m91_ebv_prob
    spearman[1+0*4L].pt05_coeff  = stats.int_obs_char_pt05_ebv_coeff
    spearman[1+0*4L].pt05_prob   = stats.int_obs_char_pt05_ebv_prob
;   print, stats.int_obs_char_m91_ebv_coeff, stats.int_obs_char_m91_ebv_prob, $
;     stats.int_obs_char_pt05_ebv_coeff, stats.int_obs_char_pt05_ebv_prob

    spearman[2+0*4L].m91_coeff   = stats.int_obs_char_m91_incl_coeff
    spearman[2+0*4L].m91_prob    = stats.int_obs_char_m91_incl_prob
    spearman[2+0*4L].pt05_coeff  = stats.int_obs_char_pt05_incl_coeff
    spearman[2+0*4L].pt05_prob   = stats.int_obs_char_pt05_incl_prob
;   print, stats.int_obs_char_m91_incl_coeff, stats.int_obs_char_m91_incl_prob, $
;     stats.int_obs_char_pt05_incl_coeff, stats.int_obs_char_pt05_incl_prob

    spearman[3+0*4L].m91_coeff   = stats.int_obs_char_m91_hasb_coeff
    spearman[3+0*4L].m91_prob    = stats.int_obs_char_m91_hasb_prob
    spearman[3+0*4L].pt05_coeff  = stats.int_obs_char_pt05_hasb_coeff
    spearman[3+0*4L].pt05_prob   = stats.int_obs_char_pt05_hasb_prob
;   print, stats.int_obs_char_m91_hasb_coeff, stats.int_obs_char_m91_hasb_prob, $
;     stats.int_obs_char_pt05_hasb_coeff, stats.int_obs_char_pt05_hasb_prob

; corrected    
    
    spearman[0+1*4L].m91_coeff   = stats.int_cor_char_m91_slope_coeff
    spearman[0+1*4L].m91_prob    = stats.int_cor_char_m91_slope_prob
    spearman[0+1*4L].pt05_coeff  = stats.int_cor_char_pt05_slope_coeff
    spearman[0+1*4L].pt05_prob   = stats.int_cor_char_pt05_slope_prob

    spearman[1+1*4L].m91_coeff   = stats.int_cor_char_m91_ebv_coeff
    spearman[1+1*4L].m91_prob    = stats.int_cor_char_m91_ebv_prob
    spearman[1+1*4L].pt05_coeff  = stats.int_cor_char_pt05_ebv_coeff
    spearman[1+1*4L].pt05_prob   = stats.int_cor_char_pt05_ebv_prob

    spearman[2+1*4L].m91_coeff   = stats.int_cor_char_m91_incl_coeff
    spearman[2+1*4L].m91_prob    = stats.int_cor_char_m91_incl_prob
    spearman[2+1*4L].pt05_coeff  = stats.int_cor_char_pt05_incl_coeff
    spearman[2+1*4L].pt05_prob   = stats.int_cor_char_pt05_incl_prob

    spearman[3+1*4L].m91_coeff   = stats.int_cor_char_m91_hasb_coeff
    spearman[3+1*4L].m91_prob    = stats.int_cor_char_m91_hasb_prob
    spearman[3+1*4L].pt05_coeff  = stats.int_cor_char_pt05_hasb_coeff
    spearman[3+1*4L].pt05_prob   = stats.int_cor_char_pt05_hasb_prob

; EW    
    
    spearman[0+2*4L].m91_coeff   = stats.int_ew_char_m91_slope_coeff
    spearman[0+2*4L].m91_prob    = stats.int_ew_char_m91_slope_prob
    spearman[0+2*4L].pt05_coeff  = stats.int_ew_char_pt05_slope_coeff
    spearman[0+2*4L].pt05_prob   = stats.int_ew_char_pt05_slope_prob

    spearman[1+2*4L].m91_coeff   = stats.int_ew_char_m91_ebv_coeff
    spearman[1+2*4L].m91_prob    = stats.int_ew_char_m91_ebv_prob
    spearman[1+2*4L].pt05_coeff  = stats.int_ew_char_pt05_ebv_coeff
    spearman[1+2*4L].pt05_prob   = stats.int_ew_char_pt05_ebv_prob

    spearman[2+2*4L].m91_coeff   = stats.int_ew_char_m91_incl_coeff
    spearman[2+2*4L].m91_prob    = stats.int_ew_char_m91_incl_prob
    spearman[2+2*4L].pt05_coeff  = stats.int_ew_char_pt05_incl_coeff
    spearman[2+2*4L].pt05_prob   = stats.int_ew_char_pt05_incl_prob

    spearman[3+2*4L].m91_coeff   = stats.int_ew_char_m91_hasb_coeff
    spearman[3+2*4L].m91_prob    = stats.int_ew_char_m91_hasb_prob
    spearman[3+2*4L].pt05_coeff  = stats.int_ew_char_pt05_hasb_coeff
    spearman[3+2*4L].pt05_prob   = stats.int_ew_char_pt05_hasb_prob

    select = tag_names(spearman)
    format = ['A0','F12.2','F12.2','F12.2','F12.2','F12.2','F12.2','A0',$
      'F12.2','F12.2','F12.2','F12.2','F12.2','F12.2']
    texcenter = ['l','c','c','c','c','c','c','c','c','c','c','c','c','c']

    supercolhead1 = '\colhead{} & \multicolumn{2}{c}{M91\tablenotemark{b}} & \colhead{} & \multicolumn{2}{c}{PT05\tablenotemark{b}} \\ '
    supercolhead2 = '\cline{2-3}\cline{5-6}'
    colhead = '\colhead{'+[ ['Variables','Coefficient','Probability','','Coefficient']+'} & ',['Probability']+'} ']

    tablenotetext = [$
      '{a}{Spearman rank correlation coefficients and the probability that the variables are \emph{uncorrelated}.  '+$
      'We compare the residuals, $\Delta\log({\rm O/H})$, for all three integrated abundance measurements (see \S\ref{sec:intoh} '+$
      'for details), against the slope of the abundance gradient, the dust reddening, $E(B-V)$, the inclination angle, '+$
      '$i$, and the H$\alpha$ surface brightness, $\Sigma({\rm H}\alpha)$.}',$
      '{b}{Oxygen abundances derived using either the \citet[M91]{mcgaugh91} or the \citet[PT05]{pilyugin05} '+$
      'abundance calibrations.  See \S\ref{sec:intoh} for details.}']

    table = html_structure_parse(spearman,tags=select,tagformats=format,tagindices=match,blank='\nodata',/keep_error)
    ntags = n_tags(table)
    tags = tag_names(table)
    niceprint, lindgen(n_tags(table)), tag_names(table)

; write out

    texfile = 'log12oh_spearman.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
       printf, lun, '\pagestyle{empty}'
    endif
;   printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\small}'
;   printf, lun, '\rotate'
    printf, lun, '\tablecaption{Residual Correlation Coefficients\label{table:spearman}\tablenotemark{a}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    printf, lun, supercolhead1
    printf, lun, supercolhead2
    niceprintf, lun, colhead
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0L, n_elements(spearman)-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else begin
             if (i lt n_elements(spearman)-1L) then suffix = ' \\ ' else suffix = ''
          endelse
          if (j eq 0L) or (j eq 3L) then begin
             if (j eq 0L) then line[j] = table[i].(j)+suffix
             if (j eq 3L) then line[j] = suffix
          endif else begin
             if (table[i].(j) gt 0.0) then line[j] = '\phs$'+table[i].(j)+'$'+suffix else $
               line[j] = '$'+table[i].(j)+'$'+suffix
          endelse
       endfor 
       printf, lun, line
       if (i eq 3L) or (i eq 7L) then printf, lun, '\cline{1-6}'
    endfor
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
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
; Table: Residual statistics
; ---------------------------------------------------------------------------    

; replicate a three-element structure, corresponding to the observed,
; corrected, and EW integrated abundances

    resid = replicate({log12oh: '', m91_mean: [-1.0,-1.0], m91_med: -1.0, empty1: '', pt05_mean: [-1.0,-1.0], pt05_med: -1.0},3)
;   resid = replicate({log12oh: '', m91_mean: [-1.0,-1.0], pt05_mean: [-1.0,-1.0], empty1: '', m91_med: -1.0, pt05_med: -1.0},3)

    resid.log12oh = ['Observed','Corrected','EWs']

    resid[0].m91_mean = stats.mean_int_obs_char_m91 & resid[0].m91_med = stats.median_int_obs_char_m91
    resid[1].m91_mean = stats.mean_int_cor_char_m91 & resid[1].m91_med = stats.median_int_cor_char_m91
    resid[2].m91_mean = stats.mean_int_ew_char_m91  & resid[2].m91_med = stats.median_int_ew_char_m91

    resid[0].pt05_mean = stats.mean_int_obs_char_pt05 & resid[0].pt05_med = stats.median_int_obs_char_pt05
    resid[1].pt05_mean = stats.mean_int_cor_char_pt05 & resid[1].pt05_med = stats.median_int_cor_char_pt05
    resid[2].pt05_mean = stats.mean_int_ew_char_pt05  & resid[2].pt05_med = stats.median_int_ew_char_pt05

    select = tag_names(resid)
    format = ['A0','F12.3','F12.3','A0','F12.3','F12.3']
    texcenter = ['l','c','c','c','c','c']

;   supercolhead1 = '\colhead{} & \multicolumn{2}{c}{Mean $\Delta\log({\rm O/H})$} & \colhead{} & '+$
;     '\multicolumn{2}{c}{Median $\Delta\log({\rm O/H})$} \\ '
;   supercolhead2 = '\cline{2-3}\cline{5-6}'
;   colhead = '\colhead{'+[ ['Method\tablenotemark{a}','M91\tablenotemark{b}','PT05\tablenotemark{b}','',$
;     'M91\tablenotemark{b}']+'} & ',['PT05\tablenotemark{b}']+'} ']

    supercolhead1 = '\colhead{} & \multicolumn{2}{c}{M91\tablenotemark{b}} & \colhead{} & '+$
      '\multicolumn{2}{c}{PT05\tablenotemark{b}} \\ '
    supercolhead2 = '\cline{2-3}\cline{5-6}'
    colhead = '\colhead{'+[ ['Method\tablenotemark{a}',$
      '$\langle\Delta\log({\rm O/H})\rangle$','$\Delta\log({\rm O/H})_{\rm median}$','',$
      '$\langle\Delta\log({\rm O/H})\rangle$']+'} & ',['$\Delta\log({\rm O/H})_{\rm median}$']+'} ']

;   colunits = '\colhead{'+[ ['(1)','(2)','(3)','','(4)']+'} & ',['(5)']+'} ']

;   tablecomments = [$
;     'Col. (1) Galaxy name; '+$
;     'Col. (2) Central abundance (at $\rho=0$) based on the derived abundance gradient; '+$
;     'Col. (3) \emph{Characteristic} abundance (at $\rho=0.4\rho_{25}$), based on the derived abundance gradient; '+$
;     'Col. (4) Abundance gradient slope; '+$
;     'Col. (5) Central oxygen abundance (at $\rho=0$); '+$
;     'Col. (6) \emph{Characteristic} oxygen abundance (at $\rho=0.4\rho_{25}$); '+$
;     'Col. (7) Slope of the abundance gradient; '+$
;     'Col. (8) Number of H~{\sc ii} regions; ',$
;     'Col. (9) H~{\sc ii}-region references.']
;   tablecomments = strjoin(tablecomments,' ')

    tablenotetext = [$
      '{a}{Integrated abundances based on the observed emission-line fluxes, the reddening-corrected '+$
      'line-fluxes, or the emission-line equivalent widths (EWs).  See \S\ref{sec:intoh} for details.}',$
      '{b}{Oxygen abundances derived using either the \citet[M91]{mcgaugh91} or the \citet[PT05]{pilyugin05} '+$
      'abundance calibrations.  See \S\ref{sec:intoh} for details.}']

    table = html_structure_parse(resid,tags=select,tagformats=format,tagindices=match,blank='\nodata',/keep_error)
    ntags = n_tags(table)
    tags = tag_names(table)
    niceprint, lindgen(n_tags(table)), tag_names(table)

; write out

    texfile = 'log12oh_residuals.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
       printf, lun, '\pagestyle{empty}'
    endif
;   printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\small}'
;   printf, lun, '\rotate'
    printf, lun, '\tablecaption{Integrated vs. Characteristic Abundance Residuals\label{table:resid}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    printf, lun, supercolhead1
    printf, lun, supercolhead2
    niceprintf, lun, colhead
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0L, n_elements(resid)-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else begin
             if (i lt n_elements(resid)-1L) then suffix = ' \\ ' else suffix = ''
          endelse
          if (j eq 0L) or (j eq 3L) then line[j] = table[i].(j)+suffix else begin ; character tags
             if (table[i].(j)[0] gt 0.0) then prefix = '\phs' else prefix = ''
             if (n_elements(table[i].(j)) eq 1L) then line[j] = prefix+'$'+table[i].(j)+'$'+suffix else $
               line[j] = prefix+'$'+strtrim(table[i].(j)[0],2)+'\pm'+strtrim(table[i].(j)[1],2)+'$'+suffix
          endelse 
       endfor
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
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
; Table 1: sample and integrated abundances
; ---------------------------------------------------------------------------    

    select = [$
      'nice_galaxy','type','m_b','r25','incl','pa','hasb','int_ebv',$
      'int_obs_log12oh_m91','int_cor_log12oh_m91','int_ew_log12oh_m91',$
      'empty1',$
      'int_obs_log12oh_pt05','int_cor_log12oh_pt05','int_ew_log12oh_pt05']

    format = [$
      'A0','A0','F12.1','F12.2','I0','I0','F12.1','F12.2',$
      'F12.2','F12.2','F12.2',$
      'A0',$
      'F12.2','F12.2','F12.2']

    texcenter = [$
      'l','l','c','c','c','c','c','c',$
      'c','c','c',$
      'c',$
      'c','c','c']

    supercolhead1 = '\multicolumn{2}{c}{} & \colhead{$M_{B}$} & \colhead{$\rho_{25}$} & '+$
      '\colhead{$i$} & \colhead{$\theta$} & \colhead{$\log\,[\Sigma({\rm H}\alpha)]$} & \colhead{$E(B-V)$} & '+$
      '\multicolumn{3}{c}{12+log(O/H)$_{\rm M91}$\tablenotemark{a}} & '+$
      '\colhead{} & \multicolumn{3}{c}{12+log(O/H)$_{\rm PT05}$\tablenotemark{a}} \\ '
    supercolhead2 = '\cline{9-11}\cline{13-15}'

    colhead = '\colhead{Galaxy} & \colhead{Type} & \colhead{(mag)} & \colhead{(arcmin)} & '+$
      '\colhead{(deg)} & \colhead{(deg)} & \colhead{(erg~s$^{-1}$~pc$^{-2}$)} & \colhead{(mag)} & '+$
      '\colhead{Observed} & \colhead{Corrected} & \colhead{EWs} & '+$
      '\colhead{} & \colhead{Observed} & \colhead{Corrected} & \colhead{EWs} \\ '
    
    colunits = '\colhead{'+[ ['(1)',$
      '(2)','(3)','(4)',$
      '(5)','(6)','(7)','(8)','(9)','(10)','(11)',$
      '',$
      '(12)','(13)']+'} & ',['(14)']+'} ']

    table = html_structure_parse(gradients,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error)
    ntags = n_tags(table)

; generate the object tags

    table.nice_galaxy = reform('\object['+strtrim(gradients.ned_galaxy,2)+']{'+strtrim(table.nice_galaxy,2)+'}',1,ngalaxy)

    tablecomments = [$
      'Col. (1) Galaxy name;',$
      'Col. (2) Morphological type from \citet{devac91};',$
      'Col. (3) Absolute $B$-band magnitude, corrected for foreground Galactic extinction '+$
      '\citep[$R_{V}=3.1$;][]{schlegel98, odonnell94}, based on the apparent magnitude and '+$
      'distance given by \citet{moustakas06a};',$
      'Col. (4) Radius of the major axis at the $B_{25}$~mag~arcsec$^{-2}$ isophote from \citet{devac91};',$
      'Col. (5) Galaxy inclination angle based on the $K_{s}$-band major-to-minor axis ratio \citep{jarrett00, jarrett03}, '+$
      'as described in \S\ref{sec:sample}, except for NGC~5194, which is from \citet{tully74};',$
      'Col. (6) Galaxy position angle, measured positive from North to East, in the $K_{s}$ band \citep{jarrett00, jarrett03}, '+$
      'as described in \S\ref{sec:sample};',$
      'Col. (7) Surface. '+$
      'Col. (8) Nebular reddening determined from the observed H$\alpha$/H$\beta$ ratio \citep{moustakas06a}, assuming '+$
      'an intrinsic ratio of $2.86$ and the \citet{odonnell94} Milky Way extinction curve;',$
      'Oxygen abundances have been computed using the observed line fluxes [Cols. (9) and (12)], '+$
      'the reddening-corrected line fluxes [Cols. (10) and (13)], and the emission-line equivalent widths (EWs) '+$
      '[Cols. (11) and (14)].']
    tablecomments = strjoin(tablecomments,' ')

    tablenotetext = [$
      '{a}{Oxygen abundances derived using either the \citet[M91]{mcgaugh91} or the \citet[PT05]{pilyugin05} '+$
      'abundance calibrations.  See \S\ref{sec:intoh} for details.}']

; write out the general properties table

    texfile = 'int_log12oh_table.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset2+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{Integrated Oxygen Abundances\label{table:intoh}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    printf, lun, supercolhead1
    printf, lun, supercolhead2
    niceprintf, lun, colhead
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0L, ngalaxy-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else begin
             if (i lt ngalaxy-1L) then suffix = ' \\ ' else suffix = ''
          endelse
          case j of
             11L: line[j] = ''+suffix ; blank tag
             else: begin
                if (n_elements(table[i].(j)) eq 1L) then begin
                   if (strcompress(table[i].(j),/remove) eq '') then $
                     line[j] = '\nodata'+suffix else line[j] = table[i].(j)+suffix
                endif else begin
                   if (strtrim(table[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else begin
                      if gradients[i].int_log12oh_lower_limit and strmatch(select[j],'*log12oh*') then $
                        line[j] = '$>'+strtrim(table[i].(j)[0],2)+'$'+suffix else $
                        line[j] = '$'+strtrim(table[i].(j)[0],2)+'\pm'+strtrim(table[i].(j)[1],2)+'$'+suffix
                   endelse
                endelse
             end
          endcase
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
; Table: HII-region abundance gradients
; ---------------------------------------------------------------------------    

    select = ['nice_galaxy',$
      'hii_m91_log12oh_central','hii_m91_log12oh_char','hii_m91_slope',$
      'empty1',$
      'hii_pt05_log12oh_central','hii_pt05_log12oh_char','hii_pt05_slope',$
      'hii_m91_nhii_used','hii_m91_texrefs']

    format = ['A0',$
      'F12.2','F12.2','F12.2',$
      'A0',$
      'F12.2','F12.2','F12.2',$
      'I0','A0']
      
    texcenter = ['l',$
      'c','c','c',$
      'c',$
      'c','c','c',$
      'c','c']

    supercolhead1 = '\colhead{} & \multicolumn{3}{c}{M91\tablenotemark{a}} & \colhead{} & '+$
      '\multicolumn{3}{c}{PT05\tablenotemark{a}} & \multicolumn{2}{c}{} \\ '
    supercolhead2 = '\cline{2-4}\cline{6-8}'
    supercolhead3 = '\colhead{} & \colhead{12+log(O/H)} & \colhead{12+log(O/H)} & \colhead{Gradient} & \colhead{} & '+$
      '\colhead{12+log(O/H)} & \colhead{12+log(O/H)} & \colhead{Gradient} & \multicolumn{2}{c}{} \\'

    colhead = '\colhead{'+[ ['Galaxy',$
      'at $\rho=0$','at $\rho=0.4\,\rho_{25}$','(dex $\rho_{25}^{-1}$)',$
      '',$
      'at $\rho=0$','at $\rho=0.4\,\rho_{25}$','(dex $\rho_{25}^{-1}$)',$
      'N(H~{\sc ii})']+'} & ',['Refs.']+'} \\ ']

    colunits = '\colhead{'+[ ['(1)',$
      '(2)','(3)','(4)',$
      '',$
      '(5)','(6)','(7)',$
      '(8)']+'} & ',['(9)']+'} ']

    tablecomments = [$
      'Col. (1) Galaxy name; '+$
      'Cols. (2) and (5) Central abundance (at radius $\rho=0$) based on the derived abundance gradient; '+$
      'Cols. (3) and (6) \emph{Characteristic} abundance (at $\rho=0.4\rho_{25}$), based on the derived abundance gradient; '+$
      'Cols. (4) and (7) Slope of the abundance gradient; '+$
      'Col. (8) Number of H~{\sc ii} regions; ',$
      'Col. (9) H~{\sc ii}-region references.']
    tablecomments = strjoin(tablecomments,' ')

    tablenotetext = [$
      '{a}{Oxygen abundances and abundance gradients derived using either the \citet[M91]{mcgaugh91} or the \citet[PT05]{pilyugin05} '+$
      'abundance calibrations.  See \S\ref{sec:intoh} for details.}']

    table = html_structure_parse(gradients,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error)
    ntags = n_tags(table)
    tags = tag_names(table)
    niceprint, lindgen(n_tags(table)), tag_names(table)

; parse the HII-region references and assign unique reference numbers

    m91_refs = strsplit(strjoin(strtrim(gradients.hii_m91_texrefs,2),','),',',/extract)
    nm91_refs = n_elements(m91_refs)

    m91_refnum = lonarr(nm91_refs)-1L
    m91_refcount = 1L
    for i = 0L, nm91_refs-1L do begin
       if (m91_refnum[i] eq -1L) then begin
          m91_match = where(m91_refs[i] eq m91_refs)
          m91_refnum[m91_match] = m91_refcount
          m91_refcount = m91_refcount + 1L
       endif
    endfor

;   niceprint, refs, refnum
    m91_urefnum = uniq(m91_refnum,sort(m91_refnum))
    m91_tablerefs = strjoin('('+string(m91_refnum[m91_urefnum],format='(I0)')+') \citet{'+$
      string(m91_refs[m91_urefnum],format='(A0)')+'}','; ')+'.'

; write out

    texfile = 'hii_log12oh_table.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
       printf, lun, '\pagestyle{empty}'
    endif
;   printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
;   printf, lun, '\rotate'
    printf, lun, '\tablecaption{\ion{H}{2}-Region Oxygen Abundance Gradients\label{table:hiioh}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    printf, lun, supercolhead1
    printf, lun, supercolhead2
    printf, lun, supercolhead3
    niceprintf, lun, colhead
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0L, ngalaxy-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else begin
             if (i lt ngalaxy-1L) then suffix = ' \\ ' else suffix = ''
          endelse
          if (j eq 0L) or (j eq 4L) then begin
             line[j] = table[i].(j)+suffix 
          endif else begin     ; character tags
             case j of
                (ntags-2L): if (table[i].(j) eq 0L) then line[j] = '\nodata'+suffix else $
                  line[j] = table[i].(j)+suffix ; N(HII) tag
                (ntags-1L): begin
                   if (table[i].hii_m91_nhii_used gt 0L) then begin
                      m91_theserefs = strsplit(strtrim(table[i].hii_m91_texrefs,2),',',/extract)
                      m91_these = cmset_op(m91_refs,'and',m91_theserefs,/index)
                      m91_thesenum = m91_refnum[m91_these[uniq(m91_refnum[m91_these],sort(m91_refnum[m91_these]))]]
                      line[j] = strjoin(string(m91_thesenum,format='(I0)'),',')+suffix
                   endif else line[j] = '\nodata'+suffix 
                end
                else: begin
                   if (strtrim(table[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else begin
                      if (median(table.(j)[0]) lt 0.0) and (table[i].(j)[0] gt 0.0) then prefix = '\phs' else prefix = '' ; slope
                      line[j] = prefix+'$'+strtrim(table[i].(j)[0],2)+'\pm'+strtrim(table[i].(j)[1],2)+'$'+suffix
                   endelse
                end
             endcase
          endelse
       endfor
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\tablerefs{'+m91_tablerefs+'}'
    printf, lun, '\end{deluxetable}'
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
