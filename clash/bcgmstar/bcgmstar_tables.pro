function get_nicecluster, sample
    return, repstr(repstr(strupcase(sample.dirname),'ABELL','Abell'),'_',' ')
end

function get_colhead, head, nobreak=nobreak
    if keyword_set(nobreak) then br = '' else br = '\\'
    nhead = n_elements(head)
    colhead = '\colhead{'+[head[0:nhead-2]+'} & ',head[nhead-1]+'} '+br]
    multi = where(strmatch(colhead,'*multicolumn*',/fold))
    if (multi[0] ne -1) then colhead[multi] = repstr(repstr(colhead[multi],'\colhead{',''),'}}','}')
return, colhead
end

pro bcgmstar_tables
; jm14may15siena - build the latex tables for my BCG/Mstar paper 

    filesuffix = ''

    paperpath = bcgmstar_path(/paper)
    sersicpath = bcgmstar_path(/sersic)
    
; by default build emulateapj format tables    
    sample = read_bcgmstar_sample(/zsort)
;   sample = sample[sort(sample.m500)]
    ncl = n_elements(sample) 

    filt = bcgmstar_filterlist(short=short)
    nfilt = n_elements(filt)

; ---------------------------------------------------------------------------    
; Table: principle Sersic fitting results split into four sub-tables

; #########################
; sub-Table 3 - nref, reref, for double Sersic clusters
    sample1 = sample[where(strtrim(sample.shortname,2) eq 'a209' or $
      strtrim(sample.shortname,2) eq 'a2261' or $
      strtrim(sample.shortname,2) eq 'rxj2248')]
    ncl1 = n_elements(sample1)
    
    texfile = paperpath+'table_doublesersic1'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    
    colhead1 = get_colhead([$
      '','','','$r_{\rm e,1}$','$r_{\rm e,2}$','',''])
    colhead2 = get_colhead([$
      'Cluster','$n_{1}$','$n_{2}$','(kpc)','(kpc)',$
      '$\chi^{2}_{\nu}$\tablenotemark{b}','$\nu$\tablenotemark{b}'],/nobreak)
    texcenter = ['l',replicate('c',6)]

    tablenotetext = [$
      '{a}{Together with Table~\ref{table:doublesersic2}, this table lists the best-fitting '+$
      'parameters for the BCGs fitted with the double-\sersic{} parametric model '+$
      'described in \S\ref{sec:sersic}.}',$
      '{b}{See Table~\ref{table:sersic1}.}']
    
; write it out    
    openw, lun, texfile, /get_lun
;   if (keyword_set(preprint) eq 0) then printf, lun, '\LongTables'
;   printf, lun, '\clearpage'
;   printf, lun, '\setlength{\tabcolsep}{0.02in}'
;   printf, lun, '\LongTables'
;   printf, lun, '\begin{turnpage}'
;   printf, lun, '\begin{landscape}'
;   printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\begin{deluxetable*}{'+strjoin(texcenter)+'}'
;   printf, lun, '\setlength{\tabcolsep}{0.005in}'
    printf, lun, '\tablecaption{Double-\sersic{} Surface-Brightness Profile Fitting '+$
      'Results\tablenotemark{a}\label{table:doublesersic1}}'
;   printf, lun, '\tableenum{\ref{fig:allsb1} \emph{continued}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    printf, lun, '}'
    printf, lun, '\startdata'
    for ic = 0, ncl1-1 do begin
       cluster = strtrim(sample1[ic].shortname,2)
       nicecluster = get_nicecluster(sample1[ic])
       if (ic eq ncl1-1) then suffix = '' else suffix = ' \\'
       sersic = read_bcgmstar_sersic(cluster,results=results)
       
       n1 = '$'+strtrim(string(results.sersic2_n1_ref,format='(F5.2)'),2)+'\pm'+$
         strtrim(string(results.sersic2_n1_ref_err,format='(F5.2)'),2)+'$'
       n2 = '$'+strtrim(string(results.sersic2_n2_ref,format='(F5.2)'),2)+'\pm'+$
         strtrim(string(results.sersic2_n2_ref_err,format='(F5.2)'),2)+'$'
       re1 = '$'+strtrim(string(results.sersic2_re1_ref,format='(F12.1)'),2)+'\pm'+$
         strtrim(string(results.sersic2_re1_ref_err>0.1,format='(F12.1)'),2)+'$'
       re2 = '$'+strtrim(string(results.sersic2_re2_ref,format='(F12.1)'),2)+'\pm'+$
         strtrim(string(results.sersic2_re2_ref_err>0.1,format='(F12.1)'),2)+'$'

       line = nicecluster+' & '+n1+' & '+n2+' & '+re1+' & '+re2+' & '+$
         '$'+strtrim(string(results.sersic2_chi2/results.sersic2_dof,format='(F5.2)'),2)+'$ & '+$
         '$'+strtrim(string(results.sersic2_dof,format='(I0)'),2)+'$'
       printf, lun, line+suffix
    endfor 
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable*}'
;   printf, lun, '\clearpage'
;   printf, lun, '\end{landscape}'
;   printf, lun, '\end{turnpage}'
    free_lun, lun

; #########################
; sub-Table 4 - double-Sersic surface brightnesses
    sample1 = sample[where(strtrim(sample.shortname,2) eq 'a209' or $
      strtrim(sample.shortname,2) eq 'a2261' or $
      strtrim(sample.shortname,2) eq 'rxj2248')]
    ncl1 = n_elements(sample1)

    texfile = paperpath+'table_doublesersic2'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    
    colhead1 = get_colhead(['','$\mu_{\rm e,1}$','$\mu_{\rm e,2}$','','$\mu_{\rm e,1}$','$\mu_{\rm e,2}$',$
      '','$\mu_{\rm e,1}$','$\mu_{\rm e,2}$'])
    colhead2 = get_colhead(['Band','\multicolumn{2}{c}{(mag arcsec$^{-2}$)}','',$
      '\multicolumn{2}{c}{(mag arcsec$^{-2}$)}','','\multicolumn{2}{c}{(mag arcsec$^{-2}$)}'],/nobreak)
    colhead3 = get_colhead(['','\multicolumn{2}{c}{Abell~209}','',$
      '\multicolumn{2}{c}{Abell~2261}','','\multicolumn{2}{c}{RXJ2248}'])
    texcenter = ['l','c','c','c','c','c','c','c','c','c']

    tablenotetext = [$
      '{a}{Together with Table~\ref{table:doublesersic1}, this table lists the best-fitting '+$
      'parameters for the BCGs fitted with the double-\sersic{} parametric model '+$
      'described in \S\ref{sec:sersic}.}']

; write it out    
    openw, lun, texfile, /get_lun
;   if (keyword_set(preprint) eq 0) then printf, lun, '\LongTables'
;   printf, lun, '\clearpage'
;   printf, lun, '\setlength{\tabcolsep}{0.02in}'
;   printf, lun, '\LongTables'
;   printf, lun, '\begin{turnpage}'
;   printf, lun, '\begin{landscape}'
;   printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\begin{deluxetable*}{'+strjoin(texcenter)+'}'
;   printf, lun, '\setlength{\tabcolsep}{0.005in}'
    printf, lun, '\tablecaption{Double-\sersic{} Surface-Brightness Profile Fitting '+$
      'Results\tablenotemark{a}\label{table:doublesersic2}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    printf, lun, '}'
    printf, lun, '\startdata'

; get the bandpasses from A209    
    junk = read_bcgmstar_sersic('a209')
    short1 = strupcase(strtrim(junk.band,2))
    nfilt1 = n_elements(short1)

;   printf, lun, '\cline{2-3}' & printf, lun, '\cline{5-6}' & printf, lun, '\cline{8-9}'
    niceprintf, lun, colhead3
    printf, lun, '\cline{2-3}' & printf, lun, '\cline{5-6}' & printf, lun, '\cline{8-9}'
    
    for ff = nfilt1-1, 0, -1 do begin
       line = short1[ff]+' & '
       
       for ic = 0, ncl1-1 do begin
          cluster = strtrim(sample1[ic].shortname,2)
          nicecluster = get_nicecluster(sample1[ic])
;         splog, cluster
;         printf, lun, '\cline{1-3}'
;         printf, lun, '\cutinhead{'+nicecluster+'}'
;         printf, lun, '\cline{1-3}'

          if (ic eq ncl1-1) then begin
             if ff eq 0 then suffix = '' else suffix = ' \\'
          endif else begin
             suffix = ' & & '
          endelse

          sersic = read_bcgmstar_sersic(cluster,results=results,band=short1[ff])
          if size(sersic,/type) ne 8 then begin
             line = line+' \nodata & \nodata & '+suffix
          endif else begin
             line = line+'$'+strtrim(string(sersic.sersic2_all_sbe1,format='(F12.2)'),2)+'\pm'+$
               strtrim(string(sersic.sersic2_all_sbe1_err,format='(F12.2)'),2)+'$'+' & $'+$
               strtrim(string(sersic.sersic2_all_sbe2,format='(F12.2)'),2)+'\pm'+$
               strtrim(string(sersic.sersic2_all_sbe2_err,format='(F12.2)'),2)+'$'+suffix
          endelse
       endfor
       printf, lun, line
       print, line
    endfor 

    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable*}'
;   printf, lun, '\clearpage'
;   printf, lun, '\end{landscape}'
;   printf, lun, '\end{turnpage}'
    free_lun, lun

; #########################
; sub-Table 1 - nref, reref, alpha, beta, etc. for single Sersic
    sample1 = sample[where(strtrim(sample.shortname,2) ne 'a209' and $
      strtrim(sample.shortname,2) ne 'a2261' and $
      strtrim(sample.shortname,2) ne 'rxj2248')]
    ncl1 = n_elements(sample1)

    texfile = paperpath+'table_sersic1'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    
    colhead1 = get_colhead(['','','$r_{\rm e}$','',''])
    colhead2 = get_colhead([$
      'Cluster','$n$','(kpc)',$
      '$\chi^{2}_{\nu}$\tablenotemark{b}','$\nu$\tablenotemark{b}'],/nobreak)
    texcenter = ['l',replicate('c',4)]
;   colhead1 = get_colhead([$
;     '','','$r_{\rm e,ref}$','','','',''])
;   colhead2 = get_colhead([$
;     'Cluster','$n_{\rm ref}$','(kpc)','$\alpha$','$\beta$',$
;     '$\chi^{2}_{\nu}$\tablenotemark{b}','$\nu$\tablenotemark{b}'],/nobreak)
;   texcenter = ['l',replicate('c',6)]

    tablenotetext = [$
      '{a}{Together with Table~\ref{table:sersic2}, this table lists the best-fitting '+$
      'parameters for the BCGs fitted with the single-\sersic{} parametric model '+$
      'described in \S\ref{sec:sersic}.}',$
      '{b}{$\chi^{2}_{\nu}$ is the reduced $\chi^{2}$ value and '+$
      '$\nu$ is the number of degrees of freedom (number of photometric data points fitted '+$
      'minus the number of free parameters).}']

;     '{c}{The parameters $\alpha$ and $\beta$ for these BCGs are not well-constrained, '+$
;     'so the wavelength-dependent parametric model is not used (i.e., $\alpha=0$ and $\beta=0$ '+$
;     'for these BCGs).}']

; write it out    
    openw, lun, texfile, /get_lun
;   if (keyword_set(preprint) eq 0) then printf, lun, '\LongTables'
;   printf, lun, '\clearpage'
;   printf, lun, '\setlength{\tabcolsep}{0.02in}'
;   printf, lun, '\LongTables'
;   printf, lun, '\begin{turnpage}'
;   printf, lun, '\begin{landscape}'
;   printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\begin{deluxetable*}{'+strjoin(texcenter)+'}'
;   printf, lun, '\setlength{\tabcolsep}{0.005in}'
    printf, lun, '\tablecaption{Single-\sersic{} Surface-Brightness Profile Fitting '+$
      'Results\tablenotemark{a}\label{table:sersic1}}'
;   printf, lun, '\tableenum{\ref{fig:allsb1} \emph{continued}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    printf, lun, '}'
    printf, lun, '\startdata'
    for ic = 0, ncl1-1 do begin
       cluster = strtrim(sample1[ic].shortname,2)
       nicecluster = get_nicecluster(sample1[ic])
       if (ic eq ncl1-1) then suffix = '' else suffix = ' \\'
       sersic = read_bcgmstar_sersic(cluster,results=results)
       
;       format_flux_error, results.sersic_n_ref, results.sersic_n_ref_err, nn, nnerr
;       format_flux_error, results.sersic_re_ref, results.sersic_re_ref_err, re, reerr
;       format_flux_error, abs(results.sersic_alpha1), results.sersic_alpha1_err, alpha, alphaerr
;       format_flux_error, abs(results.sersic_beta1), results.sersic_beta1_err, beta, betaerr
;       if results.sersic_alpha1 lt 0 then alpha = '-'+alpha
;       if results.sersic_beta1 lt 0 then beta = '-'+beta

       nn = '$'+strtrim(string(results.sersic_n_ref,format='(F5.2)'),2)+'\pm'+$
         strtrim(string(results.sersic_n_ref_err,format='(F5.2)'),2)+'$'
       re = '$'+strtrim(string(results.sersic_re_ref,format='(F12.1)'),2)+'\pm'+$
         strtrim(string(results.sersic_re_ref_err>0.1,format='(F12.1)'),2)+'$'

       if results.sersic_alpha1 gt 0 then pre = '\phs' else pre = ''
       alpha = pre+'$'+strtrim(string(results.sersic_alpha1,format='(F12.3)'),2)+'\pm'+$
         strtrim(string(results.sersic_alpha1_err,format='(F12.3)'),2)+'$'
       if results.sersic_beta1 gt 0 then pre = '\phs' else pre = ''
       beta = pre+'$'+strtrim(string(results.sersic_beta1,format='(F12.3)'),2)+'\pm'+$
         strtrim(string(results.sersic_beta1_err,format='(F12.3)'),2)+'$'

;       if results.sersic_alpha1_err eq 0.0 then begin
;          alpha = '\nodata'
;          note = '\tablenotemark{c}'
;       endif else note = ''
;       if results.sersic_beta1_err eq 0.0 then beta = '\nodata'
       
       line = nicecluster+' & '+nn+' & '+re+' & $'+$
         strtrim(string(results.sersic_chi2/results.sersic_dof,format='(F5.2)'),2)+'$ & '+$
         '$'+strtrim(string(results.sersic_dof,format='(I0)'),2)+'$'
;      line = nicecluster+note+' & '+nn+' & '+re+' & '+alpha+' & '+beta+' & '+$
;        '$'+strtrim(string(results.sersic_chi2/results.sersic_dof,format='(F5.2)'),2)+'$ & '+$
;        '$'+strtrim(string(results.sersic_dof,format='(I0)'),2)+'$'
       printf, lun, line+suffix
    endfor 
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable*}'
;   printf, lun, '\clearpage'
;   printf, lun, '\end{landscape}'
;   printf, lun, '\end{turnpage}'
    free_lun, lun

; #########################
; sub-Table 2 - single-Sersic surface brightnesses
    sample1 = sample[where(strtrim(sample.shortname,2) ne 'a209' and $
      strtrim(sample.shortname,2) ne 'a2261' and $
      strtrim(sample.shortname,2) ne 'rxj2248')]
    ncl1 = n_elements(sample1)

    texfile = paperpath+'table_sersic2'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    
    colhead1 = get_colhead([$
      '','\multicolumn{7}{c}{$\mu_{\rm e}$}'])
    colhead2 = get_colhead([$
      'Cluster','\multicolumn{7}{c}{(mag arcsec$^{-2}$)}'],/nobreak)
    texcenter = ['l',replicate('c',7)]
;    colhead1 = get_colhead([$
;      '','\multicolumn{14}{c}{$\mu_{\rm e}$ (mag arcsec$^{-2}$)}'])
;    colhead2 = get_colhead([$
;      'Cluster',$
;      'F390W','F435W','F475W','F555W','F606W','F625W','F775W','F814W','F850LP',$
;      'F105W','F110W','F125W','F140W','F160W'],/nobreak)
;    texcenter = ['l',replicate('c',14)]

    tablenotetext = [$
      '{a}{Together with Table~\ref{table:sersic1}, this table lists the best-fitting '+$
      'parameters for the BCGs fitted with the single-\sersic{} parametric model '+$
      'described in \S\ref{sec:sersic}.}']

; write it out    
    openw, lun, texfile, /get_lun
;   if (keyword_set(preprint) eq 0) then printf, lun, '\LongTables'
;   printf, lun, '\clearpage'
;   printf, lun, '\setlength{\tabcolsep}{0.02in}'
;   printf, lun, '\LongTables'
;   printf, lun, '\begin{turnpage}'
;   printf, lun, '\begin{landscape}'
;   printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\begin{deluxetable*}{'+strjoin(texcenter)+'}'
;   printf, lun, '\setlength{\tabcolsep}{0.005in}'
    printf, lun, '\tablecaption{Single-\sersic{} Surface-Brightness Profile Fitting '+$
      'Results\tablenotemark{a}\label{table:sersic2}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    printf, lun, '}'
    printf, lun, '\startdata'

; first seven filters    
    printf, lun, '\cline{1-8}'
    printf, lun, get_colhead(['','F390W','F435W','F475W','F555W','F606W','F625W','F775W'])
    printf, lun, '\cline{1-8}'
    for ic = 0, ncl1-1 do begin
       cluster = strtrim(sample1[ic].shortname,2)
       nicecluster = get_nicecluster(sample1[ic])
       sersic = read_bcgmstar_sersic(cluster,results=results)

       line = nicecluster+' & '
       for ff = 0, 6 do begin
;      for ff = 0, nfilt-1 do begin
          if ff eq 6 then suffix = ' \\' else suffix = ' & '
          this = where(short[ff] eq strtrim(sersic.band,2))
          if this[0] eq -1 then begin
             line = line+' \nodata & '
; put nodata
          endif else begin
;            format_flux_error, sersic[this].sersic_all_sbe, sersic[this].sersic_all_sbe_err, mu, muerr
             line = line+'$'+strtrim(string(sersic[this].sersic_all_sbe,format='(F12.2)'),2)+'\pm'+$
               strtrim(string(sersic[this].sersic_all_sbe_err,format='(F12.2)'),2)+'$'+suffix
          endelse
       endfor 
       printf, lun, line
    endfor 

; next seven filters    
    printf, lun, '\cline{1-8}'
    printf, lun, get_colhead(['','F814W','F850LP','F105W','F110W','F125W','F140W','F160W'])
    printf, lun, '\cline{1-8}'
    for ic = 0, ncl1-1 do begin
       cluster = strtrim(sample1[ic].shortname,2)
       nicecluster = get_nicecluster(sample1[ic])

       sersic = read_bcgmstar_sersic(cluster,results=results)

       line = nicecluster+' & '
       for ff = 7, nfilt-1 do begin
;      for ff = 0, nfilt-1 do begin
          if ff eq nfilt-1 then begin
             if (ic eq ncl1-1) then suffix = '' else suffix = ' \\'
          endif else suffix = ' & '
          this = where(short[ff] eq strtrim(sersic.band,2))
          if this[0] eq -1 then begin
             line = line+' \nodata & '
; put nodata
          endif else begin
;            format_flux_error, sersic[this].sersic_all_sbe, sersic[this].sersic_all_sbe_err, mu, muerr
             line = line+'$'+strtrim(string(sersic[this].sersic_all_sbe,format='(F12.2)'),2)+'\pm'+$
               strtrim(string(sersic[this].sersic_all_sbe_err,format='(F12.2)'),2)+'$'+suffix
          endelse
       endfor 
       printf, lun, line
    endfor 
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable*}'
;   printf, lun, '\clearpage'
;   printf, lun, '\end{landscape}'
;   printf, lun, '\end{turnpage}'
    free_lun, lun

stop    

; ---------------------------------------------------------------------------    
; Table 1: table summarizing the sample
    texfile = paperpath+'table_sample'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    
; gather the information we need
    info = replicate({cluster: '', z: '', ra: '', dec: '', m500: ''},ncl)
    info.cluster = get_nicecluster(sample)
    info.z = '$'+strtrim(string(sample.z,format='(F12.3)'),2)+'$'
    info.ra = sample.ra
    info.dec = repstr(sample.dec,'-','$-$')
    info.m500 = '$'+strtrim(string(10D^sample.m500/1D15,format='(F12.2)'),2)+'\pm'+$
      strtrim(string(sample.m500_err*10D^sample.m500*alog(10)/1D15,format='(F12.2)'),2)+'$'
    struct_print, info
    ntags = n_tags(info)

    colhead1 = get_colhead(['','','R.~A.\tablenotemark{a}','Dec.\tablenotemark{a}',$
      '$M_{\rm 500}$\tablenotemark{b}'])
    colhead2 = get_colhead(['Cluster','Redshift','(J2000)','(J2000)',$
      '$(10^{15}~M_{\odot}/h)$'],/nobreak) 
    texcenter = ['l','c','c','c','c']

    tablenotetext = [$
;     '{a}{Ancillary information for our sample of galaxy clusters.}',$
      '{a}{Cluster coordinates based on the position of the BCG.}',$
      '{b}{Total cluster mass within the $r_{500}$ radius taken from Donahue et~al. (2014, in prep.) '+$
      'for MACS2129 and Abell~1423 and from \citet{merten14a} for the other clusters.  The '+$
      '$r_{500}$ radius is the radius at which the average density of the halo is $500$ times '+$
      'the critical density of the Universe at that redshift.}']

; write it out    
    openw, lun, texfile, /get_lun
;   if (keyword_set(preprint) eq 0) then printf, lun, '\LongTables'
;   printf, lun, '\begin{landscape}'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;   printf, lun, '\tabletypesize{\small}'
    printf, lun, '\tablecaption{Cluster Sample Properties\label{table:sample}}'  
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    printf, lun, '}'
    printf, lun, '\startdata'
    for ii = 0, ncl-1 do begin
       for jj = 0, ntags-1 do begin
          if (jj eq ntags-1) then if (ii eq ncl-1) then suffix = '' else $
            suffix = ' \\' else suffix = ' & '
          printf, lun, info[ii].(jj)+suffix
       endfor
    endfor
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
;   printf, lun, '\vspace*{0.1in}'
;   printf, lun, '\clearpage'
;   printf, lun, '\end{landscape}'
    free_lun, lun

stop    
    
; ---------------------------------------------------------------------------    
; Table: prior parameters
    texfile = paperpath+'table_priors'+filesuffix+'.tex'
    splog, 'Writing '+texfile

    ngrid = 4
    info = replicate({grid: '', age: '', tau: '', Zmetal: ''},ngrid)
    info.grid = string(lindgen(ngrid)+1,format='(I2.2)')
    info.age = '$1.0-11.2$'
    info.tau = ['$0-2$','$0-2$','SSP','SSP']
    info.Zmetal = ['$0.01-0.03$','$0.03$','$0.03$','$0.01-0.03$']
    struct_print, info
    ntags = n_tags(info)

    colhead1 = get_colhead(['','Age','$\tau$',''])
    colhead2 = get_colhead(['Grid','(Gyr)','(Gyr)','Metallicity'],/nobreak)
    texcenter = ['l','c','c','c']

    tablenotetext = [$
      '{a}{Adopted prior parameters.  All grids are constructed using the FSPS \citep[v2.4;][]{'+$
      'conroy09a, conroy10b} stellar population synthesis models an}']

; write it out    
    openw, lun, texfile, /get_lun
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tablecaption{Prior Parameters\tablenotemark{a}\label{table:priors}}'  
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    printf, lun, '}'
    printf, lun, '\startdata'
    for ii = 0, ngrid-1 do begin
       for jj = 0, ntags-1 do begin
          if (jj eq ntags-1) then if (ii eq ngrid-1) then suffix = '' else $
            suffix = ' \\' else suffix = ' & '
          printf, lun, info[ii].(jj)+suffix
       endfor
    endfor
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
;   niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    free_lun, lun

stop    

; ---------------------------------------------------------------------------    
; Table: photometry
    texfile = paperpath+'table_phot'+filesuffix+'.tex'
    splog, 'Writing '+texfile

    factor = 10D^(0.4*23.9) ; conversion to microJy
    
    colhead1 = get_colhead([$
      'Radius',$
      'F390W',$
      'F435W',$
      'F475W',$
      'F555W',$
      'F606W',$
      'F625W',$
      'F775W',$
      'F814W',$
      'F850LP',$
      'F105W',$
      'F110W',$
      'F125W',$
      'F140W',$
      'F160W'])
    colhead2 = get_colhead([$
      '(kpc)',$
      '$(\mu$Jy$)$',$
      '$(\mu$Jy$)$',$
      '$(\mu$Jy$)$',$
      '$(\mu$Jy$)$',$
      '$(\mu$Jy$)$',$
      '$(\mu$Jy$)$',$
      '$(\mu$Jy$)$',$
      '$(\mu$Jy$)$',$
      '$(\mu$Jy$)$',$
      '$(\mu$Jy$)$',$
      '$(\mu$Jy$)$',$
      '$(\mu$Jy$)$',$
      '$(\mu$Jy$)$',$
      '$(\mu$Jy$)$'],/nobreak)
    texcenter = ['c',replicate('c',14)]

;    tablenotetext = [$
;;     '{a}{Ancillary information for our sample of galaxy clusters.}',$
;      '{a}{Cluster coordinates based on the position of the BCG.}',$
;      '{b}{Cluster virial mass adopted from \citet{merten14a}.}']

; write it out    
    openw, lun, texfile, /get_lun
;   if (keyword_set(preprint) eq 0) then printf, lun, '\LongTables'
;   printf, lun, '\clearpage'
    printf, lun, '\setlength{\tabcolsep}{0.02in}'
    printf, lun, '\LongTables'
;   printf, lun, '\begin{turnpage}'
;   printf, lun, '\begin{landscape}'
;   printf, lun, '\oddsidemargin=-2in'
    printf, lun, '\tablefontsize{\tiny}'
    printf, lun, '\begin{deluxetable*}{'+strjoin(texcenter)+'}'
    printf, lun, '\tablecaption{Elliptical Aperture Photometry\label{table:phot}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    printf, lun, '}'
    printf, lun, '\startdata'

    for ic = 0, ncl-1 do begin
       suffix = '\\'
       cluster = strtrim(sample[ic].shortname,2)
       nicecluster = get_nicecluster(sample[ic])
       printf, lun, '\cutinhead{'+nicecluster+'}'

       phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
       rad = phot[0].photradius_kpc
       nrad = n_elements(rad)

       for rr = 0, nrad-1 do begin
          line = '$'+strtrim(string(phot[0].photradius_kpc_in[rr],format='(G12.3)'),2)+'$--$'+$
            strtrim(string(phot[0].photradius_kpc_out[rr],format='(G12.3)'),2)+'$'
          for ff = 0, nfilt-1 do begin
             this = where(short[ff] eq strtrim(phot.band,2))
             if this[0] eq -1 then begin
                line = line+' & \nodata'
             endif else begin
                if phot[this].photradius_kpc[rr] le phot[this].rmax_kpc then begin
                   flux = phot[this].maggies[rr]*factor
                   ferr = factor/sqrt(phot[this].ivarmaggies[rr])
                   format_flux_error, flux, ferr, newflux, newferr
                   line = line+' & $'+newflux+'\pm'+newferr+'$'
                endif else line = line+' & \nodata '
             endelse
          endfor                ; close filter loop
          printf, lun, line+suffix
       endfor                   ; cloes radius loop
; add the integrated photometry
       line = 'Total'
       for ff = 0, nfilt-1 do begin
          this = where(short[ff] eq strtrim(phot.band,2))
          if this[0] eq -1 then begin
             line = line+' & \nodata'
          endif else begin
             flux = phot[this].maggies_int*factor
             ferr = factor/sqrt(phot[this].ivarmaggies_int)
             format_flux_error, flux, ferr, newflux, newferr
             line = line+' & $'+newflux+'\pm'+newferr+'$'
          endelse
       endfor                   ; close filter loop
       if ic eq ncl-1 then suffix = ''
       printf, lun, line+suffix
    endfor ; close cluster loop
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
;   niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable*}'
;   printf, lun, '\end{landscape}'
;   printf, lun, '\end{turnpage}'
;   printf, lun, '\clearpage'
    free_lun, lun


stop    

return
end



;; ---------------------------------------------------------------------------    
;; Table: Sersic fitting results
;    texfile = paperpath+'table_sersic'+filesuffix+'.tex'
;    splog, 'Writing '+texfile
;
;; build the column headings
;;   colhead1 = get_colhead([$
;;     '\multicolumn{4}{c}{F390W}','',$
;;     '\multicolumn{4}{c}{F435W}','',$
;;     '\multicolumn{4}{c}{F475W}','',$
;;     '\multicolumn{4}{c}{F555W}','',$
;;     '\multicolumn{4}{c}{F606W}','',$
;;     '\multicolumn{4}{c}{F625W}','',$
;;     '\multicolumn{4}{c}{F775W}','',$
;;     '\multicolumn{4}{c}{F814W}','',$
;;     '\multicolumn{4}{c}{F850LP}','',$
;;     '\multicolumn{4}{c}{F105W}','',$
;;     '\multicolumn{4}{c}{F110W}','',$
;;     '\multicolumn{4}{c}{F125W}','',$
;;     '\multicolumn{4}{c}{F140W}','',$
;;     '\multicolumn{4}{c}{F160W}'])
;
;;    colhead1 = get_colhead([$
;;      '','$r_{\rm min}$','$r_{\rm max}$',$
;;      '$\mu_{e}$','$r_{e}$','',$
;;      '$\mu_{e,Ser}$','$\mu_{e,dV}$','$r_{e,Ser}$','$r_{e,dV}$',''])
;;    colhead2 = get_colhead([$
;;      'Band','(kpc)','(kpc)',$
;;      '(mag~arcsec$^{-2}$)','(kpc)','$n$',$
;;      '(mag~arcsec$^{-2}$)','(mag~arcsec$^{-2}$)','(kpc)','(kpc)','$n_{Ser}$'],/nobreak)
;
;;   colhead1 = get_colhead([$
;;     '','$r_{\rm min}$','$r_{\rm max}$',$
;;     '$\mu_{e}$','$r_{e}$','','','',$
;;     '$\mu_{e,Ser}$','$\mu_{e,dV}$','$r_{e,Ser}$','$r_{e,dV}$','','',''])
;;   colhead2 = get_colhead([$
;;     'Band','(kpc)','(kpc)',$
;;     '(mag~arcsec$^{-2}$)','(kpc)','$n$','$\chi^{2}_{\nu}$','$\nu$',$
;;     '(mag~arcsec$^{-2}$)','(mag~arcsec$^{-2}$)','(kpc)','(kpc)','$n_{Ser}$','$\chi^{2}_{\nu}$','$\nu$'],/nobreak)
;    
;    colhead1 = get_colhead([$
;      '','$r_{\rm min}$','$r_{\rm max}$',$
;;     '$\mu_{e}$','$r_{e}$','','','',$
;      '$\mu_{e,Ser}$','$r_{e,Ser}$','','$\mu_{e,dV}$','$r_{e,dV}$','',''])
;    colhead2 = get_colhead([$
;      'Band','(kpc)','(kpc)',$
;;     '(mag~arcsec$^{-2}$)','(kpc)','$n$','$\chi^{2}_{\nu}$','$\nu$',$
;      '(mag~arcsec$^{-2}$)','(kpc)','$n_{Ser}$','(mag~arcsec$^{-2}$)','(kpc)','$\chi^{2}_{\nu}$','$\nu$'],/nobreak)
;    
;    texcenter = ['l',replicate('c',9)]
;;   texcenter = replicate('c',(4+1)*7)
;;   texcenter = replicate('c',(4+1)*14)
;
;;    tablenotetext = [$
;;;     '{a}{Ancillary information for our sample of galaxy clusters.}',$
;;      '{a}{Cluster coordinates based on the position of the BCG.}',$
;;      '{b}{Cluster virial mass adopted from \citet{merten14a}.}']
;
;; write it out    
;    openw, lun, texfile, /get_lun
;;   if (keyword_set(preprint) eq 0) then printf, lun, '\LongTables'
;;   printf, lun, '\clearpage'
;    printf, lun, '\LongTables'
;;   printf, lun, '\begin{turnpage}'
;;   printf, lun, '\begin{landscape}'
;    printf, lun, '\begin{deluxetable*}{'+strjoin(texcenter)+'}'
;    printf, lun, '\setlength{\tabcolsep}{0.005in}'
;    printf, lun, '\tablecaption{Surface Brightness Profile Fitting Results\label{table:sersic}}'
;    printf, lun, '\tablewidth{0pt}'
;    printf, lun, '\tabletypesize{\tiny}'
;    printf, lun, '\tablehead{'
;    niceprintf, lun, colhead1
;    niceprintf, lun, colhead2
;    printf, lun, '}'
;    printf, lun, '\startdata'
;    for ic = 0, ncl-1 do begin
;;   for ic = 0, ncl-1 do begin
;       cluster = strtrim(sample[ic].shortname,2)
;       nicecluster = get_nicecluster(sample[ic])
;       printf, lun, '\cutinhead{'+nicecluster+'}'
;
;       phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
;;      for ff = nfilt-2, nfilt-1 do begin
;;      for ff = 7, nfilt-1 do begin
;       for ff = 0, nfilt-1 do begin
;          this = where(short[ff] eq strtrim(phot.band,2))
;          if this[0] eq -1 then begin
;; put nodata
;          endif else begin
;             mu_ser = -2.5*alog10(phot[this].sersic2_sbe1)
;             mu_ser_err = 2.5*phot[this].sersic2_sbe1_err/phot[this].sersic2_sbe1/alog(10)
;             mu_dV = -2.5*alog10(phot[this].sersic2_sbe2)
;             mu_dV_err = 2.5*phot[this].sersic2_sbe2_err/phot[this].sersic2_sbe2/alog(10)
;
;             if ff eq nfilt-1 then begin
;                if (ic eq ncl-1) then suffix = '' else suffix = ' \\'
;             endif else suffix = ' \\'
;;            if ff eq nfilt-1 then if (ic eq ncl-1) then suffix = '' else $
;;              suffix = ' \\' else suffix = ' & '
;             printf, lun, strupcase(short[ff])+' & '+$
;               '$'+strtrim(string(phot[this].rmin_kpc, format='(F5.2)'),2)+'$ & '+$
;               '$'+strtrim(string(phot[this].rmax_kpc, format='(F5.2)'),2)+'$ & '+$
;;; single Sersic
;;               '$'+strtrim(string(phot[this].sersic_sbe, format='(F5.2)'),2)+'\pm'+$
;;               strtrim(string(phot[this].sersic_sbe_err, format='(F5.2)'),2)+'$ & '+$
;;               '$'+strtrim(string(phot[this].sersic_re,  format='(F5.2)'),2)+'\pm'+$
;;               strtrim(string(phot[this].sersic_re_err,  format='(F5.2)'),2)+'$ & '+$
;;               '$'+strtrim(string(phot[this].sersic_n,   format='(F5.2)'),2)+'\pm'+$
;;               strtrim(string(phot[this].sersic_n_err,   format='(F5.2)'),2)+'$ & '+$
;;;              '$'+strtrim(string((phot[this].sersic_chi2/phot[this].sersic_dof)>0.1,format='(F5.2)'),2)+'$ & '+$
;;;              '$'+strtrim(string(phot[this].sersic_dof,     format='(I0)'),2)+'$ & '+$
;; double Sersic
;               '$'+strtrim(string(mu_ser, format='(F5.2)'),2)+'\pm'+$
;               strtrim(string(mu_ser_err, format='(F5.2)'),2)+'$ & '+$
;               '$'+strtrim(string(phot[this].sersic2_re1,  format='(F5.2)'),2)+'\pm'+$
;               strtrim(string(phot[this].sersic2_re1_err,  format='(F5.2)'),2)+'$ & '+$
;               '$'+strtrim(string(phot[this].sersic2_n1,  format='(F5.2)'),2)+'\pm'+$
;               strtrim(string(phot[this].sersic2_n1_err,  format='(F5.2)'),2)+'$ & '+$
;
;               '$'+strtrim(string(mu_dV,  format='(F5.2)'),2)+'\pm'+$
;               strtrim(string(mu_dV_err,  format='(F5.2)'),2)+'$ & '+$
;               '$'+strtrim(string(phot[this].sersic2_re2,  format='(F5.2)'),2)+'\pm'+$
;               strtrim(string(phot[this].sersic2_re2_err,  format='(F5.2)'),2)+'$ & '+$
;
;               '$'+strtrim(string((phot[this].sersic2_chi2/phot[this].sersic2_dof)>0.1,format='(F5.2)'),2)+'$ & '+$
;               '$'+strtrim(string(phot[this].sersic2_dof,     format='(I0)'),2)+'$'+suffix
;          endelse
;       endfor
;    endfor 
;    
;    printf, lun, '\enddata'
;;   printf, lun, '\tablecomments{'+tablecomments+'}'
;;   niceprintf, lun, '\tablenotetext'+tablenotetext
;    printf, lun, '\end{deluxetable*}'
;;   printf, lun, '\end{landscape}'
;;   printf, lun, '\end{turnpage}'
;;   printf, lun, '\clearpage'
;    free_lun, lun
;
;
;stop    


;; ---------------------------------------------------------------------------    
;; Table: Sersic fitting results split into two tables
;    texfile = paperpath+'table_sersic'+filesuffix+'.tex'
;    splog, 'Writing '+texfile
;
;; build the column headings
;    colhead1 = get_colhead([$
;;     '$r_{\rm min}$','$r_{\rm max}$',$
;      '$n_{\rm ref}$','$r_{\rm e,ref}$','$\alpha$','$\beta$',$
;      'F390W','F435W','F475W','F555W','F606W','F625W','F775W','F814W','F850LP',$
;      'F105W','F110W','F125W','F140W','F160W','$\chi^{2}_{\nu}$','$\nu$'])
;;    colhead2 = get_colhead([$
;;      'Band','(kpc)','(kpc)',$
;;;     '(mag~arcsec$^{-2}$)','(kpc)','$n$','$\chi^{2}_{\nu}$','$\nu$',$
;;      '(mag~arcsec$^{-2}$)','(kpc)','$n_{Ser}$','(mag~arcsec$^{-2}$)','(kpc)','$\chi^{2}_{\nu}$','$\nu$'],/nobreak)
;    
;    texcenter = [replicate('c',4),replicate('c',14),'c','c']
;
;;    tablenotetext = [$
;;;     '{a}{Ancillary information for our sample of galaxy clusters.}',$
;;      '{a}{Cluster coordinates based on the position of the BCG.}',$
;;      '{b}{Cluster virial mass adopted from \citet{merten14a}.}']
;
;; write it out    
;    openw, lun, texfile, /get_lun
;;   if (keyword_set(preprint) eq 0) then printf, lun, '\LongTables'
;;   printf, lun, '\clearpage'
;    printf, lun, '\setlength{\tabcolsep}{0.02in}'
;    printf, lun, '\LongTables'
;;   printf, lun, '\begin{turnpage}'
;    printf, lun, '\begin{landscape}'
;    printf, lun, '\tabletypesize{\tiny}'
;    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;    printf, lun, '\setlength{\tabcolsep}{0.005in}'
;    printf, lun, '\tablecaption{Surface Brightness Profile Fitting Results\label{table:sersic}}'
;    printf, lun, '\tablewidth{0pt}'
;    printf, lun, '\tablehead{'
;    niceprintf, lun, colhead1
;    niceprintf, lun, colhead2
;    printf, lun, '}'
;    printf, lun, '\startdata'
;    for ic = 0, ncl-1 do begin
;;   for ic = 0, ncl-1 do begin
;       cluster = strtrim(sample[ic].shortname,2)
;       nicecluster = get_nicecluster(sample[ic])
;;      printf, lun, '\cutinhead{'+nicecluster+'}'
;       if (ic eq ncl-1) then suffix = '' else suffix = ' \\'
;
;       suf = ''
;;      suf = '-alphabetazero'
;       results = mrdfits(sersicpath+cluster+'-allsersic'+suf+'-results.fits.gz',1,/silent)
;       sersic = mrdfits(sersicpath+cluster+'-allsersic'+suf+'.fits.gz',1,/silent)
;       print, cluster, results.sersic_chi2/results.sersic2_dof, $
;         results.sersic_n_ref, results.sersic_re_ref
;
;       format_flux_error, results.sersic_n_ref, results.sersic_n_ref_err, nn, nnerr
;       format_flux_error, results.sersic_re_ref, results.sersic_re_ref_err, re, reerr
;
;       format_flux_error, abs(results.sersic_alpha1), results.sersic_alpha1_err, alpha, alphaerr
;       if results.sersic_alpha1 lt 0 then alpha = '-'+alpha
;
;       format_flux_error, abs(results.sersic_beta1), results.sersic_beta1_err, beta, betaerr
;       if results.sersic_beta1 lt 0 then beta = '-'+beta
;       
;       line = '$'+nn+'\pm'+nnerr+'$ & '+'$'+re+'\pm'+reerr+'$ & '+$
;         '$'+alpha+'\pm'+alphaerr+'$ & '+'$'+beta+'\pm'+betaerr+'$ & '
;       for ff = 0, nfilt-1 do begin
;          this = where(short[ff] eq strtrim(sersic.band,2))
;          if this[0] eq -1 then begin
;             line = line+' \nodata & '
;; put nodata
;          endif else begin
;             format_flux_error, sersic[this].sersic_all_sbe, sersic[this].sersic_all_sbe_err, mu, muerr
;             line = line+'$'+mu+'\pm'+muerr+'$ & '
;          endelse
;       endfor 
;       line = line+'$'+strtrim(string(results.sersic_chi2/results.sersic2_dof,format='(F5.2)'),2)+'$ & '+$
;         '$'+strtrim(string(results.sersic_dof,format='(I0)'),2)+'$'
;       printf, lun, line+suffix
;    endfor 
;    
;    printf, lun, '\enddata'
;;   printf, lun, '\tablecomments{'+tablecomments+'}'
;;   niceprintf, lun, '\tablenotetext'+tablenotetext
;    printf, lun, '\end{deluxetable}'
;    printf, lun, '\clearpage'
;    printf, lun, '\end{landscape}'
;;   printf, lun, '\end{turnpage}'
;    free_lun, lun
;
;    


;;; ---------------------------------------------------------------------------    
;;; Table: double-Sersic fitting results for A209 and A2261
;;    texfile = paperpath+'table_doublesersic'+filesuffix+'.tex'
;;    splog, 'Writing '+texfile
;;
;;    colhead1 = get_colhead(['','Abell~209','Abell~2261'])
;;
;;    a209 = read_bcgmstar_sersic('a209',results=a209res)
;;    a2261 = read_bcgmstar_sersic('a2261',results=a2261res)
;;
;;    texcenter = ['l','c','c']
;;;   texcenter = replicate('c',(4+1)*7)
;;;   texcenter = replicate('c',(4+1)*14)
;;
;;;    tablenotetext = [$
;;;;     '{a}{Ancillary information for our sample of galaxy clusters.}',$
;;;      '{a}{Cluster coordinates based on the position of the BCG.}',$
;;;      '{b}{Cluster virial mass adopted from \citet{merten14a}.}']
;;
;;; write it out    
;;    openw, lun, texfile, /get_lun
;;;   if (keyword_set(preprint) eq 0) then printf, lun, '\LongTables'
;;;   printf, lun, '\clearpage'
;;;   printf, lun, '\LongTables'
;;;   printf, lun, '\begin{turnpage}'
;;;   printf, lun, '\begin{landscape}'
;;    printf, lun, '\begin{deluxetable*}{'+strjoin(texcenter)+'}'
;;;   printf, lun, '\setlength{\tabcolsep}{0.005in}'
;;    printf, lun, '\tablecaption{Surface Brightness Profile Fitting Results\label{table:doublesersic}}'
;;    printf, lun, '\tablewidth{0pt}'
;;    printf, lun, '\tabletypesize{\tiny}'
;;    printf, lun, '\tablehead{'
;;    niceprintf, lun, colhead1
;;;   niceprintf, lun, colhead2
;;    printf, lun, '}'
;;    printf, lun, '\startdata'
;;
;;    for ic = 0, ncl-1 do begin
;;;   for ic = 0, ncl-1 do begin
;;       cluster = strtrim(sample[ic].shortname,2)
;;       nicecluster = get_nicecluster(sample[ic])
;;       printf, lun, '\cutinhead{'+nicecluster+'}'
;;
;;       phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
;;;      for ff = nfilt-2, nfilt-1 do begin
;;;      for ff = 7, nfilt-1 do begin
;;       for ff = 0, nfilt-1 do begin
;;          this = where(short[ff] eq strtrim(phot.band,2))
;;          if this[0] eq -1 then begin
;;; put nodata
;;          endif else begin
;;             mu_ser = -2.5*alog10(phot[this].sersic2_sbe1)
;;             mu_ser_err = 2.5*phot[this].sersic2_sbe1_err/phot[this].sersic2_sbe1/alog(10)
;;             mu_dV = -2.5*alog10(phot[this].sersic2_sbe2)
;;             mu_dV_err = 2.5*phot[this].sersic2_sbe2_err/phot[this].sersic2_sbe2/alog(10)
;;
;;             if ff eq nfilt-1 then begin
;;                if (ic eq ncl-1) then suffix = '' else suffix = ' \\'
;;             endif else suffix = ' \\'
;;;            if ff eq nfilt-1 then if (ic eq ncl-1) then suffix = '' else $
;;;              suffix = ' \\' else suffix = ' & '
;;             printf, lun, strupcase(short[ff])+' & '+$
;;               '$'+strtrim(string(phot[this].rmin_kpc, format='(F5.2)'),2)+'$ & '+$
;;               '$'+strtrim(string(phot[this].rmax_kpc, format='(F5.2)'),2)+'$ & '+$
;;;; single Sersic
;;;               '$'+strtrim(string(phot[this].sersic_sbe, format='(F5.2)'),2)+'\pm'+$
;;;               strtrim(string(phot[this].sersic_sbe_err, format='(F5.2)'),2)+'$ & '+$
;;;               '$'+strtrim(string(phot[this].sersic_re,  format='(F5.2)'),2)+'\pm'+$
;;;               strtrim(string(phot[this].sersic_re_err,  format='(F5.2)'),2)+'$ & '+$
;;;               '$'+strtrim(string(phot[this].sersic_n,   format='(F5.2)'),2)+'\pm'+$
;;;               strtrim(string(phot[this].sersic_n_err,   format='(F5.2)'),2)+'$ & '+$
;;;;              '$'+strtrim(string((phot[this].sersic_chi2/phot[this].sersic_dof)>0.1,format='(F5.2)'),2)+'$ & '+$
;;;;              '$'+strtrim(string(phot[this].sersic_dof,     format='(I0)'),2)+'$ & '+$
;;; double Sersic
;;               '$'+strtrim(string(mu_ser, format='(F5.2)'),2)+'\pm'+$
;;               strtrim(string(mu_ser_err, format='(F5.2)'),2)+'$ & '+$
;;               '$'+strtrim(string(phot[this].sersic2_re1,  format='(F5.2)'),2)+'\pm'+$
;;               strtrim(string(phot[this].sersic2_re1_err,  format='(F5.2)'),2)+'$ & '+$
;;               '$'+strtrim(string(phot[this].sersic2_n1,  format='(F5.2)'),2)+'\pm'+$
;;               strtrim(string(phot[this].sersic2_n1_err,  format='(F5.2)'),2)+'$ & '+$
;;
;;               '$'+strtrim(string(mu_dV,  format='(F5.2)'),2)+'\pm'+$
;;               strtrim(string(mu_dV_err,  format='(F5.2)'),2)+'$ & '+$
;;               '$'+strtrim(string(phot[this].sersic2_re2,  format='(F5.2)'),2)+'\pm'+$
;;               strtrim(string(phot[this].sersic2_re2_err,  format='(F5.2)'),2)+'$ & '+$
;;
;;               '$'+strtrim(string((phot[this].sersic2_chi2/phot[this].sersic2_dof)>0.1,format='(F5.2)'),2)+'$ & '+$
;;               '$'+strtrim(string(phot[this].sersic2_dof,     format='(I0)'),2)+'$'+suffix
;;          endelse
;;       endfor
;;    endfor 
;;    
;;    printf, lun, '\enddata'
;;;   printf, lun, '\tablecomments{'+tablecomments+'}'
;;;   niceprintf, lun, '\tablenotetext'+tablenotetext
;;    printf, lun, '\end{deluxetable*}'
;;;   printf, lun, '\end{landscape}'
;;;   printf, lun, '\end{turnpage}'
;;;   printf, lun, '\clearpage'
;;    free_lun, lun
;;
;;
;;stop    
;;    
;;    



;;; #########################
;;; sub-Table 4 - double-Sersic surface brightnesses
;;    sample1 = sample[where(strtrim(sample.shortname,2) eq 'a209' or $
;;      strtrim(sample.shortname,2) eq 'a2261' or $
;;      strtrim(sample.shortname,2) eq 'rxj2248')]
;;    ncl1 = n_elements(sample1)
;;
;;    texfile = paperpath+'table_doublesersic2'+filesuffix+'.tex'
;;    splog, 'Writing '+texfile
;;    
;;    colhead1 = get_colhead(['','$\mu_{\rm e,1}$','$\mu_{\rm e,2}$'])
;;    colhead2 = get_colhead(['Band','(mag arcsec$^{-2}$)','(mag arcsec$^{-2}$)'],/nobreak)
;;    texcenter = ['l','c','c']
;;
;;    tablenotetext = [$
;;      '{a}{Together with Table~\ref{table:doublesersic1}, this table lists the best-fitting '+$
;;      'parameters for the BCGs fitted with the double-\sersic{} parametric model '+$
;;      'described in \S\ref{sec:sersic}.}']
;;
;;; write it out    
;;    openw, lun, texfile, /get_lun
;;;   if (keyword_set(preprint) eq 0) then printf, lun, '\LongTables'
;;;   printf, lun, '\clearpage'
;;;   printf, lun, '\setlength{\tabcolsep}{0.02in}'
;;;   printf, lun, '\LongTables'
;;;   printf, lun, '\begin{turnpage}'
;;;   printf, lun, '\begin{landscape}'
;;;   printf, lun, '\tabletypesize{\tiny}'
;;    printf, lun, '\begin{deluxetable*}{'+strjoin(texcenter)+'}'
;;;   printf, lun, '\setlength{\tabcolsep}{0.005in}'
;;    printf, lun, '\tablecaption{Double-\sersic{} Surface-Brightness Profile Fitting '+$
;;      'Results\tablenotemark{a}\label{table:doublesersic2}}'
;;    printf, lun, '\tablewidth{0pt}'
;;    printf, lun, '\tablehead{'
;;    niceprintf, lun, colhead1
;;    niceprintf, lun, colhead2
;;    printf, lun, '}'
;;    printf, lun, '\startdata'
;;
;;    for ic = 0, ncl1-1 do begin
;;       cluster = strtrim(sample1[ic].shortname,2)
;;       nicecluster = get_nicecluster(sample1[ic])
;;;      printf, lun, '\cline{1-3}'
;;       printf, lun, '\cutinhead{'+nicecluster+'}'
;;;      printf, lun, '\cline{1-3}'
;;
;;       sersic = read_bcgmstar_sersic(cluster,results=results)
;;       short1 = strupcase(strtrim(sersic.band,2))
;;       nfilt1 = n_elements(short1)
;;
;;       for ff = nfilt1-1, 0, -1 do begin
;;          if (ic eq ncl1-1) and ff eq 0 then suffix = '' else suffix = ' \\'
;;          line = short1[ff]+' & $'+$
;;            strtrim(string(sersic[ff].sersic2_all_sbe1,format='(F12.2)'),2)+'\pm'+$
;;            strtrim(string(sersic[ff].sersic2_all_sbe1_err,format='(F12.2)'),2)+'$'+' & $'+$
;;            strtrim(string(sersic[ff].sersic2_all_sbe2,format='(F12.2)'),2)+'\pm'+$
;;            strtrim(string(sersic[ff].sersic2_all_sbe2_err,format='(F12.2)'),2)+'$'+suffix
;;          printf, lun, line
;;       endfor
;;    endfor 
;;
;;    printf, lun, '\enddata'
;;;   printf, lun, '\tablecomments{'+tablecomments+'}'
;;    niceprintf, lun, '\tablenotetext'+tablenotetext
;;    printf, lun, '\end{deluxetable*}'
;;;   printf, lun, '\clearpage'
;;;   printf, lun, '\end{landscape}'
;;;   printf, lun, '\end{turnpage}'
;;    free_lun, lun
;;
