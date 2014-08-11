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
; Table: Sersic fitting results
    texfile = paperpath+'table_sersic'+filesuffix+'.tex'
    splog, 'Writing '+texfile

; build the column headings
    colhead1 = get_colhead([$
;     '$r_{\rm min}$','$r_{\rm max}$',$
      '$n_{\rm ref}$','$r_{\rm e,ref}$','$\alpha$','$\beta$',$
      'F390W','F435W','F475W','F555W','F606W','F625W','F775W','F814W','F850LP',$
      'F105W','F110W','F125W','F140W','F160W','$\chi^{2}_{\nu}$','$\nu$'])
;    colhead2 = get_colhead([$
;      'Band','(kpc)','(kpc)',$
;;     '(mag~arcsec$^{-2}$)','(kpc)','$n$','$\chi^{2}_{\nu}$','$\nu$',$
;      '(mag~arcsec$^{-2}$)','(kpc)','$n_{Ser}$','(mag~arcsec$^{-2}$)','(kpc)','$\chi^{2}_{\nu}$','$\nu$'],/nobreak)
    
    texcenter = [replicate('c',4),replicate('c',14),'c','c']

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
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\begin{deluxetable*}{'+strjoin(texcenter)+'}'
    printf, lun, '\setlength{\tabcolsep}{0.005in}'
    printf, lun, '\tablecaption{Surface Brightness Profile Fitting Results\label{table:sersic}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    printf, lun, '}'
    printf, lun, '\startdata'
    for ic = 0, ncl-1 do begin
;   for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       nicecluster = get_nicecluster(sample[ic])
;      printf, lun, '\cutinhead{'+nicecluster+'}'
       if (ic eq ncl-1) then suffix = '' else suffix = ' \\'

       results = mrdfits(sersicpath+cluster+'-allsersic-results.fits.gz',1,/silent)
       sersic = mrdfits(sersicpath+cluster+'-allsersic.fits.gz',1,/silent)

       format_flux_error, results.sersic_n_ref, results.sersic_n_ref_err, nn, nnerr
       format_flux_error, results.sersic_n_ref, results.sersic_n_ref_err, re, reerr

       format_flux_error, abs(results.sersic_alpha1), results.sersic_alpha1_err, alpha, alphaerr
       if results.sersic_alpha1 lt 0 then alpha = '-'+alpha

       format_flux_error, abs(results.sersic_beta1), results.sersic_beta1_err, beta, betaerr
       if results.sersic_beta1 lt 0 then beta = '-'+beta
       
       line = '$'+nn+'\pm'+nnerr+'$ & '+'$'+re+'\pm'+reerr+'$ & '+$
         '$'+alpha+'\pm'+alphaerr+'$ & '+'$'+beta+'\pm'+betaerr+'$ & '
       for ff = 0, nfilt-1 do begin
          this = where(short[ff] eq strtrim(sersic.band,2))
          if this[0] eq -1 then begin
             line = line+' \nodata & '
; put nodata
          endif else begin
             format_flux_error, sersic[this].sersic_all_sbe, sersic[this].sersic_all_sbe_err, mu, muerr
             line = line+'$'+mu+'\pm'+muerr+'$ & '
          endelse
       endfor 
       line = line+'$'+strtrim(string(results.sersic_chi2/results.sersic2_dof,format='(F5.2)'),2)+'$ & '+$
         '$'+strtrim(string(results.sersic_dof,format='(I0)'),2)+'$'
       printf, lun, line+suffix
    endfor 
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
;   niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable*}'
;   printf, lun, '\end{landscape}'
;   printf, lun, '\end{turnpage}'
;   printf, lun, '\clearpage'
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

; ---------------------------------------------------------------------------    
; Table: Sersic fitting results
    texfile = paperpath+'table_sersic'+filesuffix+'.tex'
    splog, 'Writing '+texfile

; build the column headings
;   colhead1 = get_colhead([$
;     '\multicolumn{4}{c}{F390W}','',$
;     '\multicolumn{4}{c}{F435W}','',$
;     '\multicolumn{4}{c}{F475W}','',$
;     '\multicolumn{4}{c}{F555W}','',$
;     '\multicolumn{4}{c}{F606W}','',$
;     '\multicolumn{4}{c}{F625W}','',$
;     '\multicolumn{4}{c}{F775W}','',$
;     '\multicolumn{4}{c}{F814W}','',$
;     '\multicolumn{4}{c}{F850LP}','',$
;     '\multicolumn{4}{c}{F105W}','',$
;     '\multicolumn{4}{c}{F110W}','',$
;     '\multicolumn{4}{c}{F125W}','',$
;     '\multicolumn{4}{c}{F140W}','',$
;     '\multicolumn{4}{c}{F160W}'])

;    colhead1 = get_colhead([$
;      '','$r_{\rm min}$','$r_{\rm max}$',$
;      '$\mu_{e}$','$r_{e}$','',$
;      '$\mu_{e,Ser}$','$\mu_{e,dV}$','$r_{e,Ser}$','$r_{e,dV}$',''])
;    colhead2 = get_colhead([$
;      'Band','(kpc)','(kpc)',$
;      '(mag~arcsec$^{-2}$)','(kpc)','$n$',$
;      '(mag~arcsec$^{-2}$)','(mag~arcsec$^{-2}$)','(kpc)','(kpc)','$n_{Ser}$'],/nobreak)

;   colhead1 = get_colhead([$
;     '','$r_{\rm min}$','$r_{\rm max}$',$
;     '$\mu_{e}$','$r_{e}$','','','',$
;     '$\mu_{e,Ser}$','$\mu_{e,dV}$','$r_{e,Ser}$','$r_{e,dV}$','','',''])
;   colhead2 = get_colhead([$
;     'Band','(kpc)','(kpc)',$
;     '(mag~arcsec$^{-2}$)','(kpc)','$n$','$\chi^{2}_{\nu}$','$\nu$',$
;     '(mag~arcsec$^{-2}$)','(mag~arcsec$^{-2}$)','(kpc)','(kpc)','$n_{Ser}$','$\chi^{2}_{\nu}$','$\nu$'],/nobreak)
    
    colhead1 = get_colhead([$
      '','$r_{\rm min}$','$r_{\rm max}$',$
;     '$\mu_{e}$','$r_{e}$','','','',$
      '$\mu_{e,Ser}$','$r_{e,Ser}$','','$\mu_{e,dV}$','$r_{e,dV}$','',''])
    colhead2 = get_colhead([$
      'Band','(kpc)','(kpc)',$
;     '(mag~arcsec$^{-2}$)','(kpc)','$n$','$\chi^{2}_{\nu}$','$\nu$',$
      '(mag~arcsec$^{-2}$)','(kpc)','$n_{Ser}$','(mag~arcsec$^{-2}$)','(kpc)','$\chi^{2}_{\nu}$','$\nu$'],/nobreak)
    
    texcenter = ['l',replicate('c',9)]
;   texcenter = replicate('c',(4+1)*7)
;   texcenter = replicate('c',(4+1)*14)

;    tablenotetext = [$
;;     '{a}{Ancillary information for our sample of galaxy clusters.}',$
;      '{a}{Cluster coordinates based on the position of the BCG.}',$
;      '{b}{Cluster virial mass adopted from \citet{merten14a}.}']

; write it out    
    openw, lun, texfile, /get_lun
;   if (keyword_set(preprint) eq 0) then printf, lun, '\LongTables'
;   printf, lun, '\clearpage'
    printf, lun, '\LongTables'
;   printf, lun, '\begin{turnpage}'
;   printf, lun, '\begin{landscape}'
    printf, lun, '\begin{deluxetable*}{'+strjoin(texcenter)+'}'
    printf, lun, '\setlength{\tabcolsep}{0.005in}'
    printf, lun, '\tablecaption{Surface Brightness Profile Fitting Results\label{table:sersic}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    printf, lun, '}'
    printf, lun, '\startdata'
    for ic = 0, ncl-1 do begin
;   for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       nicecluster = get_nicecluster(sample[ic])
       printf, lun, '\cutinhead{'+nicecluster+'}'

       phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
;      for ff = nfilt-2, nfilt-1 do begin
;      for ff = 7, nfilt-1 do begin
       for ff = 0, nfilt-1 do begin
          this = where(short[ff] eq strtrim(phot.band,2))
          if this[0] eq -1 then begin
; put nodata
          endif else begin
             mu_ser = -2.5*alog10(phot[this].sersic2_sbe1)
             mu_ser_err = 2.5*phot[this].sersic2_sbe1_err/phot[this].sersic2_sbe1/alog(10)
             mu_dV = -2.5*alog10(phot[this].sersic2_sbe2)
             mu_dV_err = 2.5*phot[this].sersic2_sbe2_err/phot[this].sersic2_sbe2/alog(10)

             if ff eq nfilt-1 then begin
                if (ic eq ncl-1) then suffix = '' else suffix = ' \\'
             endif else suffix = ' \\'
;            if ff eq nfilt-1 then if (ic eq ncl-1) then suffix = '' else $
;              suffix = ' \\' else suffix = ' & '
             printf, lun, strupcase(short[ff])+' & '+$
               '$'+strtrim(string(phot[this].rmin_kpc, format='(F5.2)'),2)+'$ & '+$
               '$'+strtrim(string(phot[this].rmax_kpc, format='(F5.2)'),2)+'$ & '+$
;; single Sersic
;               '$'+strtrim(string(phot[this].sersic_sbe, format='(F5.2)'),2)+'\pm'+$
;               strtrim(string(phot[this].sersic_sbe_err, format='(F5.2)'),2)+'$ & '+$
;               '$'+strtrim(string(phot[this].sersic_re,  format='(F5.2)'),2)+'\pm'+$
;               strtrim(string(phot[this].sersic_re_err,  format='(F5.2)'),2)+'$ & '+$
;               '$'+strtrim(string(phot[this].sersic_n,   format='(F5.2)'),2)+'\pm'+$
;               strtrim(string(phot[this].sersic_n_err,   format='(F5.2)'),2)+'$ & '+$
;;              '$'+strtrim(string((phot[this].sersic_chi2/phot[this].sersic_dof)>0.1,format='(F5.2)'),2)+'$ & '+$
;;              '$'+strtrim(string(phot[this].sersic_dof,     format='(I0)'),2)+'$ & '+$
; double Sersic
               '$'+strtrim(string(mu_ser, format='(F5.2)'),2)+'\pm'+$
               strtrim(string(mu_ser_err, format='(F5.2)'),2)+'$ & '+$
               '$'+strtrim(string(phot[this].sersic2_re1,  format='(F5.2)'),2)+'\pm'+$
               strtrim(string(phot[this].sersic2_re1_err,  format='(F5.2)'),2)+'$ & '+$
               '$'+strtrim(string(phot[this].sersic2_n1,  format='(F5.2)'),2)+'\pm'+$
               strtrim(string(phot[this].sersic2_n1_err,  format='(F5.2)'),2)+'$ & '+$

               '$'+strtrim(string(mu_dV,  format='(F5.2)'),2)+'\pm'+$
               strtrim(string(mu_dV_err,  format='(F5.2)'),2)+'$ & '+$
               '$'+strtrim(string(phot[this].sersic2_re2,  format='(F5.2)'),2)+'\pm'+$
               strtrim(string(phot[this].sersic2_re2_err,  format='(F5.2)'),2)+'$ & '+$

               '$'+strtrim(string((phot[this].sersic2_chi2/phot[this].sersic2_dof)>0.1,format='(F5.2)'),2)+'$ & '+$
               '$'+strtrim(string(phot[this].sersic2_dof,     format='(I0)'),2)+'$'+suffix
          endelse
       endfor
    endfor 
    
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
