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
    
; by default build emulateapj format tables    
    sample = read_bcgsfhs_sample()
    sample = sample[sort(sample.mvir)]
    ncl = n_elements(sample) 

    paperpath = bcgmstar_path(/paper)

; ---------------------------------------------------------------------------    
; Table 1: table summarizing our fields
    texfile = paperpath+'table_sample'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    
; gather the information we need
    info = replicate({cluster: '', z: '', ra: '', dec: '', mvir: ''},ncl)
    info.cluster = repstr(repstr(strupcase(sample.dirname),'ABELL','Abell'),'_',' ')
    info.z = '$'+strtrim(string(sample.z,format='(F12.3)'),2)+'$'
    info.ra = sample.ra
    info.dec = sample.dec
    info.mvir = '$'+strtrim(string(sample.mvir/1D15,format='(F12.2)'),2)+'\pm'+$
      strtrim(string(sample.mvirerr/1D15,format='(F12.2)'),2)+'$'
    struct_print, info
    ntags = n_tags(info)

    colhead1 = get_colhead(['','','','','$M_{\rm vir}$\tablenotemark{b}'])
    colhead2 = get_colhead(['Cluster','z','$\alpha_{J2000}$','$\delta_{J2000}$','$(10^{15}~M_{\odot}/h)$'],/nobreak)
    texcenter = ['l','c','c','c','c']

    tablenotetext = [$
      '{a}{Ancillary information for our sample of galaxy clusters.}',$
      '{b}{Cluster virial mass adopted from \citet{merten14a} {\bf except for....}.}']

; write it out    
    openw, lun, texfile, /get_lun
;   if (keyword_set(preprint) eq 0) then printf, lun, '\LongTables'
;   printf, lun, '\begin{landscape}'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;   printf, lun, '\tabletypesize{\small}'
    printf, lun, '\tablecaption{Cluster Sample\tablenotemark{a}\label{table:sample}}'  
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

return
end

