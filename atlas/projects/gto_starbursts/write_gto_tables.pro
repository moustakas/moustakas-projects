pro write_gto_tables, final1, textest=textest, emulateapj=emulateapj
; jm06jan31uofa - based on WRITE_ATLAS_TABLES

    if keyword_set(textest) then emulateapj = 0L else emulateapj = 1L
    emulateapj = 1L
    
    version = gto_log12oh_version()

    texpath = gto_path()
    metpath = gto_path()
    if keyword_set(textest) then texpath = texpath;+'textest/'

    if (n_elements(final1) eq 0L) then final1 = mrdfits(metpath+'gto_log12oh_final_'+version+'.fits.gz',1,/silent)
    ngalaxy = n_elements(final1)

    final = replicate(create_struct('galaxy', '', 'log12oh', [0.0,0.0], 'log12oh_texref', '', 'comment', '', $
      'distance', [0.0,0.0], 'distance_texref', ''),ngalaxy)
    final.galaxy          = final1.nice_galaxy
    final.log12oh         = transpose([[final1.log12oh_final],[final1.log12oh_final_err]])
    final.distance        = transpose([[final1.distance],[final1.distance_err]])
    final.log12oh_texref  = final1.log12oh_final_ref
    final.comment         = final1.log12oh_final_comment
    final.distance_texref = final1.distance_texref

; TEMPORARY HACK BECAUSE I DON'T WANT TO RERUN WRITE_HII_REGIONS!!

    final.log12oh_texref = repstr(final.log12oh_texref,'vanzee06','vanzee06b')
    final.log12oh_texref = repstr(final.log12oh_texref,'vanbruegel85','vanbreugel85')

;   for ig = 0L, ngalaxy-1L do final[ig].texref = strjoin([final1[ig].log12oh_final_ref,final1[ig].distance_texref],',')

; ---------------------------------------------------------------------------    
; Table: final abundances
; ---------------------------------------------------------------------------    

    select = ['galaxy','log12oh','log12oh_texref','distance','distance_texref']
    format = ['A0','F12.2','A0','F12.5','A0']
    texcenter = ['l','c','c','c','c']

;   colhead1 = '\colhead{'+[ ['','','O/H','$D$\tablenotemark{b}']+'} & ',['Distance']+'} \\ ']
;   colhead2 = '\colhead{'+[ ['Galaxy','$12+\log\,({\rm O/H})$\tablenotemark{a}','Refs.','(Mpc)']+'} & ',['Refs.']+'} \\ ']
    colhead1 = '\colhead{'+[ ['Galaxy','$12+\log\,({\rm O/H})$\tablenotemark{a}','Refs.','$D$ (Mpc)\tablenotemark{b}']+'} & ',['Refs.']+'} \\ ']
    colunits = '\colhead{'+[ ['(1)','(2)','(3)','(4)']+'} & ',['(5)']+'} ']

    table = html_structure_parse(final,tags=select,tagformats=format,$
      tagindices=match,/keep_error,blank='\nodata')
    ntags = n_tags(table)
    tags = tag_names(table)

; notes    

    ohnotes = strarr(ngalaxy)
    original = where(strmatch(final.comment,'*original*',/fold) eq 1B)
    ohnotes[original] = '\tablenotemark{c}'

    hii_note = where((strmatch(final.comment,'*HII/PT05*',/fold) eq 1B) or (strmatch(final.comment,'*HII/O3N2*',/fold) eq 1B))
    ohnotes[hii_note] = '\tablenotemark{d}'
    
    mk06_note = where((strmatch(final.comment,'*MK06/PT05*',/fold) eq 1B) or (strmatch(final.comment,'*MK06/O3N2*',/fold) eq 1B))
    ohnotes[mk06_note] = '\tablenotemark{e}'
    
    tablenotetext = [$
      '{a}{Unless otherwise noted, oxygen abundances have been derived using the direct, '+$
      'electron-temperature ($T_{e}$) method applied to individual \ion{H}{2} regions \citep[e.g.,][]{skillman98}.  The electron temperature in '+$
      'the O$^{++}$ zone was derived using the [\ion{O}{3}]~$\lambda\lambda4959,5007$/[\ion{O}{3}]~$\lambda4363$ '+$
      'ratios published in the references listed.  The O$^{+}$ temperature was predicted using the relation given by '+$
      '\citet{garnett92}, and the total oxygen abundance was subsequently computed as ${\rm O/H}={\rm O}^{+}/{\rm H}+{\rm O}^{++}/{\rm H}$ '+$
      'using the electron density derived from the [\ion{S}{2}]~$\lambda\lambda6716,6731$ doublet ratio \citep{shaw95}.  '+$
      'When multiple abundance estimates were available, we adopted the mean and the standard deviation of all the measurements as '+$
      'the final metallicity and error, assuming a minimum uncertainty of $0.05$~dex in the $T_{e}$ method.}', $
      '{b}{See \S2 for details about how distances for the sample were computed.}',$
      '{c}{Emission-line fluxes for these objects were not given; therefore, we adopted the published oxygen abundances and errors.}',$
      '{d}{Unfortunately, the auroral [O~{\sc iii}]~$\lambda4363$ line in these objects was not detected; therefore, to '+$
      'estimate the oxygen abundance we use two different strong-line calibrations empirically tied to the electron temperature '+$
      'abundance scale.  For high-excitation H~{\sc ii} regions we use the \citet{pilyugin05} calibration of '+$
      '$R_{23}\equiv$([O~{\sc ii}]~$\lambda3727 + $[O~{\sc iii}]~$\lambda\lambda4959,5007)/$H$\beta$, '+$
      'and for low-excitation H~{\sc ii} regions we use the \citet{pettini04} calibration of '+$
      '([O~{\sc iii}]~$\lambda5007$/H$\beta$)/([N~{\sc ii}]~$\lambda6584$/H$\alpha$).  We adopt '+$
      '$P\equiv$[O~{\sc iii}]/([O~{\sc ii}]+[O~{\sc iii}])$=0.4$ as the boundary between low- and high-excitation '+$
      'H~{\sc ii} regions \citep{pilyugin05}, and assume a minimum uncertainty of $0.1$~dex in either method.}',$
      '{e}{For these objects no individual H~{\sc ii} regions have been observed; therefore, we applied the same methodology described in the '+$
      'preceeding footnote to the integrated emission-line fluxes published by \citet{moustakas06a}.}']

; which is calibrated against a large number of \ion{H}{2} regions with electron-temperature abundance measurements.    

    
;   tablecomments = ['']
;   tablecomments = strjoin(tablecomments,' ')

; parse the metallicity references and assign unique reference numbers

    log12ohrefs = strsplit(strjoin(strtrim(reform(table.log12oh_texref),2),','),',',/extract)
    nlog12ohrefs = n_elements(log12ohrefs)

    log12ohrefnum = lonarr(nlog12ohrefs)-1L
    refcount = 1L
    for i = 0L, nlog12ohrefs-1L do begin
       if (log12ohrefnum[i] eq -1L) then begin
          match = where(log12ohrefs[i] eq log12ohrefs)
          log12ohrefnum[match] = refcount
          refcount = refcount + 1L
       endif
    endfor

;   niceprint, log12ohrefs, log12ohrefnum
    ulog12ohrefnum = uniq(log12ohrefnum,sort(log12ohrefnum))
    tablelog12ohrefs = strjoin('('+string(log12ohrefnum[ulog12ohrefnum],format='(I0)')+') \citet{'+$
      string(log12ohrefs[ulog12ohrefnum],format='(A0)')+'}','; ')+'.'
    tablelog12ohrefs = repstr(tablelog12ohrefs,'\citet{This_paper}','this paper')

; parse the distance references and assign unique reference numbers,
; but don't start the counter over

    distrefs = strsplit(strjoin(strtrim(reform(table.distance_texref),2),','),',',/extract)
    ndistrefs = n_elements(distrefs)

    distrefnum = lonarr(ndistrefs)-1L
;   refcount = 1L
    for i = 0L, ndistrefs-1L do begin
       if (distrefnum[i] eq -1L) then begin
          match = where(distrefs[i] eq distrefs)
          distrefnum[match] = refcount
          refcount = refcount + 1L
       endif
    endfor

;   niceprint, distrefs, distrefnum
    udistrefnum = uniq(distrefnum,sort(distrefnum))
    tabledistrefs = strjoin('('+string(distrefnum[udistrefnum],format='(I0)')+') \citet{'+$
      string(distrefs[udistrefnum],format='(A0)')+'}','; ')+'.'
    tabledistrefs = repstr(repstr(tabledistrefs,'\citet{Flow_Model}','this paper'),$
      '\citet{HYPERLEDA_(Oct._2007)}','HYPERLEDA (Oct. 2007)')

; now concatenate

    tablerefs = strjoin([tablelog12ohrefs,tabledistrefs],';')
    
; write out

    texfile = 'gto_table1.tex'
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
       printf, lun, '\LongTables'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\small}'
    printf, lun, '\tablecaption{Sample Galaxies and Basic Information\label{tab:sample}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
;   niceprintf, lun, colhead2
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0L, ngalaxy-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else if (i lt ngalaxy-1L) then suffix = ' \\ ' else suffix = ''
          line[j] = table[i].(j)+suffix
          if (strmatch(tags[j],'*galaxy*',/fold) eq 0B) and (strmatch(tags[j],'*texref*',/fold) eq 0B) then begin
             if (strtrim(table[i].(j)[1],2) lt 0.0) then line[j] = '\nodata'+suffix else begin
                if (strmatch(tags[j],'*distance*',/fold) eq 1B) then begin
                   dist = table[i].(j)[0] & dist_err = table[i].(j)[1]
                   format_flux_error, dist, dist_err, newdist, newdist_err
                   line[j] = '$'+strtrim(newdist,2)+'\pm'+strtrim(newdist_err,2)+'$'+suffix
                endif else begin
                   line[j] = '$'+strtrim(table[i].(j)[0],2)+'\pm'+strtrim(table[i].(j)[1],2)+'$'+ohnotes[i]+suffix
                endelse
             endelse
          endif
          if strmatch(tags[j],'log12oh_texref',/fold) then begin
             theserefs = strsplit(strtrim(table[i].log12oh_texref,2),',',/extract)
             these = cmset_op(log12ohrefs,'and',theserefs,/index)
             thesenum = log12ohrefnum[these[uniq(log12ohrefnum[these],sort(log12ohrefnum[these]))]]
             line[j] = strjoin(string(thesenum,format='(I0)'),',')+suffix
          endif
          if strmatch(tags[j],'distance_texref',/fold) then begin
             theserefs = strsplit(strtrim(table[i].distance_texref,2),',',/extract)
             these = cmset_op(distrefs,'and',theserefs,/index)
             thesenum = distrefnum[these[uniq(distrefnum[these],sort(distrefnum[these]))]]
             line[j] = strjoin(string(thesenum,format='(I0)'),',')+suffix
          endif
       endfor
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\tablerefs{'+tablerefs+'}'
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then begin
       printf, lun, '\bibliographystyle{apj}'
       printf, lun, '\bibliography{/Users/ioannis/ay/bin/bibtex/im}'
       printf, lun, '\end{document}'
    endif
    free_lun, lun

    if keyword_set(textest) then begin
       pushd, texpath
       spawn, 'mkps '+repstr(texfile,'.tex',''), /sh
;      spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

return
end
