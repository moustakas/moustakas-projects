pro nfgs_compare_lines, atlas, nfgs, jansen, paper=paper
; jm05aug01uofa - compare EW measurements of the 11 galaxies in common
;                 between the NFGS and our survey; compare Rolf's
;                 measurements of his spectra to my measurements of
;                 his spectra, *and* my measurements of my spectra 

    paperpath = atlas_path(/papers)+'atlas/'
    
    if (n_elements(atlas) eq 0L) then atlas = read_integrated()
    if (n_elements(nfgs) eq 0L) then nfgs = read_nfgs()
    if (n_elements(jansen) eq 0L) then jansen = read_00jansen()

; which NFGS galaxies are in the atlas    

    raref = 15.0*im_hms2dec(atlas.ra)
    deref = im_hms2dec(atlas.dec)

    ra = 15.0*im_hms2dec(nfgs.ra)
    de = im_hms2dec(nfgs.dec)
    
    jansen_ra = 15.0*im_hms2dec(jansen.ra)
    jansen_de = im_hms2dec(jansen.dec)
    
    ntot = djs_angle_match(raref,deref,ra,de,dtheta=15.0/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    match = where(mindx ne -1,nmatch)
    niceprint, atlas[match].galaxy, nfgs[mindx[match]].galaxy, mdist[match]*3600.0
    splog, 'There are '+string(nmatch,format='(I0)')+' atlas galaxies in the NFGS.'

    ntot = djs_angle_match(ra[mindx[match]],de[mindx[match]],jansen_ra,jansen_de,dtheta=5.0/3600.0,$
      units='degrees',mcount=mcount,mindx=jansen_mindx,mdist=mdist)
    jansen_match = where(jansen_mindx ne -1,jansen_nmatch)
    
    matchatlas = atlas[match]
    matchnfgs = nfgs[mindx[match]]
    matchjansen = jansen[jansen_mindx[jansen_match]]

    s = sort(matchjansen.nfgs_id)
    matchjansen = matchjansen[s]
    matchnfgs = matchnfgs[s]
    matchatlas = matchatlas[s]

    niceprint, matchjansen.nfgs_id, matchjansen.galaxy, matchjansen.nedgalaxy, $
      matchatlas.galaxy, matchnfgs.galaxy, matchjansen.oii_3727_ew[0], $
      matchnfgs.oii_3727_ew[0], matchatlas.oii_3727_ew[0]
    
    print
    niceprint, matchjansen.nfgs_id, matchjansen.galaxy, matchjansen.nedgalaxy, $
      matchatlas.galaxy, matchnfgs.galaxy, matchjansen.oiii_5007_ew[0], $
      matchnfgs.oiii_5007_ew[0], matchatlas.oiii_5007_ew[0]

    print
    niceprint, matchjansen.nfgs_id, matchjansen.galaxy, matchjansen.nedgalaxy, $
      matchatlas.galaxy, matchnfgs.galaxy, matchjansen.h_beta_ew[0], $
      matchnfgs.h_beta_ew[0], matchatlas.h_beta_ew[0]
    
    print
    niceprint, matchjansen.nfgs_id, matchjansen.galaxy, matchjansen.nedgalaxy, $
      matchatlas.galaxy, matchnfgs.galaxy, matchjansen.h_alpha_ew[0], $
      matchnfgs.h_alpha_ew[0], matchatlas.h_alpha_ew[0]
    
; generate a LATEX table for the line comparison

    if keyword_set(paper) then begin
       
       openw, lun, paperpath+'nfgs_line_compare.tex', /get_lun
       printf, lun, '\begin{deluxetable}{lccccccccccc}'
       printf, lun, '\tabletypesize{\small}'
;      printf, lun, '\rotate'
       printf, lun, '\tablecolumns{12}'
       printf, lun, '\tablecaption{Emission-Line Equivalent Width Comparison\tablenotemark{a} \label{table:nfgs_compare}}'
       printf, lun, '\tablewidth{0in}'
       printf, lun, '\tablehead{'

       printf, lun, '\colhead{} & '
       printf, lun, '\multicolumn{5}{c}{NFGS} & '
       printf, lun, '\colhead{} & '
       printf, lun, '\multicolumn{5}{c}{This Paper} \\ '

       printf, lun, '\cline{2-6} \cline{8-12} \\'
       printf, lun, '\colhead{Galaxy} & '
       printf, lun, '\colhead{ID\tablenotemark{b}} & \colhead{\oii} & \colhead{\oiii} & \colhead{\hb} & \colhead{\ha} & '
       printf, lun, '\colhead{} & ' 
       printf, lun, '\colhead{ID} & \colhead{\oii} & \colhead{\oiii} & \colhead{\hb} & \colhead{\ha} \\ '

       printf, lun, '\colhead{} & \colhead{} & '
       printf, lun, '\colhead{(\AA)} & \colhead{(\AA)} & \colhead{(\AA)} & \colhead{(\AA)} & '
       printf, lun, '\colhead{} & \colhead{} & '
       printf, lun, '\colhead{(\AA)} & \colhead{(\AA)} & \colhead{(\AA)} & \colhead{(\AA)} '
       printf, lun, '}'
       printf, lun, '\startdata'
       for i = 0L, nmatch-1L do begin

          galaxy = strtrim(repstr(repstr(matchatlas[i].galaxy,'NGC','NGC '),'UGC','UGC '),2)
          nfgs_id = string(matchjansen[i].nfgs_id,format='(I3.3)')
          atlas_id = string(matchatlas[i].atlas_id,format='(I3.3)')

; table with no errors          
          
          nfgs_oii  = string(matchjansen[i].oii_3727_ew[0],format='(F5.1)')
          nfgs_oiii = string(matchjansen[i].oiii_5007_ew[0],format='(F5.1)')
          nfgs_hb   = string(matchjansen[i].h_beta_ew[0],format='(F5.1)')
          nfgs_ha   = string(matchjansen[i].h_alpha_ew[0],format='(F5.1)')
          
          atlas_oii  = string(matchatlas[i].oii_3727_ew[0],format='(F5.1)')
          atlas_oiii = string(matchatlas[i].oiii_5007_ew[0],format='(F5.1)')
          atlas_hb   = string(matchatlas[i].h_beta_ew[0],format='(F5.1)')
          atlas_ha   = string(matchatlas[i].h_alpha_ew[0],format='(F5.1)')

          printf, lun, $
            galaxy+' & '+$
            nfgs_id+' & '+'$'+nfgs_oii+'$ & $'+nfgs_oiii+'$ & $'+nfgs_hb+'$ & $'+nfgs_ha+'$ &  & '+$
            atlas_id+' & '+'$'+atlas_oii+'$ & $'+atlas_oiii+'$ & $'+atlas_hb+'$ & $'+atlas_ha+'$ \\'

; table with errors          
          
;         atlas_oii1      = matchatlas[i].oii_3727_ew[0]
;         atlas_oiii1     = matchatlas[i].oiii_5007_ew[0]
;         atlas_hb1       = matchatlas[i].h_beta_ew[0]
;         atlas_ha1       = matchatlas[i].h_alpha_ew[0]
;
;         atlas_oii1_err  = matchatlas[i].oii_3727_ew[1]
;         atlas_oiii1_err = matchatlas[i].oiii_5007_ew[1]
;         atlas_hb1_err   = matchatlas[i].h_beta_ew[1]
;         atlas_ha1_err   = matchatlas[i].h_alpha_ew[1]
;
;         format_flux_error, atlas_oii1, atlas_oii1_err, atlas_oii, atlas_oii_err
;         format_flux_error, atlas_oiii1, atlas_oiii1_err, atlas_oiii, atlas_oiii_err
;         format_flux_error, atlas_hb1, atlas_hb1_err, atlas_hb, atlas_hb_err
;         format_flux_error, atlas_ha1, atlas_ha1_err, atlas_ha, atlas_ha_err
          
;         printf, lun, $
;           nfgs_id+' & '+galaxy+' & '+$
;           '$'+nfgs_oii+'$ & $'+nfgs_oiii+'$ & $'+nfgs_hb+'$ & $'+nfgs_ha+'$ &  & '+$
;           '$'+atlas_oii+'\pm'+atlas_oii_err+'$ & '+$
;           '$'+atlas_oiii+'\pm'+atlas_oiii_err+'$ & '+$
;           '$'+atlas_hb+'\pm'+atlas_hb_err+'$ & '+$
;           '$'+atlas_ha+'\pm'+atlas_ha_err+'$ \\'

; old table - ignore          
          
;         printf, lun, $
;           string(matchjansen[i].nfgs_id,format='(I3.3)')+' & '+$
;           galaxy+' & '+$
;           string(matchjansen[i].oii_3727_ew[0],format='(F5.1)')+' & '+$
;           string(matchjansen[i].oiii_5007_ew[0],format='(F5.1)')+' & '+$
;           string(matchjansen[i].h_alpha_ew[0],format='(F5.1)')+' &  & '+$
;           '$'+string(matchatlas[i].oii_3727_ew[0],format='(F5.1)')+'\pm'+$
;           string(matchatlas[i].oii_3727_ew[1],format='(F5.1)')+'$ & '+$
;           '$'+string(matchatlas[i].oiii_5007_ew[0],format='(F5.1)')+'\pm'+$
;           string(matchatlas[i].oiii_5007_ew[1],format='(F5.1)')+'$ & '+$
;           '$'+string(matchatlas[i].h_alpha_ew[0],format='(F5.1)')+'\pm'+$
;           string(matchatlas[i].h_alpha_ew[1],format='(F5.1)')+'$ \\'

       endfor

       printf, lun, '\enddata'
       printf, lun, '\tablenotetext{a}{Comparison of integrated emission-line equivalent widths of '+$
         '\oiilam, \oiiilam, \hblam{}, and \halam{} for the galaxies in common with the NFGS.}'
       printf, lun, '\tablenotetext{b}{NFGS identification number adopted by \citet{jansen00a, jansen00b}.}'
;      printf, lun, '%\tablerefs{}'
       printf, lun, '\end{deluxetable}'
       free_lun, lun

    endif

return
end

