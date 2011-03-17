pro write_atlas_tables, atlas, int, nuc, textest=textest
; jm04mar23uofa
; jm05jul26uofa - updated
; generate all the tables for the atlas data paper

    light = 2.99792458D5        ; speed of light [km/s]
    scale = 1D15
    sigcut = 1.0

    nstub = 10L
    
    notesymbol = ['a','b','c','d','e','f','g','h','i','j',$
      'k','l','m','n','o','p','q','r','s','t','u','v','w',$
      'x','y','z']

; read the integrated and nuclear atlases

    analysis_path = atlas_path(/analysis)
    texpath = atlas_path(/papers)+'atlas/'
    mrtpath = atlas_path(/papers)+'atlas/MRT/'
    datapath = atlas_path(/papers)+'atlas/DATATABLES/'
    if keyword_set(textest) then texpath = texpath+'textest/'

    if (n_elements(atlas) eq 0L) then atlas = mrdfits(analysis_path+'atlas1d_info.fits.gz',1,/silent)
    galaxy = strtrim(atlas.galaxy,2)
    ngalaxy = n_elements(atlas)

    if (n_elements(int) eq 0L) then int = read_integrated()
    if (n_elements(nuc) eq 0L) then nuc = read_nuclear()

    int = int[sort(int.atlas_id)]
    nuc = nuc[sort(nuc.atlas_id)]
    
    nint = n_elements(int)
    nnuc = n_elements(nuc)
    
    voffset = '0.5'

; ---------------------------------------------------------------------------    
; Table 1: general properties    
; ---------------------------------------------------------------------------    

    select = [$
      'atlas_id','nice_galaxy','alt_galaxy','ra','dec',$
      'cz','ebv_mw','lit_type','d25_maj',$
      'd25_min','posangle','distance','distance_texref'$
      ]
    format = [$
      'I3.3','A0','A0','A10','A9',$
      'I0','F12.3','A0','F12.2',$
      'F12.2','I0','F12.1','A0'$
      ]
    texcenter = [$
      'c','l','l','r','r',$
      'c','l','l','c',$
      'c','c','r','c']

    colhead = '\colhead{'+[ [$
      'ID','Galaxy Name','Other Names','$\alpha_{\rm J2000.0}$','$\delta_{\rm J2000.0}$',$
      'cz','E(\bv)','Type','D$_{25}^{\rm maj}$',$
      'D$_{25}^{\rm min}$','$\theta$','D']+'} & ',['Ref.']+'} \\ ']

    colunits = '\colhead{'+[ [$
      '(1)','(2)','(3)','(4)',$
      '(5)','(6)','(7)','(8)',$
      '(9)','(10)','(11)','(12)']+'} & ',['(13)']+'} ']

    table = html_structure_parse(atlas,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata')
    ntags = n_tags(table)

; before initializing the MRT, change the distance references to
; numbers

    table[where(strmatch(table.distance_texref,'*mould*',/fold))].distance_texref = '1'
    table[where(strmatch(table.distance_texref,'*karachentsev*',/fold))].distance_texref = '2'
    table[where(strmatch(table.distance_texref,'*shapley*',/fold))].distance_texref = '3'
    table[where(strmatch(table.distance_texref,'*freedman*',/fold))].distance_texref = '4'
    table[where(strmatch(table.distance_texref,'*tully*',/fold))].distance_texref = '5'
    table[where(strmatch(table.distance_texref,'*tonry*',/fold))].distance_texref = '6'
    
    table_mrt = replicate({atlas_id: strarr(1), nice_galaxy: strarr(1), $
      alt_galaxy: strarr(1), rah: strarr(1), ram: strarr(1), ras: strarr(1), $
      ra_note: strarr(1), decsign: strarr(1), decd: strarr(1), decm: strarr(1), decs: strarr(1), $
      dec_note: strarr(1), cz: strarr(1), ebv_mw: strarr(1), lit_type: strarr(1), $
      d25_maj: strarr(1), d25_maj_note: strarr(1), d25_min: strarr(1), d25_min_note: strarr(1), $
      posangle: strarr(1), posangle_note: strarr(1), distance: strarr(1), distance_texref: strarr(1)},ngalaxy)

    table_mrt.atlas_id    = table.atlas_id
    table_mrt.nice_galaxy = table.nice_galaxy
    table_mrt.alt_galaxy  = table.alt_galaxy
    table_mrt.rah         = strmid(table.ra,0,2)
    table_mrt.ram         = strmid(table.ra,3,2)
    table_mrt.ras         = strmid(table.ra,6,4)
    table_mrt.decsign     = strmid(table.dec,0,1)
    table_mrt.decd        = strmid(table.dec,1,2)
    table_mrt.decm        = strmid(table.dec,4,2)
    table_mrt.decs        = strmid(table.dec,7,2)
    table_mrt.cz          = table.cz
    table_mrt.ebv_mw      = table.ebv_mw
    table_mrt.lit_type    = table.lit_type
    table_mrt.d25_maj     = table.d25_maj
    table_mrt.d25_min     = table.d25_min
    table_mrt.posangle    = table.posangle
    table_mrt.distance    = table.distance
    table_mrt.distance_texref = table.distance_texref

; generate the object tags

;   table.nice_galaxy = reform('\object['+strtrim(atlas.ned_galaxy,2)+']{'+strtrim(table.nice_galaxy,2)+'}',1,ngalaxy)

; identify the origin of the position measurements

    simbad = where(strmatch(atlas.position_ref,'*simbad*',/fold),nleda)
    leda = where(strmatch(atlas.position_ref,'*leda*',/fold),nleda)

    table[simbad].ra = strtrim(table[simbad].ra,2)+'\tablenotemark{a}'
    table[simbad].dec = strtrim(table[simbad].dec,2)+'\tablenotemark{a}'
    table_mrt[simbad].ra_note = 'a'
    table_mrt[simbad].dec_note = 'a'

    table[leda].ra = strtrim(table[leda].ra,2)+'\tablenotemark{b}'
    table[leda].dec = strtrim(table[leda].dec,2)+'\tablenotemark{b}'
    table_mrt[leda].ra_note = 'b'
    table_mrt[leda].dec_note = 'b'

; identify the origin of the major/minor axis measurements

    rc3 = where(strmatch(atlas.d25_origin,'*rc3*',/fold),nrc3)
    leda = where(strmatch(atlas.d25_origin,'*leda*',/fold),nleda)
    ned = where(strmatch(atlas.d25_origin,'*ned*',/fold),nned)
    ugc = where(strmatch(atlas.d25_origin,'*ugc*',/fold),nugc)
    eso = where(strmatch(atlas.d25_origin,'*eso*',/fold),neso)
    mcg = where(strmatch(atlas.d25_origin,'*mcg*',/fold),nmcg)
    twomass = where(strmatch(atlas.d25_origin,'*2mass*',/fold),ntwomass)
;   print, nrc3, nleda, nugc, nmcg, neso, nned, ntwomass

    table[leda].d25_maj = strtrim(table[leda].d25_maj,2)+'\tablenotemark{b}'
    table[leda].d25_min = strtrim(table[leda].d25_min,2)+'\tablenotemark{b}'
    table_mrt[leda].d25_maj_note = 'b'
    table_mrt[leda].d25_min_note = 'b'

    table[ned].d25_maj = strtrim(table[ned].d25_maj,2)+'\tablenotemark{c}'
    table[ned].d25_min = strtrim(table[ned].d25_min,2)+'\tablenotemark{c}'
    table_mrt[ned].d25_maj_note = 'c'
    table_mrt[ned].d25_min_note = 'c'

    table[ugc].d25_maj = strtrim(table[ugc].d25_maj,2)+'\tablenotemark{d}'
    table[ugc].d25_min = strtrim(table[ugc].d25_min,2)+'\tablenotemark{d}'
    table_mrt[ugc].d25_maj_note = 'd'
    table_mrt[ugc].d25_min_note = 'd'

    table[eso].d25_maj = strtrim(table[eso].d25_maj,2)+'\tablenotemark{e}'
    table[eso].d25_min = strtrim(table[eso].d25_min,2)+'\tablenotemark{e}'
    table_mrt[eso].d25_maj_note = 'e'
    table_mrt[eso].d25_min_note = 'e'

    table[mcg].d25_maj = strtrim(table[mcg].d25_maj,2)+'\tablenotemark{f}'
    table[mcg].d25_min = strtrim(table[mcg].d25_min,2)+'\tablenotemark{f}'
    table_mrt[mcg].d25_maj_note = 'f'
    table_mrt[mcg].d25_min_note = 'f'

; identify the origin of the position angle measurements

    rc3 = where(strmatch(atlas.posangle_origin,'*rc3*',/fold),nrc3)
    leda = where(strmatch(atlas.posangle_origin,'*leda*',/fold),nleda)
    ned = where(strmatch(atlas.posangle_origin,'*ned*',/fold),nned)
    ugc = where(strmatch(atlas.posangle_origin,'*ugc*',/fold),nugc)
    eso = where(strmatch(atlas.posangle_origin,'*eso*',/fold),neso)
    mcg = where(strmatch(atlas.posangle_origin,'*mcg*',/fold),nmcg)
    twomass_lga = where(strmatch(atlas.posangle_origin,'*2mass/lga*',/fold),ntwomass_lga)
    twomass_xsc = where(strmatch(atlas.posangle_origin,'*2mass/xsc*',/fold),ntwomass_xsc)
;   print, nrc3, nleda, nugc, nmcg, neso, nned, ntwomass_lga, ntwomass_xsc

    table[leda].posangle = strtrim(table[leda].posangle,2)+'\tablenotemark{b}'
    table_mrt[leda].posangle_note = 'b'

    table[ugc].posangle = strtrim(table[ugc].posangle,2)+'\tablenotemark{d}'
    table_mrt[ugc].posangle_note = 'd'

    table[eso].posangle = strtrim(table[eso].posangle,2)+'\tablenotemark{e}'
    table_mrt[eso].posangle_note = 'e'

    table[twomass_lga].posangle = strtrim(table[twomass_lga].posangle,2)+'\tablenotemark{g}'
    table_mrt[twomass_lga].posangle_note = 'g'

    table[twomass_xsc].posangle = strtrim(table[twomass_xsc].posangle,2)+'\tablenotemark{h}'
    table_mrt[twomass_xsc].posangle_note = 'h'

    tablecomments = [$
      'Col. (1) ID number;',$
      'Col. (2) Galaxy name;',$
      'Col. (3) Other common galaxy name or names;',$
      'Col. (4) Right ascension from NED, unless otherwise noted (J2000);',$
      'Col. (5) Declination from NED, unless otherwise noted (J2000);',$
      'Col. (6) Heliocentric redshift from NED (km~s$^{-1}$);',$
      'Col. (7) Foreground Galactic reddening from \citet{schlegel98} (mag);',$
      'Col. (8) Morphological type (see \S\ref{sec:properties});',$
      'Col. (9) Major axis diameter at the $25$~mag~arcsec$^{-2}$ isophote from the RC3, unless otherwise noted (arcmin);',$
      'Col. (10) Minor axis diameter at the $25$~mag~arcsec$^{-2}$ isophote from the RC3, unless otherwise noted (arcmin);',$
      'Col. (11) Position angle, measured positive from North to East from the RC3, unless otherwise noted (deg); ',$
      'Col. (12) Distance (Mpc); ',$
      'Col. (13) Distance reference.']
    tablecomments = strjoin(tablecomments,' ')
    tablenotetext = [$
      '{a}{Coordinates taken from SIMBAD.}',$
      '{b}{Coordinates, major/minor axis diameters, and/or position angle taken from LEDA.}',$
      '{c}{Major/minor axis diameter from NED "Basic Data".}',$
      '{d}{Major/minor axis diameter and/or position angle from the Uppsala '+$
      'General Catalog of Galaxies \citep[UGC;][]{nilson73}.}',$ 
      '{e}{Major/minor axis diameter and/or position angle from the ESO/Uppsala '+$
      'Survey of the ESO(B) Atlas \citep{lauberts82}.}',$
      '{f}{Major/minor axis diameter and/or position angle from the Morphological '+$
      'Catalogue of Galaxies \citep[MCG;][]{vorontsov62}.}',$
      '{g}{Position angle from the 2MASS Large Galaxy Atlas \citep{jarrett03}.}',$
      '{h}{Position angle from the 2MASS Extended Source Catalog \citep{jarrett00}.}']

;   tablenotetext = strjoin(tablenotetext,' ')

; write out the general properties table

    texfile = 'test.tex'
;   texfile = 'general_properties_table.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{General Properties \label{table:general_properties}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
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
          
          if (strcompress(table[i].(j),/remove) eq '') then $
            line[j] = '\nodata'+suffix else line[j] = table[i].(j)+suffix

       endfor
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\tablerefs{(1) \citet{mould00}; (2) \citet{karachentsev04}; '+$
      '(3) \citet{shapley01b}; (4) \citet{freedman01}; (5) \citet{tully96}; (6) \citet{tonry01}}'
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

stop    
    
; write out a stub of the general properties table 

    texfile = 'general_properties_table_stub.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{General Properties \label{table:general_properties}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0L, nstub-1L do begin

       line = strarr(ntags)
       for j = 0L, ntags-1L do begin

          if (j lt ntags-1L) then suffix = ' & ' else begin
             if (i lt ngalaxy-1L) then suffix = ' \\ ' else suffix = ''
          endelse

          if (strcompress(table[i].(j),/remove) eq '') then $
            line[j] = '\nodata'+suffix else line[j] = table[i].(j)+suffix

       endfor
       printf, lun, line
    endfor

    stubtext = 'Table \ref{table:general_properties} is published in its entirety in the electronic '+$
      'edition of the {\it Astrophysical Journal}.  A portion is shown '+$
      'here for guidance regarding its form and content.'
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+strjoin([tablecomments,stubtext],' ')+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\tablerefs{(1) \citet{mould00}; (2) \citet{karachentsev04}; '+$
      '(3) \citet{shapley01b}; (4) \citet{freedman01}; (5) \citet{tully96}; (6) \citet{tonry01}}'
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; write out the MRT data, labels, units, and explanation
 
    labels = ['ID','Galaxy','Altgalaxy','RAh','RAm','RAs','RANote','DE-','DEd','DEm','DEs','DENote',$
      'cz','EBV','MType','MajDiam','n_MajDiam','MinDiam','n_MinDiam','PA','n_PA','Dist','r_Dist']
    units = ['---','---','---','h','min','s','---','---','deg','arcmin','arcsec','---','km/s','mag','---',$
      'arcmin','---','arcmin','---','deg','---','Mpc','---']
    explanations = ['Unique identification number','Galaxy name','? Alternate galaxy name',$
      'Right Ascension (J2000)','Right Ascension (J2000)','Right Ascension (J2000)','? [ab] Note on Right Ascension position (1)',$
      'Sign of the Declination (J2000)','Declination (J2000)','Declination (J2000)',$
      'Declination (J2000)','? [ab] Note on Declination position (1)',$
      'Heliocentric redshift (2)','Galactic reddening (3)','Morphological type (4)','? Major axis diameter (5)',$
      '? [bcdef] Note on major axis diameter (5)','? Minor axis diameter (5)','? [bcdef] Note on minor axis diameter (5)',$
      '? Position angle (6)','? [bdegh] Note on position angle (6)','Distance (7)','Distance reference (7)']
    niceprint, labels, units, explanations
 
    openw, lun, mrtpath+'mrt_general_properties_labels.dat', /get_lun
    for i = 0L, n_elements(labels)-1L do printf, lun, labels[i]
    free_lun, lun
 
    openw, lun, mrtpath+'mrt_general_properties_units.dat', /get_lun
    for i = 0L, n_elements(units)-1L do printf, lun, units[i]
    free_lun, lun
 
    openw, lun, mrtpath+'mrt_general_properties_explanations.dat', /get_lun
    for i = 0L, n_elements(explanations)-1L do printf, lun, explanations[i]
    free_lun, lun
 
    ntags = n_tags(table_mrt)
    openw, lun, mrtpath+'mrt_general_properties_data.dat', /get_lun
    for i = 0L, ngalaxy-1L do begin
       for j = 0L, ntags-1L do begin
          if strmatch(table_mrt[i].(j),'*nodata*') then table_mrt[i].(j) = ''
          if (j eq 0L) then line = strtrim(table_mrt[i].(j),2) else $
            line = [line,strtrim(table_mrt[i].(j),2)]
       endfor
       printf, lun, strjoin(line,' & ')
    endfor
    free_lun, lun

; ---------------------------------------------------------------------------    
; Table: nuclear emission-line equivalent widths
; ---------------------------------------------------------------------------    

    select = [$
      'atlas_id','oii_3727_ew','h_delta_ew','h_gamma_ew',$
      'h_beta_ew','oiii_5007_ew','oi_6300_ew','h_alpha_ew',$
      'nii_6584_ew','sii_6716_ew','sii_6731_ew']
    format = [$
      'I3.3','F15.7','F15.7','F15.7',$
      'F15.7','F15.7','F15.7','F15.7',$
      'F15.7','F15.7','F15.7']
    texcenter = [$
      'l','c','c','c',$
      'c','c','c','c',$
      'c','c','c']
    colhead = '\colhead{'+[ [$
      'ID','[O~{\sc ii}]$~\lambda3727$','H$\delta~\lambda4101$','H$\gamma~\lambda4340$',$
      'H$\beta~\lambda4861$','[O~{\sc iii}]$~\lambda5007$','[O~{\sc i}]$~\lambda6300$','H$\alpha~\lambda6563$',$
      '[N~{\sc ii}]$~\lambda6584$','[S~{\sc ii}]$~\lambda6716$']+'} & ',['[S~{\sc ii}]$~\lambda6731$']+'}']
;   tablecomments = ''
    tablenotetext = '{a}{Nuclear rest-frame emission-line equivalent widths in Angstroms, '+$
      'corrected for underlying stellar absorption as described in \S~\ref{sec:algorithm}.  '+$
      'We give $1\sigma$ upper limits assuming two significant figures and identify them using a $<$ sign.  '+$
      'The errors only include statistical measurement uncertainties.  '+$
      'We do not give equivalent widths for the following objects because their nuclear spectra '+$
      'cannot be modeled reliably using the algorithm described in \S~\ref{sec:fitting}: '+$
      'NGC~1068 (054), NGC~1275 (070), NGC~3998 (214), NGC~4051 (223), and MRK~0315 (382).}'

    table = html_structure_parse(nuc,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error,/keep_upperlimits)
    ntags = n_tags(table)

    table_mrt = replicate({atlas_id: strarr(1)},nnuc)
    table_mrt.atlas_id = table.atlas_id
    linetags = tag_names(struct_trimtags(table,except='atlas_id'))
    for ii = 0L, n_elements(linetags)-1L do begin
       moretags_mrt = mrd_struct(linetags[ii]+['_limit','','_err'],replicate('strarr(1)',3),nnuc)
       table_mrt = struct_addtags(table_mrt,moretags_mrt)
    endfor

    table_ascii = table_mrt
    
; write out the EW table

    texfile = 'lineEW_nuc_table.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{Nuclear Emission-line Equivalent Widths\tablenotemark{a} \label{table:nuc_lineEW}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
    
    for i = 0L, nnuc-1L do begin

       line = strarr(ntags)
       for j = 0L, ntags-1L do begin

          if (j eq 0L) then line[j] = table[i].(j)+' & ' else begin ; GALAXY tag
             
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '

             if (strtrim(table[i].(j)[1],2) eq -3.0) then begin
                value = string(fix_digits(float(table[i].(j)[0]),2),format='(G0.0)')
                line[j] = '$<'+value+'$'+suffix
                table_mrt[i].(3*(j-1)+1+0) = '<'   ; offset from the ID tag
                table_mrt[i].(3*(j-1)+1+1) = value
                table_ascii[i].(3*(j-1)+1+0) = '<'   ; offset from the ID tag
                table_ascii[i].(3*(j-1)+1+1) = value
             endif

             if (strtrim(table[i].(j)[1],2) eq -2.0) then line[j] = '\nodata'+suffix

             if (strtrim(table[i].(j)[1],2) gt 0.0) then begin

                flux = table[i].(j)[0]
                ferr = table[i].(j)[1]

                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+suffix

                table_mrt[i].(3*(j-1)+1+1) = newflux
                table_mrt[i].(3*(j-1)+1+2) = newferr

                table_ascii[i].(3*(j-1)+1+0) = '+'
                table_ascii[i].(3*(j-1)+1+1) = newflux
                table_ascii[i].(3*(j-1)+1+2) = newferr

             endif

          endelse
             
       endfor
       printf, lun, line

    endfor
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    printf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; write out the ASCII table
;
;   openw, lun, datapath+repstr(texfile,'.tex','.ascii'), /get_lun
;   struct_print, table_ascii, lun=lun, /no_head;, $
;;    format=transpose([[tag_names(table_ascii[0])],[replicate('A5',n_tags(table_ascii[0]))]])
;   free_lun, lun
    
; write out the stub of the EW table

    texfile = 'lineEW_nuc_table_stub.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{Nuclear Emission-line Equivalent Widths\tablenotemark{a} \label{table:nuc_lineEW}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
    
    for i = 0L, nstub-1L do begin

       line = strarr(ntags)
       for j = 0L, ntags-1L do begin

          if (j eq 0L) then line[j] = table[i].(j)+' & ' else begin ; GALAXY tag
             
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '

             if (strtrim(table[i].(j)[1],2) eq -3.0) then $
               line[j] = '$<'+string(fix_digits(float(table[i].(j)[0]),2),format='(G0.0)')+'$'+suffix

             if (strtrim(table[i].(j)[1],2) eq -2.0) then line[j] = '\nodata'+suffix
             if (strtrim(table[i].(j)[1],2) gt 0.0) then begin

                flux = table[i].(j)[0]
                ferr = table[i].(j)[1]

                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+suffix

             endif

          endelse
             
       endfor
       printf, lun, line

    endfor
    
    stubtext = 'Table \ref{table:nuc_lineEW} is published in its entirety in the electronic '+$
      'edition of the {\it Astrophysical Journal}.  A portion is shown '+$
      'here for guidance regarding its form and content.'
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+stubtext+'}'
    printf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; write out the MRT data, labels, units, and explanation

    labels = ['ID',$
      'l_EWOII3727' ,'EWOII3727' ,'e_EWOII3727' ,'l_EWHDELTA' ,'EWHDELTA' ,'e_EWHDELTA',$
      'l_EWHGAMMA'  ,'EWHGAMMA'  ,'e_EWHGAMMA'  ,'l_EWHBETA'  ,'EWHBETA'  ,'e_EWHBETA',$
      'l_EWOIII5007','EWOIII5007','e_EWOIII5007','l_EWOI6300' ,'EWOI6300' ,'e_EWOI6300',$
      'l_EWHALPHA'  ,'EWHALPHA'  ,'e_EWHALPHA'  ,'l_EWNII6584','EWNII6584','e_EWNII6584',$
      'l_EWSII6716' ,'EWSII6716' ,'e_EWSII6716' ,'l_EWSII6731','EWSII6731','e_EWSII6731']
    units = ['---','---','0.1nm','0.1nm','---','0.1nm','0.1nm',$
      '---','0.1nm','0.1nm','---','0.1nm','0.1nm',$
      '---','0.1nm','0.1nm','---','0.1nm','0.1nm',$
      '---','0.1nm','0.1nm','---','0.1nm','0.1nm',$
      '---','0.1nm','0.1nm','---','0.1nm','0.1nm']
    explanations = ['Unique identification number (1)',$
      '? Upper limit (3sigma) flag on EW([O II] 3727)', '? [O II] 3727 equivalent width', '? [O II] 3727 equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW(H-DELTA)',     '? H-DELTA equivalent width',     '? H-DELTA equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW(H-GAMMA)',     '? H-GAMMA equivalent width',     '? H-GAMMA equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW(H-BETA)',      '? H-BETA equivalent width',      '? H-BETA equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW([O III] 5007)','? [O III] 5007 equivalent width','? [O III] 5007 equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW([O I] 6300)',  '? [O I] 6300 equivalent width',  '? [O I] 6300 equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW(H-ALPHA)',     '? H-ALPHA equivalent width',     '? H-ALPHA equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW([N II] 6584)', '? [N II] 6584 equivalent width', '? [N II] 6584 equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW([S II] 6716)', '? [S II] 6716 equivalent width', '? [S II] 6716 equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW([S II] 6731)', '? [S II] 6731 equivalent width', '? [S II] 6731 equivalent width uncertainty (1sigma)']
;   niceprint, labels, units, explanations

    openw, lun, mrtpath+'mrt_lineEW_nuc_labels.dat', /get_lun
    for i = 0L, n_elements(labels)-1L do printf, lun, labels[i]
    free_lun, lun

    openw, lun, mrtpath+'mrt_lineEW_nuc_units.dat', /get_lun
    for i = 0L, n_elements(units)-1L do printf, lun, units[i]
    free_lun, lun

    openw, lun, mrtpath+'mrt_lineEW_nuc_explanations.dat', /get_lun
    for i = 0L, n_elements(explanations)-1L do printf, lun, explanations[i]
    free_lun, lun

    ntags = n_tags(table_mrt)
    openw, lun, mrtpath+'mrt_lineEW_nuc_data.dat', /get_lun
    for i = 0L, nnuc-1L do begin
       for j = 0L, ntags-1L do begin
          if strmatch(table_mrt[i].(j),'*nodata*') then table_mrt[i].(j) = ''
          if (j eq 0L) then line = strtrim(table_mrt[i].(j),2) else $
            line = [line,strtrim(table_mrt[i].(j),2)]
       endfor
       printf, lun, strjoin(line,' & ')
    endfor
    free_lun, lun

; ---------------------------------------------------------------------------    
; Table: integrated emission-line equivalent widths
; ---------------------------------------------------------------------------    

    select = [$
      'atlas_id','oii_3727_ew','h_delta_ew','h_gamma_ew',$
      'h_beta_ew','oiii_5007_ew','oi_6300_ew','h_alpha_ew',$
      'nii_6584_ew','sii_6716_ew','sii_6731_ew']

    format = [$
      'I3.3','F15.7','F15.7','F15.7',$
      'F15.7','F15.7','F15.7','F15.7',$
      'F15.7','F15.7','F15.7']
    texcenter = [$
      'l','c','c','c',$
      'c','c','c','c',$
      'c','c','c']
    colhead = '\colhead{'+[ [$
      'ID','[O~{\sc ii}]$~\lambda3727$','H$\delta~\lambda4101$','H$\gamma~\lambda4340$',$
      'H$\beta~\lambda4861$','[O~{\sc iii}]$~\lambda5007$','[O~{\sc i}]$~\lambda6300$','H$\alpha~\lambda6563$',$
      '[N~{\sc ii}]$~\lambda6584$','[S~{\sc ii}]$~\lambda6716$']+'} & ',['[S~{\sc ii}]$~\lambda6731$']+'}']
;   tablecomments = ''
    tablenotetext = '{a}{Integrated rest-frame emission-line equivalent widths in Angstroms, '+$
      'corrected for underlying stellar absorption as described in \S~\ref{sec:algorithm}.  '+$
      'We give $1\sigma$ upper limits assuming two significant figures and identify them using a $<$ sign.  '+$
      'The errors only include statistical measurement uncertainties.  '+$
      'We do not give equivalent widths for the following objects because their integrated spectra '+$
      'cannot be modeled reliably using the algorithm described in \S~\ref{sec:fitting}: '+$
      'NGC~1275 (070), IRAS~05189-2524 (087), UGC~08058 (279), NGC~7469 (381),  MRK~0315 (382), NGC~7674 (400), and ARP~182 (401).}'

    table = html_structure_parse(int,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error,/keep_upperlimits)
    ntags = n_tags(table)

    table_mrt = replicate({atlas_id: strarr(1)},nint)
    table_mrt.atlas_id = table.atlas_id
    linetags = tag_names(struct_trimtags(table,except='atlas_id'))
    for ii = 0L, n_elements(linetags)-1L do begin
       moretags_mrt = mrd_struct(linetags[ii]+['_limit','','_err'],replicate('strarr(1)',3),nint)
       table_mrt = struct_addtags(table_mrt,moretags_mrt)
    endfor
    
; write out the EW table

    texfile = 'lineEW_int_table.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{Integrated Emission-line Equivalent Widths\tablenotemark{a} \label{table:int_lineEW}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
    
    for i = 0L, nint-1L do begin

       line = strarr(ntags)
       for j = 0L, ntags-1L do begin

          if (j eq 0L) then line[j] = table[i].(j)+' & ' else begin ; GALAXY tag
             
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '

; how should we treat upper limits?  how many decimal points?             
             
             if (strtrim(table[i].(j)[1],2) eq -3.0) then begin
                value = string(fix_digits(float(table[i].(j)[0]),2),format='(G0.0)')
                line[j] = '$<'+value+'$'+suffix
                table_mrt[i].(3*(j-1)+1+0) = '<'   ; offset from the ID tag
                table_mrt[i].(3*(j-1)+1+1) = value
             endif

             if (strtrim(table[i].(j)[1],2) eq -2.0) then line[j] = '\nodata'+suffix

             if (strtrim(table[i].(j)[1],2) gt 0.0) then begin

                flux = table[i].(j)[0]
                ferr = table[i].(j)[1]

                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+suffix

                table_mrt[i].(3*(j-1)+1+1) = newflux
                table_mrt[i].(3*(j-1)+1+2) = newferr

             endif

          endelse
             
       endfor
       printf, lun, line

    endfor
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    printf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; write out a stub of the EW table

    texfile = 'lineEW_int_table_stub.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{Integrated Emission-line Equivalent Widths\tablenotemark{a} \label{table:int_lineEW}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
    
    for i = 0L, nstub-1L do begin

       line = strarr(ntags)
       for j = 0L, ntags-1L do begin

          if (j eq 0L) then line[j] = table[i].(j)+' & ' else begin ; GALAXY tag
             
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '

; how should we treat upper limits?  how many decimal points?             
             
             if (strtrim(table[i].(j)[1],2) eq -3.0) then $
               line[j] = '$<'+string(fix_digits(float(table[i].(j)[0]),2),format='(G0.0)')+'$'+suffix

             if (strtrim(table[i].(j)[1],2) eq -2.0) then line[j] = '\nodata'+suffix

             if (strtrim(table[i].(j)[1],2) gt 0.0) then begin

                flux = table[i].(j)[0]
                ferr = table[i].(j)[1]

                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+suffix

             endif

          endelse
             
       endfor
       printf, lun, line

    endfor
    
    stubtext = 'Table \ref{table:int_lineEW} is published in its entirety in the electronic '+$
      'edition of the {\it Astrophysical Journal}.  A portion is shown '+$
      'here for guidance regarding its form and content.'
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+stubtext+'}'
    printf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; write out the MRT data, labels, units, and explanation

    labels = ['ID',$
      'l_EWOII3727' ,'EWOII3727' ,'e_EWOII3727' ,'l_EWHDELTA' ,'EWHDELTA' ,'e_EWHDELTA',$
      'l_EWHGAMMA'  ,'EWHGAMMA'  ,'e_EWHGAMMA'  ,'l_EWHBETA'  ,'EWHBETA'  ,'e_EWHBETA',$
      'l_EWOIII5007','EWOIII5007','e_EWOIII5007','l_EWOI6300' ,'EWOI6300' ,'e_EWOI6300',$
      'l_EWHALPHA'  ,'EWHALPHA'  ,'e_EWHALPHA'  ,'l_EWNII6584','EWNII6584','e_EWNII6584',$
      'l_EWSII6716' ,'EWSII6716' ,'e_EWSII6716' ,'l_EWSII6731','EWSII6731','e_EWSII6731']
    units = ['---','---','0.1nm','0.1nm','---','0.1nm','0.1nm',$
      '---','0.1nm','0.1nm','---','0.1nm','0.1nm',$
      '---','0.1nm','0.1nm','---','0.1nm','0.1nm',$
      '---','0.1nm','0.1nm','---','0.1nm','0.1nm',$
      '---','0.1nm','0.1nm','---','0.1nm','0.1nm']
    explanations = ['Unique identification number (1)',$
      '? Upper limit (3sigma) flag on EW([O II] 3727)', '? [O II] 3727 equivalent width', '? [O II] 3727 equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW(H-DELTA)',     '? H-DELTA equivalent width',     '? H-DELTA equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW(H-GAMMA)',     '? H-GAMMA equivalent width',     '? H-GAMMA equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW(H-BETA)',      '? H-BETA equivalent width',      '? H-BETA equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW([O III] 5007)','? [O III] 5007 equivalent width','? [O III] 5007 equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW([O I] 6300)',  '? [O I] 6300 equivalent width',  '? [O I] 6300 equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW(H-ALPHA)',     '? H-ALPHA equivalent width',     '? H-ALPHA equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW([N II] 6584)', '? [N II] 6584 equivalent width', '? [N II] 6584 equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW([S II] 6716)', '? [S II] 6716 equivalent width', '? [S II] 6716 equivalent width uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on EW([S II] 6731)', '? [S II] 6731 equivalent width', '? [S II] 6731 equivalent width uncertainty (1sigma)']
;   niceprint, labels, units, explanations

    openw, lun, mrtpath+'mrt_lineEW_int_labels.dat', /get_lun
    for i = 0L, n_elements(labels)-1L do printf, lun, labels[i]
    free_lun, lun

    openw, lun, mrtpath+'mrt_lineEW_int_units.dat', /get_lun
    for i = 0L, n_elements(units)-1L do printf, lun, units[i]
    free_lun, lun

    openw, lun, mrtpath+'mrt_lineEW_int_explanations.dat', /get_lun
    for i = 0L, n_elements(explanations)-1L do printf, lun, explanations[i]
    free_lun, lun

    ntags = n_tags(table_mrt)
    openw, lun, mrtpath+'mrt_lineEW_int_data.dat', /get_lun
    for i = 0L, nint-1L do begin
       for j = 0L, ntags-1L do begin
          if strmatch(table_mrt[i].(j),'*nodata*') then table_mrt[i].(j) = ''
          if (j eq 0L) then line = strtrim(table_mrt[i].(j),2) else $
            line = [line,strtrim(table_mrt[i].(j),2)]
       endfor
       printf, lun, strjoin(line,' & ')
    endfor
    free_lun, lun

; ---------------------------------------------------------------------------    
; Table: nuclear emission-line fluxes
; ---------------------------------------------------------------------------    

    select = [$
      'atlas_id','oii_3727','h_delta','h_gamma',$
      'h_beta','oiii_5007','oi_6300','h_alpha',$
      'nii_6584','sii_6716','sii_6731']

; multiply the line fluxes and errors by SCALE    
    
    nucselect = struct_trimtags(nuc,select=select)
    for i = 0L, n_elements(nuc)-1L do begin 
       for iline = 1L, n_tags(nucselect)-1L do begin ; offset from ATLAS_ID
          flux = nucselect[i].(iline)[0]
          ferr = nucselect[i].(iline)[1]
          if (ferr eq -3.0) then flux = flux*scale ; upper limit
          if (ferr gt 0.0) then begin
             flux = flux*scale
             ferr = ferr*scale
             if (flux/ferr lt sigcut) then ferr = -3.0 ; here!!
          endif
          nucselect[i].(iline)[0] = flux
          nucselect[i].(iline)[1] = ferr
       endfor
    endfor

    format = [$
      'I3.3','F15.7','F15.7','F15.7',$
      'F15.7','F15.7','F15.7','F15.7',$
      'F15.7','F15.7','F15.7']
    texcenter = [$
      'l','c','c','c',$
      'c','c','c','c',$
      'c','c','c']
    colhead = '\colhead{'+[ [$
      'ID','[O~{\sc ii}]$~\lambda3727$','H$\delta~\lambda4101$','H$\gamma~\lambda4340$',$
      'H$\beta~\lambda4861$','[O~{\sc iii}]$~\lambda5007$','[O~{\sc i}]$~\lambda6300$','H$\alpha~\lambda6563$',$
      '[N~{\sc ii}]$~\lambda6584$','[S~{\sc ii}]$~\lambda6716$']+'} & ',['[S~{\sc ii}]$~\lambda6731$']+'}']
;   tablecomments = ''
    tablenotetext = '{a}{Nuclear emission-line fluxes in units of $10^{-15}$~erg~s$^{-1}$~cm$^{-2}$.  '+$
      'We give $1\sigma$ upper limits assuming two significant figures and identify them using a $<$ sign.  '+$
      'These flux measurements have been '+$
      'corrected for foreground Galactic extinction using the reddening values in Table~\ref{table:general_properties} '+$
      'and the \citet{odonnell94} Milky Way extinction curve assuming $R_{\rm V}=3.1$, as well as '+$
      'underlying stellar absorption as '+$
      'described in \S~\ref{sec:algorithm}.  The errors only include statistical measurement uncertainties.  '+$
      'We do not give fluxes for the following objects because their nuclear spectra '+$
      'cannot be modeled reliably using the algorithm described in \S~\ref{sec:fitting}: '+$
      'NGC~1068 (054), NGC~1275 (070), NGC~3998 (214), NGC~4051 (223), and MRK~0315 (382).}'

    table = html_structure_parse(nucselect,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error,/keep_upperlimits)
    ntags = n_tags(table)

    table_mrt = replicate({atlas_id: strarr(1)},nnuc)
    table_mrt.atlas_id = table.atlas_id
    linetags = tag_names(struct_trimtags(table,except='atlas_id'))
    for ii = 0L, n_elements(linetags)-1L do begin
       moretags_mrt = mrd_struct(linetags[ii]+['_limit','','_err'],replicate('strarr(1)',3),nnuc)
       table_mrt = struct_addtags(table_mrt,moretags_mrt)
    endfor
    
; write out the line flux table

    texfile = 'lineflux_nuc_table.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{Nuclear Emission-line Fluxes\tablenotemark{a} \label{table:nuc_lineflux}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
    
    for i = 0L, nnuc-1L do begin

       line = strarr(ntags)
       for j = 0L, ntags-1L do begin

          if (j eq 0L) then line[j] = table[i].(j)+' & ' else begin ; GALAXY tag
             
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '

; how should we treat upper limits?  how many decimal points?             
             
             if (strtrim(table[i].(j)[1],2) eq -3.0) then begin
                value = string(fix_digits(float(table[i].(j)[0]),2),format='(G0.0)')
                line[j] = '$<'+value+'$'+suffix
                table_mrt[i].(3*(j-1)+1+0) = '<'   ; offset from the ID tag
                table_mrt[i].(3*(j-1)+1+1) = value
             endif

             if (strtrim(table[i].(j)[1],2) eq -2.0) then line[j] = '\nodata'+suffix

             if (strtrim(table[i].(j)[1],2) gt 0.0) then begin

                flux = table[i].(j)[0]
                ferr = table[i].(j)[1]

                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+suffix

                table_mrt[i].(3*(j-1)+1+1) = newflux
                table_mrt[i].(3*(j-1)+1+2) = newferr

             endif

          endelse
             
       endfor
       printf, lun, line

    endfor
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    printf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; write out the stub of the line flux table

    texfile = 'lineflux_nuc_table_stub.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{Nuclear Emission-line Fluxes\tablenotemark{a} \label{table:nuc_lineflux}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
    
    for i = 0L, nstub-1L do begin

       line = strarr(ntags)
       for j = 0L, ntags-1L do begin

          if (j eq 0L) then line[j] = table[i].(j)+' & ' else begin ; GALAXY tag
             
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '

; how should we treat upper limits?  how many decimal points?             
             
             if (strtrim(table[i].(j)[1],2) eq -3.0) then $
               line[j] = '$<'+string(fix_digits(float(table[i].(j)[0]),2),format='(G0.0)')+'$'+suffix

             if (strtrim(table[i].(j)[1],2) eq -2.0) then line[j] = '\nodata'+suffix

             if (strtrim(table[i].(j)[1],2) gt 0.0) then begin

                flux = table[i].(j)[0]
                ferr = table[i].(j)[1]

                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+suffix

             endif

          endelse
             
       endfor
       printf, lun, line

    endfor
    
    stubtext = 'Table \ref{table:nuc_lineflux} is published in its entirety in the electronic '+$
      'edition of the {\it Astrophysical Journal}.  A portion is shown '+$
      'here for guidance regarding its form and content.'
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+stubtext+'}'
    printf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; write out the MRT data, labels, units, and explanation

    labels = ['ID',$
      'l_OII3727','OII3727','e_OII3727','l_HDELTA','HDELTA','e_HDELTA',$
      'l_HGAMMA','HGAMMA','e_HGAMMA','l_HBETA','HBETA','e_HBETA',$
      'l_OIII5007','OIII5007','e_OIII5007','l_OI6300','OI6300','e_OI6300',$
      'l_HALPHA','HALPHA','e_HALPHA','l_NII6584','NII6584','e_NII6584',$
      'l_SII6716','SII6716','e_SII6716','l_SII6731','SII6731','e_SII6731']
    units = ['---','---','10-15mW/m2','10-15mW/m2','---','10-15mW/m2','10-15mW/m2',$
      '---','10-15mW/m2','10-15mW/m2','---','10-15mW/m2','10-15mW/m2',$
      '---','10-15mW/m2','10-15mW/m2','---','10-15mW/m2','10-15mW/m2',$
      '---','10-15mW/m2','10-15mW/m2','---','10-15mW/m2','10-15mW/m2',$
      '---','10-15mW/m2','10-15mW/m2','---','10-15mW/m2','10-15mW/m2']
    explanations = ['Unique identification number (1)',$
      '? Upper limit (3sigma) flag on [O II] 3727', '? [O II] 3727 flux', '? [O II] 3727 flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on H-DELTA',     '? H-DELTA flux',     '? H-DELTA flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on H-GAMMA',     '? H-GAMMA flux',     '? H-GAMMA flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on H-BETA',      '? H-BETA flux',      '? H-BETA flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on [O III] 5007','? [O III] 5007 flux','? [O III] 5007 flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on [O I] 6300',  '? [O I] 6300 flux',  '? [O I] 6300 flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on H-ALPHA',     '? H-ALPHA flux',     '? H-ALPHA flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on [N II] 6584', '? [N II] 6584 flux', '? [N II] 6584 flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on [S II] 6716', '? [S II] 6716 flux', '? [S II] 6716 flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on [S II] 6731', '? [S II] 6731 flux', '? [S II] 6731 flux uncertainty (1sigma)']
;   niceprint, labels, units, explanations

    openw, lun, mrtpath+'mrt_lineflux_nuc_labels.dat', /get_lun
    for i = 0L, n_elements(labels)-1L do printf, lun, labels[i]
    free_lun, lun

    openw, lun, mrtpath+'mrt_lineflux_nuc_units.dat', /get_lun
    for i = 0L, n_elements(units)-1L do printf, lun, units[i]
    free_lun, lun

    openw, lun, mrtpath+'mrt_lineflux_nuc_explanations.dat', /get_lun
    for i = 0L, n_elements(explanations)-1L do printf, lun, explanations[i]
    free_lun, lun

    ntags = n_tags(table_mrt)
    openw, lun, mrtpath+'mrt_lineflux_nuc_data.dat', /get_lun
    for i = 0L, nnuc-1L do begin
       for j = 0L, ntags-1L do begin
          if strmatch(table_mrt[i].(j),'*nodata*') then table_mrt[i].(j) = ''
          if (j eq 0L) then line = strtrim(table_mrt[i].(j),2) else $
            line = [line,strtrim(table_mrt[i].(j),2)]
       endfor
       printf, lun, strjoin(line,' & ')
    endfor
    free_lun, lun

; ---------------------------------------------------------------------------    
; Table: integrated emission-line fluxes
; ---------------------------------------------------------------------------    

    select = [$
      'atlas_id','oii_3727','h_delta','h_gamma',$
      'h_beta','oiii_5007','oi_6300','h_alpha',$
      'nii_6584','sii_6716','sii_6731']

; multiply the line fluxes and errors by SCALE    
    
    intselect = struct_trimtags(int,select=select)
    for i = 0L, n_elements(int)-1L do begin 
       for iline = 1L, n_tags(intselect)-1L do begin ; offset from ATLAS_ID
          flux = intselect[i].(iline)[0]
          ferr = intselect[i].(iline)[1]
          if (ferr eq -3.0) then flux = flux*scale ; upper limit
          if (ferr gt 0.0) then begin
             flux = flux*scale
             ferr = ferr*scale
             if (flux/ferr lt sigcut) then ferr = -3.0 ; here!!
          endif
          intselect[i].(iline)[0] = flux
          intselect[i].(iline)[1] = ferr
       endfor
    endfor

    format = [$
      'I3.3','F15.7','F15.7','F15.7',$
      'F15.7','F15.7','F15.7','F15.7',$
      'F15.7','F15.7','F15.7']
    texcenter = [$
      'l','c','c','c',$
      'c','c','c','c',$
      'c','c','c']
    colhead = '\colhead{'+[ [$
      'ID','[O~{\sc ii}]$~\lambda3727$','H$\delta~\lambda4101$','H$\gamma~\lambda4340$',$
      'H$\beta~\lambda4861$','[O~{\sc iii}]$~\lambda5007$','[O~{\sc i}]$~\lambda6300$','H$\alpha~\lambda6563$',$
      '[N~{\sc ii}]$~\lambda6584$','[S~{\sc ii}]$~\lambda6716$']+'} & ',['[S~{\sc ii}]$~\lambda6731$']+'}']
;   tablecomments = ''
    tablenotetext = '{a}{Integrated emission-line fluxes in units of $10^{-15}$~erg~s$^{-1}$~cm$^{-2}$.  '+$
      'We give $1\sigma$ upper limits assuming two significant figures and identify them using a $<$ sign.  '+$
      'These flux measurements have been '+$
      'corrected for foreground Galactic extinction using the reddening values in Table~\ref{table:general_properties} '+$
      'and the \citet{odonnell94} Milky Way extinction curve assuming $R_{\rm V}=3.1$, as well as '+$
      'underlying stellar absorption as '+$
      'described in \S~\ref{sec:algorithm}.  The errors only include statistical measurement uncertainties.  '+$
      'We do not give fluxes for the following objects because their integrated spectra '+$
      'cannot be modeled reliably using the algorithm described in \S~\ref{sec:fitting}: '+$
      'NGC~1275 (070), IRAS~05189-2524 (087), UGC~08058 (279), NGC~7469 (381),  MRK~0315 (382), NGC~7674 (400), and ARP~182 (401).}'
    table = html_structure_parse(intselect,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error,/keep_upperlimits)
    ntags = n_tags(table)

    table_mrt = replicate({atlas_id: strarr(1)},nint)
    table_mrt.atlas_id = table.atlas_id
    linetags = tag_names(struct_trimtags(table,except='atlas_id'))
    for ii = 0L, n_elements(linetags)-1L do begin
       moretags_mrt = mrd_struct(linetags[ii]+['_limit','','_err'],replicate('strarr(1)',3),nint)
       table_mrt = struct_addtags(table_mrt,moretags_mrt)
    endfor
    
; write out the line flux table

    texfile = 'lineflux_int_table.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{Integrated Emission-line Fluxes\tablenotemark{a} \label{table:int_lineflux}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
    
    for i = 0L, nint-1L do begin

       line = strarr(ntags)
       for j = 0L, ntags-1L do begin

          if (j eq 0L) then line[j] = table[i].(j)+' & ' else begin ; GALAXY tag
             
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '

; how should we treat upper limits?  how many decimal points?             
             
             if (strtrim(table[i].(j)[1],2) eq -3.0) then begin
                value = string(fix_digits(float(table[i].(j)[0]),2),format='(G0.0)')
                line[j] = '$<'+value+'$'+suffix
                table_mrt[i].(3*(j-1)+1+0) = '<'   ; offset from the ID tag
                table_mrt[i].(3*(j-1)+1+1) = value
             endif

             if (strtrim(table[i].(j)[1],2) eq -2.0) then line[j] = '\nodata'+suffix

             if (strtrim(table[i].(j)[1],2) gt 0.0) then begin

                flux = table[i].(j)[0]
                ferr = table[i].(j)[1]

                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+suffix

                table_mrt[i].(3*(j-1)+1+1) = newflux
                table_mrt[i].(3*(j-1)+1+2) = newferr

             endif

          endelse
             
       endfor
       printf, lun, line

    endfor
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    printf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; write out a stub of the line flux table

    texfile = 'lineflux_int_table_stub.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{Integrated Emission-line Fluxes\tablenotemark{a} \label{table:int_lineflux}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
;   niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
    
    for i = 0L, nstub-1L do begin

       line = strarr(ntags)
       for j = 0L, ntags-1L do begin

          if (j eq 0L) then line[j] = table[i].(j)+' & ' else begin ; GALAXY tag
             
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '

; how should we treat upper limits?  how many decimal points?             
             
             if (strtrim(table[i].(j)[1],2) eq -3.0) then $
               line[j] = '$<'+string(fix_digits(float(table[i].(j)[0]),2),format='(G0.0)')+'$'+suffix

             if (strtrim(table[i].(j)[1],2) eq -2.0) then line[j] = '\nodata'+suffix

             if (strtrim(table[i].(j)[1],2) gt 0.0) then begin

                flux = table[i].(j)[0]
                ferr = table[i].(j)[1]

                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+suffix

             endif

          endelse
             
       endfor
       printf, lun, line

    endfor
    
    stubtext = 'Table \ref{table:int_lineflux} is published in its entirety in the electronic '+$
      'edition of the {\it Astrophysical Journal}.  A portion is shown '+$
      'here for guidance regarding its form and content.'
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+stubtext+'}'
    printf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; write out the MRT data, labels, units, and explanation

    labels = ['ID',$
      'l_OII3727','OII3727','e_OII3727','l_HDELTA','HDELTA','e_HDELTA',$
      'l_HGAMMA','HGAMMA','e_HGAMMA','l_HBETA','HBETA','e_HBETA',$
      'l_OIII5007','OIII5007','e_OIII5007','l_OI6300','OI6300','e_OI6300',$
      'l_HALPHA','HALPHA','e_HALPHA','l_NII6584','NII6584','e_NII6584',$
      'l_SII6716','SII6716','e_SII6716','l_SII6731','SII6731','e_SII6731']
    units = ['---','---','10-15mW/m2','10-15mW/m2','---','10-15mW/m2','10-15mW/m2',$
      '---','10-15mW/m2','10-15mW/m2','---','10-15mW/m2','10-15mW/m2',$
      '---','10-15mW/m2','10-15mW/m2','---','10-15mW/m2','10-15mW/m2',$
      '---','10-15mW/m2','10-15mW/m2','---','10-15mW/m2','10-15mW/m2',$
      '---','10-15mW/m2','10-15mW/m2','---','10-15mW/m2','10-15mW/m2']
    explanations = ['Unique identification number (1)',$
      '? Upper limit (3sigma) flag on [O II] 3727', '? [O II] 3727 flux', '? [O II] 3727 flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on H-DELTA',     '? H-DELTA flux',     '? H-DELTA flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on H-GAMMA',     '? H-GAMMA flux',     '? H-GAMMA flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on H-BETA',      '? H-BETA flux',      '? H-BETA flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on [O III] 5007','? [O III] 5007 flux','? [O III] 5007 flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on [O I] 6300',  '? [O I] 6300 flux',  '? [O I] 6300 flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on H-ALPHA',     '? H-ALPHA flux',     '? H-ALPHA flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on [N II] 6584', '? [N II] 6584 flux', '? [N II] 6584 flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on [S II] 6716', '? [S II] 6716 flux', '? [S II] 6716 flux uncertainty (1sigma)',$
      '? Upper limit (3sigma) flag on [S II] 6731', '? [S II] 6731 flux', '? [S II] 6731 flux uncertainty (1sigma)']
;   niceprint, labels, units, explanations

    openw, lun, mrtpath+'mrt_lineflux_int_labels.dat', /get_lun
    for i = 0L, n_elements(labels)-1L do printf, lun, labels[i]
    free_lun, lun

    openw, lun, mrtpath+'mrt_lineflux_int_units.dat', /get_lun
    for i = 0L, n_elements(units)-1L do printf, lun, units[i]
    free_lun, lun

    openw, lun, mrtpath+'mrt_lineflux_int_explanations.dat', /get_lun
    for i = 0L, n_elements(explanations)-1L do printf, lun, explanations[i]
    free_lun, lun

    ntags = n_tags(table_mrt)
    openw, lun, mrtpath+'mrt_lineflux_int_data.dat', /get_lun
    for i = 0L, nint-1L do begin
       for j = 0L, ntags-1L do begin
          if strmatch(table_mrt[i].(j),'*nodata*') then table_mrt[i].(j) = ''
          if (j eq 0L) then line = strtrim(table_mrt[i].(j),2) else $
            line = [line,strtrim(table_mrt[i].(j),2)]
       endfor
       printf, lun, strjoin(line,' & ')
    endfor
    free_lun, lun

; ---------------------------------------------------------------------------    
; Table 2: photometric properties
; ---------------------------------------------------------------------------    
    
; TODO: perhaps add the FIR/IR flux and/or luminosity; perhaps add the
; FIR or IR/optical luminosity ratio
    
    select = [$
      'atlas_id','rc3_u','rc3_u_err','rc3_b',$
      'rc3_b_err','rc3_v','rc3_v_err','twomass_j',$
      'twomass_j_err','twomass_h','twomass_h_err','twomass_ks',$
      'twomass_ks_err','iras_12','iras_12_err','iras_25',$
      'iras_25_err','iras_60','iras_60_err','iras_100',$
      'iras_100_err']
    format = [$
      'I3.3','F12.3','F12.3','F12.3',$
      'F12.3','F12.3','F12.3','F12.3',$
      'F12.3','F12.3','F12.3','F12.3',$
      'F12.3','F12.3','F12.3','F12.3',$
      'F12.3','F12.3','F12.3','F12.3',$
      'F12.3']
    texcenter = ['l','c','c','c','c','c','c','c','c','c','c']
    note = '' ; '^{a}'
    colhead = '\colhead{'+[ [$
      'ID','U','B','V','J','H','K$_{\rm s}$','$S_{\nu}(12~\micron)'+note+'$',$
      '$S_{\nu}(25~\micron)'+note+'$','$S_{\nu}(60~\micron)'+note+'$']+$
      '} & ',['$S_{\nu}(100~\micron)'+note+'$']+'} \\ ']
    colunits = '\colhead{'+[ [$
      '(1)','(2)','(3)','(4)',$
      '(5)','(6)','(7)','(8)',$
      '(9)','(10)']+'} & ',['(11)']+'} ']
;   colunits = '\colhead{'+[ [$
;     ' ','(mag)','(mag)','(mag)','(mag)','(mag)','(mag)','(Jy)','(Jy)','(Jy)']+'} & ',['(Jy)']+'} ']
    tablecomments = [$
      'Col. (1) ID number;',$
      'Col. (2) U magnitude;',$
      'Col. (3) B magnitude;',$
      'Col. (4) V magnitude;',$
      'Col. (5) 2MASS total J magnitude;',$
      'Col. (6) 2MASS total H magnitude;',$
      'Col. (7) 2MASS total K$_{\mathrm s}$ magnitude;',$
      'Col. (8) IRAS $12~\micron$ flux density (Jy);',$
      'Col. (9) IRAS $25~\micron$ flux density (Jy);',$
      'Col. (10) IRAS $60~\micron$ flux density (Jy);',$
      'Col. (11) IRAS $100~\micron$ flux density (Jy).']
    tablecomments = strjoin(tablecomments,' ')

    table = html_structure_parse(atlas,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata')
    ntags = n_tags(table)

    ncols = 1 + (ntags-1L)/2L ; <-- NOT GENERAL IF TAGS ARE ADDED!

    moretags_mrt = replicate({atlas_id: strarr(1), atlas_id_note: strarr(1)},ngalaxy)
    table_mrt = struct_addtags(moretags_mrt,struct_trimtags(table,except=['atlas_id']))
    table_mrt.atlas_id = table.atlas_id
    
; M51 has a problem with its IRAS photometry, so add a tablenotemark    
    
    m51 = where(strmatch(strcompress(atlas.galaxy,/remove),'*ngc5194*',/fold) eq 1B,nm51)
    tablenotetext = '{a}{\citet{soifer89} do not provide measurement errors for the IRAS fluxes of NGC~5194, '+$
      'therefore we assume a fixed $25\%$ uncertainty at every wavelength.}'
    table[m51].atlas_id = table[m51].atlas_id+'\tablenotemark{a}'
    table_mrt[m51].atlas_id_note = 'a'
    
; write out the photometric properties table

    texfile = 'photometric_properties_table.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{Photometric Properties \label{table:photometric_properties}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
    
    for i = 0L, ngalaxy-1L do begin
       line = strarr(ncols)
       for j = 0L, ncols-1L do begin
          if (j lt ncols-1L) then suffix = ' & ' else suffix = ' \\ '

          if (j eq 0L) then line[j] = table[i].(j)+suffix else begin ; ID number

             flux = table[i].(2*j-1L)
             ferr = table[i].(2*j)

; well-measured flux and error             
             
             if (strmatch(flux,'*nodata*') eq 0B) and (strmatch(ferr,'*nodata*') eq 0B) then begin
                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+suffix
             endif 

; well-measured flux and but no error; THIS CASE MUST BE FIXED (NGC5194)
             
             if (strmatch(flux,'*nodata*') eq 0B) and (strmatch(ferr,'*nodata*') eq 1B) then begin
                line[j] = '$'+flux+'$'+suffix
                print, 'Hey! '+atlas[i].galaxy
             endif

; no flux measurement (assume no error)
             
             if (strmatch(flux,'*nodata*') eq 1B) then line[j] = '\nodata'+suffix

; flux upper limit (no error)
             
             if (strmatch(flux,'*-*') eq 1B) then line[j] = '$<'+string(-float(flux),format='(G0.0)')+'$'+suffix

          endelse

       endfor
          
       printf, lun, line

    endfor
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
    printf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; write out a stub of the photometric properties table

    texfile = 'photometric_properties_table_stub.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{Photometric Properties \label{table:photometric_properties}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
    
    for i = 0L, nstub-1L do begin
       line = strarr(ncols)
       for j = 0L, ncols-1L do begin
          if (j lt ncols-1L) then suffix = ' & ' else suffix = ' \\ '

          if (j eq 0L) then line[j] = table[i].(j)+suffix else begin ; ID number

             flux = table[i].(2*j-1L)
             ferr = table[i].(2*j)

; well-measured flux and error             
             
             if (strmatch(flux,'*nodata*') eq 0B) and (strmatch(ferr,'*nodata*') eq 0B) then begin
                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+suffix
             endif 

; well-measured flux and but no error; THIS CASE MUST BE FIXED (NGC5194)
             
             if (strmatch(flux,'*nodata*') eq 0B) and (strmatch(ferr,'*nodata*') eq 1B) then begin
                line[j] = '$'+flux+'$'+suffix
                print, 'Hey! '+atlas[i].galaxy
             endif

; no flux measurement (assume no error)
             
             if (strmatch(flux,'*nodata*') eq 1B) then line[j] = '\nodata'+suffix

; flux upper limit (no error)
             
             if (strmatch(flux,'*-*') eq 1B) then line[j] = '$<'+string(-float(flux),format='(G0.0)')+'$'+suffix

          endelse

       endfor
          
       printf, lun, line

    endfor
    
    stubtext = 'Table \ref{table:photometric_properties} is published in its entirety in the electronic '+$
      'edition of the {\it Astrophysical Journal}.  A portion is shown '+$
      'here for guidance regarding its form and content.'
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+strjoin([tablecomments,stubtext],' ')+'}'
    printf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; write out the MRT data, labels, units, and explanation

    labels = ['ID','n_ID','Umag','e_Umag','Bmag','e_Bmag','Vmag','e_Vmag',$
      'Jmag','e_Jmag','Hmag','e_Hmag','Ksmag','e_Ksmag','iras12','e_iras12',$
      'iras25','e_iras25','iras60','e_iras60','iras100','e_iras100']
    units = ['---','---','mag','mag','mag','mag','mag','mag','mag','mag',$
      'mag','mag','mag','mag','Jy','Jy','Jy','Jy','Jy','Jy','Jy','Jy']
    explanations = ['Unique identification number','? [a] Note on ID (1)',$
      '? U magnitude (2)','? U magnitude error (2)','? B magnitude (2)','? B magnitude error (2)',$
      '? V magnitude (2)','? V magnitude error (2)','? J magnitude (3)','? J magnitude error (3)',$
      '? H magnitude (3)','? H magnitude error (3)','? Ks magnitude (3)','? Ks magnitude error (3)',$
      '? 12-micron flux (4)','? Error in the 12-micron flux (4)','? 25-micron flux (4)','? Error in the 25-micron flux (4)',$
      '? 60-micron flux (4)','? Error in the 60-micron flux (4)','? 100-micron flux (4)','? Error in the 100-micron flux (4)']
    niceprint, labels, units, explanations

    openw, lun, mrtpath+'mrt_photometric_properties_labels.dat', /get_lun
    for i = 0L, n_elements(labels)-1L do printf, lun, labels[i]
    free_lun, lun

    openw, lun, mrtpath+'mrt_photometric_properties_units.dat', /get_lun
    for i = 0L, n_elements(units)-1L do printf, lun, units[i]
    free_lun, lun

    openw, lun, mrtpath+'mrt_photometric_properties_explanations.dat', /get_lun
    for i = 0L, n_elements(explanations)-1L do printf, lun, explanations[i]
    free_lun, lun

    ntags = n_tags(table_mrt)
    openw, lun, mrtpath+'mrt_photometric_properties_data.dat', /get_lun
    for i = 0L, ngalaxy-1L do begin
       for j = 0L, ntags-1L do begin
          if strmatch(table_mrt[i].(j),'*nodata*') then table_mrt[i].(j) = ''
          if (j eq 0L) then line = strtrim(table_mrt[i].(j),2) else $
            line = [line,strtrim(table_mrt[i].(j),2)]
       endfor
       printf, lun, strjoin(line,' & ')
    endfor
    free_lun, lun

; ---------------------------------------------------------------------------    
; Table 3: journal of observations
; ---------------------------------------------------------------------------    

    select = [$
      'atlas_id','nice_galaxy','drift_scan','drift_ap','drift_posangle',$
      'drift_oexptime','drift_photflag','drift_absphoterr','drift_comments']
    format = [$
      'I3.3','A0','I0','I0','I0',$
      'I0','A0','F15.3','A0']
    texcenter = [$
      'l','l','c','c','c',$
      'c','c','c','l']
    colhead = '\colhead{'+[ [$
      'ID','Galaxy name','$\Delta_{\rm scan}$','Ap.','$\theta_{\rm slit}$',$
      '$t$','Clear','$\Delta$m']+'} & ',['Remarks\tablenotemark{a}']+'} \\ ']
    colunits = '\colhead{'+[ [$
      '(1)','(2)','(3)','(4)',$
      '(5)','(6)','(7)','(8)']+'} & ',['(9)']+'} ']
    tablecomments = [$
      'Col. (1) ID number;',$
      'Col. (2) Galaxy name;',$
      'Col. (3) Drift scan length perpendicular to the slit (arcsec);',$
      'Col. (4) Extraction aperture along the slit (arcsec);',$
      'Col. (5) Slit position angle measured positive from North to East (deg);',$
      'Col. (6) Total exposure time (seconds);',$
      'Col. (7) Flag indicating clear (Y) or non-photometric (N) observing conditions;',$
      'Col. (8) Absolute spectrophotometric uncertainty ($1\sigma$) '+$
      'based only on the scatter in the observed sensitivity function (mag);',$
      'Col. (9) Remarks regarding the object or the spectral extraction.']
    tablecomments = strjoin(tablecomments,' ')
    tablenotetext = '{a}{In this footnote we clarify some of the remarks that appear in column (9): '+$
      '{\em stellar contamination}: indicates foreground stellar contamination '+$
      'that could not be subtracted; {\em foreground star(s) excluded}: a smaller extraction aperture '+$
      'was adopted to avoid one or more foreground stars; {\em scan avoids star(s)}: '+$
      'the drift scan center and length were selected to '+$
      'avoid nearby bright stars.}'

    table = html_structure_parse(atlas,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata')
    ntags = n_tags(table)

    moretags_mrt = replicate({atlas_id: strarr(1), atlas_id_note: strarr(1)},ngalaxy)
    table_mrt = struct_addtags(moretags_mrt,struct_trimtags(table,except=['atlas_id']))
    table_mrt.atlas_id = table.atlas_id

; annotate particular galaxies    

    annotategals = ['ngc1614','ngc2782','ngc3079','ngc3448','ngc3985','ngc4096',$
      'ngc4194','ngc4389','ngc4676','ngc4676a','ngc4676b','ic0883','ngc6090',$
      'ugc12490','mrk0331']
    nannotate = n_elements(annotategals)
    notesymbols = ['b','c','d','e','f','g','h','i','j',$
      'k','l','m','n','o','p']
    tablenotes = [ $
      ['{Variance-weighted average of three '+$
      'indistinguishable spectra extracted using $80\arcsec\times25\arcsec$ (PA=$180^{\circ}$), '+$
      '$80\arcsec\times60\arcsec$ (PA=$180^{\circ}$), and '+$
      '$80\arcsec\times60\arcsec$ (PA=$90^{\circ}$) apertures.}'], $
      ['{Variance-weighted average of two '+$
      'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
      'and extracted using $110\arcsec\times60\arcsec$ and $110\arcsec\times100\arcsec$ apertures.}'], $
      ['{Variance-weighted average of two '+$
      'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
      'and extracted using $90\arcsec\times360\arcsec$ and $100\arcsec\times360\arcsec$ apertures.}'], $
      ['{Variance-weighted average of two '+$
      'indistinguishable spectra extracted using $160\arcsec\times60\arcsec$ (PA=$70^{\circ}$) and '+$
      '$140\arcsec\times90\arcsec$ (PA=$90^{\circ}$) apertures.}'], $
      ['{Variance-weighted average of two '+$
      'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
      'and extracted using $80\arcsec\times40\arcsec$ and $80\arcsec\times50\arcsec$ apertures.}'], $
      ['{Variance-weighted average of two '+$
      'indistinguishable spectra taken at $114^{\circ}$ and $108^{\circ}$.}'], $
      ['{Variance-weighted average of two '+$
      'indistinguishable spectra extracted using $125\arcsec\times30\arcsec$ (PA=$165^{\circ}$) and '+$
      '$60\arcsec\times60\arcsec$ (PA=$90^{\circ}$) apertures.}'], $
      ['{Variance-weighted average of three '+$
      'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
      'and extracted using $140\arcsec\times60\arcsec$ and two $140\arcsec\times75\arcsec$ apertures.}'], $
      ['{Variance-weighted average of two '+$
      'indistinguishable spectra extracted using $75\arcsec\times150\arcsec$ (PA=$90^{\circ}$) and '+$
      '$95\arcsec\times155\arcsec$ (PA=$180^{\circ}$) apertures.}'], $
      ['{Variance-weighted average of two '+$
      'indistinguishable spectra extracted using $40\arcsec\times150\arcsec$ (PA=$90^{\circ}$) and '+$
      '$45\arcsec\times155\arcsec$ (PA=$180^{\circ}$) apertures.}'], $
      ['{Variance-weighted average of two '+$
      'indistinguishable spectra extracted using $35\arcsec\times150\arcsec$ (PA=$90^{\circ}$) and '+$
      '$50\arcsec\times155\arcsec$ (PA=$180^{\circ}$) apertures.}'], $
      ['{Variance-weighted average of three '+$
      'indistinguishable spectra extracted using $90\arcsec\times20\arcsec$ (PA=$130^{\circ}$), '+$
      '$90\arcsec\times20\arcsec$ (PA=$130^{\circ}$), and '+$
      '$50\arcsec\times30\arcsec$ (PA=$90^{\circ}$) apertures.}'], $
      ['{Variance-weighted average of two '+$
      'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
      'and extracted using $45\arcsec\times10\arcsec$ and $45\arcsec\times20\arcsec$ apertures.}'], $
      ['{Variance-weighted average of two '+$
      'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
      'and extracted using $50\arcsec\times30\arcsec$ and $50\arcsec\times40\arcsec$ apertures.}'], $
      ['{Variance-weighted average of two '+$
      'indistinguishable spectra taken at $90^{\circ}$ and $153^{\circ}$.}'] ]

    for ii = 0L, nannotate-1L do begin
       galmatch = where(strmatch(galaxy,annotategals[ii],/fold),ngalmatch)
       if (ngalmatch eq 0L) then message, 'Problem!'
       tablenotetext = [tablenotetext,'{'+notesymbols[ii]+'}'+tablenotes[*,ii]]
       table[galmatch].atlas_id = strtrim(table[galmatch].atlas_id,2)+'\tablenotemark{'+notesymbols[ii]+'}'
       table_mrt[galmatch].atlas_id_note = notesymbols[ii]
    endfor
    
; write out the journal of observations table

    texfile = 'observing_journal_table.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
;   printf, lun, '\rotate'
    printf, lun, '\tablecaption{Summary of Integrated Spectrophotometric Observations\label{table:journal}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
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
          line[j] = table[i].(j)+suffix
       endfor 
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
    for inote = 0L, n_elements(tablenotetext)-1L do printf, lun, '\tablenotetext'+tablenotetext[inote] ; offset from the blank string
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; write out the stub of the journal of observations table

    texfile = 'observing_journal_table_stub.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
;   printf, lun, '\rotate'
    printf, lun, '\tablecaption{Summary of Integrated Spectrophotometric Observations\label{table:journal}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0L, nstub-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else begin
             if (i lt ngalaxy-1L) then suffix = ' \\ ' else suffix = ''
          endelse
          line[j] = table[i].(j)+suffix
       endfor 
       printf, lun, line
    endfor
    
    stubtext = 'Table \ref{table:journal} is published in its entirety in the electronic '+$
      'edition of the {\it Astrophysical Journal}.  A portion is shown '+$
      'here for guidance regarding its form and content.'
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+strjoin([tablecomments,stubtext],' ')+'}'
    for inote = 0L, n_elements(tablenotetext)-1L do printf, lun, '\tablenotetext'+tablenotetext[inote] ; offset from the blank string
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; write out the MRT data, labels, units, and explanation

    labels = ['ID','n_ID','Galaxy','Scan','Aperture','PA','Exptime','Photflag','Absphoterr','Remarks']
    units = ['---','---','---','arcsec','arcsec','deg','s','---','mag','---']
    explanations = ['Unique identification number','? [bp] Note on ID (1)','Galaxy name','Drift scan length',$
      'Extraction aperture','Position angle','Exposure time','Flag indicating clear or non-photometric conditions',$
      'Absolute spectrophotometric uncertainty','? Remarks (2)']
;   niceprint, labels, units, explanations

    openw, lun, mrtpath+'mrt_observing_journal_labels.dat', /get_lun
    for i = 0L, n_elements(labels)-1L do printf, lun, labels[i]
    free_lun, lun

    openw, lun, mrtpath+'mrt_observing_journal_units.dat', /get_lun
    for i = 0L, n_elements(units)-1L do printf, lun, units[i]
    free_lun, lun

    openw, lun, mrtpath+'mrt_observing_journal_explanations.dat', /get_lun
    for i = 0L, n_elements(explanations)-1L do printf, lun, explanations[i]
    free_lun, lun

    ntags = n_tags(table_mrt)
    openw, lun, mrtpath+'mrt_observing_journal_data.dat', /get_lun
    for i = 0L, ngalaxy-1L do begin
       for j = 0L, ntags-1L do begin
          if strmatch(table_mrt[i].(j),'*nodata*') then table_mrt[i].(j) = ''
          if (j eq 0L) then line = strtrim(table_mrt[i].(j),2) else $
            line = [line,strtrim(table_mrt[i].(j),2)]
       endfor
       printf, lun, strjoin(line,' & ')
    endfor
    free_lun, lun

; ###########################################################################    
; integrated & nuclear spectra    
; ###########################################################################    
    
; ---------------------------------------------------------------------------    
; Table: integrated continuum and absorption-line indices
; ---------------------------------------------------------------------------    

    select = [$
      'atlas_id','D4000_narrow','C41_50','lick_hb',$
      'lick_mgb','lick_fe','lick_hd_a','babs_h_alpha_ew',$
      'babs_h_beta_ew','babs_h_gamma_ew']
    format = [$
      'I3.3','F15.7','F15.7','F15.7',$
      'F15.7','F15.7','F15.7','F15.7',$
      'F15.7','F15.7']
    texcenter = [$
      'l','c','c','c',$
      'c','c','c','c',$
      'c','c']
    colhead = '\colhead{'+[ [$
      'ID','D$_{\rm n}(4000)$','41-50','H$\beta$',$
      'Mg~$b$','$\langle {\rm Fe}\rangle$','H$\delta_{\rm A}$','EW(H$\alpha$)',$
      'EW(H$\beta$)']+'} & ',['EW(H$\gamma$)']+'} \\']
    colunits = '\colhead{'+[ [$
      '(1)','(2)','(3)','(4)',$
      '(5)','(6)','(7)','(8)',$
      '(9)']+'} & ',['(10)']+'} ']
    tablecomments = [$
      'Col. (1) ID number;',$
      'Col. (2) $4000$-\AA{} break \citep{balogh99};',$
      'Col. (3) Continuum color index \citep{kenn92a};',$
      'Col. (4) H$\beta$ Lick index in Angstroms \citep{trager98};',$
      'Col. (5) Mg~$b$ Lick index in Angstroms \citep{trager98};',$
      'Col. (6) Average Fe absorption-line strength in Angstroms \citep{thomas03};',$
      'Col. (7) H$\delta_{\rm A}$ absorption-line index \citep{worthey97};',$
      'Col. (8) Equivalent width of the H$\alpha$ Balmer absorption-line correction in Angstroms;',$
      'Col. (9) Equivalent width of the H$\beta$ Balmer absorption-line correction in Angstroms;',$
      'Col. (10) Equivalent width of the H$\gamma$ Balmer absorption-line correction in Angstroms.']
    tablecomments = strjoin(tablecomments,' ')
    tablenotetext = '{a}{Integrated continuum and absorption-line indices.}'

    table = html_structure_parse(int,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error,/keep_upperlimits)
    ntags = n_tags(table)

; write out the EW table

    texfile = 'indices_int_table.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{Integrated Continuum and Absorption-Line Indices\tablenotemark{a} \label{table:int_indices}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
    
    for i = 0L, nint-1L do begin

       line = strarr(ntags)
       for j = 0L, ntags-1L do begin

          if (j eq 0L) then line[j] = table[i].(j)+' & ' else begin ; GALAXY tag
             
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '

             if (strtrim(table[i].(j)[1],2) eq -3.0) then $
               line[j] = '$<'+string(fix_digits(float(table[i].(j)[0]),4),format='(G0.0)')+'$'+suffix

             if (strtrim(table[i].(j)[1],2) eq -2.0) then line[j] = '\nodata'+suffix
             if (strtrim(table[i].(j)[1],2) gt 0.0) then begin

                flux = table[i].(j)[0]
                ferr = table[i].(j)[1]

                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+suffix

             endif
          endelse
       endfor
       printf, lun, line

    endfor
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
    printf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; ---------------------------------------------------------------------------    
; Table: nuclear continuum and absorption-line indices
; ---------------------------------------------------------------------------    

    select = [$
      'atlas_id','D4000_narrow','C41_50','lick_hb',$
      'lick_mgb','lick_fe','lick_hd_a','babs_h_alpha_ew',$
      'babs_h_beta_ew','babs_h_gamma_ew']
    format = [$
      'I3.3','F15.7','F15.7','F15.7',$
      'F15.7','F15.7','F15.7','F15.7',$
      'F15.7','F15.7']
    texcenter = [$
      'l','c','c','c',$
      'c','c','c','c',$
      'c','c']
    colhead = '\colhead{'+[ [$
      'ID','D$_{\rm n}(4000)$','41-50','H$\beta$',$
      'Mg~$b$','$\langle {\rm Fe}\rangle$','H$\delta_{\rm A}$','EW(H$\alpha$)',$
      'EW(H$\beta$)']+'} & ',['EW(H$\gamma$)']+'} \\']
    colunits = '\colhead{'+[ [$
      '(1)','(2)','(3)','(4)',$
      '(5)','(6)','(7)','(8)',$
      '(9)']+'} & ',['(10)']+'} ']
    tablecomments = [$
      'Col. (1) ID number;',$
      'Col. (2) $4000$-\AA{} break \citep{balogh99};',$
      'Col. (3) Continuum color index \citep{kenn92a};',$
      'Col. (4) H$\beta$ Lick index in Angstroms \citep{trager98};',$
      'Col. (5) Mg~$b$ Lick index in Angstroms \citep{trager98};',$
      'Col. (6) Average Fe absorption-line strength in Angstroms \citep{thomas03};',$
      'Col. (7) H$\delta_{\rm A}$ absorption-line index \citep{worthey97};',$
      'Col. (8) Equivalent width of the H$\alpha$ Balmer absorption-line correction in Angstroms;',$
      'Col. (9) Equivalent width of the H$\beta$ Balmer absorption-line correction in Angstroms;',$
      'Col. (10) Equivalent width of the H$\gamma$ Balmer absorption-line correction in Angstroms.']
    tablecomments = strjoin(tablecomments,' ')
    tablenotetext = '{a}{Nuclear continuum and absorption-line indices.}'

    table = html_structure_parse(nuc,tags=select,tagformats=format,$
      tagindices=match,blank='\nodata',/keep_error,/keep_upperlimits)
    ntags = n_tags(table)

; write out the EW table

    texfile = 'indices_nuc_table.tex'
    splog, 'Writing '+texpath+texfile+'.'
    openw, lun, texpath+texfile, /get_lun
    if keyword_set(textest) then begin
       printf, lun, '\documentclass[10pt,preprint]{aastex}'
       printf, lun, '\begin{document}'
    endif
    printf, lun, '\voffset'+voffset+'in'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\rotate'
    printf, lun, '\tablecaption{Nuclear Continuum and Absorption-Line Indices\tablenotemark{a} \label{table:nuc_indices}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'
    
    for i = 0L, nnuc-1L do begin

       line = strarr(ntags)
       for j = 0L, ntags-1L do begin

          if (j eq 0L) then line[j] = table[i].(j)+' & ' else begin ; GALAXY tag
             
             if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '

             if (strtrim(table[i].(j)[1],2) eq -3.0) then $
               line[j] = '$<'+string(fix_digits(float(table[i].(j)[0]),4),format='(G0.0)')+'$'+suffix

             if (strtrim(table[i].(j)[1],2) eq -2.0) then line[j] = '\nodata'+suffix
             if (strtrim(table[i].(j)[1],2) gt 0.0) then begin

                flux = table[i].(j)[0]
                ferr = table[i].(j)[1]

                format_flux_error, flux, ferr, newflux, newferr
                line[j] = '$'+newflux+'\pm'+newferr+'$'+suffix

             endif
          endelse
       endfor
       printf, lun, line

    endfor
    
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
    printf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    if keyword_set(textest) then printf, lun, '\end{document}'
    free_lun, lun
    if keyword_set(textest) then begin
       pushd, texpath
       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
       popd
    endif

; OLD

;;;; write out the journal of observations table
;;;
;;;    texfile = 'observing_journal_table.tex'
;;;    splog, 'Writing '+texpath+texfile+'.'
;;;    openw, lun, texpath+texfile, /get_lun
;;;    if keyword_set(textest) then begin
;;;       printf, lun, '\documentclass[10pt,preprint]{aastex}'
;;;       printf, lun, '\begin{document}'
;;;    endif
;;;    printf, lun, '\voffset'+voffset+'in'
;;;    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;;;    printf, lun, '\tabletypesize{\tiny}'
;;;;   printf, lun, '\rotate'
;;;    printf, lun, '\tablecaption{Summary of Integrated Spectrophotometric Observations\label{table:journal}}'
;;;    printf, lun, '\tablewidth{0pt}'
;;;    printf, lun, '\tablehead{'
;;;    niceprintf, lun, colhead
;;;    niceprintf, lun, colunits
;;;    printf, lun, '}'
;;;    printf, lun, '\startdata'
;;;
;;;    notecounter = 1L ; offset from the caption note
;;;    for i = 0L, ngalaxy-1L do begin
;;;       line = strarr(ntags)
;;;       for j = 0L, ntags-1L do begin
;;;          note = ''
;;;          if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '
;;;          if (j eq 0L) then begin
;;;             if (strmatch(galaxy[i],'ic0883',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of three '+$
;;;                  'indistinguishable spectra extracted using $90\arcsec\times20\arcsec$ (PA=$130^{\circ}$), '+$
;;;                  '$90\arcsec\times20\arcsec$ (PA=$130^{\circ}$), and '+$
;;;                  '$50\arcsec\times30\arcsec$ (PA=$90^{\circ}$) apertures.}']
;;;                note = '\tablenotemark{'+note+'}' ; this has to happen after TABLENOTETEXT
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'mrk0331',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra taken at $90^{\circ}$ and $153^{\circ}$.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc1614',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of three '+$
;;;                  'indistinguishable spectra extracted using $80\arcsec\times25\arcsec$ (PA=$180^{\circ}$), '+$
;;;                  '$80\arcsec\times60\arcsec$ (PA=$180^{\circ}$), and '+$
;;;                  '$80\arcsec\times60\arcsec$ (PA=$90^{\circ}$) apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc2782',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
;;;                  'and extracted using $110\arcsec\times60\arcsec$ and $110\arcsec\times100\arcsec$ apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc3079',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
;;;                  'and extracted using $90\arcsec\times360\arcsec$ and $100\arcsec\times360\arcsec$ apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc3448',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra extracted using $160\arcsec\times60\arcsec$ (PA=$70^{\circ}$) and '+$
;;;                  '$140\arcsec\times90\arcsec$ (PA=$90^{\circ}$) apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc3985',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
;;;                  'and extracted using $80\arcsec\times40\arcsec$ and $80\arcsec\times50\arcsec$ apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc4096',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra taken at $114^{\circ}$ and $108^{\circ}$.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc4194',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra extracted using $125\arcsec\times30\arcsec$ (PA=$165^{\circ}$) and '+$
;;;                  '$60\arcsec\times60\arcsec$ (PA=$90^{\circ}$) apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc4389',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of three '+$
;;;                  'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
;;;                  'and extracted using $140\arcsec\times60\arcsec$ and two $140\arcsec\times75\arcsec$ apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc4676',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra extracted using $75\arcsec\times150\arcsec$ (PA=$90^{\circ}$) and '+$
;;;                  '$95\arcsec\times155\arcsec$ (PA=$180^{\circ}$) apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc4676a',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra extracted using $40\arcsec\times150\arcsec$ (PA=$90^{\circ}$) and '+$
;;;                  '$45\arcsec\times155\arcsec$ (PA=$180^{\circ}$) apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc4676b',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra extracted using $35\arcsec\times150\arcsec$ (PA=$90^{\circ}$) and '+$
;;;                  '$50\arcsec\times155\arcsec$ (PA=$180^{\circ}$) apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc6090',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
;;;                  'and extracted using $45\arcsec\times10\arcsec$ and $45\arcsec\times20\arcsec$ apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ugc12490',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
;;;                  'and extracted using $50\arcsec\times30\arcsec$ and $50\arcsec\times40\arcsec$ apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;          endif
;;;          line[j] = table[i].(j)+note+suffix
;;;       endfor
;;;       printf, lun, line
;;;    endfor
;;;    
;;;    printf, lun, '\enddata'
;;;    printf, lun, '\tablecomments{'+tablecomments+'}'
;;;    for inote = 0L, n_elements(tablenotetext)-1L do printf, lun, '\tablenotetext'+tablenotetext[inote] ; offset from the blank string
;;;    printf, lun, '\end{deluxetable}'
;;;    if keyword_set(textest) then printf, lun, '\end{document}'
;;;    free_lun, lun
;;;    if keyword_set(textest) then begin
;;;       pushd, texpath
;;;       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
;;;       popd
;;;    endif
;;;
;;;; write out the stub of the journal of observations table
;;;
;;;    texfile = 'observing_journal_table_stub.tex'
;;;    splog, 'Writing '+texpath+texfile+'.'
;;;    openw, lun, texpath+texfile, /get_lun
;;;    if keyword_set(textest) then begin
;;;       printf, lun, '\documentclass[10pt,preprint]{aastex}'
;;;       printf, lun, '\begin{document}'
;;;    endif
;;;    printf, lun, '\voffset'+voffset+'in'
;;;    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;;;    printf, lun, '\tabletypesize{\tiny}'
;;;;   printf, lun, '\rotate'
;;;    printf, lun, '\tablecaption{Summary of Integrated Spectrophotometric Observations\label{table:journal}}'
;;;    printf, lun, '\tablewidth{0pt}'
;;;    printf, lun, '\tablehead{'
;;;    niceprintf, lun, colhead
;;;    niceprintf, lun, colunits
;;;    printf, lun, '}'
;;;    printf, lun, '\startdata'
;;;
;;;    notecounter = 1L ; offset from the caption note
;;;    for i = 0L, nstub-1L do begin
;;;       line = strarr(ntags)
;;;       for j = 0L, ntags-1L do begin
;;;          note = ''
;;;          if (j lt ntags-1L) then suffix = ' & ' else suffix = ' \\ '
;;;          if (j eq 0L) then begin
;;;             if (strmatch(galaxy[i],'ic0883',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of three '+$
;;;                  'indistinguishable spectra extracted using $90\arcsec\times20\arcsec$ (PA=$130^{\circ}$), '+$
;;;                  '$90\arcsec\times20\arcsec$ (PA=$130^{\circ}$), and '+$
;;;                  '$50\arcsec\times30\arcsec$ (PA=$90^{\circ}$) apertures.}']
;;;                note = '\tablenotemark{'+note+'}' ; this has to happen after TABLENOTETEXT
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'mrk0331',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra taken at $90^{\circ}$ and $153^{\circ}$.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc1614',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of three '+$
;;;                  'indistinguishable spectra extracted using $80\arcsec\times25\arcsec$ (PA=$180^{\circ}$), '+$
;;;                  '$80\arcsec\times60\arcsec$ (PA=$180^{\circ}$), and '+$
;;;                  '$80\arcsec\times60\arcsec$ (PA=$90^{\circ}$) apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc2782',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
;;;                  'and extracted using $110\arcsec\times60\arcsec$ and $110\arcsec\times100\arcsec$ apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc3079',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
;;;                  'and extracted using $90\arcsec\times360\arcsec$ and $100\arcsec\times360\arcsec$ apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc3448',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra extracted using $160\arcsec\times60\arcsec$ (PA=$70^{\circ}$) and '+$
;;;                  '$140\arcsec\times90\arcsec$ (PA=$90^{\circ}$) apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc3985',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
;;;                  'and extracted using $80\arcsec\times40\arcsec$ and $80\arcsec\times50\arcsec$ apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc4096',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra taken at $114^{\circ}$ and $108^{\circ}$.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc4194',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra extracted using $125\arcsec\times30\arcsec$ (PA=$165^{\circ}$) and '+$
;;;                  '$60\arcsec\times60\arcsec$ (PA=$90^{\circ}$) apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc4389',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of three '+$
;;;                  'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
;;;                  'and extracted using $140\arcsec\times60\arcsec$ and two $140\arcsec\times75\arcsec$ apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc4676',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra extracted using $75\arcsec\times150\arcsec$ (PA=$90^{\circ}$) and '+$
;;;                  '$95\arcsec\times155\arcsec$ (PA=$180^{\circ}$) apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc4676a',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra extracted using $40\arcsec\times150\arcsec$ (PA=$90^{\circ}$) and '+$
;;;                  '$45\arcsec\times155\arcsec$ (PA=$180^{\circ}$) apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc4676b',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra extracted using $35\arcsec\times150\arcsec$ (PA=$90^{\circ}$) and '+$
;;;                  '$50\arcsec\times155\arcsec$ (PA=$180^{\circ}$) apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ngc6090',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
;;;                  'and extracted using $45\arcsec\times10\arcsec$ and $45\arcsec\times20\arcsec$ apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;             if (strmatch(galaxy[i],'ugc12490',/fold)) then begin
;;;                note = '{'+notesymbol[notecounter]+'}'
;;;                tablenotetext = [tablenotetext,note+'{Variance-weighted average of two '+$
;;;                  'indistinguishable spectra obtained at a $90^{\circ}$ slit position angle '+$
;;;                  'and extracted using $50\arcsec\times30\arcsec$ and $50\arcsec\times40\arcsec$ apertures.}']
;;;                note = '\tablenotemark{'+note+'}'
;;;                notecounter = notecounter + 1L
;;;             endif
;;;          endif
;;;          line[j] = table[i].(j)+note+suffix
;;;       endfor
;;;       printf, lun, line
;;;    endfor
;;;    
;;;    stubtext = 'Table \ref{table:journal} is published in its entirety in the electronic '+$
;;;      'edition of the {\it Astrophysical Journal}.  A portion is shown '+$
;;;      'here for guidance regarding its form and content.'
;;;    
;;;    printf, lun, '\enddata'
;;;    printf, lun, '\tablecomments{'+strjoin([tablecomments,stubtext],' ')+'}'
;;;    for inote = 0L, n_elements(tablenotetext)-1L do printf, lun, '\tablenotetext'+tablenotetext[inote] ; offset from the blank string
;;;    printf, lun, '\end{deluxetable}'
;;;    if keyword_set(textest) then printf, lun, '\end{document}'
;;;    free_lun, lun
;;;    if keyword_set(textest) then begin
;;;       pushd, texpath
;;;       spawn, ['latex '+texfile+' ; dvips -o '+repstr(texfile,'.tex','.ps')+' '+repstr(texfile,'.tex','')]
;;;       popd
;;;    endif
    
return
end
