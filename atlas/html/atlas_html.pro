;+
; NAME:
;       ATLAS_HTML
;
; PURPOSE:
;       Generate web page visualizations for ATLAS.
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       This routine assumes that all the appropriate directories
;       already exist.
;
; PROCEDURES USED:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 April 05
;-

pro atlas_populate_frames, atlas, suffix=suffix, http=http

    webpath = atlas_path(/web)
    thumbpath = webpath+'thumb/'
    masspath = webpath+'mass/'

    httphome = '"'+http+'"'
    baseref = http+'research/spectral_atlas/'
    httpstyle = baseref+'atlas.css'
    
    if (n_elements(suffix) eq 0L) then suffix = '' else $
      if (strcompress(suffix,/remove) ne '') then suffix = '_'+suffix

    galaxy = strtrim(strlowcase(atlas.galaxy),2)
    nedgalaxy = repstr(strtrim(atlas.ned_galaxy,2),'+','%2B')
    altgalaxy = strtrim(atlas.alt_galaxy,2)
;   nicegalaxy = strtrim(atlas.nice_galaxy,2)
    ngalaxy = n_elements(galaxy)

    nicegalaxy = strtrim(atlas.nice_galaxy,2)
    altgalaxy = strtrim(atlas.alt_galaxy,2)
;   nicegalaxy = strarr(ngalaxy)
;   for i = 0L, ngalaxy-1L do if (altgalaxy[i] eq '') then $
;     nicegalaxy[i] = strtrim(atlas[i].nice_galaxy,2) else $
;     nicegalaxy[i] = strtrim(atlas[i].nice_galaxy,2)+' / '+altgalaxy[i]

    htmlfile = 'html/'+galaxy+suffix+'.html'

    specpng = '../png/'+galaxy+'_spec.png'
    specnucpng = '../png/'+galaxy+'_nuclear.png'
    thumbpng = '../thumb/'+galaxy+'_thumb.png'

; -------------------------
; WEBATLAS_SUFFIX.HTML
; -------------------------

    mainhtml = 'webatlas'+suffix+'.html'
    tablehtml = 'webatlas_table'+suffix+'.html'
    datahtml = 'webatlas_data'+suffix+'.html'

    pushd, webpath

    splog, 'Writing '+mainhtml+'.'
    openw, lun, mainhtml, /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
    printf, lun, '<html>'
    printf, lun, '<head>'
;   printf, lun, '<base href="'+webpath+'" />'
    printf, lun, '<link rel="stylesheet" type="text/css" href="'+httpstyle+'" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '<title>An Integrated Spectrophotometric Survey of Nearby Star-Forming Galaxies</title>'
    printf, lun, '</head>'
    printf, lun, '<frameset cols="25%,*">'
    printf, lun, '<frame src="'+tablehtml+'"></frame>'
    printf, lun, '<frame src="'+datahtml+'" name="data"></frame>'
    printf, lun, '</frameset>'
    printf, lun, '</html>'
    free_lun, lun

; -------------------------
; WEBATLAS_TABLE_SUFFIX.HTML
; -------------------------

    digits = string(strlen(string(ngalaxy,format='(I0)')),format='(I0)')
    index = string(lindgen(ngalaxy)+1,format='(I'+digits+'.'+digits+')')
    
    splog, 'Writing '+tablehtml+'.'
    openw, lun, tablehtml, /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
    printf, lun, '<html>'
    printf, lun, '<head>'
;   printf, lun, '<base href="'+webpath+'" />'
    printf, lun, '<link rel="stylesheet" type="text/css" href="'+httpstyle+'" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '</head>'
    printf, lun, '<body>'
    printf, lun, '<h1>Sample</h1>'
;   printf, lun, '<p class="menu">'
;   printf, lun, '<span class="left"><a href="'+baseref+'index.html" target="_top">Home</a></span>'
;   printf, lun, '<br /></p>'
     printf, lun, '<table border="2" width="100%">'
    printf, lun, '<tbody>'
;   for i = 0, ngalaxy-1L do $
;     printf, lun, '<tr><td class="atlastable"><a href="'+htmlfile[i]+$
;       '" target="data">'+nicegalaxy[i]+'</a></td></tr>'
    for i = 0, ngalaxy-1L do $
      printf, lun, '<tr><td class="atlastable"><a href="'+htmlfile[i]+$
        '" target="data">'+index[i]+' '+nicegalaxy[i]+'</a></td></tr>'
    printf, lun, '</tbody>'
    printf, lun, '</table>'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

; -------------------------
; WEBATLAS_DATA_SUFFIX.HTML
; -------------------------

    splog, 'Writing '+datahtml+'.'
    openw, lun, datahtml, /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
    printf, lun, '<html>'
    printf, lun, '<head>'
;   printf, lun, '<base href="'+webpath+'" />'
    printf, lun, '<link rel="stylesheet" type="text/css" href="'+httpstyle+'" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '</head>'
    printf, lun, '<body>'
    printf, lun, '<h1>Visualizations</h1>'
    printf, lun, '<br />'
    printf, lun, '<p>Please select a galaxy.</p>'
    printf, lun, '<br />'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

; ---------------------------------------------------------------------------
; Generate individual HTML files for each object
; ---------------------------------------------------------------------------

; galaxy properties

    galtags = [$
      'ra',         $
      'dec',        $
      'cz',         $
      'distance',   $
      'rc3_type',   $
      'd25_maj',    $
      'd25_min',    $
      'posangle',   $
      'inclination' $
      ]
    galformats = [$
      'A0',   $
      'A0',   $
      'I0',   $
      'F12.1',$
      'A0',   $
      'F12.1',$
      'F12.1',$
      'I0',   $
      'I0'    $
      ]
    gallabels = [$
      'RA',  $
      'DEC', $
      'cz',  $
      'D',   $
      'Type',$
      'D25', $
      'd25', $
      'PA',  $
      'incl' $
      ]
    galunits = [$
      '[J2000]', $
      '[J2000]', $
      '[km/s]',  $
      '[Mpc]',   $
      '[RC3]',   $
      '[arcmin]',$
      '[arcmin]',$
      '[degree]',$
      '[degree]' $
      ]

    galdata = html_structure_parse(atlas,tags=galtags,tagformats=galformats,tagindices=galindices)
    gallabels = gallabels[galindices] & galunits = galunits[galindices]
    ngaltags = n_tags(galdata)

; photometric properties

    phototags = [$
      'rc3_m_b',     $
      'twomass_m_ks',$
      'rc3_b',       $
      'rc3_v',       $
      'twomass_j',   $
      'twomass_ks',  $
      'ir_lum',      $
      'L_ir_L_b'     $
      ]
    photoformats = [$
      'F12.1',$
      'F12.1',$
      'F12.1',$
      'F12.1',$
      'F12.1',$
      'F12.1',$
      'F12.2',$
      'F12.1' $
      ]
    photolabels = [$
      'M(B)',      $
      'M(Ks)',     $
      'B',         $
      'V',         $
      'J',         $
      'Ks',        $
      'L(IR)',     $
      'L(IR)/L(B)' $
      ]
    photounits = [$
      '[mag]',   $
      '[mag]',   $
      '[mag]',   $
      '[mag]',   $
      '[mag]',   $
      '[mag]',   $
      '[Lsun]',  $
      '&nbsp'    $
      ]

    photodata = html_structure_parse(atlas,tags=phototags,tagformats=photoformats,tagindices=photoindices)
    photolabels = photolabels[photoindices] & photounits = photounits[photoindices]
    nphototags = n_tags(photodata)

; integrated spectroscopic observing parameters; the last tag has to
; be COMMENTS!

    spectags = [$
      'drift_date',     $
      'drift_strap',    $
      'drift_posangle', $
      'drift_comments'  $
      ]
    specformats = [$
      'A0',$
      'A0',$
      'I0',$
      'A0' $
      ]
    speclabels = [  $
      'Observed',   $
      'Aperture',   $
      'PA',         $
      'Comments'    $
      ]
    specunits = [$
      '&nbsp',            $
      '[arcsec x arcsec]',$
      '[degree]',         $
      '&nbsp'             $
      ]

    specdata = html_structure_parse(atlas,tags=spectags,tagformats=specformats,tagindices=specindices)

    speclabels = speclabels[specindices]
    specunits = specunits[specindices]
    nspectags = n_tags(specdata)

; nuclear spectroscopic observing parameters; the last tag has to
; be COMMENTS!

    nucspectags = [$
      'nuclear_date',     $
      'nuclear_strap',    $
      'nuclear_posangle', $
      'nuclear_comments'  $
      ]

    nucspecdata = html_structure_parse(atlas,tags=nucspectags,$
      tagformats=specformats,tagindices=specindices)

; make the individual pages    
    
    ra = repstr(atlas.ra,':','%3A')
    dec = repstr(repstr(atlas.dec,':','%3A'),'+','%2B') ; negative declinations are not accounted for here!

    splog, 'Generating the individual ATLAS webpages.'
    for i = 0L, ngalaxy-1L do begin

       if (altgalaxy[i] ne '') then alt = ' = '+altgalaxy[i] else alt = ''
       
       openw, lun, htmlfile[i], /get_lun
       printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
       printf, lun, '<html>'
       printf, lun, '<head>'
       printf, lun, '<link rel="stylesheet" type="text/css" href="'+httpstyle+'" />'
       printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
       printf, lun, '</head>'
       printf, lun, '<body>'
;      printf, lun, '<h1>'+nicegalaxy[i]+'</a></h1>'
       printf, lun, '<h2><a href="http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?objname='+nedgalaxy[i]+$
         '&extend=no&out_csys=Equatorial&out_equinox=B2000.0&obj_sort=RA+or+Longitude&zv=z&zv_breaker=10000.0" target="_top">'+$
         nicegalaxy[i]+'</a><span class="gold">'+alt+'</span></h2>'
;      printf, lun, '<h2><a target="_top" href="http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?search_type=Near+'+$
;        'Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon='+ra[i]+'&lat='+dec[i]+'&radius=1&out_csys=Equatorial'+$
;        '&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp='+$
;        'YES&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY">'+$
;        nicegalaxy[i]+'</a></h2>'
       printf, lun, '<div class="menu right">'
       if i eq ngalaxy-1L then nextfile = '../'+htmlfile[0] else nextfile = '../'+htmlfile[i+1L]
       printf, lun, '<a href="'+nextfile+'">Next</a><br />'
       if i eq 0L then backfile = '../'+htmlfile[ngalaxy-1L] else backfile = '../'+htmlfile[i-1L]
       printf, lun, '<a href="'+backfile+'">Back</a>'
       printf, lun, '</div>'

       printf, lun, '<div class="menu left">'
       printf, lun, '<a href="../index.html" target="_top">Home</a>'
       printf, lun, '</div><br clear="bottom"/>';<br />'

;      printf, lun, '<div class="menu right">Page '+string(i+1,format='(I0)')+'</div><br />'

; image + spectral visualizations
       
       printf, lun, '<table width="100%" align="center" border="0" cellspacing="0" cellpadding="0">' ; open table
       printf, lun, '<tbody>'

       printf, lun, '<tr>'
       printf, lun, '<td><a class="nolinkcolor" href="'+thumbpng[i]+'" target="_top">'+$
         '<img width="99%" src="'+thumbpng[i]+'" alt="'+nicegalaxy[i]+' Thumbnail" /></a></td>'
       printf, lun, '<td width="55%"><a class="nolinkcolor" href="'+specpng[i]+'" target="_top">'+$
         '<img width="99%" src="'+specpng[i]+'" alt="'+nicegalaxy[i]+' Thumbnail" /></a></td>'
       printf, lun, '</tr>'

       printf, lun, '</tbody>' 
       printf, lun, '</table>' ; close table

       printf, lun, '<br clear="bottom"/>'
       printf, lun, '<p class="smaller">The dashed purple ellipse shows the 25 mag per sq. arcsec '+$
         'surface brightness isophote.  The red box delineates our integrated spectroscopic aperture.</p>'

; galaxy properties table
       
       printf, lun, '<table width="100%" border="2" cellspacing="0" cellpadding="2" class="props">'
       printf, lun, '<caption>General Properties</caption>'
       printf, lun, '<thead>'
       printf, lun, '<tr>'
       for itag = 0L, ngaltags-1L do printf, lun, '<th>'+gallabels[itag]+'</th>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       for itag = 0L, ngaltags-1L do printf, lun, '<th><span class="smaller">'+galunits[itag]+'</span></th>'
       printf, lun, '</tr>'
       printf, lun, '</thead>'
       printf, lun, '<tbody>'
       printf, lun, '<tr>'
       for itag = 0L, ngaltags-1L do printf, lun, '<td>'+(galdata[i].(itag))[0]+'</td>'
       printf, lun, '</tr>'
       printf, lun, '</tbody>'
       printf, lun, '</table>' 
       printf, lun, '<br clear="bottom"/>' ; close the properties table
       
; photometric properties table
       
       printf, lun, '<table width="100%" border="2" cellspacing="0" cellpadding="2" class="props">'
       printf, lun, '<caption>Photometry</caption>'
       printf, lun, '<thead>'
       printf, lun, '<tr>'
       for itag = 0L, nphototags-1L do printf, lun, '<th>'+photolabels[itag]+'</th>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       for itag = 0L, nphototags-1L do printf, lun, '<th><span class="smaller">'+photounits[itag]+'</span></th>'
       printf, lun, '</tr>'
       printf, lun, '</thead>'
       printf, lun, '<tbody>'
       printf, lun, '<tr>'
       for itag = 0L, nphototags-1L do printf, lun, '<td>'+(photodata[i].(itag))[0]+'</td>'
       printf, lun, '</tr>'
       printf, lun, '</tbody>'
       printf, lun, '</table>' 
       printf, lun, '<br clear="bottom"/>' ; close the properties table
       
; spectroscopic observing parameters
       
       printf, lun, '<table width="100%" border="2" cellspacing="0" cellpadding="2" class="props">'
       printf, lun, '<caption>Optical Spectroscopy</caption>'
       printf, lun, '<thead>'
       printf, lun, '<tr>'
       printf, lun, '<th>Spectrum</th>'
       for itag = 0L, nspectags-1L do printf, lun, '<th>'+speclabels[itag]+'</th>'
       printf, lun, '<th colspan="2">Download</th>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       printf, lun, '<th>&nbsp</th>'
       for itag = 0L, nspectags-1L do printf, lun, '<th><span class="smaller">'+specunits[itag]+'</span></th>'
       printf, lun, '<th colspan="2">1D</th>'
       printf, lun, '</tr>'
       printf, lun, '</thead>'
 
       printf, lun, '<tbody>'
; integrated
       printf, lun, '<tr>'
       printf, lun, '<td>Integrated</td>'       
       if (atlas[i].drift eq 0L) then begin
          for itag = 0L, nspectags-2L do printf, lun, '<td>&nbsp</td>' ; no spectroscopic data at all
          printf, lun, '<td>unavailable</td>'                ; "comments" field
          printf, lun, '<td>&nbsp</td><td>&nbsp</td>'                  ; no ascii/fits/fits2d spectrum
       endif else begin
          for itag = 0L, nspectags-1L do printf, lun, '<td>'+(specdata[i].(itag))[0]+'</td>'
          asciifile1 = strtrim(atlas[i].drift_asciifile,2)
          if file_test(webpath+'ascii/'+asciifile1,/regular) then $
            printf, lun, '<td><a href="../ascii/'+asciifile1+'">ASCII</a></td>' else $
            printf, lun, '<td>&nbsp</td>' ; no ascii spectrum
          fitsfile1 = strtrim(atlas[i].drift_file,2)
          if file_test(webpath+'fits/'+fitsfile1,/regular) then $
            printf, lun, '<td><a href="../fits/'+fitsfile1+'">FITS</a></td>' else $
            printf, lun, '<td>&nbsp</td>' ; no fits spectrum
       endelse
       printf, lun, '</tr>'
; nuclear
       printf, lun, '<tr>'
       printf, lun, '<td>Nuclear</td>'       
       if (atlas[i].nuclear eq 0L) then begin
          for itag = 0L, nspectags-2L do printf, lun, '<td>&nbsp</td>' ; no spectroscopic data at all
          printf, lun, '<td>unavailable</td>'                ; "comments" field
          printf, lun, '<td>&nbsp</td><td>&nbsp</td>'                  ; no ascii/fits/fits2d spectrum
       endif else begin
          for itag = 0L, nspectags-1L do printf, lun, '<td>'+(nucspecdata[i].(itag))[0]+'</td>'
          asciifile1 = strtrim(atlas[i].nuclear_asciifile,2)
          if file_test(webpath+'ascii/'+asciifile1,/regular) then $
            printf, lun, '<td><a href="../ascii/'+asciifile1+'">ASCII</a></td>' else $
            printf, lun, '<td>&nbsp</td>' ; no ascii spectrum
          fitsfile1 = strtrim(atlas[i].nuclear_file,2)
          if file_test(webpath+'fits/'+fitsfile1,/regular) then $
            printf, lun, '<td><a href="../fits/'+fitsfile1+'">FITS</a></td>' else $
            printf, lun, '<td>&nbsp</td>' ; no fits spectrum
       endelse
       printf, lun, '</tr>'

       printf, lun, '</tbody>'
       printf, lun, '</table>' 
       printf, lun, '<br clear="bottom"/>' ; close the spectroscopy table
       
; finish the page       
       
       printf, lun, '<div class="menu right">'
       if i eq ngalaxy-1L then nextfile = '../'+htmlfile[0] else nextfile = '../'+htmlfile[i+1L]
       printf, lun, '<a href="'+nextfile+'">Next</a><br />'
       if i eq 0L then backfile = '../'+htmlfile[ngalaxy-1L] else backfile = '../'+htmlfile[i-1L]
       printf, lun, '<a href="'+backfile+'">Back</a>'
       printf, lun, '</div>'

       printf, lun, '<div class="menu right">Page '+string(i+1,format='(I0)')+'</div><br />'

       printf, lun, '<div class="menu left">'
       printf, lun, '<a href="../index.html" target="_top">Home</a><br />'
       printf, lun, '</div><br /><br clear="bottom"/>'

       printf, lun, '<p class="email">Last update '+im_today()+'.  Email <a href="mailto:john.moustakas@gmail.com">'+$
         'John Moustakas</a> with questions, comments, or to report errors.</p>'

       printf, lun, '</body>'
       printf, lun, '</html>'
       free_lun, lun
       
    endfor
    
    popd

return
end
    
pro atlas_html, atlas, integrated, nuclear, make_png=make_png, make_thumbs=make_thumbs, $
  make_links=make_links, make_tarballs=make_tarballs, debug=debug, _extra=extra

; datapaths    
    
    datapath = atlas_path(/analysis)
    dsspath = atlas_path(/dss)
    webpath = atlas_path(/web)
    specfitpath = atlas_path(/specfit)
    specfitpath_data = atlas_path(/spectral_atlas)
    spec1dpath = atlas_path(/atlas1d)
    ascii1dpath = atlas_path(/ascii)
    htmlpath = webpath+'html/'
    pngpath = webpath+'png/'
    thumbpath = webpath+'thumb/'
    asciipath = webpath+'ascii/'
    fitspath = webpath+'fits/'

    http = 'http://sdss.physics.nyu.edu/ioannis/'
    httphome = '"'+http+'"'
    baseref = http+'research/spectral_atlas/'
    httpstyle = baseref+'atlas.css'
    
    version = atlas_version(/specfit)

    if (n_elements(atlas) eq 0L) then begin
       splog, 'Reading the spectral fitting results.'
       atlas = atlas_read_info()
       atlas = atlas[sort(atlas.galaxy)] ; alphabetize
    endif
    if (n_elements(integrated) eq 0L) then begin
       integrated = read_integrated(silent=0)
    endif
    if (n_elements(nuclear) eq 0L) then begin
       nuclear = read_nuclear(silent=0)
    endif
    ngalaxy = n_elements(atlas)

    galaxy = strtrim(strlowcase(atlas.galaxy),2)

    htmlfile = 'html/'+galaxy+'.html'
    specpng = galaxy+'_spec.png'
    specps = galaxy+'_spec.ps'
    thumbfits = galaxy+'.fits.gz'
    thumbpng = galaxy+'_thumb.png'

; ---------------------------------------------------------------------------
; define subsamples and generate a WEBATLAS.HTML, WEBATLAS_TABLE.HTML
; and WEBATLAS_DATA.HTML file for each subsample   
; ---------------------------------------------------------------------------

; full sample    
    
    atlas_populate_frames, atlas, suffix=suffix, http=http

;; Ursa Major Cluster sub-sample
;
;    umajor = index_umajor(atlas)
;    if (umajor[0] ne -1L) then $
;      atlas_populate_frames, atlas[umajor[sort(atlas[umajor].galaxy)]], suffix='umajor', http=http
;
;; 11 Mpc sub-sample
;
;    mpc11 = index_11mpc(atlas)
;    if (mpc11[0] ne -1L) then $
;      atlas_populate_frames, atlas[mpc11[sort(atlas[mpc11].galaxy)]], suffix='11mpc', http=http
;
;; Markarian sub-sample
;
;    mrk = index_markarian(atlas)
;    if (mrk[0] ne -1L) then $
;      atlas_populate_frames, atlas[mrk[sort(atlas[mrk].galaxy)]], suffix='mrk', http=http
;
;; IR-selected sub-sample
;
;    ir = index_irselected(atlas)
;    if (ir[0] ne -1L) then $
;      atlas_populate_frames, atlas[ir[sort(atlas[ir].galaxy)]], suffix='ir', http=http
;
;; GTO starbursts sub-sample
;
;    gto = index_gto_starbursts(atlas)
;    if (gto[0] ne -1L) then $
;      atlas_populate_frames, atlas[gto[sort(atlas[gto].galaxy)]], suffix='gto', http=http
;
;; HCN survey sub-sample
;
;    hcn = index_hcn_sample(atlas,data=atlas1)
;    if (hcn[0] ne -1L) then $
;      atlas_populate_frames, atlas1[hcn[sort(atlas1[hcn].galaxy)]], suffix='hcn', http=http

; -------------------------
; INDEX.HTML
; -------------------------

    pushd, webpath
    
    splog, 'Generating the ATLAS home page.'
    openw, lun, 'index.html', /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<link rel="stylesheet" type="text/css" href="'+httpstyle+'" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '<title>An Integrated Spectrophotometric Survey of Nearby Star-Forming Galaxies</title>'
    printf, lun, '</head>'
    printf, lun, '<body>'
    printf, lun, '<br />'
    printf, lun, '<h1 class="smaller">Integrated Spectrophotometric Atlas of Nearby Star-Forming Galaxies</h1>'
    printf, lun, '<table width="100%" align="center" border="0" cellspacing="0" cellpadding="0">' ; open table
    printf, lun, '<tbody>'
    printf, lun, '<tr>'
    printf, lun, '<td width="40%" class="normal">'+$
      '<a class="nolinkcolor" href="'+baseref+'atlas_poster.png">'+$
      '<img width="100%" src="'+baseref+'atlas_poster.png" alt="Poster" /></td>'
    printf, lun, '<td class="normal">'
    printf, lun, '<div class="border">'
    printf, lun, '<p>We have obtained integrated and nuclear optical spectrophotometry for a sample of '
    printf, lun, '417 nearby galaxies, targeting a diverse range of galaxy types, including'
    printf, lun, 'starbursts, peculiar galaxies, interacting/merging systems, and dusty, infrared-luminous galaxies.'
    printf, lun, 'The following webpages present one- and two-dimensional visualizations of the data, and various'
    printf, lun, 'ancillary observations.  You can also download all the spectra as a single tarball, or'
    printf, lun, 'the (FITS or ASCII) spectrum of an individual galaxy.</p>'
;   printf, lun, '<p>Welcome to the home page of the integrated optical spectrophotometric '+$
;     'atlas of nearby star-forming galaxies.</p>'
    printf, lun, '<ul>'
    printf, lun, '<li><a href="'+baseref+'webatlas.html">Visualizations</a></li>'
    printf, lun, '<li><a href="'+baseref+'atlasdownload.html">Download</a></li>'
    printf, lun, '<li><a href="'+baseref+'atlaspublications.html">Publications</a></li>'
;   printf, lun, '<li><a href="'+baseref+'webatlas_umajor.html">ursa major cluster sub-sample</a></li>'
;   printf, lun, '<li><a href="'+baseref+'webatlas_11mpc.html">11 Mpc sub-sample</a></li>'
;   printf, lun, '<li><a href="'+baseref+'webatlas_mrk.html">markarian sub-sample</a></li>'
;   printf, lun, '<li><a href="'+baseref+'webatlas_ir.html">infrared sub-sample</a></li>'
;   printf, lun, '<li><a href="'+baseref+'webatlas_gto.html">GTO Starbursts</a></li>'
;   printf, lun, '<li><a href="'+baseref+'webatlas_hcn.html">HCN Survey</a></li>'
    printf, lun, '</ul>'
    printf, lun, '</td>'
    printf, lun, '</tr>'
    printf, lun, '</tbody>'
    printf, lun, '</table>'
    printf, lun, '</div>'
    printf, lun, '<br clear="bottom"/>'

    printf, lun, '<p class="email">Last update '+im_today()+'.  Email <a href="mailto:john.moustakas@gmail.com">'+$
      'John Moustakas</a> with questions, comments, or to report errors.</p>'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

    popd

;   printf, lun, 'Our observations consist of spatially'+$
;     ' integrated, high signal-to-noise spectrophotometry for 351'+$
;     ' galaxies between 3600 and 6900 A at a FWHM resolution of'+$
;     ' 8.0 A.  We also present unpublished integrated spectrophotometry of'+$
;     ' 22 merging galaxies from the Ph.~D. thesis by Anne Turner in addition to'+$
;     ' improved emission-line flux and equivalent width measurements for 32 star-forming'+$
;     ' galaxies from the original Kennicutt spectral atlas.  We develop a population-synthesis'+$
;     ' fitting technique to achieve accurate continuum subtraction and stellar absorption'+$
;     ' corrections to the nebular emission lines.  These web pages'+$
;     ' provide visualizations for individual galaxies, the final, fully-calibrated'+$
;     ' spectra in electronic format, and tabulated emission-line fluxes, equivalent widths,'+$
;     ' and continuum index measurements.  These data are suitable for a'+$
;     ' variety of applications in galaxy formation and evolution, both as a'+$
;     ' reference sample for more distant galaxies and for improving our'+$
;     ' understanding of nearby star-forming galactic systems.</p>'

; -------------------------
; ATLASPUBLICATIONS.HTML
; -------------------------

    pushd, webpath
    
    splog, 'Generating the publications page.'
    openw, lun, 'atlaspublications.html', /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<link rel="stylesheet" type="text/css" href="'+httpstyle+'" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '<title>Spectral Atlas: Publications</title>'
    printf, lun, '</head>'
    printf, lun, '<body>'
    printf, lun, '<h1>Publications</h1>'
    printf, lun, '<br />'
    printf, lun, '<div class="box download">'
    printf, lun, '<ul class="third">'
    printf, lun, '<li><a href="http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2006ApJS..164...81M&amp ;db_key=AST&amp;data_type=HTML&amp;format=&amp;high=452a3b59ab21525">"An Integrated Spectrophotometric Survey of Nearby Star-forming Galaxies", Moustakas, J. & Kennicutt, Robert C., ApJS, 164, 81</a></li>'
    printf, lun, '<br />'
    printf, lun, '<li><a href="http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2006ApJ...642..775M&amp;db_key=AST&amp;data_type=HTML&amp;format=&amp;high=452a3b59ab21525">"Optical Star Formation Rate Indicators", Moustakas, J., Kennicutt, Robert C., & Tremonti, C., ApJ, 642, 775</a></li>'
    printf, lun, '<br />'
    printf, lun, '<li><a href="http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2006ApJ...651..155M&amp;db_key=AST&amp;data_type=HTML&amp;format=&amp;high=452a3b59ab21525">"Integrated Nebular Abundances of Disk Galaxies", Moustakas, J. & Kennicutt, Robert C., ApJ, 651, 155</a></li>'
    printf, lun, '</ul>'
    printf, lun, '</div>'
    printf, lun, '<br clear="bottom"/>'

    printf, lun, '<p class="email">Last update '+im_today()+'.  Email <a href="mailto:john.moustakas@gmail.com">'+$
      'John Moustakas</a> with questions, comments, or to report errors.</p>'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

    popd

; -------------------------
; ATLASDOWNLOAD.HTML
; -------------------------

    integrated_fitdate = integrated[0].fitdate
    nuclear_fitdate = nuclear[0].fitdate

    pushd, webpath
    
    splog, 'Generating the download page.'
    openw, lun, 'atlasdownload.html', /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<link rel="stylesheet" type="text/css" href="'+httpstyle+'" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '<title>Spectral Atlas: Download</title>'
    printf, lun, '</head>'
    printf, lun, '<body>'
    printf, lun, '<h1>Download</h1>'
    printf, lun, '<br />'
    printf, lun, '<div class="box download">'
    printf, lun, '<h2>Spectrophotometry</h2>'
    printf, lun, '<ul>'
    printf, lun, '<li><a href="atlas_integrated_spectra_ascii.tar.gz">'+$
      'Integrated ASCII spectra</a><span class="smaller"> (gzipped tarball) ['+im_today()+']</span></li>'
    printf, lun, '<li><a href="atlas_integrated_spectra_fits.tar.gz">'+$
      'Integrated FITS spectra</a><span class="smaller"> (gzipped tarball) ['+im_today()+']</span></li>'
    printf, lun, '<br />'
    printf, lun, '<li><a href="atlas_nuclear_spectra_ascii.tar.gz">'+$
      'Nuclear ASCII spectra</a><span class="smaller"> (gzipped tarball) ['+im_today()+']</span></li>'
    printf, lun, '<li><a href="atlas_nuclear_spectra_fits.tar.gz">'+$
      'Nuclear FITS spectra</a><span class="smaller"> (gzipped tarball) ['+im_today()+']</span></li>'
    printf, lun, '</ul>'
;   printf, lun, '<br />'
    printf, lun, '<h2>Spectral Fitting Results</h2>'
    printf, lun, '<ul>'

    printf, lun, '<li><a href="atlas_integrated_speclinefit_'+version+'.fits.gz">Integrated Spectral Fits</a>'+$
      '<span class="smaller"> (<a href="atlas_integrated_speclinefit_'+version+'_nodust.fits.gz">Reddening-corrected</a>)</span>'+$
      '<span class="smaller"> (<a href="atlas_integrated_speclinefit_'+version+'.ps.gz">Postscript</a>) ['+$
      integrated_fitdate+']</span></li>'
    printf, lun, '<li><a href="atlas_nuclear_speclinefit_'+version+'.fits.gz">Nuclear Spectral Fits</a>'+$
      '<span class="smaller"> (<a href="atlas_nuclear_speclinefit_'+version+'_nodust.fits.gz">Reddening-corrected</a>)</span>'+$
      '<span class="smaller"> (<a href="atlas_nuclear_speclinefit_'+version+'.ps.gz">Postscript</a>) ['+$
      nuclear_fitdate+']</span></li>'

;   printf, lun, '<li><a href="atlas_integrated_speclinefit_'+version+'.fits.gz">Integrated Spectral Fits</a><span class="smaller">'
;   printf, lun, '(gzipped binary FITS table) ['+integrated_fitdate+']</span></li>'
;   printf, lun, '<li><a href="atlas_integrated_speclinefit_'+version+'.ps.gz">Quality Assurance Plots</a><span class="smaller">'
;   printf, lun, '(gzipped postscript) ['+integrated_fitdate+']</span></li>'
;   printf, lun, '<br />'
;   printf, lun, '<li><a href="atlas_nuclear_speclinefit_'+version+'.fits.gz">Nuclear Spectral Fits</a><span class="smaller">'
;   printf, lun, '(gzipped binary FITS table) ['+nuclear_fitdate+']</span></li>'
;   printf, lun, '<li><a href="atlas_nuclear_speclinefit_'+version+'.ps.gz">Quality Assurance Plots</a><span class="smaller">'
;   printf, lun, '(gzipped postscript) ['+nuclear_fitdate+']</span></li>'

;   printf, lun, '<li><a href="README.results">Fitting Results README</a></li>'
    printf, lun, '</ul>'
    printf, lun, '</div>'
    printf, lun, '<br clear="bottom"/>'

    printf, lun, '<p class="email">Last update '+im_today()+'.  Email <a href="mailto:john.moustakas@gmail.com">'+$
      'John Moustakas</a> with questions, comments, or to report errors.</p>'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

    popd

; -------------------------
; ATLASINTERNAL.HTML
; -------------------------

    pushd, webpath
    
    splog, 'Generating the internal page.'
    openw, lun, 'atlasinternal.html', /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<link rel="stylesheet" type="text/css" href="'+httpstyle+'" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '<title>Spectral Atlas: Internal</title>'
    printf, lun, '</head>'
    printf, lun, '<body>'
    printf, lun, '<h1>Analysis</h1>'
    printf, lun, '<br />'
    printf, lun, '<div class="box download">'
    printf, lun, '<ul>'
    printf, lun, '<li><a href="analysis/properties.html">Sample Properties</a></li>'
    printf, lun, '<li><a href="analysis/sfrs.html">Optical Star-Formation Rate Indicators</a></li>'
    printf, lun, '<li><a href="analysis/ir_sfrs.html">Infrared Star-Formation Rates</a></li>'
    printf, lun, '<li><a href="analysis/abundances.html">Galaxy Abundance Diagnostics</a></li>'
    printf, lun, '<li><a href="analysis/zindicators.html">HII-Region Abundance Diagnostics</a></li>'
    printf, lun, '<li><a href="analysis/zintegrated.html">Integrated Nebular Abundances of Disk Galaxies</a></li>'
    printf, lun, '<br />'
    printf, lun, '<li><a href="analysis/tremonti_compare.html">Comparison with Christy</a></li>'
    printf, lun, '<li><a href="analysis/rob.html">Plots for Rob</a></li>'
    printf, lun, '</ul>'
    printf, lun, '</div>'
;   printf, lun, '<h1>Spectroscopic Sub-Samples</h1>'
;   printf, lun, '<br />'
;   printf, lun, '<div class="box">'
;   printf, lun, '<table width="100%" style="font-size: 125%;" border="1">'
;   printf, lun, '<tbody>'
;   printf, lun, '<tr>'
;   printf, lun, '<td>Spectral Atlas</td>'
;   printf, lun, '<td><span class="smaller"><a href="subsample_atlas.fits.gz">FITS</a></span></td>'
;   printf, lun, '<td><span class="smaller"><a href="subsample_atlas.dat">ASCII</a></span></td>'
;   printf, lun, '<td><span class="smaller"><a href="subsample_atlas.html">HTML</span></a></td>'
;   printf, lun, '</tr>'
;   printf, lun, '<tr>'
;   printf, lun, '<td>NFGS</td>'
;   printf, lun, '<td><span class="smaller"><a href="subsample_nfgs.fits.gz">FITS</a></span></td>'
;   printf, lun, '<td><span class="smaller"><a href="subsample_nfgs.dat">ASCII</a></span></td>'
;   printf, lun, '<td><span class="smaller"><a href="subsample_nfgs.html">HTML</span></a></td>'
;   printf, lun, '</tr>'
;   printf, lun, '<tr>'
;   printf, lun, '<td>Ursa Major Cluster</td>'
;   printf, lun, '<td><span class="smaller"><a href="subsample_ursamajor.fits.gz">FITS</a></span></td>'
;   printf, lun, '<td><span class="smaller"><a href="subsample_ursamajor.dat">ASCII</a></span></td>'
;   printf, lun, '<td><span class="smaller"><a href="subsample_ursamajor.html">HTML</span></a></td>'
;   printf, lun, '</tr>'
;   printf, lun, '<tr>'
;   printf, lun, '<td>11 Mpc Survey</td>'
;   printf, lun, '<td><span class="smaller"><a href="subsample_11mpc.fits.gz">FITS</a></span></td>'
;   printf, lun, '<td><span class="smaller"><a href="subsample_11mpc.dat">ASCII</a></span></td>'
;   printf, lun, '<td><span class="smaller"><a href="subsample_11mpc.html">HTML</span></a></td>'
;   printf, lun, '</tr>'
;   printf, lun, '</tbody>'
;   printf, lun, '</table>'
;   printf, lun, '</div>'
;   printf, lun, '<br />'

; these are the iSPEC2d visualizations, the links to which have been
; removed to save disk space in PUBLIC_HTML; the code below can be
; easily restored (jm05aug03uofa)
    
;   printf, lun, '<h1>iSPEC2d Visualizations</h1>'
;   printf, lun, '<br />'
;   printf, lun, '<div class="box reductions">'
;   printf, lun, '<br />'
;   printf, lun, '<table class="props" width="100%" style="font-size: 100%;" border="1">'
;   printf, lun, '<caption>Integrated & Nuclear Spectrophotometry</caption>'
;   printf, lun, '<tbody>'
;   printf, lun, '<tr>'
;   printf, lun, '<td><a href="redux/98mar.html">1998 March</a><br /><span class="smaller"> [March 21]</span></td>'
;   printf, lun, '<td><a href="redux/98apr.html">1998 April</a><br /><span class="smaller"> [April 27 - May 2]</span></td>'
;   printf, lun, '</tr>'
;   printf, lun, '<tr>'
;   printf, lun, '<td><a href="redux/98jun.html">1998 June</a><br /><span class="smaller"> [June 24 - 26]</span></td>'
;   printf, lun, '<td><a href="redux/98oct.html">1998 October</a><br /><span class="smaller"> [October 15 - 19]</span></td>'
;   printf, lun, '</tr>'
;   printf, lun, '<tr>'
;   printf, lun, '<td><a href="redux/99apr.html">1999 April</a><br /><span class="smaller"> [April 19; 21 - 23]</span></td>'
;   printf, lun, '<td><a href="redux/99may.html">1999 May</a><br /><span class="smaller"> [May 8 - 11]</span></td>'
;   printf, lun, '</tr>'
;   printf, lun, '<tr>'
;   printf, lun, '<td><a href="redux/99nov.html">1999 November</a><br /><span class="smaller"> [November 2 - 4; 29]</span></td>'
;   printf, lun, '<td><a href="redux/00apr.html">2000 April</a><br /><span class="smaller"> [March 29 - April 3]</span></td>'
;   printf, lun, '</tr>'
;   printf, lun, '<tr>'
;   printf, lun, '<td><a href="redux/02apr.html">2002 April</a><br /><span class="smaller"> [April 13 - 15]</span></td>'
;   printf, lun, '<td><a href="redux/03may.html">2003 May</a><br /><span class="smaller"> [May 26 - 27]</span></td>'
;   printf, lun, '</tr>'
;   printf, lun, '</tbody>'
;   printf, lun, '</table>'
;   printf, lun, '<br />'
;   printf, lun, '<table class="props" width="100%" style="font-size: 100%;" border="1">'
;   printf, lun, '<caption>Turner Mergers</caption>'
;   printf, lun, '<tbody>'
;   printf, lun, '<tr>'
;   printf, lun, '<td><a href="redux/94nov.html">1994 November</a><br /><span class="smaller"> [November 28]</span></td>'
;   printf, lun, '<td><a href="redux/95mar.html">1995 March</a><br /><span class="smaller"> [March 25 - 26]</span></td>'
;   printf, lun, '</tr>'
;   printf, lun, '<tr>'
;   printf, lun, '<td><a href="redux/95oct.html">1995 October</a><br /><span class="smaller"> [October 17 - 18]</span></td>'
;   printf, lun, '<td><a href="redux/96apr.html">1996 April</a><br /><span class="smaller"> [April 22]</span></td>'
;   printf, lun, '</tr>'
;   printf, lun, '<tr>'
;   printf, lun, '<td><a href="redux/97apr.html">1997 April</a><br /><span class="smaller"> [April 10]</span></td>'
;   printf, lun, '<td>&nbsp</td>'
;   printf, lun, '</tr>'
;   printf, lun, '</tbody>'
;   printf, lun, '</table>'
;   printf, lun, '<br />'
;   printf, lun, '<table class="props" width="100%" style="font-size: 100%;" border="1">'
;   printf, lun, '<caption>Spitzer/SINGS Spectrophotometry</caption>'
;   printf, lun, '<tbody>'
;   printf, lun, '<tr>'
;   printf, lun, '<td><a href="redux/01nov.html">2001 November</a><br /><span class="smaller"> [November 10 - 13]</span></td>'
;   printf, lun, '<td><a href="redux/01dec.html">2001 December</a><br /><span class="smaller"> [December 20 - 24]</span></td>'
;   printf, lun, '</tr>'
;   printf, lun, '<tr>'
;   printf, lun, '<td><a href="redux/02feb.html">2002 February</a><br /><span class="smaller"> [February 7 - 10]</span></td>'
;   printf, lun, '<td><a href="redux/02may.html">2002 May</a><br /><span class="smaller"> [May 11 - 13]</span></td>'
;   printf, lun, '</tr>'
;   printf, lun, '</tbody>'
;   printf, lun, '</table>'
;   printf, lun, '</div>'
;   printf, lun, '<br clear="bottom"/>'

    printf, lun, '<p class="email">Last update '+im_today()+'.  Email <a href="mailto:john.moustakas@gmail.com">'+$
      'John Moustakas</a> with questions, comments, or to report errors.</p>'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

    popd

; ---------------------------------------------------------------------------
; Generate tarballs of all the FITS and ASCII spectra
; ---------------------------------------------------------------------------

    if keyword_set(make_tarballs) then begin

       splog, 'Generating tarballs.'

; ##########       
; integrated       
; ##########       
       
       spec1dlist = strjoin([atlas.drift_file],' ')
       ascii1dlist = strjoin([atlas.drift_asciifile],' ')

       pushd, spec1dpath & spawn, ['tar cvzf '+webpath+'atlas_integrated_spectra_fits.tar.gz '+spec1dlist], /sh & popd
       pushd, ascii1dpath & spawn, ['tar cvzf '+webpath+'atlas_integrated_spectra_ascii.tar.gz '+ascii1dlist], /sh & popd

; ##########       
; nuclear
; ##########       
       
       spec1dlist = strjoin([atlas.nuclear_file],' ')
       ascii1dlist = strjoin([atlas.nuclear_asciifile],' ')
       
       pushd, spec1dpath & spawn, ['tar cvzf '+webpath+'atlas_nuclear_spectra_fits.tar.gz '+spec1dlist], /sh & popd
       pushd, ascii1dpath & spawn, ['tar cvzf '+webpath+'atlas_nuclear_spectra_ascii.tar.gz '+ascii1dlist], /sh & popd

    endif
       
; ---------------------------------------------------------------------------
; Generate symbolic links to the FITS and ASCII spectra
; ---------------------------------------------------------------------------

    if keyword_set(make_links) then begin

       splog, 'Generating symbolic links.'

; integrated spectral fitting results       

       integrated_mjdstr = integrated[0].mjdstr
       psfile_integrated = specfitpath+integrated_mjdstr+'_integrated_atlas_specfit.ps.gz'

       spawn, ['/bin/cp -fp '+specfitpath_data+'integrated_atlas_speclinefit_'+version+'.fits.gz'+$
         ' '+webpath+'atlas_integrated_speclinefit_'+version+'.fits.gz'], /sh
       spawn, ['/bin/cp -fp '+specfitpath_data+'integrated_atlas_speclinefit_'+version+'_nodust.fits.gz'+$
         ' '+webpath+'atlas_integrated_speclinefit_'+version+'_nodust.fits.gz'], /sh
       spawn, ['/bin/cp -fp '+psfile_integrated+' '+webpath+'atlas_integrated_speclinefit_'+version+'.ps.gz'], /sh

; nuclear spectral fitting results       
       
       nuclear_mjdstr = nuclear[0].mjdstr
       psfile_nuclear = specfitpath+nuclear_mjdstr+'_nuclear_atlas_specfit.ps.gz'
 
       spawn, ['/bin/cp -fp '+specfitpath_data+'nuclear_atlas_speclinefit_'+version+'.fits.gz'+$
         ' '+webpath+'atlas_nuclear_speclinefit_'+version+'.fits.gz'], /sh
       spawn, ['/bin/cp -fp '+specfitpath_data+'nuclear_atlas_speclinefit_'+version+'_nodust.fits.gz'+$
         ' '+webpath+'atlas_nuclear_speclinefit_'+version+'_nodust.fits.gz'], /sh
       spawn, ['/bin/cp -fp '+psfile_nuclear+' '+webpath+'atlas_nuclear_speclinefit_'+version+'.ps.gz'], /sh

; ##########       
; integrated
; ##########       
       
; tarballs of all the spectra and links to the individual spectra       
       
;      spawn, ['/bin/cp -fp '+spec1dpath+'atlas_integrated_spectra_fits.tar.gz '+$
;        webpath+'atlas_integrated_spectra_fits.tar.gz'], /sh
;      spawn, ['/bin/cp -fp '+ascii1dpath+'atlas_integrated_spectra_ascii.tar.gz '+$
;        webpath+'atlas_integrated_spectra_ascii.tar.gz'], /sh

       for k = 0L, n_elements(integrated)-1L do begin

          fitsfile = strtrim(integrated[k].drift_file,2)
          asciifile = strtrim(integrated[k].drift_asciifile,2)

          spawn, ['/bin/cp -fp '+spec1dpath+fitsfile+' '+fitspath+fitsfile], /sh
          spawn, ['/bin/cp -fp '+ascii1dpath+asciifile+' '+asciipath+asciifile], /sh
          
       endfor

; ##########       
; nuclear
; ##########       
       
;      spawn, ['/bin/cp -fp '+spec1dpath+'atlas_nuclear_spectra_fits.tar.gz '+webpath+'atlas_nuclear_spectra_fits.tar.gz'], /sh
;      spawn, ['/bin/cp -fp '+ascii1dpath+'atlas_nuclear_spectra_ascii.tar.gz '+webpath+'atlas_nuclear_spectra_ascii.tar.gz'], /sh

       for k = 0L, n_elements(nuclear)-1L do begin

          fitsfile = strtrim(nuclear[k].drift_file,2)
          asciifile = strtrim(nuclear[k].drift_asciifile,2)

          spawn, ['/bin/cp -fp '+spec1dpath+fitsfile+' '+fitspath+fitsfile], /sh
          spawn, ['/bin/cp -fp '+ascii1dpath+asciifile+' '+asciipath+asciifile], /sh
          
       endfor

    endif
       
; ---------------------------------------------------------------------------
; Generate PNG files of all the spectra and spectral fits.
; ---------------------------------------------------------------------------

    if keyword_set(make_png) then begin

       splog, 'Generating spectral PNG files.'

       scale = 1E15
       ytitle = 'Flux [10^{-15} '+flam_units()+']'
       xtitle = 'Rest Wavelength [\AA]'

       postthick = 3.0
       charsize = 1.3

       sigrej = 5.0
       smoothbox = 3.0
       
       stime0 = systime(1)

       for k = 0L, ngalaxy-1L do begin

          if (not keyword_set(debug)) then $
            print, format='("Object ",I0,"/",I0,".",A1,$)', k+1, ngalaxy, string(13b)

          if (not keyword_set(debug)) then set_plot, 'Z'

          atlas_display_spectrum, atlas[k], lcharsize=1.2, labeltype=3L, /with_nuclear

          if keyword_set(debug) then begin
             print, strupcase(galaxy[k])
             cc = get_kbrd(1) 
          endif else begin
             img = tvrd()
             tvlct, r, g, b, /get
             write_png, pngpath+specpng[k], img, r, g, b
             set_plot, 'X'
          endelse

       endfor 

       splog, format='("Total time = ",G0," minutes.")', (systime(1)-stime0)/60.0
    
    endif
    
; ---------------------------------------------------------------------------
; Generate PNG files of the visualizations
; ---------------------------------------------------------------------------

    if keyword_set(make_thumbs) then begin

       splog, 'Generating image PNG files.'

       postthick = 3.0
       postthick2 = 2.0

       stime0 = systime(1)

;      for k = 306, 308 do begin
       for k = 0L, ngalaxy-1L do begin

          if (not keyword_set(debug)) then $
            print, format='("Object ",I0,"/",I0,".",A1,$)', k+1, ngalaxy, string(13b)

          if (not keyword_set(debug)) then set_plot, 'Z'

          atlas_display_image, atlas[k], imagepath=dsspath, imposition=pos, $
            lcharsize=1.5, pcharsize=pcharsize, /preserve_aspect, $
            pspath=pspath, labeltype=1L, _extra=extra, postscript=postscript, $
            barlabelcolor='black'

; capture the image and write out             

          if keyword_set(debug) then begin
             print, strupcase(galaxy[k])
             cc = get_kbrd(1) 
          endif else begin
             x0 = fix((convert_coord(pos[0:1],/normal,/to_device))[0])
             nx = fix((convert_coord(pos[2:3],/normal,/to_device))[0])-x0
             y0 = fix((convert_coord(pos[0:1],/normal,/to_device))[1])
             ny = fix((convert_coord(pos[2:3],/normal,/to_device))[1])-y0

             img = tvrd(x0+1L,y0+1L-1L,nx-1L,ny-1L)
             
             tvlct, r, g, b, /get
             write_png, thumbpath+thumbpng[k], img, r, g, b
             set_plot, 'X'
          endelse

          delvarx, pos

       endfor 

       splog, format='("Total time = ",G0," minutes.")', (systime(1)-stime0)/60.0
    
    endif

return
end
