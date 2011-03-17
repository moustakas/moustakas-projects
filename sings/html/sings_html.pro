;+
; NAME:
;       SINGS_HTML
;
; PURPOSE:
;       Generate web page visualizations for SINGS.
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
;       jm05jul27uofa - updated
;-

pro sings_populate_frames, sings, irs, sings_drift56, sings_drift20, $
  sings_nuclear, suffix=suffix, http=http

    webpath = sings_path(/web)
    thumbpath = webpath+'thumb/'
    masspath = webpath+'mass/'

    httphome = '"'+http+'"'
    baseref = http+'research/sings/'
    httpstyle = baseref+'sings.css'

    if (n_elements(suffix) eq 0L) then suffix = '' else $
      if (strcompress(suffix,/remove) ne '') then suffix = '_'+suffix

    galaxy = strtrim(strlowcase(sings.galaxy),2)
    nedgalaxy = repstr(strtrim(sings.ned_galaxy,2),'+','%2B')
    nicegalaxy = strtrim(sings.nice_galaxy,2)
    ngalaxy = n_elements(galaxy)

    htmlfile = 'html/'+galaxy+suffix+'.html'

    specpng = '../png/'+galaxy+'_spec.png'
    spec56png = '../png/'+galaxy+'_drift56.png'
    spec20png = '../png/'+galaxy+'_drift20.png'
    specnucpng = '../png/'+galaxy+'_nuclear.png'
    thumbpng = '../thumb/'+galaxy+'_thumb.png'

    specfitsname = '../fits/'+galaxy+'.fits'
    specasciiname = '../ascii/'+galaxy+'.ascii'

; -------------------------
; WEBSINGS_SUFFIX.HTML
; -------------------------

    mainhtml = 'websings'+suffix+'.html'
    tablehtml = 'websings_table'+suffix+'.html'
    datahtml = 'websings_data'+suffix+'.html'

    pushd, webpath

    splog, 'Writing '+mainhtml+'.'
    openw, lun, mainhtml, /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
    printf, lun, '<html>'
    printf, lun, '<head>'
;   printf, lun, '<base href="'+webpath+'" />'
    printf, lun, '<link rel="stylesheet" type="text/css" href="'+httpstyle+'" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '<title>SINGS Spectra</title>'
    printf, lun, '</head>'
    printf, lun, '<frameset cols="20%,*">'
    printf, lun, '<frame src="'+tablehtml+'"></frame>'
    printf, lun, '<frame src="'+datahtml+'" name="data"></frame>'
    printf, lun, '</frameset>'
    printf, lun, '</html>'
    free_lun, lun

; -------------------------
; WEBSINGS_TABLE_SUFFIX.HTML
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
;   printf, lun, '<br />'
    printf, lun, '<h1>SINGS</h1>'
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
; WEBSINGS_DATA_SUFFIX.HTML
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
    printf, lun, '<h1>SINGS</h1>'
    printf, lun, '<br />'
    printf, lun, '<p>Please select a galaxy.</p>'
;   printf, lun, '<p>Matched-aperture optical spectroscopy in the 3600-7000 Angstrom spectral range has been '+$
;     'obtained using the Steward Observatory 2.3-meter and the CTIO 1.5-meter telescopes.  '+$
;     'See <a href="'+$
;     'http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2003PASP..115..928K&amp ;db_key=AST&amp;high=419b7ae4ee20795>'+$
;     'Kennicutt et al. 2003</a> for details.  Please select an object from the menu on the left.</p>'
    printf, lun, '<br />'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

; ---------------------------------------------------------------------------
; Generate individual HTML files for each object
; ---------------------------------------------------------------------------

; galaxy properties

    galtags = [$
      'ra', $
      'dec', $
      'sings_type',$
      'd25_maj',$
      'd25_min',$
      'pa',     $
      'distance',$
      'rc3_m_b',$
      'twomass_m_ks'$
      ]
    galformats = [$
      'A0',$
      'A0',$
      'A0',$
      'F12.1',$
      'F12.1',$
      'I0',   $
      'F12.1',$
      'F12.1',$
      'F12.1' $
      ]
    gallabels = [$
      'RA',$
      'DEC',$
      'Type',$
      'D25',$
      'd25',$
      'PA',$
      'D',$
      'M(B)',$
      'M(Ks)'$
      ]
    galunits = [$
      '[J2000]',$
      '[J2000]',$
      '&nbsp',$
      '[arcmin]',$
      '[arcmin]',$
      '[degree]',$
      '[Mpc]',$
      '[mag]',$
      '[mag]'$
      ]

    galdata = html_structure_parse(sings,tags=galtags,tagformats=galformats,tagindices=galindices)
    gallabels = gallabels[galindices] & galunits = galunits[galindices]
    ngaltags = n_tags(galdata)

; IRS observing parameters; generate an empty structure for objects
; that have not been observed with IRS

;   irsmatched = im_empty_structure(irs,ncopies=ngalaxy)
;   irsmatched.galaxy = sings.galaxy
;   indx = cmset_op(strtrim(irsmatched.galaxy,2),'AND',strtrim(irs.nice_galaxy,2),/index)
;   
;   junk = irsmatched[indx] & copy_struct, irs, junk & irsmatched[indx] = junk

    irstags = [$
      'date' ,$
      'ra',   $
      'dec',  $
      'step' ,$
      'pa'   ,$
      'length_2x',$
      'width'$
      ]
    irsformats = [$
      'A0',$
      'A0',$
      'A0',$
      'A0',$
      'I0',$
      'F12.1',$
      'F12.1' $
      ]
    irslabels = [   $
      'Observed',   $
      'RA',         $
      'DEC',        $
      'Step',       $
      'PA',         $
      'Length (2x)',$
      'Width'       $
      ]
    irsunits = [$
      '&nbsp',$
      '[J2000]',$
      '[J2000]',$
      '&nbsp',$
      '[degree]',$
      '[arcsec]',$
      '[arcsec]' $
      ]

    mirs = irs
    mirs.pa = im_angle_format(mirs.pa) ; NOTE!
    mirs.width = mirs.width*60.0
    mirs.length_2x = mirs.length_2x*60.0
    
    irsdata = html_structure_parse(mirs,tags=irstags,tagformats=irsformats,tagindices=irsindices)
    irslabels = irslabels[irsindices] & irsunits = irsunits[irsindices]
    nirstags = n_tags(irsdata)

; optical spectroscopic observing parameters; the last tag has to be
; COMMENTS! 

    spectags56 = [$
      'drift56_date',     $
      'drift56_strap',       $
      'drift56_posangle',       $
      'drift56_comments'  $
      ]
    spectags20 = [$
      'drift20_date',     $
      'drift20_strap',       $
      'drift20_posangle',       $
      'drift20_comments'  $
      ]
    spectagsnuc = [$
      'nuclear_date',     $
      'nuclear_strap',       $
      'nuclear_posangle',       $
      'nuclear_comments'  $
      ]
    specformats = [$
      'A0',$
      'A0',$
      'I0',$
      'A0' $
      ]
    speclabels = [  $
      'Observed',   $
      'Aperture',      $
      'PA',         $
      'Comments'    $
      ]
    specunits = [$
      '&nbsp',$
      '[arcsec x arcsec]',$
      '[degree]',$
      '&nbsp' $
      ]

    specdata56 = html_structure_parse(sings,tags=spectags56,tagformats=specformats,tagindices=specindices)
    specdata20 = html_structure_parse(sings,tags=spectags20,tagformats=specformats,tagindices=specindices)
    specdatanuc = html_structure_parse(sings,tags=spectagsnuc,tagformats=specformats,tagindices=specindices)

    specdata56.drift56_date = reform(strtrim(sings.drift56_observat,2),1,n_elements(sings))+'/'+specdata56.drift56_date
    specdata20.drift20_date = reform(strtrim(sings.drift20_observat,2),1,n_elements(sings))+'/'+specdata20.drift20_date
    specdatanuc.nuclear_date = reform(strtrim(sings.nuclear_observat,2),1,n_elements(sings))+'/'+specdatanuc.nuclear_date

    w56 = where(sings.drift56 eq 0L,nw56)
    if (nw56 ne 0L) then begin
       specdata56[w56].drift56_date = '&nbsp'
       specdata56[w56].drift56_posangle = '&nbsp'
    endif
    w20 = where(sings.drift20 eq 0L,nw20)
    if (nw20 ne 0L) then begin
       specdata20[w20].drift20_date = '&nbsp'
       specdata20[w20].drift20_posangle = '&nbsp'
    endif
    wnuc = where(sings.nuclear eq 0L,nnuc)
    if (nnuc ne 0L) then begin
       specdatanuc[wnuc].nuclear_date = '&nbsp'
       specdatanuc[wnuc].nuclear_posangle = '&nbsp'
    endif

    speclabels = speclabels[specindices] & specunits = specunits[specindices]
    nspectags = n_tags(specdata56)

; make the individual pages    
    
    ra = repstr(sings.ra,':','%3A')
    dec = repstr(repstr(sings.dec,':','%3A'),'+','%2B') ; negative declinations are not accounted for here!
    
    splog, 'Generating the individual SINGS webpages.'
    for i = 0L, ngalaxy-1L do begin

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
         nicegalaxy[i]+'</a></h2>'
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
       printf, lun, '<a href="../index.html" target="_top">Home</a>';<br />'
       printf, lun, '</div><br clear="bottom"/>'

; image + spectral visualizations
       
       printf, lun, '<table width="100%" align="center" border="0" cellspacing="0" cellpadding="0">' ; open table
       printf, lun, '<tbody>'

       printf, lun, '<tr>'
       printf, lun, '<td><a class="nolinkcolor" href="'+thumbpng[i]+'" target="_top">'+$
         '<img width="99%" src="'+thumbpng[i]+'" alt="'+nicegalaxy[i]+' Thumbnail" /></a></td>'
       printf, lun, '<td width="50%"><a class="nolinkcolor" href="'+specpng[i]+'" target="_top">'+$
         '<img width="99%" src="'+specpng[i]+'" alt="'+nicegalaxy[i]+' Thumbnail" /></a></td>'
       printf, lun, '</tr>'

       printf, lun, '</tbody>' 
       printf, lun, '</table>' ; close table

       printf, lun, '<br clear="bottom"/>'
       printf, lun, '<p class="smaller">The yellow outline shows the configuration of the IRS LL/SL '+$
         'map with 2x coverage, when available.  The red and blue rectangles show the '+$
         'apertures corresponding to the radial strip and circumnuclear spectra, respectively.</p>'

; galaxy properties table
       
       printf, lun, '<table width="100%" border="2" cellspacing="0" cellpadding="2" class="props">'
       printf, lun, '<caption>Galaxy Properties</caption>'
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
       
; IRS observing parameters
       
       match, strlowcase(galaxy[i]), strlowcase(strtrim(irs.nice_galaxy,2)), irsindx1, irsindx2

       printf, lun, '<table width="100%" border="2" cellspacing="0" cellpadding="2" class="props">'
       printf, lun, '<caption>IRS LL/SL Spectroscopy</caption>'
       printf, lun, '<thead>'
       printf, lun, '<tr>'
       
       for itag = 0L, nirstags-1L do printf, lun, '<th>'+irslabels[itag]+'</th>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       for itag = 0L, nirstags-1L do printf, lun, '<th><span class="smaller">'+irsunits[itag]+'</span></th>'
       printf, lun, '</tr>'
       printf, lun, '</thead>'
       printf, lun, '<tbody>'
       printf, lun, '<tr>'
       if (irsindx1[0] eq -1L) then $
         for itag = 0L, nirstags-1L do printf, lun, '<td>&nbsp</td>' else $
         for itag = 0L, nirstags-1L do printf, lun, '<td>'+(irsdata[irsindx2].(itag))[0]+'</td>'
       printf, lun, '</tr>'
       printf, lun, '</tbody>'
       printf, lun, '</table>' 
       printf, lun, '<br clear="bottom"/>' ; close the IRS table

; optical spectroscopic observing parameters
       
       match, strlowcase(galaxy[i]), strlowcase(strtrim(sings.galaxy,2)), specindx1, specindx2

       printf, lun, '<table width="100%" border="2" cellspacing="0" cellpadding="2" class="props">'
       printf, lun, '<caption>Optical Spectroscopy</caption>'
       printf, lun, '<thead>'
       printf, lun, '<tr>'
       printf, lun, '<th>Spectrum</th>'
       for itag = 0L, nspectags-1L do printf, lun, '<th>'+speclabels[itag]+'</th>'
       printf, lun, '<th colspan="3">Download</th>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       printf, lun, '<th>&nbsp</th>'
       for itag = 0L, nspectags-1L do printf, lun, '<th><span class="smaller">'+specunits[itag]+'</span></th>'
       printf, lun, '<th colspan="2">1D</th><th>2D</th>'
       printf, lun, '</tr>'
       printf, lun, '</thead>'

       printf, lun, '<tbody>'
; drift56
       printf, lun, '<tr>'
       printf, lun, '<td>Radial Strip</td>'       
       if (specindx1[0] eq -1L) then begin
;         for itag = 0L, nspectags-1L do printf, lun, '<td>&nbsp</td>' ; no spectroscopic data at all
          for itag = 0L, nspectags-2L do printf, lun, '<td>&nbsp</td>' ; no spectroscopic data at all
          printf, lun, '<td>unavailable</td>'                ; "comments" field
          printf, lun, '<td>&nbsp</td><td>&nbsp</td><td>&nbsp</td>'    ; no ascii/fits/fits2d spectrum
       endif else begin
          for itag = 0L, nspectags-1L do printf, lun, '<td>'+(specdata56[specindx2].(itag))[0]+'</td>'
          asciifile1 = strtrim(sings[specindx2].drift56_asciifile,2)
          if file_test(webpath+'ascii/'+asciifile1,/regular) then $
            printf, lun, '<td><a href="../ascii/'+asciifile1+'">ASCII</a></td>' else $
            printf, lun, '<td>&nbsp</td>' ; no ascii spectrum
          fitsfile1 = strtrim(sings[specindx2].drift56_file,2)
          if file_test(webpath+'fits/'+fitsfile1,/regular) then $

            printf, lun, '<td><a href="../fits/'+fitsfile1+'">FITS</a></td>' else $
            printf, lun, '<td>&nbsp</td>' ; no fits spectrum
          fits2dfile1 = strtrim(sings[specindx2].drift56_2dfile,2)
          if file_test(webpath+'fits/'+fits2dfile1,/regular) then $
            printf, lun, '<td><a href="../fits/'+fits2dfile1+'">FITS</a></td>' else $
            printf, lun, '<td>&nbsp</td>' ; no 2D fits spectrum
       endelse
       printf, lun, '</tr>'
; drift20
       printf, lun, '<tr>'
       printf, lun, '<td>Circumnuclear</td>'       
       if (specindx1[0] eq -1L) then begin
;         for itag = 0L, nspectags-1L do printf, lun, '<td>&nbsp</td>'
          for itag = 0L, nspectags-2L do printf, lun, '<td>&nbsp</td>' ; no spectroscopic data at all
          printf, lun, '<td>unavailable</td>'                ; "comments" field
          printf, lun, '<td>&nbsp</td><td>&nbsp</td><td>&nbsp</td>'
       endif else begin
          for itag = 0L, nspectags-1L do printf, lun, '<td>'+(specdata20[specindx2].(itag))[0]+'</td>'
          asciifile1 = strtrim(sings[specindx2].drift20_asciifile,2)
          if file_test(webpath+'ascii/'+asciifile1,/regular) then $
            printf, lun, '<td><a href="../ascii/'+asciifile1+'">ASCII</a></td>' else $
            printf, lun, '<td>&nbsp</td>' ; no ascii spectrum
          fitsfile1 = strtrim(sings[specindx2].drift20_file,2)
          if file_test(webpath+'fits/'+fitsfile1,/regular) then $
            printf, lun, '<td><a href="../fits/'+fitsfile1+'">FITS</a></td>' else $
            printf, lun, '<td>&nbsp</td>' ; no fits spectrum
          fits2dfile1 = strtrim(sings[specindx2].drift20_2dfile,2)
          if file_test(webpath+'fits/'+fits2dfile1,/regular) then $
            printf, lun, '<td><a href="../fits/'+fits2dfile1+'">FITS</a></td>' else $
            printf, lun, '<td>&nbsp</td>' ; no 2D fits spectrum
       endelse
       printf, lun, '</tr>'
; nuclear
       printf, lun, '<tr>'
       printf, lun, '<td>Nuclear</td>'       
       if (specindx1[0] eq -1L) then begin
;         for itag = 0L, nspectags-1L do printf, lun, '<td>&nbsp</td>' 
          for itag = 0L, nspectags-2L do printf, lun, '<td>&nbsp</td>' ; no spectroscopic data at all
          printf, lun, '<td>unavailable</td>'                ; "comments" field
          printf, lun, '<td>&nbsp</td><td>&nbsp</td><td>&nbsp</td>'
       endif else begin
          for itag = 0L, nspectags-1L do printf, lun, '<td>'+(specdatanuc[specindx2].(itag))[0]+'</td>'
          asciifile1 = strtrim(sings[specindx2].nuclear_asciifile,2)
          if file_test(webpath+'ascii/'+asciifile1,/regular) then $
            printf, lun, '<td><a href="../ascii/'+asciifile1+'">ASCII</a></td>' else $
            printf, lun, '<td>&nbsp</td>' ; no ascii spectrum
          fitsfile1 = strtrim(sings[specindx2].nuclear_file,2)
          if file_test(webpath+'fits/'+fitsfile1,/regular) then $
            printf, lun, '<td><a href="../fits/'+fitsfile1+'">FITS</a></td>' else $
            printf, lun, '<td>&nbsp</td>' ; no fits spectrum
          fits2dfile1 = strtrim(sings[specindx2].nuclear_2dfile,2)
          if file_test(webpath+'fits/'+fits2dfile1,/regular) then $
            printf, lun, '<td><a href="../fits/'+fits2dfile1+'">FITS</a></td>' else $
            printf, lun, '<td>&nbsp</td>' ; no 2D fits spectrum
       endelse
       printf, lun, '</tr>'
       printf, lun, '</tbody>'
       printf, lun, '</table>' 
       printf, lun, '<br clear="bottom"/>' ; close the spectroscopy table
       
;; download the data
;       
;       printf, lun, '<div class="box download">'
;       printf, lun, '<h2>Download</h2>'
;       printf, lun, '<ul>'
;       printf, lun, '<li><a href="'+specfitsname[i]+'">FITS  Spectrum</a></li>'
;       printf, lun, '<li><a href="'+specasciiname[i]+'">ASCII Spectrum</a></li>'
;       printf, lun, '</ul>'
;       printf, lun, '</div>'
       
; finish the page       
       
       printf, lun, '<div class="menu right">'
       if i eq ngalaxy-1L then nextfile = '../'+htmlfile[0] else nextfile = '../'+htmlfile[i+1L]
       printf, lun, '<a href="'+nextfile+'">Next</a><br />'
       if i eq 0L then backfile = '../'+htmlfile[ngalaxy-1L] else backfile = '../'+htmlfile[i-1L]
       printf, lun, '<a href="'+backfile+'">Back</a>'
       printf, lun, '</div>'

       printf, lun, '<div class="menu left">'
       printf, lun, '<a href="../index.html" target="_top">Home</a><br />'
       printf, lun, '</div><br /><br clear="bottom"/>'

       printf, lun, '<div class="menu right">Page '+string(i+1,format='(I0)')+'</div><br />'
       
       printf, lun, '<p class="email">Last update '+im_today()+'.  Email <a href="mailto:john.moustakas@gmail.com">'+$
         'John Moustakas</a> with questions, comments, or to report errors.</p>'

       printf, lun, '</body>'
       printf, lun, '</html>'
       free_lun, lun
       
    endfor
    
    popd

return
end
    
pro sings_html, sings, sings_drift56, sings_drift20, sings_nuclear, make_png=make_png, $
  make_thumbs=make_thumbs, make_links=make_links, make_tarballs=make_tarballs, $
  debug=debug, rasort=rasort

; datapaths    
    
    datapath = sings_path(/analysis)
    obspath = sings_path(/observing)
    dsspath = sings_path(/dss)
    webpath = sings_path(/web)
    spec1dpath = sings_path(/spec1d)
;   spec2dpath = sings_path(/spec2d)+'rectified/' ; <-- NOTE!!
    ascii1dpath = sings_path(/ascii)
    specfitpath = sings_path(/specfit)

    htmlpath = webpath+'html/'
    pngpath = webpath+'png/'
    thumbpath = webpath+'thumb/'
    asciipath = webpath+'ascii/'
    fitspath = webpath+'fits/'

    http = 'http://sdss.physics.nyu.edu/ioannis/'
    httphome = '"'+http+'"'
    baseref = http+'research/sings/'
    httpstyle = baseref+'sings.css'

    rasort = 1L ; NOTE!
    
    version = sings_version(/specfit)
    info_version = sings_version(/ancillary)

    if (n_elements(sings) eq 0L) then begin
       splog, 'Reading the SINGS information structure.'
       sings = mrdfits(datapath+'sings_info_'+info_version+'.fits.gz',1,/silent)
       if keyword_set(rasort) then sings = sings[sort(sings.ra)] else $ ; ra-sorted
         sings = sings[sort(sings.galaxy)] ; alphabetize; alphabetize
    endif
    if (n_elements(sings_drift56) eq 0L) then begin
       splog, 'Reading the DRIFT56 spectral fitting results.'
       sings_drift56 = read_sings(/drift56)
    endif
    if (n_elements(sings_drift20) eq 0L) then begin
       splog, 'Reading the DRIFT20 spectral fitting results.'
       sings_drift20 = read_sings(/drift20)
    endif
    if (n_elements(sings_nuclear) eq 0L) then begin
       splog, 'Reading the NUCLEAR spectral fitting results.'
       sings_nuclear = read_sings(/nuclear)
    endif

    if (file_test(obspath+'irs_ll_params.txt',/regular) eq 0L) then begin
       splog, 'IRS file '+obspath+'irs_ll_params.txt not found.'
       return
    endif
    irs = rsex(obspath+'irs_ll_params.txt')
    ngalaxy = n_elements(sings)

    galaxy = strtrim(strlowcase(sings.galaxy),2)

    htmlfile = 'html/'+galaxy+'.html'

    specpng = galaxy+'_spec.png'
    specps = galaxy+'_spec.ps'
    spec56png = galaxy+'_drift56.png'
    spec56ps = galaxy+'_drift56.ps'
    spec20png = galaxy+'_drift20.png'
    spec20ps = galaxy+'_drift20.ps'
    specnucpng = galaxy+'_nuclear.png'
    specnucps = galaxy+'_nuclear.ps'

    thumbfits = galaxy+'.fits.gz'
    thumbpng = galaxy+'_thumb.png'

; ---------------------------------------------------------------------------
; define subsamples and generate a WEBSINGS.HTML, WEBSINGS_TABLE.HTML
; and WEBSINGS_DATA.HTML file for each subsample   
; ---------------------------------------------------------------------------

; full sample    
    
    sings_populate_frames, sings, irs, sings_drift56, sings_drift20, $
      sings_nuclear, suffix=suffix, http=http

; -------------------------
; INDEX.HTML
; -------------------------

    pushd, webpath
    
    splog, 'Generating the SINGS home page.'
    openw, lun, 'index.html', /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<link rel="stylesheet" type="text/css" href="'+httpstyle+'" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '<title>Spitzer/SINGS Optical Spectroscopy</title>'
    printf, lun, '</head>'
    printf, lun, '<body>'
    printf, lun, '<br /><br />'
    printf, lun, '<table width="100%" align="center" border="0" cellspacing="0" cellpadding="0">' ; open table
    printf, lun, '<tbody>'
    printf, lun, '<tr>'
    printf, lun, '<td width="40%" class="normal">'+$
      '<a class="nolinkcolor" href="http://www.spitzer.caltech.edu/Media/releases/ssc2004-12/index.shtml">'+$
      '<img width="100%" src="'+baseref+'ngc7331_sings.jpg" alt="NGC7331" /></a></td>'
    printf, lun, '<td class="normal">'
    printf, lun, '<div class="border">'
    printf, lun, '<p>Welcome to the Spitzer Infrared Nearby Galaxies Survey '+$
      '(SINGS) optical spectroscopy web pages.</p>'
    printf, lun, '<ul>'
    printf, lun, '<li><a href="http://sings.stsci.edu">SINGS Legacy Home Page</a></li>'
    printf, lun, '<li><a href="'+baseref+'websings.html">Optical Spectroscopy Visualizations</a></li>'
    printf, lun, '<li><a href="'+baseref+'websings_download.html">Download</a></li>'

;   printf, lun, '<li><span class="gold">Download</span></li>'
    printf, lun, '</ul>'
;   printf, lun, '<p>Download:</p>'
;   printf, lun, '<ul>'
;   printf, lun, '<li><a href="fits/sings_optspec1d_ascii.tar.gz"><span class="smaller">One-dimensional ASCII spectra</span></a></li>'
;   printf, lun, '<li><a href="fits/sings_optspec1d.tar.gz"><span class="smaller">One-dimensional FITS spectra</span></a></li>'
;   printf, lun, '<li><a href="fits/sings_optspec2d.tar.gz"><span class="smaller">Two-dimensional FITS spectra</span></a></li>'
;   printf, lun, '</ul>'
    printf, lun, '</div>'
    printf, lun, '</td>'
    printf, lun, '</tr>'
;   printf, lun, '<div class="widebox">'
;   printf, lun, '<p>Welcome to the Spitzer/SINGS optical spectroscopy web pages.</p>'
;   printf, lun, '<ul>'
;   printf, lun, '<li><a href="'+baseref+'websings.html">Visualizations</a></li>'
;   printf, lun, '</ul>'
;   printf, lun, '<br />'
;   printf, lun, '</div>' ; close the row
    printf, lun, '</tbody>'
    printf, lun, '</table>'
    printf, lun, '<br clear="bottom"/>'

    printf, lun, '<p class="email">Last update '+im_today()+'.  Email <a href="mailto:john.moustakas@gmail.com">'+$
      'John Moustakas</a> with questions, comments, or to report errors.</p>'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

    popd

; -------------------------
; WEBSINGS_DOWNLOAD.HTML
; -------------------------

    drift56_fitdate = sings_version(/specfit)
    drift20_fitdate = sings_version(/specfit)
    nuclear_fitdate = sings_version(/specfit)
;   drift56_fitdate = sings_drift56[0].fitdate
;   drift20_fitdate = sings_drift20[0].fitdate
;   nuclear_fitdate = sings_nuclear[0].fitdate
    
    pushd, webpath
    
    splog, 'Generating the download page.'
    openw, lun, 'websings_download.html', /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<link rel="stylesheet" type="text/css" href="'+httpstyle+'" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '<title>SINGS Optical Spectroscopy: Download</title>'
    printf, lun, '</head>'
    printf, lun, '<body>'
    printf, lun, '<h1>Download</h1>'
    printf, lun, '<br />'
    printf, lun, '<div class="box download">'
    printf, lun, '<h2>Spectrophotometry</h2>'
    printf, lun, '<ul>'
    printf, lun, '<li><a href="sings_optspec1d_ascii.tar.gz">'+$
      'One-dimensional ASCII spectra</a><span class="smaller"> (gzipped tarball) ['+im_today()+']</span></li>'
    printf, lun, '<li><a href="sings_optspec1d.tar.gz">'+$
      'One-dimensional FITS spectra</a><span class="smaller"> (gzipped tarball) ['+im_today()+']</span></li>'
;   printf, lun, '<li><a href="sings_optspec2d.tar.gz">'+$
;     'Two-dimensional FITS spectra</a><span class="smaller"> (gzipped tarball) ['+im_today()+']</span></li>'
    printf, lun, '</ul>'
    printf, lun, '<h2>Spectral Fitting Results</h2>'
    printf, lun, '<ul>'

    printf, lun, '<li><a href="sings_drift56_speclinefit_'+version+'.fits.gz">Radial Strip Fitting Results</a>'+$
;     '<span class="smaller"> (<a href="sings_drift56_speclinefit_'+version+'_nodust.fits.gz">Reddening-corrected</a>)</span>'+$
      '<span class="smaller"> (<a href="sings_drift56_specfit_'+version+'.ps.gz">Postscript</a>) ['+drift56_fitdate+']</span></li>'

    printf, lun, '<li><a href="sings_drift20_speclinefit_'+version+'.fits.gz">Circumnuclear Fitting Results</a>'+$
;     '<span class="smaller"> (<a href="sings_drift20_speclinefit_'+version+'_nodust.fits.gz">Reddening-corrected</a>)</span>'+$
      '<span class="smaller"> (<a href="sings_drift20_specfit_'+version+'.ps.gz">Postscript</a>) ['+drift20_fitdate+']</span></li>'

    printf, lun, '<li><a href="sings_nuclear_speclinefit_'+version+'.fits.gz">Nuclear Fitting Results</a>'+$
;     '<span class="smaller"> (<a href="sings_nuclear_speclinefit_'+version+'_nodust.fits.gz">Reddening-corrected</a>)</span>'+$
      '<span class="smaller"> (<a href="sings_nuclear_specfit_'+version+'.ps.gz">Postscript</a>) ['+nuclear_fitdate+']</span></li>'
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

; ---------------------------------------------------------------------------
; Generate tarballs of all the FITS and ASCII spectra
; ---------------------------------------------------------------------------

    if keyword_set(make_tarballs) then begin

       splog, 'Generating tarballs.'

; only tar up the 2D spectra that correspond to the 1D spectra
       
       spec1dlist = strjoin([sings.drift56_file,sings.drift20_file,sings.nuclear_file],' ')
       spec2dlist = strjoin([sings.drift56_2dfile,sings.drift20_2dfile,sings.nuclear_2dfile],' ')
       ascii1dlist = strjoin([sings.drift56_asciifile,sings.drift20_asciifile,sings.nuclear_asciifile],' ')
       
       pushd, spec1dpath & spawn, ['tar cvzf '+webpath+'sings_optspec1d.tar.gz '+spec1dlist], /sh & popd
;      pushd, spec2dpath & spawn, ['tar cvzf '+webpath+'sings_optspec2d.tar.gz '+spec2dlist], /sh & popd
       pushd, ascii1dpath & spawn, ['tar cvzf '+webpath+'sings_optspec1d_ascii.tar.gz '+ascii1dlist], /sh & popd

    endif
       
; ---------------------------------------------------------------------------
; Generate symbolic links to the FITS and ASCII spectra
; ---------------------------------------------------------------------------

    if keyword_set(make_links) then begin

       splog, 'Generating symbolic links.'

; spectral fitting results

;      drift56_fitdate = sings_drift56[0].fitdate
;      drift20_fitdate = sings_drift20[0].fitdate
;      nuclear_fitdate = sings_nuclear[0].fitdate
;      drift56_mjdstr = sings_drift56[0].mjdstr
;      drift20_mjdstr = sings_drift20[0].mjdstr
;      nuclear_mjdstr = sings_nuclear[0].mjdstr

; HACK!       
       drift56_mjdstr = max(strmid(file_basename(file_search(specfitpath+$
         version+'/*drift56*_specfit.ps.gz')),0,5))
       drift20_mjdstr = max(strmid(file_basename(file_search(specfitpath+$
         version+'/*drift20*_specfit.ps.gz')),0,5))
       nuclear_mjdstr = max(strmid(file_basename(file_search(specfitpath+$
         version+'/*nuclear*_specfit.ps.gz')),0,5))

       psfile_drift56 = specfitpath+version+'/'+drift56_mjdstr+'_sings_drift56_all_specfit.ps.gz'
       psfile_drift20 = specfitpath+version+'/'+drift20_mjdstr+'_sings_drift20_all_specfit.ps.gz'
       psfile_nuclear = specfitpath+version+'/'+nuclear_mjdstr+'_sings_nuclear_all_specfit.ps.gz'

; instead of symbolic links make copies so that the webpages are
; self-contained        
       
       spawn, ['/bin/cp -fp '+specfitpath+'sings_drift56_speclinefit_'+version+'.fits.gz'+$
         ' '+webpath+'sings_drift56_speclinefit_'+version+'.fits.gz'], /sh
;      spawn, ['/bin/cp -fp '+specfitpath+'sings_drift56_speclinefit_'+version+'_nodust.fits.gz'+$
;        ' '+webpath+'sings_drift56_speclinefit_'+version+'_nodust.fits.gz'], /sh
       spawn, ['/bin/cp -fp '+psfile_drift56+' '+webpath+'sings_drift56_specfit_'+version+'.ps.gz'], /sh

       spawn, ['/bin/cp -fp '+specfitpath+'sings_drift20_speclinefit_'+version+'.fits.gz'+$
         ' '+webpath+'sings_drift20_speclinefit_'+version+'.fits.gz'], /sh
;      spawn, ['/bin/cp -fp '+specfitpath+'sings_drift20_speclinefit_'+version+'_nodust.fits.gz'+$
;        ' '+webpath+'sings_drift20_speclinefit_'+version+'_nodust.fits.gz'], /sh
       spawn, ['/bin/cp -fp '+psfile_drift20+' '+webpath+'sings_drift20_specfit_'+version+'.ps.gz'], /sh

       spawn, ['/bin/cp -fp '+specfitpath+'sings_nuclear_speclinefit_'+version+'.fits.gz'+$
         ' '+webpath+'sings_nuclear_speclinefit_'+version+'.fits.gz'], /sh
;      spawn, ['/bin/cp -fp '+specfitpath+'sings_nuclear_speclinefit_'+version+'_nodust.fits.gz'+$
;        ' '+webpath+'sings_nuclear_speclinefit_'+version+'_nodust.fits.gz'], /sh
       spawn, ['/bin/cp -fp '+psfile_nuclear+' '+webpath+'sings_nuclear_specfit_'+version+'.ps.gz'], /sh

; individual spectra       
       
;      spawn, ['/bin/cp -fp '+spec1dpath+'sings_optspec1d.tar.gz '+fitspath+'sings_optspec1d.tar.gz'], /sh
;      spawn, ['/bin/cp -fp '+spec2dpath+'sings_optspec2d.tar.gz '+fitspath+'sings_optspec2d.tar.gz'], /sh
;      spawn, ['/bin/cp -fp '+ascii1dpath+'sings_optspec1d_ascii.tar.gz '+fitspath+'sings_optspec1d_ascii.tar.gz'], /sh
       
       for k = 0L, n_elements(sings)-1L do begin

          if sings[k].nuclear then begin
             fitsfile = strtrim(sings[k].nuclear_file,2)
             fits2dfile = strtrim(sings[k].nuclear_2dfile,2)
             asciifile = strtrim(sings[k].nuclear_asciifile,2)
             spawn, ['/bin/cp -fp '+spec1dpath+fitsfile+' '+fitspath+fitsfile], /sh
;            spawn, ['/bin/cp -fp '+spec2dpath+fits2dfile+' '+fitspath+fits2dfile], /sh
             spawn, ['/bin/cp -fp '+ascii1dpath+asciifile+' '+asciipath+asciifile], /sh
          endif

          if sings[k].drift20 then begin
             fitsfile = strtrim(sings[k].drift20_file,2)
             fits2dfile = strtrim(sings[k].drift20_2dfile,2)
             asciifile = strtrim(sings[k].drift20_asciifile,2)
             spawn, ['/bin/cp -fp '+spec1dpath+fitsfile+' '+fitspath+fitsfile], /sh
;            spawn, ['/bin/cp -fp '+spec2dpath+fits2dfile+' '+fitspath+fits2dfile], /sh
             spawn, ['/bin/cp -fp '+ascii1dpath+asciifile+' '+asciipath+asciifile], /sh
          endif

          if sings[k].drift56 then begin
             fitsfile = strtrim(sings[k].drift56_file,2)
             fits2dfile = strtrim(sings[k].drift56_2dfile,2)
             asciifile = strtrim(sings[k].drift56_asciifile,2)
             spawn, ['/bin/cp -fp '+spec1dpath+fitsfile+' '+fitspath+fitsfile], /sh
;            spawn, ['/bin/cp -fp '+spec2dpath+fits2dfile+' '+fitspath+fits2dfile], /sh
             spawn, ['/bin/cp -fp '+ascii1dpath+asciifile+' '+asciipath+asciifile], /sh
          endif

       endfor 

;       spawn, ['ln -fs '+specfitpath+'sings_drift56_speclinefit.fits.gz'+$
;         ' '+webpath+'sings_drift56_speclinefit.fits.gz'], /sh
;       spawn, ['ln -fs '+specfitpath+'sings_drift56_speclinefit_nodust.fits.gz'+$
;         ' '+webpath+'sings_drift56_speclinefit_nodust.fits.gz'], /sh
;       spawn, ['ln -fs '+psfile_drift56+' '+webpath+'sings_drift56_specfit.ps.gz'], /sh
;
;       spawn, ['ln -fs '+specfitpath+'sings_drift20_speclinefit.fits.gz'+$
;         ' '+webpath+'sings_drift20_speclinefit.fits.gz'], /sh
;       spawn, ['ln -fs '+specfitpath+'sings_drift20_speclinefit_nodust.fits.gz'+$
;         ' '+webpath+'sings_drift20_speclinefit_nodust.fits.gz'], /sh
;       spawn, ['ln -fs '+psfile_drift20+' '+webpath+'sings_drift20_specfit.ps.gz'], /sh
;
;       spawn, ['ln -fs '+specfitpath+'sings_nuclear_speclinefit.fits.gz'+$
;         ' '+webpath+'sings_nuclear_speclinefit.fits.gz'], /sh
;       spawn, ['ln -fs '+specfitpath+'sings_nuclear_speclinefit_nodust.fits.gz'+$
;         ' '+webpath+'sings_nuclear_speclinefit_nodust.fits.gz'], /sh
;       spawn, ['ln -fs '+psfile_nuclear+' '+webpath+'sings_nuclear_specfit.ps.gz'], /sh
;
;; individual spectra       
;       
;       spawn, ['ln -fs '+spec1dpath+'sings_optspec1d.tar.gz '+fitspath+'sings_optspec1d.tar.gz'], /sh
;       spawn, ['ln -fs '+spec2dpath+'sings_optspec2d.tar.gz '+fitspath+'sings_optspec2d.tar.gz'], /sh
;       spawn, ['ln -fs '+ascii1dpath+'sings_optspec1d_ascii.tar.gz '+fitspath+'sings_optspec1d_ascii.tar.gz'], /sh
;       
;       for k = 0L, n_elements(sings)-1L do begin
;
;          if sings[k].nuclear then begin
;             fitsfile = strtrim(sings[k].nuclear_file,2)
;             fits2dfile = strtrim(sings[k].nuclear_2dfile,2)
;             asciifile = strtrim(sings[k].nuclear_asciifile,2)
;             spawn, ['ln -fs '+spec1dpath+fitsfile+' '+fitspath+fitsfile], /sh
;             spawn, ['ln -fs '+spec2dpath+fits2dfile+' '+fitspath+fits2dfile], /sh
;             spawn, ['ln -fs '+ascii1dpath+asciifile+' '+asciipath+asciifile], /sh
;          endif
;
;          if sings[k].drift20 then begin
;             fitsfile = strtrim(sings[k].drift20_file,2)
;             fits2dfile = strtrim(sings[k].drift20_2dfile,2)
;             asciifile = strtrim(sings[k].drift20_asciifile,2)
;             spawn, ['ln -fs '+spec1dpath+fitsfile+' '+fitspath+fitsfile], /sh
;             spawn, ['ln -fs '+spec2dpath+fits2dfile+' '+fitspath+fits2dfile], /sh
;             spawn, ['ln -fs '+ascii1dpath+asciifile+' '+asciipath+asciifile], /sh
;          endif
;
;          if sings[k].drift56 then begin
;             fitsfile = strtrim(sings[k].drift56_file,2)
;             fits2dfile = strtrim(sings[k].drift56_2dfile,2)
;             asciifile = strtrim(sings[k].drift56_asciifile,2)
;             spawn, ['ln -fs '+spec1dpath+fitsfile+' '+fitspath+fitsfile], /sh
;             spawn, ['ln -fs '+spec2dpath+fits2dfile+' '+fitspath+fits2dfile], /sh
;             spawn, ['ln -fs '+ascii1dpath+asciifile+' '+asciipath+asciifile], /sh
;          endif
;
;       endfor 

    endif
       
; ---------------------------------------------------------------------------
; Generate PNG files of all the spectra and spectral fits.
; ---------------------------------------------------------------------------

    if keyword_set(make_png) then begin

       splog, 'Generating spectral PNG files.'

       scale = 1E15
       ytitle = 'Flux (10^{-15} '+flam_units()+')'
       xtitle = 'Rest Wavelength (\AA)'

       postthick = 2.0
       charsize = 1.3

       sigrej = 5.0
       smoothbox = 3.0
       
       stime0 = systime(1)

;      for k = 3L, ngalaxy-1L do begin
       for k = 0L, ngalaxy-1L do begin

          if (not keyword_set(debug)) then $
            print, format='("Object ",I0,"/",I0,".",A1,$)', k+1, ngalaxy, string(13b)

          match, strlowcase(galaxy[k]), strlowcase(strtrim(sings.galaxy,2)), indx1, indx2
          if (indx1[0] eq -1L) then $
            specinfo = {drift56: 0L, drift20: 0L, nuclear: 0L} else $
            specinfo = sings[indx2]

; generate a 3-panel visualization of the three spectra

          if (not keyword_set(debug)) then set_plot, 'Z'
          
          polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')
          pagemaker, nx=1, ny=3, xspace=0.0, yspace=0.0, xmargin=[1.2,0.2], $
            ymargin=[0.3,1.4], position=pos, /normal

          if specinfo.drift56 then begin

             specdata = rd1dspec(strtrim(specinfo.drift56_file,2),datapath=spec1dpath,/silent)
             wave = specdata.wave
             flux = scale*specdata.spec
             smoothflux = smooth(flux,smoothbox)

             refwave = wave
             
             stats = im_stats(flux,sigrej=sigrej)
;            stats = im_stats(smoothflux,sigrej=sigrej)
;            yrange = [stats.min,stats.maxrej]
;            yrange = [stats.minrej,stats.maxrej]
             yrange = [stats.minrej,stats.max]
             xrange = minmax(wave)
             
             plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=2.0, xtickname=replicate(' ',10), $
               ythick=2.0, xtitle='', ytitle='', charsize=charsize, charthick=2.0, $
               yrange=yrange, xrange=xrange, thick=2.0, color=djs_icolor('white'), position=pos[*,0], $
               yminor=3
             djs_oplot, wave, flux, ps=10, thick=0.8, color='red'

             label = ['Radial Strip']
             legend, textoidl(label), /left, /top, box=0, charsize=charsize, $
               charthick=postthick, color=djs_icolor('white')

          endif else begin

             plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=2.0, xtickname=replicate(' ',10), $
               ythick=2.0, xtitle='', ytitle='', charsize=charsize, charthick=2.0, $
               position=pos[*,0], ytickname=replicate(' ',10), yticklen=1E-10, xticklen=1E-10
             xyouts, (pos[2,0]-pos[0,0])/2.0+pos[0,0], (pos[3,0]-pos[1,0])/2.0+pos[1,0], $
               'Spectrum Unavailable', align=0.5, charthick=postthick, /normal
             
          endelse

          if specinfo.drift20 then begin

             specdata = rd1dspec(strtrim(specinfo.drift20_file,2),datapath=spec1dpath,/silent)
             wave = specdata.wave
             flux = scale*specdata.spec

             if (specinfo.drift56 eq 0L) then refwave = wave else begin
                flux = interpol(flux,wave,refwave)
                wave = refwave
             endelse
             
             smoothflux = smooth(flux,smoothbox)
             
             stats = im_stats(flux,sigrej=sigrej)
;               stats = im_stats(smoothflux,sigrej=sigrej)
;               yrange = [stats.min,stats.maxrej]
;               yrange = [stats.minrej,stats.maxrej]
             yrange = [stats.minrej,stats.max]
             xrange = minmax(wave)
             
             plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=2.0, xtickname=replicate(' ',10), $
               ythick=2.0, xtitle='', ytitle='', charsize=charsize, charthick=2.0, /noerase, $
               yrange=yrange, xrange=xrange, thick=2.0, color=djs_icolor('white'), position=pos[*,1], $
               yminor=3
             djs_oplot, wave, flux, ps=10, thick=0.8, color='blue'

             label = ['Circumnuclear']
             legend, textoidl(label), /left, /top, box=0, charsize=charsize, $
               charthick=postthick, color=djs_icolor('white')

          endif else begin

             plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=2.0, xtickname=replicate(' ',10), $
               ythick=2.0, xtitle='', ytitle='', charsize=charsize, charthick=2.0, /noerase, $
               position=pos[*,1], ytickname=replicate(' ',10), yticklen=1E-10, xticklen=1E-10
             xyouts, (pos[2,1]-pos[0,1])/2.0+pos[0,1], (pos[3,1]-pos[1,1])/2.0+pos[1,1], $
               'Spectrum Unavailable', align=0.5, charthick=postthick, /normal
             
          endelse
          
          if specinfo.nuclear then begin

             specdata = rd1dspec(strtrim(specinfo.nuclear_file,2),datapath=spec1dpath,/silent)
             wave = specdata.wave
             flux = scale*specdata.spec

             if (specinfo.drift56 eq 0L) and (specinfo.drift20 eq 0L) then refwave = wave else begin
                flux = interpol(flux,wave,refwave)
                wave = refwave
             endelse
             
             smoothflux = smooth(flux,smoothbox)
             
             stats = im_stats(flux,sigrej=sigrej)
;            stats = im_stats(smoothflux,sigrej=sigrej)
;            yrange = [stats.min,stats.maxrej]
;            yrange = [stats.minrej,stats.maxrej]
             yrange = [stats.minrej,stats.max]
             xrange = minmax(wave)
             
             plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=2.0, $
               ythick=2.0, xtitle=textoidl(xtitle), ytitle='', charsize=charsize, charthick=2.0, /noerase, $
               yrange=yrange, xrange=xrange, thick=2.0, color=djs_icolor('white'), position=pos[*,2], $
               yminor=3
             djs_oplot, wave, flux, ps=10, thick=0.8, color='dark green'

             label = ['Nuclear']
             legend, textoidl(label), /left, /top, box=0, charsize=charsize, $
               charthick=postthick, color=djs_icolor('white')

          endif else begin

             if (specinfo.drift56 eq 0L) and (specinfo.drift20 eq 0L) then begin
                
                plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=2.0, xtickname=replicate(' ',10), $
                  ythick=2.0, xtitle='', ytitle='', charsize=charsize, charthick=2.0, /noerase, $
                  position=pos[*,2], ytickname=replicate(' ',10), yticklen=1E-10, xticklen=1E-10
                
             endif else begin

                plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=2.0, xrange=xrange, $
                  ythick=2.0, xtitle=textoidl(xtitle), ytitle='', charsize=charsize, charthick=2.0, /noerase, $
                  position=pos[*,2], ytickname=replicate(' ',10), yticklen=1E-10, xticklen=1E-10
                
             endelse

             xyouts, (pos[2,2]-pos[0,2])/2.0+pos[0,2], (pos[3,2]-pos[1,2])/2.0+pos[1,2], $
               'Spectrum Unavailable', align=0.5, charthick=postthick, /normal

          endelse
          
; ytitle

          xyouts, 0.3*pos[0,0], (pos[3,0]-pos[1,2])/2.0+pos[1,2], textoidl(ytitle), $
            align=0.5, orientation=90, charsize=charsize, charthick=postthick, /normal

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

       postthick = 2.3

       stime0 = systime(1)

;      for k = 5L, 25L do begin
       for k = 0L, ngalaxy-1L do begin

          print, format='("Object ",I0,"/",I0,".",A1,$)', k+1, ngalaxy, string(13b)

          if (not keyword_set(debug)) then set_plot, 'Z'
          
          sings_display_image, sings[k], imagepath=dsspath, imposition=pos, $
            lcharsize=1.2, pcharsize=pcharsize, /preserve_aspect, /noscanbox, $
            pspath=pspath, labeltype=1L, _extra=extra, postscript=postscript, $
            astr=astr, imageinfo=imageinfo, /arcminlabel, barlabelcolor='black';, /norc3box

          match, strlowcase(galaxy[k]), strlowcase(strtrim(irs.nice_galaxy,2)), indx1, indx2
          if (indx1[0] ne -1L) then begin
             irsinfo = irs[indx2]
             irs_ra = irsinfo.ra & irs_dec = irsinfo.dec
          endif else begin
             irsinfo = -1L
             irs_ra = sings[k].ra & irs_dec = sings[k].dec
          endelse

; visualize the IRS radial strip aperture             
          
          if (size(irsinfo,/type) eq 8L) then $
            im_oplot_box, irsinfo.width, irsinfo.length_2x, im_angle_format(irsinfo.pa), $
              line=0, color=djs_icolor('yellow'), thick=postthick

; visualize the spectroscopic apertures

          match, strlowcase(galaxy[k]), strlowcase(strtrim(sings.galaxy,2)), indx1, indx2
          if (indx1[0] eq -1L) then $
            specinfo = {drift56: 0L, drift20: 0L, nuclear: 0L} else $
            specinfo = sings[indx2]

          if specinfo.drift56 then begin

             gsssadxy, astr, 15.0*im_hms2dec(sings[k].ra), im_hms2dec(sings[k].dec), x56, y56
;            gsssadxy, astr, 15.0*im_hms2dec(specinfo.drift56_ra), im_hms2dec(specinfo.drift56_dec), x56, y56
             x56 = (x56 - imageinfo.xcen)*imageinfo.xpixscale
             y56 = (y56 - imageinfo.ycen)*imageinfo.ypixscale

             im_oplot_box, specinfo.drift56_scan/60.0, specinfo.drift56_ap/60.0, $
               specinfo.drift56_posangle, /noplot, corners=corners

             indx = [0,1,2,3,0]
             djs_oplot, [corners[0,indx],corners[0,indx+1]]+x56, $
               [corners[1,indx],corners[1,indx+1]]+y56, line=0, $
               color=djs_icolor('red'), thick=postthick
             
          endif

          if specinfo.drift20 then begin

             gsssadxy, astr, 15.0*im_hms2dec(sings[k].ra), im_hms2dec(sings[k].dec), x20, y20
;            gsssadxy, astr, 15.0*im_hms2dec(specinfo.drift20_ra), im_hms2dec(specinfo.drift20_dec), x20, y20
             x20 = (x20 - imageinfo.xcen)*imageinfo.xpixscale
             y20 = (y20 - imageinfo.ycen)*imageinfo.ypixscale

             im_oplot_box, specinfo.drift20_scan/60.0, specinfo.drift20_ap/60.0, $
               specinfo.drift20_posangle, /noplot, corners=corners

             indx = [0,1,2,3,0]
             djs_oplot, [corners[0,indx],corners[0,indx+1]]+x20, $
               [corners[1,indx],corners[1,indx+1]]+y20, line=0, $
               color=djs_icolor('blue'), thick=postthick
             
;            im_oplot_box, specinfo.drift20_scan/60.0, specinfo.drift20_ap/60.0, $
;              specinfo.drift20_posangle, line=0, color=djs_icolor('blue'), thick=postthick

          endif
          
; capture the image and write out             

          x0 = fix((convert_coord(pos[0:1],/normal,/to_device))[0])
          nx = fix((convert_coord(pos[2:3],/normal,/to_device))[0])-x0
          y0 = fix((convert_coord(pos[0:1],/normal,/to_device))[1])
          ny = fix((convert_coord(pos[2:3],/normal,/to_device))[1])-y0
          img = tvrd(x0+1L,y0+1L-1L,nx-1L,ny-1L)
          
          tvlct, r, g, b, /get
          write_png, thumbpath+thumbpng[k], img, r, g, b

          delvarx, pos
          
          if keyword_set(debug) then cc = get_kbrd(1) else set_plot, 'X'

       endfor 

       splog, format='("Total time = ",G0," minutes.")', (systime(1)-stime0)/60.0
    
    endif

return
end
