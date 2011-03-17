pro html_kenn92atlas, webpath, baseref, httphome, httpstyle, nicegalaxy, htmlfile
; jm04may04uofa

    if (n_params() ne 6L) then begin
       print, 'Syntax - html_kenn92atlas, webpath, baseref, httphome, httpstyle'
       return
    endif
    
    pushd, webpath

; ---------------------------------------------------------------------------
; KENN92ATLAS.HTML
; ---------------------------------------------------------------------------

    splog, 'Writing KENN92ATLAS.HTML.'
    openw, lun, 'kenn92atlas.html',/get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Frameset//EN" '+httphome+'>' 
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<base href='+baseref+' />'
    printf, lun, '<link rel="stylesheet" type="text/css" href='+httpstyle+' />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '<title>K92 Spectral Atlas Visualizations</title>'
    printf, lun, '</head>'
    printf, lun, '<frameset cols="20%,*">'
    printf, lun, '<frame src="kenn92table.html"></frame>'
    printf, lun, '<frame src="kenn92intro.html" name="galaxy"></frame>'
    printf, lun, '</frameset>'
    printf, lun, '</html>'
    free_lun, lun

; ---------------------------------------------------------------------------
; KENN92TABLE.HTML
; ---------------------------------------------------------------------------

    splog, 'Writing KENN92TABLE.HTML.'

    splog, 'Writing the kenn92 galaxy list.'
    openw, lun, 'kenn92table.html', /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>' 
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<base href='+baseref+' />'
    printf, lun, '<link rel="stylesheet" type="text/css" href='+httpstyle+' />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '</head>'
    printf, lun, '<body>'
    printf, lun, '<h1>Sample</h1>'
    printf, lun, '<p class="menu">'
    printf, lun, '<span class="left"><a href="../index.html" target="_top">Home</a></span>'
    printf, lun, '<span class="right"><a href="kenn92intro.html" target="galaxy">Intro</a></span>'
    printf, lun, '<br /></p>'
    printf, lun, '<table border="2" width="100%">'
    printf, lun, '<tbody>'
    for i = 0, n_elements(htmlfile)-1L do $
      printf, lun, '<tr><td class="atlastable"><a href="'+htmlfile[i]+$
        '" target="galaxy">'+nicegalaxy[i]+'</a></td></tr>'
    printf, lun, '</tbody>'
    printf, lun, '</table>'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

; ---------------------------------------------------------------------------
; KENN92INTRO.HTML
; ---------------------------------------------------------------------------

    splog, 'Writing KENN92INTRO.HTML.'
    openw, lun, 'kenn92intro.html', /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<base href='+baseref+' />'
    printf, lun, '<link rel="stylesheet" type="text/css" href='+httpstyle+' />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '<style type="text/css">'
    printf, lun, 'body {'
    printf, lun, '   background-image: url(kenn92_poster.png);'
;   printf, lun, '   background-repeat: no-repeat;'
    printf, lun, '}'
    printf, lun, '</style>'
    printf, lun, '</head>'
    printf, lun, '<body>'
    printf, lun, '<h1>K92 Spectral Atlas Visualizations</h1>'
    printf, lun, '<br />'
;   printf, lun, '<p class="centered"><a href="kenn92_poster.ps.gz"><img width="99%" '+$
;     'src="kenn92_poster.png" alt="K92 Atlas Poster"></a></p>'
    printf, lun, '<p class="centered">Welcome to the web page visualizations for the K92 spectral atlas.'
    printf, lun, 'Please select an object from the list on the left.</p>'
    printf, lun, '<br />'
;   printf, lun, '<p class="email">Last update '+im_today()+'.  Email <a href="mailto:jmoustakas@as.arizona.edu">'+$
;     'J. Moustakas</a> with questions or comments.</p>'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun
    
    popd    

return
end
    
