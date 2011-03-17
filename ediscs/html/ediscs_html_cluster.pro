pro ediscs_html_cluster, cluster
; jm04feb17uofa

    if n_elements(cluster) eq 0L then begin
       print, 'Syntax - '
       return
    endif

    doit = execute('datapath = ediscs_path(/'+repstr(cluster,'-','_')+')')
    rootpath = ediscs_path()
    webpath = ediscs_path(/html)
    clpath = webpath+cluster+'/'
    htmlpath = clpath+'html/'
    pngpath = clpath+'png/'
    
    splog, 'Generating the web page for '+cluster+'.'

; make a soft link for the cluster image

    spawn, ['ln -fs '+rootpath+cluster+'/'+cluster+'_color.jpg '+clpath+cluster+'_color.jpg'], /sh
    
; grab FITS file headers

    forage = ediscs_headgrab('*fits.gz',datapath=datapath)
    galaxy = forage.galaxy
    ngalaxy = n_elements(galaxy)

    htmlfile = 'html/'+galaxy+'.html'
    
;   pngname = '../png/'+galaxy+'.png'
;   psname = '../png/'+galaxy+'.ps.gz'
    specnamepng = '../png/'+galaxy+'_spec.png'
    specnameps = '../png/'+galaxy+'_spec.ps.gz'
    spec1namepng = '../png/'+galaxy+'_spec1.png'
    spec1nameps = '../png/'+galaxy+'_spec1.ps.gz'
    spec2namepng = '../png/'+galaxy+'_spec2.png'
    spec2nameps = '../png/'+galaxy+'_spec2.ps.gz'
    spec3namepng = '../png/'+galaxy+'_spec3.png'
    spec3nameps = '../png/'+galaxy+'_spec3.ps.gz'
    spec4namepng = '../png/'+galaxy+'_spec4.png'
    spec4nameps = '../png/'+galaxy+'_spec4.ps.gz'

; read the photometric catalog for this cluster, crop to the relevant
; galaxies, and select specific structure tags (note that not all
; clusters will have the same structure tags)

    splog, 'Reading the photometric catalog for '+cluster+'.'
    cat = ediscs_readcat(cluster)

    doit = match_string(galaxy,cat.old_ediscsid,/exact,index=keep)
    ccat = cat[keep]
;   niceprint, galaxy, ccat.old_ediscsid

    tags = ['ediscsid','mautb','mautv','mauti','mautj','mautk']
    tformats = ['A21','F5.2','F5.2','F5.2','F5.2','F5.2']
    doit = match_string(tags,tag_names(ccat),/exact,index=match,/silent)
    good = where(match ne -1L)
    tags = tags[good]
    tformats = tformats[good]

;   subcat = struct_trimtags(ccat,select=tags)

; setup the index page with frames

    basepath = clpath
    pushd, clpath
    
    splog, 'Writing the FRAMES page for '+cluster+'.'
    openw, lun, cluster+'.html',/get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Frameset//EN" "http://cerebus.as.arizona.edu/~ioannis">' 
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<base href="'+basepath+'" />'
    printf, lun, '<link rel="stylesheet" type="text/css" href="http://cerebus.as.arizona.edu/~ioannis/styles/ediscs.css" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '<title>EDisCS: '+cluster+'</title>'
    printf, lun, '</head>'
    printf, lun, '<frameset cols="22%,*">'
    printf, lun, '<frame src="'+cluster+'_table.html"></frame>'
    printf, lun, '<frame src="'+cluster+'_data.html" name="data"></frame>'
    printf, lun, '</frameset>'
    printf, lun, '</html>'
    free_lun, lun

; setup the table frame

    splog, 'Writing the galaxy list for '+cluster+'.'
    openw, lun, cluster+'_table.html', /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://cerebus.as.arizona.edu/~ioannis">' 
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<base href="'+basepath+'" />'
    printf, lun, '<link rel="stylesheet" type="text/css" href="http://cerebus.as.arizona.edu/~ioannis/styles/ediscs.css" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '</head>'
    printf, lun, '<body>'
    printf, lun, '<h1>Sample</h1>'
    printf, lun, '<p class="menu">'
    printf, lun, '<span class="left"><a href="../index.html" target="_top">Home</a></span>'
    printf, lun, '<br /></p>'
    printf, lun, '<table border="2" width="100%">'
    printf, lun, '<tbody>'
    for i = 0, ngalaxy-1L do $
      printf, lun, '<tr><td class="ediscstable"><a href="'+htmlfile[i]+$
        '" target="data">'+galaxy[i]+'</a></td></tr>'
    printf, lun, '</tbody>'
    printf, lun, '</table>'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

; initialize the data page

    splog, 'Writing the data page for '+cluster+'.'
    openw, lun, cluster+'_data.html', /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://cerebus.as.arizona.edu/~ioannis">'
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<base href="'+basepath+'" />'
    printf, lun, '<link rel="stylesheet" type="text/css" href="http://cerebus.as.arizona.edu/~ioannis/styles/ediscs.css" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '</head>'
    printf, lun, '<body>'
    printf, lun, '<h1>'+cluster+'</h1>'
    printf, lun, '<br />'
    printf, lun, '<img width="90%" src="'+cluster+'_color.jpg'+' />'
;   printf, lun, '<p class="centered">Please select an object from the alphabetical list on the left.</p>'
    printf, lun, '<br />'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

; initialize some data structures

; photometric properties

    tags = ['mautb','mautv','mauti','mautj','mautk']
    tformats = ['F5.2','F5.2','F5.2','F5.2','F5.2']
    pplabels = ['B','V','I','J','Ks']
    ppunits = ['mag','mag','mag','mag','mag']
    
    doit = match_string(tags,tag_names(ccat),/exact,index=match,/silent)
    good = where(match ne -1L) & tags = tags[good] & tformats = tformats[good]
    pplabels = pplabels[good] & ppunits = ppunits[good]
    
    pp = struct_trimtags(ccat,select=tags,format=tformats)
    npptags = n_tags(pp)

    for itag = 0L, npptags-1L do begin
       prop = strtrim(reform(pp.(itag)),2)
       no = where((prop gt 90.0) or (prop lt -90),nno,comp=go,ncomp=ngo)
       if nno ne 0L then prop[no] = '-'
       if ngo ne 0L then prop[go] = prop[go]
       pp.(itag) = strtrim(reform(prop,1,ngalaxy),2)
    endfor
    
; write the individual web pages

    splog, 'Writing the individual web pages for '+cluster+'.'
    for i = 0L, ngalaxy-1L do begin

       openw, lun, htmlfile[i], /get_lun
       printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional'+$
         '//EN" "http://cerebus.as.arizona.edu/~ioannis">' 
       printf, lun, '<html>'
       printf, lun, '<head>'
       printf, lun, '<link rel="stylesheet" type="text/css" href="http://cerebus.as.arizona.edu/~ioannis/styles/ediscs.css" />'
       printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
;      printf, lun, '<title>'+galaxy[i]+'</title>'
       printf, lun, '</head>'
       printf, lun, '<body>'
       printf, lun, '<h1>'+galaxy[i]+'</a></h1>'

       printf, lun, '<div class="menu right">'
       if i eq ngalaxy-1L then nextfile = '../'+htmlfile[0] else nextfile = '../'+htmlfile[i+1L]
       printf, lun, '<a href="'+nextfile+'">Next</a><br />'
       if i eq 0L then backfile = '../'+htmlfile[ngalaxy-1L] else backfile = '../'+htmlfile[i-1L]
       printf, lun, '<a href="'+backfile+'">Back</a>'
       printf, lun, '</div>'

       printf, lun, '<div class="menu left">'
       printf, lun, '<a href="../../index.html" target="_top">Home</a><br />'
       printf, lun, '</div><br /><br clear="bottom"/>'

; spectrum

       printf, lun, '<table width="100%" border="0" cellspacing="1" cellpadding="1">'
       printf, lun, '<tbody>'

       printf, lun, '<tr><td width="100%" rowspan="2" valign="top"><a style="color: silver;" href="'+$
         specnameps[i]+'"><img width="90%" src="'+specnamepng[i]+'" alt="'+galaxy[i]+' Thumbnail" /></a></td></tr>'

;      printf, lun, '<tr><td width="100%" align="center" valign="center">'
;      printf, lun, '   <table style="font-size: 150%;" border="1">'
;      printf, lun, '   <tr>'
;      printf, lun, '      <td>Observed Spectrum</td><td><a href="'+spec1name[i]+'">png</a> <a href="'+spec1nameps[i]+'">ps</a></td>'
;      printf, lun, '      <td>Continuum Model</td><td><a href="'+spec2name[i]+'">png</a> <a href="'+spec2nameps[i]+'">ps</a></td>'
;      printf, lun, '   </tr>'
;      printf, lun, '   <tr>'
;      printf, lun, '      <td>Emission Line Model</td><td><a href="'+spec3name[i]+'">png</a> <a href="'+spec3nameps[i]+'">ps</a></td>'
;      printf, lun, '      <td>Composite Model</td><td><a href="'+spec4name[i]+'">png</a> <a href="'+spec4nameps[i]+'">ps</a></td>'
;      printf, lun, '   </tr>'
;      printf, lun, '   </table>'
;      printf, lun, '</td></tr>'

       printf, lun, '</tbody>'
       printf, lun, '</table>'
       printf, lun, '<br clear="bottom"/>'
       printf, lun, '<br />'
       
; photometric properties table       
       
       printf, lun, '<div class="tablecenter">'
       printf, lun, '<table width="90%" border="2" cellspacing="0" cellpadding="2" class="props">'
       printf, lun, '<caption>Photometric Properties</caption>'
       printf, lun, '<thead>'
       printf, lun, '<tr>'
       for itag = 0L, npptags-1L do printf, lun, '<th>'+pplabels[itag]+'</th>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       for itag = 0L, npptags-1L do printf, lun, '<th><span class="smaller">'+ppunits[itag]+'</span></th>'
       printf, lun, '</tr>'
       printf, lun, '</thead>'
       printf, lun, '<tbody>'
       printf, lun, '<tr>'
       for itag = 0L, npptags-1L do printf, lun, '<td>'+pp[i].(itag)+'</td>'
       printf, lun, '</tr>'
       printf, lun, '</tbody>'
       printf, lun, '</table>'
       printf, lun, '</div>'
       printf, lun, '<br /><br />'

       printf, lun, '<p class="menu right">'
       if i eq ngalaxy-1L then nextfile = '../'+htmlfile[0] else nextfile = '../'+htmlfile[i+1L]
       printf, lun, '<a href="'+nextfile+'">Next</a><br />'
       if i eq 0L then backfile = '../'+htmlfile[ngalaxy-1L] else backfile = '../'+htmlfile[i-1L]
       printf, lun, '<a href="'+backfile+'">Back</a>'
       printf, lun, '</p>'
 
       printf, lun, '<div class="menu left">'
       printf, lun, '<a href="../../index.html" target="_top">Home</a><br />'
       printf, lun, '</div><br /><br clear="bottom"/>'

       printf, lun, '<p class="email">Last update '+im_today()+'.  Email <a href="mailto:jmoustakas@as.arizona.edu">'+$
         'J. Moustakas</a> with questions or comments.</p>'

       printf, lun, '</body>'
       printf, lun, '</html>'
       free_lun, lun
       
    endfor
    
    popd
    
return
end
