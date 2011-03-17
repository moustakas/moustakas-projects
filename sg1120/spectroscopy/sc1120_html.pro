;+
; NAME:
;       SC1120_HTML
;
; PURPOSE:
;       Generate web page visualizations for the SC1120 spectra. 
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
;       MRDFITS(), MATCH_STRING(), VISUALIZE_ATLAS,
;       VISUALIZE_APERTURES, 
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jan 28, U of A - written
;-

pro sc1120_populate_frames, sc1120, suffix=suffix

    webpath = sc1120_path(/web)
    thumbpath = webpath+'thumb/'

    baseref = 'http://cerebus.as.arizona.edu/~ioannis/research/sc1120/'
    httphome = '"http://cerebus.as.arizona.edu/~ioannis"'
    httpstyle = 'http://cerebus.as.arizona.edu/~ioannis/styles/sc1120.css'

    if (n_elements(suffix) eq 0L) then suffix = '' else $
      if (strcompress(suffix,/remove) ne '') then suffix = '_'+suffix

    galaxy = strtrim(sc1120.galaxy,2)
    nicegalaxy = repstr(galaxy,'_',' ')
    ngalaxy = n_elements(galaxy)

    htmlfile = 'html/'+galaxy+suffix+'.html'
    specpng = '../png/'+galaxy+'_spec.png'
    specps = '../png/'+galaxy+'_spec.ps.gz'

; -------------------------
; WEBSC1120_SUFFIX.HTML
; -------------------------

    mainhtml = 'websc1120'+suffix+'.html'
    tablehtml = 'websc1120_table'+suffix+'.html'
    datahtml = 'websc1120_data'+suffix+'.html'

    pushd, webpath

    splog, 'Writing '+mainhtml+'.'
    openw, lun, mainhtml, /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
    printf, lun, '<html>'
    printf, lun, '<head>'
;   printf, lun, '<base href="'+webpath+'" />'
    printf, lun, '<link rel="stylesheet" type="text/css" href="'+httpstyle+'" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '<title>SC1120 Spectra</title>'
    printf, lun, '</head>'
    printf, lun, '<frameset cols="30%,*">'
    printf, lun, '<frame src="'+tablehtml+'"></frame>'
    printf, lun, '<frame src="'+datahtml+'" name="data"></frame>'
    printf, lun, '</frameset>'
    printf, lun, '</html>'
    free_lun, lun

; -------------------------
; WEBSC1120_TABLE_SUFFIX.HTML
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
    printf, lun, '<p class="menu">'
    printf, lun, '<span class="left"><a href="'+baseref+'index.html" target="_top">Home</a></span>'
    printf, lun, '<br /></p>'
    printf, lun, '<table border="2" width="100%">'
    printf, lun, '<tbody>'
;   for i = 0, ngalaxy-1L do begin
;      printf, lun, '<tr>'
;      printf, lun, '<td class="atlastable"><a href="'+htmlfile[i]+$
;        '" target="data">'+index[i]+'</a></td>'
;      printf, lun, '<td class="atlastable"><a href="'+htmlfile[i]+$
;        '" target="data">'+nicegalaxy[i]+'</a></td>'
;      printf, lun, '</tr>'
;   endfor
    for i = 0, ngalaxy-1L do $
      printf, lun, '<tr><td class="atlastable"><a href="'+htmlfile[i]+$
        '" target="data">'+index[i]+' '+nicegalaxy[i]+'</a></td></tr>'
    printf, lun, '</tbody>'
    printf, lun, '</table>'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

; -------------------------
; WEBSC1120_DATA_SUFFIX.HTML
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
    printf, lun, '<h1>SC1120</h1>'
    printf, lun, '<br />'
    printf, lun, '<p>Select an object from the left to begin.</p>'
    printf, lun, '<br />'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

; ---------------------------------------------------------------------------
; Generate individual HTML files for each object
; ---------------------------------------------------------------------------

; initialize various data structures

    tags = [$
      'GALAXY',$
      'MASK',$
      'QUADRANT',$
      'SLIT',$
      'Q',$
      'H_ALPHA_EW',$
      'H_ALPHA',$
      'H_ALPHA_LUM',$
      'SFR_H_ALPHA'$
      ]
    tformats = [$
      'A0',$
      'A0',$
      'A0',$
      'A0',$
      'I0',$
      'F12.1',$
      'E12.2',$
      'F12.1',$
      'F12.1' $
      ]
    labels = [$
      'Galaxy',$
      'Mask',$
      'Quadrant',$
      'Slit',$
      'Q',$
      'EW(Ha)',$
      'f(Ha)',$
      'log L(Ha)',$
      'log (SFR)'$
      ]
    units = [$
      '',$
      '',$
      '',$
      '',$
      '[Quality]',$
      '[\AA]',$
      '[erg/s/cm2]',$
      '[Lsun]',$
      '[Msun/yr]'$
      ]

    doit = match_string(tags,tag_names(sc1120),/exact,index=match,/silent)
    good = where(match ne -1L) & tags = tags[good] & tformats = tformats[good]
    labels = labels[good] & units = units[good]
    
    rdata = struct_trimtags(sc1120,select=tags)
    data = struct_trimtags(sc1120,select=tags,format=tformats)
    ntags = n_tags(data)

    for itag = 0L, ntags-1L do begin
       property = reform(strtrim(reform(data.(itag)),2))
       rproperty = rdata.(itag)
       dims = size(data.(itag),/dimension)
       if (dims[0] eq 2L) then begin
          no = where((property[1,*] lt 0.0),nno,comp=go,ncomp=ngo) ; flag on the error
          property[1,*] = ''
          if (nno ne 0L) then property[0,no] = '&nbsp'
          data.(itag) = strtrim(reform(property,2,ngalaxy),2)
       endif else begin
          if (size(rproperty,/type) eq 7L) then begin
             data.(itag) = strtrim(reform(property,1,ngalaxy),2)
          endif else begin
             no = where((property lt -900.0),nno,comp=go,ncomp=ngo)
             if (nno ne 0L) then property[no] = '&nbsp'
             data.(itag) = strtrim(reform(property,1,ngalaxy),2)
          endelse
       endelse
    endfor

    ra = repstr(sc1120.ra,':','%3A')
    dec = repstr(repstr(sc1120.dec,':','%3A'),'+','%2B') ; negative declinations are not accounted for here!
    
    splog, 'Generating the individual SC1120 web pages.'
    for i = 0L, ngalaxy-1L do begin

       openw, lun, htmlfile[i], /get_lun
       printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
       printf, lun, '<html>'
       printf, lun, '<head>'
       printf, lun, '<link rel="stylesheet" type="text/css" href="'+httpstyle+'" />'
       printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
       printf, lun, '</head>'
       printf, lun, '<body>'
;      printf, lun, '<h2>'+nicegalaxy[i]+'</h2>'
       printf, lun, '<h2><a target="_top" href="http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?search_type=Near+'+$
         'Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon='+ra[i]+'&lat='+dec[i]+'&radius=1&out_csys=Equatorial'+$
         '&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp='+$
         'YES&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY">'+$
         nicegalaxy[i]+'</a></h2>'

       printf, lun, '<div class="menu right">'
       if i eq ngalaxy-1L then nextfile = '../'+htmlfile[0] else nextfile = '../'+htmlfile[i+1L]
       printf, lun, '<a href="'+nextfile+'">Next</a><br />'
       if i eq 0L then backfile = '../'+htmlfile[ngalaxy-1L] else backfile = '../'+htmlfile[i-1L]
       printf, lun, '<a href="'+backfile+'">Back</a>'
       printf, lun, '</div>'

       printf, lun, '<div class="menu left">'
       printf, lun, '<a href="../index.html" target="_top">Home</a>';<br />'
       printf, lun, '</div><br clear="bottom"/>'

       printf, lun, '<table align="center" border="0" cellspacing="0" cellpadding="0">'
       printf, lun, '<tbody>'

; spectrum       
       
       printf, lun, '<tr><td width="100%" valign="center"><a href="'+specpng[i]+'" target="_top"><img src="'+specpng[i]+$
         '" alt="'+nicegalaxy[i]+' Thumbnail" /></a></td></tr>'

; data properties table
       
       printf, lun, '<tr>'
 
       printf, lun, '<table width="100%" border="2" cellspacing="0" cellpadding="2" class="props">'
       printf, lun, '<thead>'
       printf, lun, '<tr>'
       for itag = 0L, ntags-1L do printf, lun, '<th>'+labels[itag]+'</th>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       for itag = 0L, ntags-1L do printf, lun, '<th><span class="smaller">'+units[itag]+'</span></th>'
       printf, lun, '</tr>'
       printf, lun, '</thead>'
       printf, lun, '<tbody>'
       printf, lun, '<tr>'
       for itag = 0L, ntags-1L do printf, lun, '<td>'+(data[i].(itag))[0]+'</td>'
       printf, lun, '</tr>'
       printf, lun, '</tbody>'
       printf, lun, '</table>' ; close the properties table
 
       printf, lun, '</tr>' ; close the properties row
       printf, lun, '</tbody>'
       printf, lun, '</table>'
       printf, lun, '<br clear="bottom"/>'
       
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
       
       printf, lun, '<p class="email">Last update '+im_today()+'.  Email <a href="mailto:jmoustakas@as.arizona.edu">'+$
         'J. Moustakas</a> with questions or comments.</p>'

       printf, lun, '</body>'
       printf, lun, '</html>'
       free_lun, lun
       
    endfor
    
    popd

return
end
    
pro sc1120_html, sc1120, make_png=make_png, make_thumbs=make_thumbs

; datapaths    
    
    datapath = sc1120_path(/spec1d)
    webpath = sc1120_path(/web)
    htmlpath = webpath+'html/'
    pngpath = webpath+'png/'
    thumbpath = webpath+'thumb/'

    baseref = 'http://cerebus.as.arizona.edu/~ioannis/research/sc1120/'
    httphome = '"http://cerebus.as.arizona.edu/~ioannis"'
    httpstyle = 'http://cerebus.as.arizona.edu/~ioannis/styles/sc1120.css'
    
    splog, 'Reading the spectral fitting results.'
    if (n_elements(sc1120) eq 0L) then sc1120 = read_sc1120(/silent)
;   if (n_elements(sc1120) eq 0L) then sc1120 = sc1120_forage(file_search('*.fits'),datapath=datapath)
    ngalaxy = n_elements(sc1120)

    galaxy = strtrim(sc1120.galaxy,2)
    nicegalaxy = repstr(galaxy,'_',' ')
    xscale = dangular(sc1120.z_obj,/kpc)/206265.0D ; [kpc/arcsec]
 
    htmlfile = 'html/'+galaxy+'.html'
    specnamepng = galaxy+'_spec.png'
    specnameps = galaxy+'_spec.ps'
 
; ---------------------------------------------------------------------------
; define subsamples and generate a WEBSC1120.HTML, WEBSC1120_TABLE.HTML
; and WEBSC1120_DATA.HTML file for each subsample   
; ---------------------------------------------------------------------------

; full sample    
    
    sc1120_populate_frames, sc1120, suffix=suffix

; star-formation sample    
    
    suffix = 'SF'

    sigcut = 5.0
    SFsample = where($
      (sc1120.h_alpha[0]/sc1120.h_alpha[1] gt sigcut)$
;     (sc1120.oii_3727[0]/sc1120.oii_3727[1] gt sigcut)$
      ,ngalaxy_SF)

    if (ngalaxy_SF ne 0L) then sc1120_populate_frames, sc1120[SFsample], suffix=suffix

; -------------------------
; INDEX.HTML
; -------------------------

    pushd, webpath
    
    splog, 'Generating the SC1120 home page.'
    openw, lun, 'index.html', /get_lun
    printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" '+httphome+'>'
    printf, lun, '<html>'
    printf, lun, '<head>'
    printf, lun, '<link rel="stylesheet" type="text/css" href="'+httpstyle+'" />'
    printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
    printf, lun, '<title>SC1120</title>'
    printf, lun, '</head>'
    printf, lun, '<body>'
    printf, lun, '<h1>SC1120</h1>'
    printf, lun, '<div class="widebox">'
    printf, lun, '<p>Spectroscopic Visualizations</p>'
    printf, lun, '<ul>'
    printf, lun, '<li><a href="'+baseref+'websc1120.html">Full Sample ('+$
      string(ngalaxy,format='(I0)')+')</a></li>'
    printf, lun, '<li><a href="'+baseref+'websc1120_SF.html">Emission-Line Sample ('+$
      string(ngalaxy_sf,format='(I0)')+')</a></li>'
    printf, lun, '</ul>'
    printf, lun, '<p>Standard-Star Reductions</p>'
    printf, lun, '<ul>'
    printf, lun, '<li><a href="'+baseref+'redux_Q1.html">Quadrant 1</a></li>'
    printf, lun, '<li><a href="'+baseref+'redux_Q2.html">Quadrant 2</a></li>'
    printf, lun, '<li><a href="'+baseref+'redux_Q3.html">Quadrant 3</a></li>'
    printf, lun, '<li><a href="'+baseref+'redux_Q4.html">Quadrant 4</a></li>'
    printf, lun, '</ul>'
    printf, lun, '<br />'
    printf, lun, '</div>'

    printf, lun, '<p class="email">Last update '+im_today()+'.  Email <a href="mailto:jmoustakas@as.arizona.edu">'+$
      'J. Moustakas</a> with questions or comments.</p>'
    printf, lun, '</body>'
    printf, lun, '</html>'
    free_lun, lun

    popd

; ---------------------------------------------------------------------------
; Generate PNG files of all the spectra and spectral fits.
; ---------------------------------------------------------------------------
 
    splog, 'Generating spectral PNG files.'
    if keyword_set(make_png) then begin
 
       scale = 1E17
       stime0 = systime(1)
 
       for k = 0L, ngalaxy-1L do begin
 
          print, format='("Object ",I4,"/",I4,".",A1,$)', k+1, ngalaxy, string(13b)
 
          ytitle = 'Flux [10^{-17} '+flam_units()+']'
          xtitle = 'Rest Wavelength [\AA]'
 
          specdata = read_sc1120_specfit(galaxy[k],/silent)
          
          wave = reform(specdata[*,0])
          flux = reform(specdata[*,1])
          cflux = reform(specdata[*,2])
          eflux = reform(specdata[*,3]) > 0.0 ; NOTE!
          modelflux = cflux + eflux
 
          smoothflux = scale*smooth(flux,5)
          smoothmodelflux = scale*smooth(modelflux,5)
          smoothcflux = scale*smooth(cflux,5)
          
          stats = im_stats(smoothflux,sigrej=8.0)
          yrange = [stats.minrej,stats.maxrej]
          yrange[0] = yrange[0]<0
          xrange = minmax(wave)
          
; postscript output - white background          
          
          dfpsplot, pngpath+specnameps[k], /landscape, /color
 
          djs_plot, wave, smoothflux, xsty=3, ysty=3, ps=10, xthick=8.0, $
            ythick=8.0, xtitle=xtitle, ytitle=ytitle, charsize=1.5, charthick=8.0, $
            yrange=yrange, thick=2.0, color='grey'
          djs_oplot, wave, smoothmodelflux, ps=10, color='red', thick=3.0
          djs_oplot, wave, smoothcflux, ps=10, color='blue', thick=3.0
          legend, [nicegalaxy[k],'z = '+string(sc1120[k].z_obj,format='(F6.4)')], /left, /top, $
            box=0, charsize=1.5, charthick=8.0
 
          dfpsclose
          spawn, ['gzip -f '+pngpath+specnameps[k]], /sh

; PNG output - black background          
 
          set_plot, 'Z'
          
          polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')
 
          plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=2.0, $
            ythick=2.0, xtitle=xtitle, ytitle=textoidl(ytitle), charsize=1.5, charthick=2.0, $
            yrange=yrange, xrange=xrange, thick=2.0, color=djs_icolor('white')
          djs_oplot, wave, smoothflux, ps=10, thick=0.5, color='grey'
          djs_oplot, wave, smoothmodelflux, ps=10, color='red', thick=1.0
          djs_oplot, wave, smoothcflux, ps=10, color='blue', thick=1.0
          legend, [nicegalaxy[k],'z = '+string(sc1120[k].z_obj,format='(F6.4)')], /left, /top, $
            box=0, charsize=1.5, charthick=1.5, color=djs_icolor('white')
 
          img = tvrd()
          
          tvlct, r, g, b, /get
          write_png, pngpath+specnamepng[k], img, r, g, b
          
          set_plot, 'X'
 
;         splog, 'Writing and compressing '+specnameps[k]+' and '+specnamepng[k]+'.'
;         spawn, ['convert -geometry 512 -rotate 270 '+pngpath+specnameps[k]+' '+pngpath+specnamepng[k]], /sh
 
; ###########################################################################          
          
;          ytitle = 'Flux [10^{-17} '+flam_units()+']'
;          xtitle = 'Rest Wavelength [\AA]'
; 
;          specdata = rd1dspec(sc1120[k].file,datapath=datapath,/silent)
;          
;          wave = specdata.wave / (1+sc1120[k].z)
;          flux = specdata.spec * (1+sc1120[k].z)
; 
;          smoothflux = scale*smooth(flux,3)
;          
;          stats = im_stats(smoothflux,sigrej=8.0)
;          yrange = [stats.minrej,stats.maxrej]
;          yrange[0] = yrange[0]<0
;          xrange = minmax(wave)
;          
;; postscript output - white background          
;          
;          dfpsplot, pngpath+specnameps[k], /landscape, /color
; 
;          djs_plot, wave, smoothflux, xsty=3, ysty=3, ps=10, xthick=8.0, $
;            ythick=8.0, xtitle=xtitle, ytitle=ytitle, charsize=1.5, charthick=8.0, $
;            yrange=yrange, thick=2.0, color='green'
;          legend, [nicegalaxy[k],'z = '+string(sc1120[k].z,format='(F6.4)')], /left, /top, $
;            box=0, charsize=1.5, charthick=8.0
; 
;          dfpsclose
;          spawn, ['gzip -f '+pngpath+specnameps[k]], /sh
;          
;; PNG output - black background          
; 
;          set_plot, 'Z'
;          
;          polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')
; 
;          plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=2.0, $
;            ythick=2.0, xtitle=xtitle, ytitle=textoidl(ytitle), charsize=1.5, charthick=2.0, $
;            yrange=yrange, xrange=xrange, thick=2.0, color=djs_icolor('white')
;          djs_oplot, wave, smoothflux, ps=10, thick=0.5, color='white'
;          legend, [nicegalaxy[k],'z = '+string(sc1120[k].z,format='(F6.4)')], /left, /top, $
;            box=0, charsize=1.5, charthick=1.5, color=djs_icolor('white')
; 
;          img = tvrd()
;          
;          tvlct, r, g, b, /get
;          write_png, pngpath+specnamepng[k], img, r, g, b
;          
;          set_plot, 'X'
;
;          icleanup, specdata
          
       endfor 
 
       splog, format='("Total time = ",G0," minutes.")', (systime(1)-stime0)/60.0
    
    endif
    
; ---------------------------------------------------------------------------
; Generate PNG files of all the image cutouts.
; ---------------------------------------------------------------------------
;
;   splog, 'Generating image PNG files.'
;   if keyword_set(make_thumbs) then begin
;
;      topvalue = 255L
;      minvalue = -10L
;
;      stime0 = systime(1)
;
;      for k = 0L, ngalaxy-1L do begin
;
;         print, format='("Object ",I4,"/",I4,".",A1,$)', k+1, ngalaxy, string(13b)
;
;         info = file_info(thumbpath+Bwthumbnamefits[k])
;         if (info.size lt 1E3) then begin ; no image cutout available
;
;         endif else begin
;
;            image = mrdfits(info.name,0,h,/silent)
;            pltscl = sxpar(h,'PIXSCAL1') ; [arcsec/pixel]
;            
;            imsize = size(image,/dimension)
;            xsize = imsize[0] & xcen = xsize/2.0
;            ysize = imsize[1] & ycen = ysize/2.0
;
;            xaxis = (findgen(xsize)-xcen)*pltscl ; centered on the image [arcsec]
;            yaxis = (findgen(ysize)-ycen)*pltscl ; centered on the image [arcsec]
;
;            scale = 10.0/pltscl ; 10 arcseconds in pixels
;            
; display the image
;            
;            set_plot, 'Z'
;
;;             logimage = alog10(image)
;;             stats = im_stats(logimage,sigrej=5.0)
;;;            stats = im_stats(logimage[xcen-10:xcen+10,ycen-10:ycen+10],sigrej=5.0)
;;
;;;            maximg = median(logimage[xcen-10:xcen+10,ycen-10:ycen+10])
;;             maximg = stats.median+5*stats.sigma_rej
;;             minimg = stats.median-1*stats.sigma_rej
;;             
;;             img = bytscl(logimage,min=minimg,max=maximg)
;
;            img = imgscl(image,/log)
;
;            topvalue = 255L
;            minvalue = -50 ; byte(median(img)-3*djsig(img))
;            
;            img = bytscl(topvalue-img,min=minvalue,top=topvalue)
;
;            window, 0, xs=xsize*3, ys=ysize*3
;            imdisp, img, /erase, position=[0,0,1,1], /normal
;            tvimage, img, /erase, /keep_aspect, xrange=[0,xsize], yrange=[0,ysize], $
;              position=[0,0,1,1], /normal
;            plotimage, img, /noaxes, /preserve_aspect, position=pos, /normal
;            tvcircle, 1.5/pltscl, xsize/2.0, ysize/2.0, /data, color=djs_icolor('dark red'), $
;              thick=2.0
;            legend, 'Bw', /left, /top, box=0, charsize=4.0, charthick=3.5, $
;              textcolor=djs_icolor('dark red'), margin=0, position=[0.22,0.86], /normal
;
;            x0 = fix((convert_coord(pos[0:1],/normal,/to_device))[0])
;            nx = fix((convert_coord(pos[2:3],/normal,/to_device))[0])-X0
;            y0 = fix((convert_coord(pos[0:1],/normal,/to_device))[1])
;            ny = fix((convert_coord(pos[2:3],/normal,/to_device))[1])-y0
;            
;            img = tvrd(x0,y0,nx,ny)
;            
;            tvlct, r, g, b, /get
;            write_png, thumbpath+Bwthumbnamepng[k], img, r, g, b
;
;            set_plot, 'X'
;
;            delvarx, pos
;            
;         endelse
;

;      endfor 
;
;      splog, format='("Total time = ",G0," minutes.")', (systime(1)-stime0)/60.0
;   
;   endif

return
end
