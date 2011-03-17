;+
; NAME:
;       SINGS_SPECTRUM_MOSAIC
;
; PURPOSE:
;       Mosaic all the integrated spectra nuclear spectra with the DSS
;       visualizations.  
;
; CALLING SEQUENCE:
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
; PROCEDURES USED:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;
;-

pro sings_spectrum_mosaic, sings, postscript=postscript, make_png=make_png, $
  nperpage=nperpage, single=single
; jm03oct16uofa
; jm05jul29uofa - updated    
; generate the multi-panel figures for the sings data paper; fit 20
; spectra per page, resulting in 18 pages

; jm03dec03uofa - rebin the spectra to lower resolution to reduce the
;                 image size
; jm05oct27uofa - generate a single postscript page per object 
    
; initialize path names

    paperpath = sings_path(/papers)
    
    if keyword_set(single) and keyword_set(make_png) then message, 'Not supported yet.'
    
    if keyword_set(make_png) then begin
       pspath = sings_path(/web)
       postscript = 0L
       labelcolor = 'white'
       speccolor = 'white'
       no_specfit = 1L
       postthick = 1.5
       plotthick = 1.0
       lcharsize = 0.5
    endif else begin
       if keyword_set(single) then pspath = paperpath+'sings/SUBMIT/' else $
         pspath = paperpath+'sings/FIG_SINGS/'
       speccolor = ''
       labelcolor = ''
       no_specfit = 1L
       if keyword_set(single) then begin
          lcharsize = 1.0
          pcharsize = 1.2
       endif else begin
          lcharsize = 0.6
          pcharsize = 1.5
       endelse
       if keyword_set(postscript) then begin
          postthick = 3.0
          plotthick = 3.0
       endif else begin
          im_window, 0, xratio=0.8, yratio=0.7
          postthick = 1.0
          plotthick = 1.0
       endelse
    endelse

; restore all the fitting results

    if (n_elements(sings) eq 0L) then sings = sings_read_info()
    
    nicegalaxy = strtrim(sings.nice_galaxy,2)
    id = sings.sings_id
    galaxy = strcompress(sings.galaxy,/remove)
    ngalaxy = n_elements(galaxy)

    if keyword_set(single) then begin

       nperpage = 1L
       npage = ngalaxy

       ncols = 2L
       nrows = 1L

       width = [3.0,6.0]
       height = 3.0
       
       xspace = 0.0
       yspace = 0.0

       xmargin = [0.2,0.2]
       ymargin = [0.2,0.6]

       psname = 'f8_'+string(lindgen(npage)+1L,format='(I0)')+'.eps'

    endif else begin

; plotting variables

       if (n_elements(nperpage) eq 0L) then nperpage = 10.0 ; number of objects per page
       npage = ceil(ngalaxy/float(nperpage))
       
       ncols = 4.0              ; 2 objects, each with an image and a spectrum
       nrows = nperpage/(ncols/2.0)

       width = [1.5,3.0,1.5,3.0]
       height = replicate(1.5,nrows)

       xspace = [0.0,0.3,0.0]
       yspace = 0.0
       xmargin = [0.2,0.2]
       ymargin = [0.2,0.4]

;       pagemaker, nx=ncols, ny=nrows, xspace=0.0, yspace=0.0, $
;         xmargin=xmargin, ymargin=ymargin, width=width, /normal, $
;         height=height, position=pos, xpage=xpage, ypage=ypage, $
;         /landscape

       psname = 'sings_figure_'+string(lindgen(npage),format='(I2.2)')+'.eps'

    endelse
       
    xpage = total(width)+total(xmargin)+total(xspace)
    ypage = total(height)+total(ymargin)+total(yspace)

    pngname = repstr(psname,'.eps','.png')

    for j = 0L, npage-1L do begin

       xspace1 = xspace & yspace1 = yspace
       xmargin1 = xmargin & ymargin1 = ymargin
       width1 = width & height1 = height
       xpage1 = xpage & ypage1 = ypage

       if keyword_set(postscript) then begin
          splog, 'Writing '+pspath+psname[j]+'.'
          arm_plotconfig, /landscape, nx=ncols, ny=nrows, xmargin=xmargin1, $
            ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
            height=height1, coord=pos, xpage=ypage1, ypage=xpage1, $
            psfile=pspath+psname[j], /writeover, bw=1L
;         dfpsplot, pspath+psname[j], /landscape, xsize=xpage, ysize=ypage, /color
       endif else begin
          if keyword_set(make_png) then set_plot, 'Z'
          arm_plotconfig, /landscape, nx=ncols, ny=nrows, xmargin=xmargin1, $
            ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
            height=height1, coord=pos, xpage=ypage1, ypage=xpage1, bw=0
       endelse

; define a common wavelength vector for this subset of objects 

       i1 = j*nperpage
       i2 = (i1+(nperpage-1L))<(ngalaxy-1L)
       
       minwave = max(sings[i2:i2].drift56_minwave/(1+sings[i1:i2].z))
       maxwave = min(sings[i2:i2].drift56_maxwave/(1+sings[i1:i2].z))
       dwave = 2.75
       interpwave = findgen((maxwave-minwave)/dwave+1)*dwave+minwave

; visualize each object       
       
       for k = 0L, nperpage-1L do begin
          
          indx = j*nperpage+k

          if (k mod (ncols/2.0)) eq 0L then begin
             ytickname = replicate(' ',10)
             ytitle = ''
          endif else begin
             ytickname = replicate(' ',10)
;            delvarx, ytickname
             ytitle = ''
;            ytitle = 'Normalized Flux'
          endelse

          if (j eq (npage-1L)) then cutit = (ngalaxy mod nperpage)-(ncols/2.0) else cutit = nperpage-(ncols/2.0)
          if (k lt cutit) then begin
             xtickname = replicate(' ',10)
             xtitle = ''
          endif else begin
             delvarx, xtickname
             xtitle = 'Rest Wavelength [\AA]'
          endelse
          
          if (indx[0] lt ngalaxy) then begin

             sings_display_image, sings[indx], imagepath=dsspath, imposition=pos[*,2*k], $
               lcharsize=lcharsize, pcharsize=pcharsize, /preserve_aspect, /nolabelbar, $;/nobar, $
               pspath=pspath, labeltype=3L, _extra=extra, spacing=0, $ ;postscript=postscript, $
               noerase=(k ne 0L), xtickname=xtickname, ytickname=ytickname, $
               postthick=postthick

; generate a 3-panel visualization of the three spectra

             match, strlowcase(galaxy[k]), strlowcase(strtrim(sings.galaxy,2)), indx1, indx2
             if (indx1[0] eq -1L) then $
               specinfo = {drift56: 0L, drift20: 0L, nuclear: 0L} else $
                 specinfo = sings[indx2]

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
                  ythick=2.0, xtitle=xtitle, ytitle='', charsize=charsize, charthick=2.0, /noerase, $
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
                     ythick=2.0, xtitle=xtitle, ytitle='', charsize=charsize, charthick=2.0, /noerase, $
                     position=pos[*,2], ytickname=replicate(' ',10), yticklen=1E-10, xticklen=1E-10
                   
                endelse

                xyouts, (pos[2,2]-pos[0,2])/2.0+pos[0,2], (pos[3,2]-pos[1,2])/2.0+pos[1,2], $
                  'Spectrum Unavailable', align=0.5, charthick=postthick, /normal

             endelse
             
; ytitle

             xyouts, 0.3*pos[0,0], (pos[3,0]-pos[1,2])/2.0+pos[1,2], textoidl(ytitle), $
               align=0.5, orientation=90, charsize=charsize, charthick=postthick, /normal

             
             
;            sings_display_spectrum, sings[indx], lcharsize=1.0, labeltype=0L, $
;              position=pos[*,2*k+1], xtickname=xtickname, ytickname=ytickname, $
;              _extra=extra, /noerase, xtitle=xtitle, ytitle=ytitle, $
;              no_specfit=no_specfit, interpwave=interpwave, speccolor=speccolor, $
;              labelcolor=labelcolor, postthick=postthick, plotthick=plotthick, $
;              plotcolor='', pcharsize=pcharsize

          endif

       endfor

       if keyword_set(postscript) then dfpsclose else begin
          if keyword_set(make_png) then begin
             img = tvrd()
             tvlct, r, g, b, /get
             write_png, pspath+pngname[j], img, r, g, b
             set_plot, 'X'
          endif else cc = get_kbrd(1)
       endelse

    endfor

return
end    
