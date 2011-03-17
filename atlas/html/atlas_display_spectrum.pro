;+
; NAME:
;       ATLAS_DISPLAY_SPECTRUM
;
; PURPOSE:
;       Display an atlas spectrum and the best-fitting model.  
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
;       J. Moustakas, 2005 May 16, U of A
;       jm05jul25uofa - added WITH_NUCLEAR keyword; added NSMOOTH
;                       optional input
;-

pro atlas_display_spectrum, atlas, allatlas=allatlas, postscript=postscript, pcharsize=pcharsize, $
  lcharsize=lcharsize, labeltype=labeltype, nsmooth=nsmooth, $
  with_nuclear=with_nuclear, position=position, _extra=extra, $
  no_specfit=no_specfit, interpwave=interpwave, speccolor=speccolor, $
  labelcolor=labelcolor, plotcolor=plotcolor, postthick=postthick, plotthick1=plotthick1 , $
  plotthick2=plotthick2, nuclear=nuclear, offsetscale=offsetscale, only_continuum=only_continuum, $
  setyrange=setyrange, factorymax=factorymax, specfitseparate=specfitseparate

    atlas1dpath = atlas_path(/atlas1d)

    natlas = n_elements(atlas)
    if (natlas eq 0L) then begin
       print, 'Syntax - atlas_display_spectrum, atlas'
       return
    endif

    if (n_elements(pcharsize) eq 0L) then pcharsize = 1.5
    if (n_elements(lcharsize) eq 0L) then lcharsize = 1.2
    if (n_elements(labeltype) eq 0L) then labeltype = 1L
    if (n_elements(nsmooth) eq 0L) then nsmooth = 1L
    if (n_elements(speccolor) eq 0L) then speccolor = 'grey'
    if (n_elements(labelcolor) eq 0L) then labelcolor = 'white'
    if (n_elements(plotcolor) eq 0L) then plotcolor = 'white'
    if (n_elements(offsetscale) eq 0L) then offsetscale = 1.0
    if (n_elements(setyrange) eq 0L) then setyrange = 1L
    if (n_elements(factorymax) eq 0L) then factorymax = 1.2

; call this routine recursively

    if (natlas gt 1L) then begin
       for i = 0L, natlas-1L do begin
          atlas_display_spectrum, atlas[i], allatlas=allatlas, postscript=postscript, pcharsize=pcharsize, $
            lcharsize=lcharsize, labeltype=labeltype, nsmooth=nsmooth, $
            with_nuclear=with_nuclear, position=position, _extra=extra, $
            no_specfit=no_specfit, interpwave=interpwave, speccolor=speccolor, $
            labelcolor=labelcolor, plotcolor=plotcolor, postthick=postthick, plotthick1=plotthick1, $
            plotthick2=plotthick2, nuclear=nuclear, offsetscale=offsetscale, only_continuum=only_continuum, $
            setyrange=setyrange, factorymax=factorymax, specfitseparate=specfitseparate
          if (not keyword_set(postscript)) then cc = get_kbrd(1)
       endfor
       return
    endif

    galaxy = strtrim(atlas.galaxy,2)
    nicegalaxy = strtrim(atlas.nice_galaxy,2)
    driftfile = strtrim(atlas.drift_file,2)

; two thicknesses are required depending on whether or not a polyfill
; black background is placed before displaying this image (compare,
; for example, ATLAS_HTML, and INTEGRATED_ABUNDANCES)

    if (n_elements(postthick) eq 0L) then if keyword_set(postscript) then $
      postthick = 5.0 else postthick = 2.0
    if (n_elements(plotthick1) eq 0L) then if keyword_set(postscript) then $
      plotthick1 = 5.0 else plotthick1 = 1.0
    if (n_elements(plotthick2) eq 0L) then if keyword_set(postscript) then $
      plotthick2 = 5.0 else plotthick2 = 1.0
    
; display the spectra

    if (n_elements(extra) ne 0) then begin
       if (tag_exist(extra,'xtickname')) then xtickname = extra.xtickname
       if (tag_exist(extra,'ytickname')) then ytickname = extra.ytickname
    endif
    
    scale = 1E15
    ytitle = 'Flux [10^{-15} '+flam_units()+']'
    xtitle = 'Rest Wavelength [\AA]'

    if keyword_set(with_nuclear) and (n_elements(position) eq 0L) then begin
       pagemaker, nx=1, ny=2, xspace=0.0, yspace=0.0, xmargin=[1.2,0.2], $
         ymargin=[0.3,1.4], position=pos, /normal
       xtickname = replicate(' ',10)
       xtitle2 = ''
       ytitle2 = ''
    endif else begin
       if (n_elements(position) eq 0L) then $
         pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, xmargin=[1.2,0.2], $
           ymargin=[0.3,1.4], position=pos, /normal else $
         pos = position
       xtitle2 = xtitle
       ytitle2 = ytitle
    endelse
    
    specdata = read_atlas_specfit(galaxy,specdata=allatlas,nuclear=nuclear,/silent)
    if (specdata[0] eq -1L) then begin
       specdata = rd1dspec(driftfile,datapath=atlas1dpath,/silent)
       wave = specdata.wave
       flux = specdata.spec
       cflux = flux*0.0-1.0
       eflux = flux*0.0-1.0
    endif else begin
       wave = reform(specdata[*,0])
       flux = reform(specdata[*,1])
       cflux = reform(specdata[*,2])
       eflux = reform(specdata[*,3]) > 0.0 ; NOTE!
    endelse

    modelflux = cflux + eflux

    if (n_elements(interpwave) ne 0L) then begin
       linterp, wave, flux, interpwave, flux, missing='NaN'
       linterp, wave, cflux, interpwave, cflux, missing='NaN'
       linterp, wave, eflux, interpwave, eflux, missing='NaN'
       linterp, wave, modelflux, interpwave, modelflux, missing='NaN'
;      flux = interpol(flux,wave,interpwave)
;      cflux = interpol(cflux,wave,interpwave)
;      eflux = interpol(eflux,wave,interpwave)
;      modelflux = interpol(modelflux,wave,interpwave)
       wave = interpwave
    endif
    
    smoothflux = scale*smooth(flux,nsmooth)
    smoothmodelflux = scale*smooth(modelflux,nsmooth)
    smoothcflux = scale*smooth(cflux,nsmooth)
    xrange = minmax(wave)

    case setyrange of
       1L: begin
          stats = im_stats(smoothflux,sigrej=8.0)
          yrange = [stats.minrej,stats.maxrej]
          yrange[0] = (yrange[0]<(min(smoothcflux)*offsetscale))>0
       end
       2L: begin
          stats = im_stats(scale*cflux,sigrej=8.0)
          yrange = [stats.minrej,factorymax*stats.maxrej]
          yrange[0] = (yrange[0]<(min(scale*cflux)*offsetscale))>0
       end
       3L: begin
          yrange = minmax(scale*modelflux)
       end
    endcase

    plot, [0], [0], /nodata, xthick=postthick, xsty=3, ysty=3, $
      ythick=postthick, xtitle=textoidl(xtitle2), ytitle=textoidl(ytitle2), yminor=4, $
      charsize=pcharsize, charthick=postthick, yrange=yrange, xrange=xrange, $
      thick=postthick, color=djs_icolor(plotcolor), _extra=extra, $
      position=pos[*,0], xtickname=xtickname, ytickname=ytickname

    if keyword_set(no_specfit) then begin
       djs_oplot, wave, smoothflux, ps=10, color=speccolor, thick=plotthick1
    endif else begin
       djs_oplot, wave, smoothflux, ps=10, color=speccolor, thick=plotthick1
       if keyword_set(only_continuum) then begin
          djs_oplot, wave, smoothcflux*offsetscale, ps=10, color='red', thick=plotthick2 
       endif else begin
          if keyword_set(specfitseparate) then begin
             djs_oplot, wave, smoothmodelflux*offsetscale, ps=10, color='green', thick=plotthick2
             djs_oplot, wave, smoothcflux*offsetscale, ps=10, color='red', thick=plotthick2
          endif else begin
             djs_oplot, wave, smoothmodelflux*offsetscale, ps=10, color='red', thick=plotthick2
          endelse
       endelse
    endelse

    case labeltype of
       0L: label = ''
       1L: label = nicegalaxy
       2L: label = [nicegalaxy,strtrim(atlas.lit_type,2)]
       3L: label = nicegalaxy+' - Integrated'
       else: label = nicegalaxy
    endcase
    
    legend, textoidl(label), /left, /top, box=0, charsize=lcharsize, $
      charthick=postthick, textcolor=djs_icolor(labelcolor), $
      /normal, _extra=extra
    
    if keyword_set(with_nuclear) then begin

       specdata = read_atlas_specfit(galaxy,/silent,/nuclear)

       if (specdata[0] eq -1L) then begin

          plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=postthick, $
            ythick=postthick, xtitle=textoidl(xtitle), ytitle='', charsize=pcharsize, $
            charthick=postthick, /noerase, position=pos[*,1], xrange=xrange, $
            ytickname=replicate(' ',10), yticklen=1E-10, color=djs_icolor(plotcolor)
          xyouts, (pos[2,1]-pos[0,1])/2.0+pos[0,1], (pos[3,1]-pos[1,1])/2.0+pos[1,1], $
            'Nuclear Spectrum Unavailable', align=0.5, charthick=postthick, /normal

       endif else begin

          wave = reform(specdata[*,0])
          flux = reform(specdata[*,1])
          cflux = reform(specdata[*,2])
          eflux = reform(specdata[*,3]) > 0.0 ; NOTE!
          modelflux = cflux + eflux

          if (n_elements(interpwave) ne 0L) then begin
             linterp, wave, flux, interpwave, flux, missing='NaN'
             linterp, wave, cflux, interpwave, cflux, missing='NaN'
             linterp, wave, eflux, interpwave, eflux, missing='NaN'
             linterp, wave, modelflux, interpwave, modelflux, missing='NaN'
;            flux = interpol(flux,wave,interpwave)
;            cflux = interpol(cflux,wave,interpwave)
;            eflux = interpol(eflux,wave,interpwave)
;            modelflux = interpol(modelflux,wave,interpwave)
             wave = interpwave
          endif
          
          smoothflux = scale*smooth(flux,nsmooth)
          smoothmodelflux = scale*smooth(modelflux,nsmooth)
          smoothcflux = scale*smooth(cflux,nsmooth)
          
          stats = im_stats(smoothflux,sigrej=8.0)
          yrange = [stats.minrej,stats.maxrej]
          yrange[0] = yrange[0]>0

          plot, [0], [0], /nodata, xthick=postthick, xsty=3, ysty=3, $
            ythick=postthick, xtitle=textoidl(xtitle), ytitle='', /noerase, yminor=4, $
            charsize=pcharsize, charthick=postthick, yrange=yrange, xrange=xrange, $
            thick=postthick, color=djs_icolor(plotcolor), _extra=extra, $
            position=pos[*,1]
          djs_oplot, wave, smoothflux, ps=10, color=speccolor, thick=plotthick1 ; thick=(postthick-2L)>1L
          djs_oplot, wave, smoothmodelflux, ps=10, color='red', thick=plotthick2 ; thick=postthick
;         djs_oplot, wave, smoothcflux, ps=10, color='yellow', thick=postthick

          case labeltype of
             0L: label = ''
             1L: label = 'Nuclear'
             2L: label = 'Nuclear'
             3L: label = 'Nuclear'
             else: label = nicegalaxy
          endcase
          
          legend, textoidl(label), /left, /top, box=0, charsize=lcharsize, $
            charthick=postthick, textcolor=djs_icolor(labelcolor), $
            /normal, _extra=extra
          
       endelse

       xyouts, 0.3*pos[0,0], (pos[3,0]-pos[1,1])/2.0+pos[1,1], textoidl(ytitle), $
         align=0.5, orientation=90, charsize=pcharsize, charthick=postthick, /normal
       
    endif
    
return
end
    
