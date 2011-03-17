pro atlas_example_specfit, atlas, atlasinfo, postscript=postscript, blackwhite=blackwhite
; jm05aug02uofa - generate a couple examples of our continuum fitting
;                 and stellar absorption correction procedure
; jm05oct21uofa - make it a multi-page figure

    if (n_elements(atlas) eq 0L) then atlas = read_integrated()
    if (n_elements(atlasinfo) eq 0L) then atlasinfo = atlas_read_info()

;   indx = speclinefit_locate(atlas,['NGC2415']) 
    indx = speclinefit_locate(atlas,['NGC0034','NGC4651','NGC3893','NGC0337'])
;   indx = speclinefit_locate(atlas,['NGC3893','NGC2415'])
    ngalaxy = n_elements(indx)

    nx = 2L
    ny = 4L
    
    width1 = [4.5,1.85]
    height1 = [1.85,1.85,1.85,1.85]

;   xmargin1 = [0.0,0.0]
;   ymargin1 = [0.0,0.0]
    xmargin1 = [0.8,0.3]
    ymargin1 = [0.3,0.8]

    xspace1 = 0.1
    yspace1 = 0.0
    
    xpage1 = total(width1)+total(xmargin1)+total(xspace1)
    ypage1 = total(height1)+total(ymargin1)+total(yspace1)

;   pagemaker, nx=2, ny=1, xspace=xspace, yspace=yspace, xmargin=xmargin, $
;     ymargin=ymargin, width=width, height=height, position=pos, $
;     xpage=xpage, ypage=ypage, /normal

    if keyword_set(blackwhite) then begin
       psname = 'example_specfit_blackwhite.eps'
    endif else begin
       psname = 'example_specfit.eps'
    endelse
    
    pspath = atlas_path(/papers)+'atlas/FIG_ATLAS/'
    if keyword_set(postscript) then begin
       postthick = 5.0
       plotthick1 = 2.0
       plotthick2 = 4.5
       arm_plotconfig, landscape=0, nx=nx, ny=ny, xmargin=xmargin1, $
         ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
         height=height1, coord=pos, xpage=xpage1, ypage=ypage1, $
         psfile=pspath+psname, /writeover, bw=blackwhite
;      arm_plotconfig, /landscape, nx=nx, ny=ny, xmargin=xmargin1, $
;        ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
;        height=height1, coord=pos, xpage=ypage1, ypage=xpage1, $
;        psfile=pspath+psname, /writeover, color=color
;      dfpsplot, pspath+psname, /color, xsize=xpage1, ysize=ypage1, encapsulated=1
;      pagemaker, nx=nx, ny=ny, xmargin=xmargin1, ymargin=ymargin1, xspace=xspace1, $
;        yspace=yspace1, width=width1, height=height1, position=pos, xpage=ypage1, $
;        ypage=xpage1, /normal, /landscape
    endif else begin
       im_window, 0, xratio=0.7, yratio=0.7
       postthick = 2.0
       plotthick1 = 1.0
       plotthick2 = 1.0
       arm_plotconfig, /landscape, nx=nx, ny=ny, xmargin=xmargin1, $
         ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
         height=height1, coord=pos, xpage=ypage1, ypage=xpage1, $
         /writeover, bw=blackwhite
    endelse

    lcharsize = 1.2             ; 0.7
    pcharsize = 2.0
    
    width = 120.0
    wave0 = 4861.0
    dwave = 2.75
    interpwave = findgen((2*width)/dwave+1)*dwave+wave0-width

    ytitle = ''
;   ytitle = 'Relative Flux'
    ytitle2 = 'Relative Flux [arbitrary units]'
    
    for i = 0L, ngalaxy-1L do begin

       if (i lt ngalaxy-1L) then begin
          xtickname = replicate(' ',10)
          xtitle = ''
       endif else begin
          delvarx, xtickname
          delvarx, xtitle
       endelse

       if (i eq 0L) then factorymax = 1.2 else factorymax = 1.15
       
; full spectrum plus continuum fit and then the zoom into Balmer line     
    
       atlas_display_spectrum, atlas[indx[i]], lcharsize=lcharsize, labeltype=0L, $
         position=pos[*,2*i], pcharsize=pcharsize, _extra=extra, labelcolor='', $
         postscript=postscript, postthick=postthick, plotthick=plotthick1, $
         spacing=0.8, speccolor='', plotcolor='', /only_continuum, offsetscale=0.75, $
         ytickname=replicate(' ',10), ytitle=ytitle, noerase=(i gt 0L), $
         xtickname=xtickname, xtitle=xtitle, setyrange=2L, factorymax=1.3
       xyouts, 3850.0, !y.crange[1]*0.85, strtrim(atlas[indx[i]].nice_galaxy,2), align=0.0, $
         /data, charsize=lcharsize, charthick=postthick
       
       atlas_display_spectrum, atlas[indx[i]], lcharsize=1.0, labeltype=0L, $
         position=pos[*,2*i+1], ytickname=replicate(' ',10), /noerase, $
         interpwave=interpwave, ytitle='', pcharsize=pcharsize, $
         postscript=postscript, postthick=postthick, plotthick=plotthick2, $
         speccolor='', plotcolor='', /only_continuum, offsetscale=0.9, $
         xtickname=xtickname, xtitle=xtitle, xtickinterval=100.0, setyrange=2L, $
         factorymax=factorymax

;      if (i eq 0L) then xyouts, pos[0,0]*0.8, (pos[3,0]-pos[1,2*ngalaxy-1])/2.0+pos[1,2*ngalaxy-1], $
;        ytitle2, align=0.5, charsize=pcharsize*0.6, charthick=postthick, /normal, $
;        orientation=90
       if (i eq 0L) then begin
          xyouts, pos[0,0]*0.8, (pos[3,0]-pos[1,2])/2.0+pos[1,2], $
            ytitle2, align=0.5, charsize=pcharsize*0.5, charthick=postthick, /normal, $
            orientation=90
          xyouts, pos[0,0]*0.8, (pos[3,4]-pos[1,6])/2.0+pos[1,6], $
            ytitle2, align=0.5, charsize=pcharsize*0.5, charthick=postthick, /normal, $
            orientation=90
       endif

; give the Balmer-line EW's    

       print, atlas[indx[i]].h_beta_ew[0]-atlas[indx[i]].babs_h_beta_ew[0]
       print, atlas[indx[i]].h_beta_ew
       print, atlas[indx[i]].babs_h_beta_ew
       print

       label = ['EW(H\beta)_{em} = '+strtrim(string(atlas[indx[i]].h_beta_ew[0],$
         format='(F12.1)'),2),'EW(H\beta)_{abs} = '+strtrim(string(atlas[indx[i]].babs_h_beta_ew[0],$
         format='(F12.1)'),2)]+' '+angstrom()
;      legend, textoidl(label), /left, /bottom, box=0, charsize=0.5, $
;        charthick=postthick, /normal, spacing=0
       if (i eq 0L) then begin
          legend, textoidl('H\beta'), /left, /top, box=0, charsize=lcharsize, $
            charthick=postthick, /normal, spacing=0
       endif
    
    endfor

    if keyword_set(postscript) then dfpsclose
    cleanplot, /silent

stop    

; code used to pick out a couple examples:    
    
    bd = atlas.h_beta_ew[0]-atlas.babs_h_beta_ew[0]
    indx = where((bd gt 4.5) and (bd lt 5.2),nindx)
    infomatch = lonarr(nindx)
    for iindx = 0L, nindx-1L do infomatch[iindx] = where(atlasinfo.atlas_id eq atlas[indx[iindx]].atlas_id)
    niceprint, atlas[indx].galaxy, atlasinfo[infomatch].galaxy, bd[indx]
    atlas_html, atlasinfo[infomatch]
    atlas_display_spectrum, atlasinfo[infomatch], interpwave=findgen(500)+4611.0

    indx = where((atlas.h_beta_ew[0] gt 5.0) and (atlas.h_beta_ew[0] lt 10.0),nindx)
    infomatch = lonarr(nindx)
    for iindx = 0L, nindx-1L do infomatch[iindx] = where(atlasinfo.atlas_id eq atlas[indx[iindx]].atlas_id)
    niceprint, atlas[indx].galaxy, atlasinfo[infomatch].galaxy, atlas[indx].h_beta_ew[0]
    
    atlas_display_spectrum, atlasinfo[infomatch];, interpwave=findgen(500)+4611.0
    
return
end
    
