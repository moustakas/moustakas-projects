pro atlas_example_specfit_single, atlas, atlasinfo, postscript=postscript, blackwhite=blackwhite
; jm05aug02uofa - generate a couple examples of our continuum fitting
;                 and stellar absorption correction procedure - one object

    if (n_elements(atlas) eq 0L) then atlas = read_integrated()
    if (n_elements(atlasinfo) eq 0L) then atlasinfo = atlas_read_info()

;   indx = speclinefit_locate(atlas,['NGC2415']) 
    indx = speclinefit_locate(atlas,['NGC3893']) 
;   indx = speclinefit_locate(atlas,['NGC3893','NGC2415']) 

    width1 = [6.0,3.0]
    height1 = 3.0

    xmargin1 = [1.1,0.2]
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
       bw = 1L
    endif else begin
       psname = 'example_specfit.eps'
       bw = 0L
    endelse
    
    pspath = atlas_path(/papers)+'atlas/FIG_ATLAS/'
    if keyword_set(postscript) then begin
       postthick = 5.0
       plotthick1 = 2.0
       plotthick2 = 4.5
       arm_plotconfig, /landscape, nx=2, ny=1, xmargin=xmargin1, $
         ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
         height=height1, coord=pos, xpage=ypage1, ypage=xpage1, $
         psfile=pspath+psname, /writeover, bw=bw
;      dfpsplot, pspath+'example_specfit.eps', /color, xsize=xpage, ysize=ypage, /encapsulated
    endif else begin
       im_window, 0, xratio=0.8, yratio=0.6
       postthick = 1.0
       plotthick1 = 1.0
       plotthick2 = 1.0
       arm_plotconfig, /landscape, nx=2, ny=1, xmargin=xmargin1, $
         ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
         height=height1, coord=pos, xpage=ypage1, ypage=xpage1, $
         /writeover, bw=bw
    endelse

    lcharsize = 1.2             ; 0.7
    pcharsize = 1.4             ; 0.7
    
    width = 120.0
    wave0 = 4861.0
    dwave = 2.75
    interpwave = findgen((2*width)/dwave+1)*dwave+wave0-width
    
; full spectrum plus continuum fit and then the zoom into Balmer line     
    
    atlas_display_spectrum, atlas[indx], lcharsize=lcharsize, labeltype=0L, $
      position=pos[*,0], pcharsize=pcharsize, _extra=extra, labelcolor='', $
      postscript=postscript, postthick=postthick, plotthick=plotthick1, $
      spacing=0.8, speccolor='', plotcolor='', /only_continuum, offsetscale=0.8, $
      ytickname=replicate(' ',10), ytitle='Relative Flux [arbitrary units]'
    xyouts, 3800.0, 380.0, strtrim(atlas[indx].nice_galaxy,2), align=0.0, $
      /data, charsize=lcharsize, charthick=postthick

    atlas_display_spectrum, atlas[indx], lcharsize=1.0, labeltype=0L, $
      position=pos[*,1], ytickname=replicate(' ',10), /noerase, $
      interpwave=interpwave, ytitle='', pcharsize=pcharsize, $
      postscript=postscript, postthick=postthick, plotthick=plotthick2, $
      speccolor='', plotcolor='', /only_continuum, offsetscale=0.9

; give the Balmer-line EW's    

    print, atlas[indx].h_beta_ew
    print, atlas[indx].babs_h_beta_ew

    label = ['EW(H\beta)_{em} = '+strtrim(string(atlas[indx].h_beta_ew[0],$
      format='(F12.1)'),2),'EW(H\beta)_{abs} = '+strtrim(string(atlas[indx].babs_h_beta_ew[0],$
      format='(F12.1)'),2)]+' '+angstrom()
;   legend, textoidl(label), /left, /bottom, box=0, charsize=0.5, $
;     charthick=postthick, /normal, spacing=0
    legend, textoidl('H\beta'), /left, /top, box=0, charsize=lcharsize, $
      charthick=postthick, /normal, spacing=0
    
    if keyword_set(postscript) then dfpsclose

stop    

; code used to pick out a couple examples:    
    
;   bd = atlas.h_beta_ew[0]-atlas.babs_h_beta_ew[0]
;   indx = where((bd gt 4.5) and (bd lt 5.2))
;
;   nindx = n_elements(indx)
;   
;   infomatch = lonarr(nindx)
;   for iindx = 0L, nindx-1L do infomatch[iindx] = where(atlasinfo.atlas_id eq atlas[indx[iindx]].atlas_id)
;
;   niceprint, atlas[indx].galaxy, atlasinfo[infomatch].galaxy, bd[indx]
;   atlas_html, atlasinfo[infomatch]
;
;   atlas_display_spectrum, atlasinfo[infomatch], interpwave=findgen(500)+4611.0

return
end
    
