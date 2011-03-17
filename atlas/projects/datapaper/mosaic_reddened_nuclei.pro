pro mosaic_reddened_nuclei, atlas1, postscript=postscript, make_png=make_png
; jm05aug23uofa - generate a mosaic of all the highly reddened nuclei
;                 in the spectral atlas
    
; initialize path names

    if keyword_set(make_png) then begin
       pspath = atlas_path(/web)
       postscript = 0L
       labelcolor = 'white'
       speccolor = 'white'
       no_specfit = 1L
       postthick = 1.5
       plotthick = 1.0
       lcharsize = 0.5
    endif else begin
       pspath = '/home/ioannis/temp/'
;      pspath = atlas_path(/papers)+'atlas/FIG_ATLAS/'
       speccolor = ''
       labelcolor = ''
       no_specfit = 1L
       lcharsize = 0.6
       if keyword_set(postscript) then begin
          postthick = 5.0
          plotthick = 5.0
       endif else begin
          im_window, 0, xratio=0.8, yratio=0.7
          postthick = 2.0
          plotthick = 1.0
       endelse
    endelse

; restore all the fitting results

    if (n_elements(atlas1) eq 0L) then atlas1 = read_nuclear()

    srt = reverse(sort(atlas1.continuum_ebv[0]))
    atlas = atlas1[srt]

    keep = where(atlas.continuum_ebv[0] gt 0.5)
    atlas = atlas[keep]
    
    nicegalaxy = strtrim(atlas.nice_galaxy,2)
    id = atlas.atlas_id
    galaxy = strcompress(atlas.galaxy,/remove)
    ngalaxy = n_elements(galaxy)

; plotting variables

    nperpage = 10.0              ; number of objects per page
    npage = ceil(ngalaxy/float(nperpage))
    
    xmargin = [0.2,0.2] & ymargin = [0.2,0.4]

    ncols = 4.0 ; 2 objects, each with an image and a spectrum
    nrows = nperpage/(ncols/2.0)

    width = [1.5,3.0,1.5,3.0]
    height = replicate(1.5,nrows)

    xspace = [0.0,0.3,0.0] & yspace = 0.0
    
    xpage = total(width)+total(xmargin)+total(xspace)
    ypage = total(height)+total(ymargin)+total(yspace)

;    pagemaker, nx=ncols, ny=nrows, xspace=0.0, yspace=0.0, $
;      xmargin=xmargin, ymargin=ymargin, width=width, /normal, $
;      height=height, position=pos, xpage=xpage, ypage=ypage, $
;      /landscape

    psname = 'reddened_nuclei_'+string(lindgen(npage),format='(I2.2)')+'.ps'
    pngname = repstr(psname,'.ps','.png')

;   for j = 41, npage-1L do begin
;   for j = 4L, 4L do begin
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
            psfile=pspath+psname[j], /writeover
;         dfpsplot, pspath+psname[j], /landscape, xsize=xpage, ysize=ypage, /color
       endif else begin
          if keyword_set(make_png) then set_plot, 'Z'
          arm_plotconfig, /landscape, nx=ncols, ny=nrows, xmargin=xmargin1, $
            ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
            height=height1, coord=pos, xpage=ypage1, ypage=xpage1
       endelse

; define a common wavelength vector for this subset of objects 

       i1 = j*nperpage
       i2 = (i1+(nperpage-1L))<(ngalaxy-1L)
       
       minwave = max(atlas[i2:i2].drift_minwave/(1+atlas[i1:i2].z))
       maxwave = min(atlas[i2:i2].drift_maxwave/(1+atlas[i1:i2].z))
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
             xtitle = 'Rest Wavelength ['+angstrom()+']'
          endelse
          
          if (indx[0] lt ngalaxy) then begin
             atlas_display_image, atlas[indx], imagepath=dsspath, imposition=pos[*,2*k], $
               lcharsize=lcharsize, pcharsize=pcharsize, /preserve_aspect, /nobar, /nolabelbar, $
               pspath=pspath, labeltype=1L, _extra=extra, spacing=0, $ ;postscript=postscript, $
               noerase=(k ne 0L), xtickname=xtickname, ytickname=ytickname, $
               postthick=postthick, /norc3box, /noscanbox
             atlas_display_spectrum, atlas[indx], lcharsize=1.0, labeltype=0L, $
               position=pos[*,2*k+1], xtickname=xtickname, ytickname=ytickname, $
               _extra=extra, /noerase, xtitle=xtitle, ytitle=ytitle, $
               no_specfit=no_specfit, interpwave=interpwave, speccolor=speccolor, $
               labelcolor=labelcolor, postthick=postthick, plotthick=plotthick, $
               /nuclear

             ebv = strtrim(string(atlas[indx].continuum_ebv[0],format='(F12.2)'),2)
             age = strtrim(string(atlas[indx].continuum_age_b/1E3,format='(F12.2)'),2)
             label = ['E(B-V) = '+ebv+' mag','Age = '+age+' Gyr']
             
             legend, textoidl(label), /left, /top, box=0, charsize=1.2, $
               charthick=postthick, /normal, _extra=extra, spacing=0
             
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
