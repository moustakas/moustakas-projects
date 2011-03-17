pro ages_spectrum_mosaic, ages, agesancillary, postscript=postscript, make_png=make_png, $
  nperpage=nperpage, single=single, prefix=prefix, pspath=pspath, mypng=mypng, _extra=extra
; jm05nov18uofa - written based on ATLAS_SPECTRUM_MOSAIC
    
; initialize path names

    if keyword_set(single) and keyword_set(make_png) then message, 'Not supported yet.'

    if (n_elements(prefix) eq 0L) then prefix = 'ages_figure'
    
    if keyword_set(make_png) then begin
       if (n_elements(pspath) eq 0L) then pspath = cwd() ; ages_path(/papers)+'MOSAIC/'
;      pspath = ages_path(/web)
       postscript = 0L
       labelcolor = 'white'
       speccolor = 'white'
       no_specfit = 1L
       postthick = 1.5
       plotthick = 1.0
       lcharsize = 0.5
    endif else begin
       if (n_elements(pspath) eq 0L) then pspath = cwd() ; ages_path(/papers)+'MOSAIC/'
       if keyword_set(mypng) then begin
          plotcolor = 'white' 
          speccolor = 'white'
          labelcolor = 'white'
       endif else begin
          plotcolor = 'black'
          speccolor = 'black'
          labelcolor = 'black'
       endelse
       no_specfit = 0L ; 1L
       if keyword_set(single) then begin
          lcharsize = 1.0
          pcharsize = 1.2
       endif else begin
          lcharsize = 0.9
          pcharsize = 1.1
       endelse
       if keyword_set(postscript) then begin
          if keyword_set(mypng) then begin
             postthick = 5.0
             plotthick = 4.0
          endif else begin
             postthick = 4.0
             plotthick = 2.0
          endelse
       endif else begin
          im_window, 0, xratio=0.8, yratio=0.7
          postthick = 1.0
          plotthick = 1.0
       endelse
    endelse

; restore all the fitting results

    if (n_elements(ages) eq 0L) then ages = read_ages_mz_sample(agesancillary=agesancillary)

    galaxy = strcompress(ages.galaxy,/remove)
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

       psname = 'f8_'+string(lindgen(npage)+1L,format='(I0)')+'.ps'

    endif else begin

; plotting variables

       if (n_elements(nperpage) eq 0L) then nperpage = 15.0 ; number of objects per page
       npage = ceil(ngalaxy/float(nperpage))
       
       ncols = 3.0
       nrows = nperpage/(ncols)

       width = [3.3,3.3,3.3]
       height = replicate(1.5,nrows)

       xspace = [0.15,0.15]
       yspace = 0.0
       xmargin = [0.6,0.3]
       ymargin = [0.3,0.6]

;       pagemaker, nx=ncols, ny=nrows, xspace=0.0, yspace=0.0, $
;         xmargin=xmargin, ymargin=ymargin, width=width, /normal, $
;         height=height, position=pos, xpage=xpage, ypage=ypage, $
;         /landscape

       psname = prefix+'_'+string(lindgen(npage),format='(I3.3)')+'.ps'

    endelse
       
    xpage = total(width)+total(xmargin)+total(xspace)
    ypage = total(height)+total(ymargin)+total(yspace)

    pngname = repstr(repstr(psname,'.ps','.png'),'.eps','.png')

;   for j = 9L, npage-1L do begin
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
            psfile=pspath+psname[j], /writeover, bw=0L
          cleanplot, /silent
;         dfpsplot, pspath+psname[j], /landscape, xsize=xpage, ysize=ypage, /color
       endif else begin
          if keyword_set(make_png) then set_plot, 'Z'
          arm_plotconfig, /landscape, nx=ncols, ny=nrows, xmargin=xmargin1, $
            ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
            height=height1, coord=pos, xpage=ypage1, ypage=xpage1, bw=0
          cleanplot, /silent
       endelse

; define a common wavelength vector for this subset of objects 

       i1 = j*nperpage
       i2 = (i1+(nperpage-1L))<(ngalaxy-1L)

       minwave = 3650.0 & maxwave = 5150.0
;      minwave = max(ages[i2:i2].wavemin/(1+ages[i1:i2].z_obj))
;      maxwave = min(ages[i2:i2].wavemax/(1+ages[i1:i2].z_obj))
       dwave = 1.6
       interpwave = findgen((maxwave-minwave)/dwave+1)*dwave+minwave

       if keyword_set(mypng) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')
       
; visualize each object       
       
       for k = 0L, nperpage-1L do begin
          
          indx = j*nperpage+k
;         print, indx

          if (k mod (ncols)) eq 0L then begin   ; label every column
;         if (k mod (2*ncols)) eq 0L then begin ; label every second column
             delvarx, ytickname
             ytitle = 'Relative Flux'
             ytickname = replicate(' ',10)
          endif else begin
             ytickname = replicate(' ',10)
             ytitle = ''
          endelse

          if (j eq (npage-1L)) then cutit = (ngalaxy mod nperpage)-(ncols) else cutit = nperpage-(ncols)
          if (k lt cutit) then begin
             xtickname = replicate(' ',10)
             xtitle = ''
          endif else begin
             delvarx, xtickname
             xtitle = 'Rest Wavelength (\AA)'
          endelse

          if (indx[0] lt ngalaxy) then begin

             if keyword_set(mypng) then begin
                noerase = 1L 
             endif else noerase = (k ne 0L)

             ages_display_spectrum, ages[indx], agesancillary[indx], lcharsize=lcharsize, labeltype=2L, $
               position=pos[*,k], xtickname=xtickname, ytickname=ytickname, $
               _extra=extra, xtitle=xtitle, ytitle=ytitle, $
               no_specfit=no_specfit, interpwave=interpwave, speccolor=speccolor, $
               labelcolor=labelcolor, postthick=postthick, plotthick=plotthick, $
               plotcolor=plotcolor, pcharsize=pcharsize, setyrange=4L, noerase=noerase, /silent
             labelinfo1 = ['M_{B} = '+string(agesancillary[indx].m_b,format='(F5.1)'),$
               '12+log(O/H) = '+string(ages[indx].zstrong_ew_12oh_kk04,format='(F4.2)')]
             legend, textoidl(labelinfo1), /left, /top, box=0, charsize=lcharsize, $
               charthick=postthick, /data, margin=0, position=[3775.0,!y.crange[1]*0.9], $
               textcolor=djs_icolor(labelcolor)
          endif

       endfor 

       if keyword_set(postscript) then begin
          dfpsclose 
          if keyword_set(mypng) then spawn, ['convert '+pspath+psname[j]+' '+pspath+repstr(psname[j],'.ps','.png')], /sh
       endif else begin
          if keyword_set(make_png) then begin
             img = tvrd()
             tvlct, r, g, b, /get
             write_png, pspath+pngname[j], img, r, g, b
             set_plot, 'X'
          endif else cc = get_kbrd(1)
       endelse

       cleanplot, /silent

    endfor 

return
end    
