pro zindicators_lineplot, x, y, xerr, yerr, xgal=xgal, ygal=ygal, xerrgal=xerrgal, $
  yerrgal=yerrgal, plottype=plottype, legendtype=legendtype, postscript=postscript, $
  xfrac=xfrac, yfrac=yfrac, hiipsize=hiipsize, hiicolor=hiicolor, _extra=extra, $
  errorleft=errorleft, xreverse=xreverse, yreverse=yreverse, talk=talk, postthick=postthick
; jm04jan08uofa
; jm04dec09uofa - optional input - hiipsize

    if keyword_set(talk) then talkcolor = 'white' else talkcolor = 'black'

    nindx = n_elements(x)
    ngal = n_elements(xgal)

    hiisym = 108
    hiifill = 1L
    if (n_elements(hiicolor) eq 0L) then hiicolor = 'grey'
    if (n_elements(hiipsize) eq 0L) then hiipsize = 0.6

    galsym = 106
    galfill = 0L
    galcolor = 'red' ; 'blue'
    galpsize = 1.0

    label = 'HII Region'
    lsym = hiisym
    lfill = hiifill
    color = hiicolor
    
    if (ngal gt 0L) then begin

       label = [label,'Galaxy']
       lsym = [lsym,galsym]
       lfill = [lfill,galfill]
       color = [color,galcolor]

       bigxerr = [xerr,xerrgal]
       bigyerr = [yerr,yerrgal]

    endif else begin

       bigxerr = xerr
       bigyerr = yerr
       
    endelse

    xerravg = djs_median(bigxerr) ; djs_mean(bigxerr)
    yerravg = djs_median(bigyerr) ; djs_mean(bigyerr)

    if size(extra,/type) ne 8L then extra = {dummy: ''}
    if (n_elements(postthick) eq 0L) then postthick = 5.0
    
    if tag_exist(extra,'XSTYLE') then xstyle=extra.xstyle else xstyle = 3
    if tag_exist(extra,'YSTYLE') then ystyle=extra.ystyle else ystyle = 3
    if tag_exist(extra,'PSIZE') then psize=extra.psize else psize = 1.0
    if tag_exist(extra,'LCHARSIZE') then lcharsize=extra.lcharsize else lcharsize = 1.5
    if tag_exist(extra,'CHARSIZE') then charsize=extra.charsize else charsize = 2.0
    if tag_exist(extra,'XRANGE') then xrange=extra.xrange else xrange = [0.0,1.0]
    if tag_exist(extra,'YRANGE') then yrange=extra.yrange else yrange = [0.0,1.0]
    if tag_exist(extra,'NOERASE') then noerase=extra.noerase else noerase = 0

    if tag_exist(extra,'XTITLE') then begin
       xtitle = textoidl(extra.xtitle)
       extra = struct_trimtags(extra,except='XTITLE')
    endif else xtitle = ''
    if tag_exist(extra,'YTITLE') then begin
       ytitle = textoidl(extra.ytitle) 
       extra = struct_trimtags(extra,except='YTITLE')
    endif else ytitle = ''

    if keyword_set(postscript) then charthick = postthick else charthick = 2.0
    if keyword_set(postscript) then lcharthick = postthick else lcharthick = 2.0
    if keyword_set(postscript) then xthick = postthick else xthick = 2.0
    if keyword_set(postscript) then ythick = postthick else ythick = 2.0
    if keyword_set(postscript) then errthick = postthick else errthick = 2.0

    if keyword_set(talk) then begin
;      if (noerase ne 1L) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')
       plot, [0], [0], charsize=charsize, charthick=charthick, xthick=xthick, $
         ythick=ythick, xstyle=xstyle, ystyle=ystyle, color=djs_icolor('white'), $
         xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), _extra=extra, /noerase
    endif else begin
       if (not keyword_set(overplot)) then $
         plot, [0], [0], charsize=charsize, charthick=charthick, xthick=xthick, $
           ythick=ythick, xstyle=xstyle, ystyle=ystyle, xtitle=textoidl(xtitle), $
           ytitle=textoidl(ytitle), _extra=extra
    endelse

;   if (not keyword_set(overplot)) then $
;     djs_plot, [0], [0], charsize=charsize, charthick=charthick, xthick=xthick, $
;       ythick=ythick, xstyle=xstyle, ystyle=ystyle, _extra=extra

    if n_elements(plottype) eq 0L then plottype = 1L

    if (plottype eq 3L) then begin ; points, average error bar

       if (n_elements(xfrac) eq 0L) then xfrac = 2.0
       if (n_elements(yfrac) eq 0L) then yfrac = 2.0

       xoff = (xfrac*xerravg) > (0.075*(!x.crange[1]-!x.crange[0]))
       yoff = (yfrac*yerravg) > (0.075*(!y.crange[1]-!y.crange[0]))

       if keyword_set(errorleft) then begin
          if keyword_set(xreverse) then xc = !x.crange[0]-xoff else xc = !x.crange[0]+xoff
       endif else begin
          if keyword_set(xreverse) then xc = !x.crange[1]+xoff else xc = !x.crange[1]-xoff
       endelse

       if keyword_set(yreverse) then yc = !y.crange[0]-yoff else yc = !y.crange[0]+yoff

       if keyword_set(postscript) then $
         oploterror, xc, yc, xerravg, yerravg, ps=3, errthick=errthick, $ ; /nohat, $
           color=djs_icolor(talkcolor), errcolor=djs_icolor(talkcolor) else $
         oploterror, xc, yc, xerravg, yerravg, ps=3, errthick=errthick;, /nohat

       plottype = 1L

    endif
    
    case plottype of
       1L: begin ; points, no error bars
          im_symbols, hiisym, fill=hiifill, psize=hiipsize, color=djs_icolor(hiicolor)
          djs_oplot, x, y, psym=8
          if (ngal ne 0L) then begin
             im_symbols, galsym, fill=galfill, psize=galpsize, color=djs_icolor(galcolor)
             djs_oplot, xgal, ygal, psym=8
          endif
       end
       2L: begin ; points, error bars
          oploterror, x, y, xerr, yerr, psym=3, errthick=errthick, $ ; /nohat, $
            color=djs_icolor(hiicolor), errcolor=djs_icolor(hiicolor), errstyle=0
          if (ngal ne 0L) then begin
             oploterror, xgal, ygal, xerrgal, yerrgal, psym=3, $ ; /nohat, 
               errthick=errthick, color=djs_icolor(galcolor), errcolor=djs_icolor(galcolor), errstyle=0
          endif
       end 
    endcase
    
    if not keyword_set(legendtype) then legendtype = 0L

    case legendtype of
       0L: ; no legend
       1L: begin
          im_legend, textoidl(label), psym=lsym, charsize=1.5, charthick=lcharthick, $
            box=0, color=djs_icolor(color), symsize=1.3, fill=lfill, spacing=1.7, _extra=extra
       end
    endcase
    
    if tag_exist(extra,'XTITLE') and tag_exist(extra,'YTITLE') then begin
;      print, extra.xtitle+' vs '+extra.ytitle+': '+strn(nindx,format='(I0)')+' data points: '
;      print, '['+string([!x.crange,minmax(x)],format='(4G0.0)')+'] ['+$
;        string([!y.crange,minmax(y)],format='(4F0.0)')+']'
    endif

    xout = where((x gt !x.crange[1]) or (x lt !x.crange[0]),nxout)
    yout = where((y gt !y.crange[1]) or (y lt !y.crange[0]),nyout)
    if (nxout+nyout gt 0L) then $
      print, 'There are ('+strn(nxout)+', '+strn(nyout)+') points outside the (x,y) data range.'

return
end    
