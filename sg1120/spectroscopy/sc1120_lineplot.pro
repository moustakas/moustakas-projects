pro sc1120_lineplot, x, y, xerr, yerr, xregion=xregion, $
  yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, $
  xlocal=xlocal, ylocal=ylocal, xerrlocal=xerrlocal, yerrlocal=yerrlocal, $
  xdgss=xdgss, ydgss=ydgss, xerrdgss=xerrdgss, yerrdgss=yerrdgss, $
  xtkrs=xtkrs, ytkrs=ytkrs, xerrtkrs=xerrtkrs, yerrtkrs=yerrtkrs, $
  xcfrs=xcfrs, ycfrs=ycfrs, xerrcfrs=xerrcfrs, yerrcfrs=yerrcfrs, $
  xliang=xliang, yliang=yliang, xerrliang=xerrliang, yerrliang=yerrliang, $
  xmaier=xmaier, ymaier=ymaier, xerrmaier=xerrmaier, yerrmaier=yerrmaier, $
  xkz99=xkz99, ykz99=ykz99, xerrkz99=xerrkz99, yerrkz99=yerrkz99, $
  hiilabel=hiilabel, plottype=plottype, legendtype=legendtype, $
  postscript=postscript, talk=talk, sc1120psize=sc1120psize, _extra=extra, $
  xreverse=xreverse, yreverse=yreverse
; jm04jan8uofa

    nsc1120 = n_elements(x)
    nlocal = n_elements(xlocal)
    ndgss = n_elements(xdgss)
    ntkrs = n_elements(xtkrs)
    ncfrs = n_elements(xcfrs)
    nliang = n_elements(xliang)
    nmaier = n_elements(xmaier)
    nkz99 = n_elements(xkz99)
    nregion = n_elements(xregion)

    sc1120sym = 106
    localsym = 108
    dgsssym = 105
    tkrssym = 115
    cfrssym = 122
    liangsym = 102
    maiersym = 104
    kz99sym = 121
    hiisym = 108

    sc1120fill = 1
    localfill = 0
    dgssfill = 1
    tkrsfill = 1
    cfrsfill = 1
    liangfill = 1
    maierfill = 1
    kz99fill = 1
    hiifill = 1

    if (n_elements(sc1120psize) eq 0L) then sc1120psize = 1.1
    localpsize = 0.6
    dgsspsize = 1.3
    tkrspsize = 1.6
    cfrspsize = 1.5
    liangpsize = 1.0
    maierpsize = 1.0
    kz99psize = 1.3
    hiipsize = 0.4

    sc1120color = 'dark green'
    localcolor = 'light blue'
    dgsscolor = 'blue'
    tkrscolor = 'dark red'
    cfrscolor = 'dark green'
    liangcolor = 'green'
    maiercolor = 'purple'
    kz99color = 'magenta'
    hiicolor = 'purple'
    
; figure out the plot legend    
    
    label = ''
    lsym = [-1L]
    lfill = [0L]
    color = ''
    
    if (nsc1120 gt 0L) then begin
       label = [label,'SC1120']
       lsym = [lsym,sc1120sym]
       color = [color,sc1120color]
       lfill = [lfill,sc1120fill]
    endif
    if (nlocal gt 0L) then begin
       label = [label,'Local Galaxy']
       lsym = [lsym,localsym]
       color = [color,localcolor]
       lfill = [lfill,localfill]
    endif
    if (ndgss gt 0L) then begin
       label = [label,'DGSS']
       lsym = [lsym,dgsssym]
       color = [color,dgsscolor]
       lfill = [lfill,dgssfill]
    endif
    if (ntkrs gt 0L) then begin
       label = [label,'TKRS']
       lsym = [lsym,tkrssym]
       color = [color,tkrscolor]
       lfill = [lfill,tkrsfill]
    endif
    if (ncfrs gt 0L) then begin
       label = [label,'CFRS']
       lsym = [lsym,cfrssym]
       color = [color,cfrscolor]
       lfill = [lfill,cfrsfill]
    endif
    if (nliang gt 0L) then begin
       label = [label,'LIANG']
       lsym = [lsym,liangsym]
       color = [color,liangcolor]
       lfill = [lfill,liangfill]
    endif
    if (nmaier gt 0L) then begin
       label = [label,'MAIER']
       lsym = [lsym,maiersym]
       color = [color,maiercolor]
       lfill = [lfill,maierfill]
    endif
    if (nkz99 gt 0L) then begin
       label = [label,'KZ99']
       lsym = [lsym,kz99sym]
       color = [color,kz99color]
       lfill = [lfill,kz99fill]
    endif
    if (nregion gt 0L) then begin
       label = [label,'HII Region']
       lsym = [lsym,hiisym]
       color = [color,hiicolor]
       lfill = [lfill,hiifill]
    endif

    label = label[1:n_elements(label)-1]
    lsym = lsym[1:n_elements(lsym)-1]
    color = color[1:n_elements(color)-1]
    lfill = lfill[1:n_elements(lfill)-1]
       
    xerravg = mean(xerr,/nan) ; djs_median(xerr)
    yerravg = mean(yerr,/nan) ; djs_median(yerr)

    if size(extra,/type) ne 8L then extra = {dummy: ''}
    postthick = 8.0
    
    if tag_exist(extra,'XSTYLE') then xstyle=extra.xstyle else xstyle = 3
    if tag_exist(extra,'YSTYLE') then ystyle=extra.ystyle else ystyle = 3
    if tag_exist(extra,'PSIZE') then psize=extra.psize else psize = 1.0
    if tag_exist(extra,'LCHARSIZE') then lcharsize=extra.lcharsize else lcharsize = 1.5
    if tag_exist(extra,'CHARSIZE') then charsize=extra.charsize else charsize = 2.0
    if tag_exist(extra,'XRANGE') then xrange=extra.xrange else xrange = [0.0,1.0]
    if tag_exist(extra,'YRANGE') then yrange=extra.yrange else yrange = [0.0,1.0]
    if tag_exist(extra,'NOERASE') then noerase=extra.noerase else noerase = 0

    if tag_exist(extra,'XTITLE') then begin
       xtitle=textoidl(extra.xtitle)
       extra = struct_trimtags(extra,except='XTITLE')
    endif else xtitle = ''
    if tag_exist(extra,'YTITLE') then begin
       ytitle=textoidl(extra.ytitle) 
       extra = struct_trimtags(extra,except='YTITLE')
    endif else ytitle = ''

    if keyword_set(postscript) then charthick = postthick else charthick = 2.0
    if keyword_set(postscript) then lcharthick = postthick else lcharthick = 2.0
    if keyword_set(postscript) then xthick = postthick else xthick = 2.0
    if keyword_set(postscript) then ythick = postthick else ythick = 2.0
    if keyword_set(postscript) then errthick = postthick else errthick = 2.0

    if keyword_set(talk) then begin
       if (noerase ne 1L) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')
       plot, [0], [0], charsize=charsize, charthick=charthick, xthick=xthick, $
         ythick=ythick, xstyle=xstyle, ystyle=ystyle, color=djs_icolor(talkcolor), $
         xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), _extra=extra, /noerase
    endif else begin
       plot, [0], [0], charsize=charsize, charthick=charthick, xthick=xthick, $
         ythick=ythick, xstyle=xstyle, ystyle=ystyle, xtitle=textoidl(xtitle), $
         ytitle=textoidl(ytitle), _extra=extra
    endelse

;;; display the data
;;
;;    bin = 0.025
;;    
;;    if keyword_set(xreverse) then ximage = (xrange[0] - x) / bin else ximage = (x - xrange[0]) / bin
;;    if keyword_set(yreverse) then yimage = (yrange[0] - y) / bin else yimage = (y - yrange[0]) / bin
;;    
;;    nx = fix(abs(xrange[1]-xrange[0])/bin+1)
;;    ny = fix(abs(yrange[1]-yrange[0])/bin+1)
;;    image = fltarr(nx,ny)
;;
;;    populate_image, image, ximage, yimage, assign='cic'
;;
;;    topvalue = 255L
;;    minvalue = -100L
;;
;;    logimage = alog10(image+1)
;;    img = bytscl(logimage,min=min(logimage),max=max(logimage),top=topvalue)
;;
;;    if keyword_set(postscript) then img = bytscl(topvalue-img,min=minvalue,top=topvalue)
;;
;;    plotimage, img, charsize=charsize, charthick=charthick, xthick=xthick, $
;;      ythick=ythick, xstyle=xstyle, ystyle=ystyle, imgxrange=xrange, imgyrange=yrange, $
;;      xtitle=xtitle, ytitle=ytitle, _extra=extra

    if n_elements(plottype) eq 0L then plottype = 1L

    if (plottype eq 3L) then begin ; points, average error bar

       xoff = (2.0*xerravg) > (0.075*(!x.crange[1]-!x.crange[0]))
       yoff = (2.0*yerravg) > (0.075*(!y.crange[1]-!y.crange[0]))

       xc = !x.crange[1]-xoff
       yc = !y.crange[0]+yoff

       oploterror, xc, yc, xerravg, yerravg, ps=3, errthick=errthick, /nohat, $
         color=djs_icolor(talkcolor), errcolor=djs_icolor(talkcolor)
       plottype = 1L

    endif
    
    case plottype of
       1L: begin ; points, no error bars
          if (nlocal ne 0L) then begin
             im_symbols, localsym, fill=localfill, psize=localpsize, color=djs_icolor(localcolor)
             if (nlocal eq 1L) then plots, xlocal, ylocal, psym=8 else djs_oplot, xlocal, ylocal, psym=8
          endif
          if (ndgss ne 0L) then begin
             im_symbols, dgsssym, fill=dgssfill, psize=dgsspsize, color=djs_icolor(dgsscolor)
             if (ndgss eq 1L) then plots, xdgss, ydgss, psym=8 else djs_oplot, xdgss, ydgss, psym=8
          endif
          if (ntkrs ne 0L) then begin
             im_symbols, tkrssym, fill=tkrsfill, psize=tkrspsize, color=djs_icolor(tkrscolor)
             if (ntkrs eq 1L) then plots, xtkrs, ytkrs, psym=8 else djs_oplot, xtkrs, ytkrs, psym=8
          endif
          if (ncfrs ne 0L) then begin
             im_symbols, cfrssym, fill=cfrsfill, psize=cfrspsize, color=djs_icolor(cfrscolor)
             if (ncfrs eq 1L) then plots, xcfrs, ycfrs, psym=8 else djs_oplot, xcfrs, ycfrs, psym=8
          endif
          if (nliang ne 0L) then begin
             im_symbols, liangsym, fill=liangfill, psize=liangpsize, color=djs_icolor(liangcolor)
             if (nliang eq 1L) then plots, xliang, yliang, psym=8 else djs_oplot, xliang, yliang, psym=8
          endif
          if (nmaier ne 0L) then begin
             im_symbols, maiersym, fill=maierfill, psize=maierpsize, color=djs_icolor(maiercolor)
             if (nmaier eq 1L) then plots, xmaier, ymaier, psym=8 else djs_oplot, xmaier, ymaier, psym=8
          endif
          if (nkz99 ne 0L) then begin
             im_symbols, kz99sym, fill=kz99fill, psize=kz99psize, color=djs_icolor(kz99color)
             if (nkz99 eq 1L) then plots, xkz99, ykz99, psym=8 else djs_oplot, xkz99, ykz99, psym=8
          endif
          if (nregion ne 0L) then begin
             im_symbols, hiisym, fill=hiifill, psize=hiipsize, color=djs_icolor(hiicolor)
             djs_oplot, xregion, yregion, psym=8
          endif
          if (nsc1120 ne 0L) then begin
             im_symbols, sc1120sym, fill=sc1120fill, psize=sc1120psize, color=djs_icolor(sc1120color)
             djs_oplot, x, y, psym=8
          endif
       end
       2L: begin ; points, error bars
          if (nsc1120 ne 0L) then oploterror, x, y, xerr, yerr, psym=3, $
            /nohat, errthick=errthick, color=djs_icolor(sc1120color), errcolor=djs_icolor(sc1120color)
          if (nlocal ne 0L) then oploterror, xlocal, ylocal, xerrlocal, yerrlocal, psym=3, $
            /nohat, errthick=errthick, color=djs_icolor(localcolor), errcolor=djs_icolor(localcolor)
          if (ndgss ne 0L) then oploterror, xdgss, ydgss, xerrdgss, yerrdgss, psym=3, $
            /nohat, errthick=errthick, color=djs_icolor(dgsscolor), errcolor=djs_icolor(dgsscolor)
          if (ntkrs ne 0L) then oploterror, xtkrs, ytkrs, xerrtkrs, yerrtkrs, psym=3, $
            /nohat, errthick=errthick, color=djs_icolor(tkrscolor), errcolor=djs_icolor(tkrscolor)
          if (ncfrs ne 0L) then oploterror, xcfrs, ycfrs, xerrcfrs, yerrcfrs, psym=3, $
            /nohat, errthick=errthick, color=djs_icolor(cfrscolor), errcolor=djs_icolor(cfrscolor)
          if (nliang ne 0L) then oploterror, xliang, yliang, xerrliang, yerrliang, psym=3, $
            /nohat, errthick=errthick, color=djs_icolor(liangcolor), errcolor=djs_icolor(liangcolor)
          if (nmaier ne 0L) then oploterror, xmaier, ymaier, xerrmaier, yerrmaier, psym=3, $
            /nohat, errthick=errthick, color=djs_icolor(maiercolor), errcolor=djs_icolor(maiercolor)
          if (nkz99 ne 0L) then oploterror, xkz99, ykz99, xerrkz99, yerrkz99, psym=3, $
            /nohat, errthick=errthick, color=djs_icolor(kz99color), errcolor=djs_icolor(kz99color)
          if (nregion ne 0L) then oploterror, xregion, yregion, xerrregion, yerrregion, psym=3, $
            /nohat, errthick=errthick, color=djs_icolor(hiicolor), errcolor=djs_icolor(hiicolor)
       end
    endcase

    if not keyword_set(legendtype) then legendtype = 0L

    case legendtype of
       0L: ; no legend
       1L: begin
          im_legend, textoidl(label), psym=lsym, charsize=lcharsize, charthick=lcharthick, $
            box=0, color=djs_icolor(color), symsize=1.3, fill=lfill, spacing=1.7, _extra=extra;, /clear
       end
    endcase
    
    if tag_exist(extra,'XTITLE') and tag_exist(extra,'YTITLE') then begin
       print, extra.xtitle+' vs '+extra.ytitle+': '+string(nsc1120,format='(I0)')+' data points: '
    endif

    xout = where((x gt !x.crange[1]) or (x lt !x.crange[0]),nxout)
    yout = where((y gt !y.crange[1]) or (y lt !y.crange[0]),nyout)
    if (nxout+nyout gt 0L) then $
      print, 'There are ('+strn(nxout)+', '+strn(nyout)+') points outside the (x,y) data range.'

return
end    
