pro atlas1d_lineplot, xatlas, yatlas, xerratlas, yerratlas, xregion=xregion, $
  yregion=yregion, xerrregion=xerrregion, yerrregion=yerrregion, xnfgs=xnfgs, $
  ynfgs=ynfgs, xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs, plottype=plottype, $
  legendtype=legendtype, postscript=postscript, _extra=extra, xfrac=xfrac, $
  yfrac=yfrac, overplot=overplot, nodata=nodata, xreverse=xreverse, yreverse=yreverse, $
  errorleft=errorleft, atlaspsize=atlaspsize, nfgspsize=nfgspsize, $
  hiicolor=hiicolor, atlascolor=atlascolor, talk=talk, atlasfill=atlasfill, $
  nfgsfill=nfgsfill, blackwhite=blackwhite, symthick=symthick1, coloroffset=coloroffset, $
  hiifill=hiifill, hiipsize=hiipsize, atlassym=atlassym, postthick=postthick
; jm04jun10uofa

    loadct, 0, /silent ; keep this to prevent FSC_COLOR from taking over all the color indices
    
    if keyword_set(talk) then talkcolor = 'white' else talkcolor = 'black'
    
    natlas = n_elements(xatlas)
    nnfgs = n_elements(xnfgs)
    nregion = n_elements(xregion)

    if (n_elements(atlassym) eq 0L) then atlassym = 106
    nfgssym = 105
    hiisym = 108

    if (n_elements(atlasfill) eq 0L) then atlasfill = 1L
    if (n_elements(nfgsfill) eq 0L) then nfgsfill = 1L
    if (n_elements(hiifill) eq 0L) then hiifill = 1L

    if (n_elements(coloroffset) eq 0L) then coloroffset = 0L

    if keyword_set(blackwhite) then begin
       if (n_elements(atlascolor) eq 0L) then atlascolor = 'dark gray' ; 'gray'
       if (n_elements(nfgscolor) eq 0L) then nfgscolor = 'black'
       if (n_elements(hiicolor) eq 0L) then hiicolor = 'gray' ; 'dark gray'
       if (n_elements(atlaspsize) eq 0L) then atlaspsize = 0.8
       if (n_elements(nfgspsize) eq 0L) then nfgspsize = 1.0
       if (n_elements(hiipsize) eq 0L) then hiipsize = 0.5
    endif else begin
       if (n_elements(atlascolor) eq 0L) then if keyword_set(talk) then atlascolor = 'sky blue' else atlascolor = 'blue'
       if (n_elements(nfgscolor) eq 0L) then nfgscolor = 'red'
       if (n_elements(hiicolor) eq 0L) then hiicolor = 'purple'
       if (n_elements(atlaspsize) eq 0L) then if keyword_set(talk) then atlaspsize = 1.1 else atlaspsize = 0.8
       if (n_elements(nfgspsize) eq 0L) then if keyword_set(talk) then nfgspsize = 1.4 else nfgspsize = 1.0
       if (n_elements(hiipsize) eq 0L) then hiipsize = 0.5
    endelse

;   if (natlas ne 0L) then splog, 'ATLAS: '+string(natlas,format='(I0)')+' data points.'
;   if (nnfgs ne 0L) then splog, 'NFGS: '+string(nnfgs,format='(I0)')+' data points.'
;   if (nregion ne 0L) then splog, 'HII Regions: '+string(nregion,format='(I0)')+' data points.'
    
; figure out the plot legend    
    
    label = ''
    lsym = [-1L]
    lfill = [0L]
    color = ''

    bigxerr = xerratlas
    bigyerr = yerratlas
    
    if (natlas gt 0L) then begin
       label = [label,'K/M ATLAS']
       lsym = [lsym,atlassym]
       color = [color,atlascolor]
       lfill = [lfill,atlasfill]
    endif
    if (nnfgs gt 0L) then begin
       label = [label,'NFGS']
       lsym = [lsym,nfgssym]
       color = [color,nfgscolor]
       lfill = [lfill,nfgsfill]
       bigxerr = [bigxerr,xerrnfgs]
       bigyerr = [bigyerr,yerrnfgs]
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
    
    xerravg = djs_mean(bigxerr) ; djs_median(xerr)
    yerravg = djs_mean(bigyerr) ; djs_median(yerr)

    if size(extra,/type) ne 8L then extra = {dummy: ''}
    if (n_elements(postthick) eq 0L) then postthick = 5.0
    
    if tag_exist(extra,'XSTYLE') then xstyle=extra.xstyle else xstyle = 3
    if tag_exist(extra,'YSTYLE') then ystyle=extra.ystyle else ystyle = 3
    if tag_exist(extra,'PSIZE') then psize=extra.psize else psize = 1.0

    if tag_exist(extra,'LCHARSIZE') then lcharsize=extra.lcharsize else lcharsize = 1.5
    if tag_exist(extra,'CHARSIZE') then charsize=extra.charsize else charsize = 2.0
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
    
    if (n_elements(symthick1) eq 0L) then if keyword_set(postscript) then symthick = 5L else $
      symthick = 2.0 else symthick = symthick1
    if keyword_set(postscript) then thick = postthick else thick = 2.0
    if keyword_set(postscript) then charthick = postthick else charthick = 2.0
    if keyword_set(postscript) then lcharthick = postthick else lcharthick = 2.0
    if keyword_set(postscript) then xthick = postthick else xthick = 2.0
    if keyword_set(postscript) then ythick = postthick else ythick = 2.0
    if keyword_set(postscript) then if keyword_set(talk) then errthick = 10.0 else errthick = postthick else errthick = 2.0

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

    if (not keyword_set(nodata)) then begin
       
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

;         oploterror, xc, yc, xerravg, yerravg, ps=3, errthick=errthick, /nohat

          plottype = 1L

       endif
       
       case plottype of
          1L: begin             ; points, no error bars
             if (nregion ne 0L) then begin
                col = fsc_color(hiicolor,!d.table_size-2)
;               col = djs_icolor(hiicolor)
                im_symbols, hiisym, fill=hiifill, psize=hiipsize, color=col, thick=symthick
                djs_oplot, xregion, yregion, psym=8
             endif
             if (natlas ne 0L) then begin
                col = fsc_color(atlascolor,!d.table_size-3-coloroffset)
;               col = djs_icolor(atlascolor)
                im_symbols, atlassym, fill=atlasfill, psize=atlaspsize, color=col, thick=symthick
                djs_oplot, xatlas, yatlas, psym=8
             endif
             if (nnfgs ne 0L) then begin
                col = fsc_color(nfgscolor,!d.table_size-4)
;               col = djs_icolor(nfgscolor)
                im_symbols, nfgssym, fill=nfgsfill, psize=nfgspsize, color=col, thick=symthick
                if (nnfgs eq 1L) then plots, xnfgs, ynfgs, psym=8 else djs_oplot, xnfgs, ynfgs, psym=8
             endif
          end
          2L: begin             ; points, error bars
             if (nregion ne 0L) then oploterror, xregion, yregion, xerrregion, yerrregion, psym=3, $
               /nohat, errthick=errthick, color=djs_icolor(hiicolor), errcolor=djs_icolor(hiicolor)
             if (natlas ne 0L) then oploterror, xatlas, yatlas, xerratlas, yerratlas, psym=3, $
               /nohat, errthick=errthick, color=djs_icolor(atlascolor), errcolor=djs_icolor(atlascolor)
             if (nnfgs ne 0L) then oploterror, xnfgs, ynfgs, xerrnfgs, yerrnfgs, psym=3, $
               /nohat, errthick=errthick, color=djs_icolor(nfgscolor), errcolor=djs_icolor(nfgscolor)
          end
          4L: begin             ; points, error bars
             if (nregion ne 0L) then begin
                col = fsc_color(hiicolor,!d.table_size-2)
                im_symbols, hiisym, fill=hiifill, psize=hiipsize, color=col, thick=symthick
                oploterror, xregion, yregion, xerrregion, yerrregion, psym=8, $
                  errthick=errthick, color=col, errcolor=col
             endif
             if (natlas ne 0L) then begin
                col = fsc_color(atlascolor,!d.table_size-3-coloroffset)
                im_symbols, atlassym, fill=atlasfill, psize=atlaspsize, color=col, thick=symthick
                oploterror, xatlas, yatlas, xerratlas, yerratlas, psym=8, $
                  errthick=errthick, color=col, errcolor=col
             endif
             if (nnfgs ne 0L) then begin
                col = fsc_color(nfgscolor,!d.table_size-4)
                im_symbols, nfgssym, fill=nfgsfill, psize=nfgspsize, color=col, thick=symthick
                oploterror, xnfgs, ynfgs, xerrnfgs, yerrnfgs, psym=8, $
                  errthick=errthick, color=col, errcolor=col
             endif
          end
       endcase

       if not keyword_set(legendtype) then legendtype = 0L

       if tag_exist(extra,'POSITION') then extra = struct_trimtags(extra,except='POSITION')

       case legendtype of
          0L:                   ; no legend
          1L: begin
             im_legend, textoidl(label), psym=lsym, charsize=lcharsize, charthick=lcharthick, $
               box=0, color=djs_icolor(color), symsize=1.3, fill=lfill, spacing=1.7, thick=thick, $
               _extra=extra, clear=keyword_set(postscript)
;         legend, ['('+strn(nindx)+')'], /right, /top, box=0, charsize=1.6, charthick=5.0, /clear
          end
       endcase
       
       if tag_exist(extra,'XTITLE') and tag_exist(extra,'YTITLE') then begin
;      print, extra.xtitle+' vs '+extra.ytitle+': '+strn(natlas,format='(I0)')+' data points: '
;      print, '['+string([!x.crange,minmax(x)],format='(4G0.0)')+'] ['+$
;        string([!y.crange,minmax(y)],format='(4F0.0)')+']'
       endif

       xout = where((xatlas gt !x.crange[1]) or (xatlas lt !x.crange[0]),nxout)
       yout = where((yatlas gt !y.crange[1]) or (yatlas lt !y.crange[0]),nyout)
;      if (nxout+nyout gt 0L) then $
;        print, 'There are ('+strn(nxout)+', '+strn(nyout)+') points outside the (x,y) data range.'

    endif
       
return
end    
