pro kenn92_poster, postscript=postscript
; jm04may03uofa
; put the whole kenn92 mergers on one page

; initialize path names
    
    webpath = atlas_path(/web)+'kenn92/'

; restore all the fitting results

    atlas = read_kenn92(/silent)

    srt = reverse(sort(atlas.h_alpha_ew[0])) ; reverse sort by EW(Ha)
    atlas = atlas[srt]
    
    nicegalaxy = atlas.nice_galaxy
    galaxy = strcompress(atlas.galaxy,/remove)
    ngalaxy = n_elements(galaxy)

    if not keyword_set(postscript) then window, 0, xs=600*10.0/8.0, ys=600

    rebinpix = 600L
    bin1 = lindgen(rebinpix/2-1)*2
    bin2 = bin1+3L
    nsubpix = n_elements(bin1)
   
    nx = 7L
    ny = ceil(ngalaxy/float(nx))
    
    pagemaker, nx=nx, ny=ny, xspace=0.1, yspace=0.1, $;, xpage=40.0, ypage=32.0
      xmargin=0, ymargin=0, position=pos, /normal, /landscape

    psname = 'kenn92_poster.ps'
    pngname = repstr(psname,'.ps','.png')
    im_openclose, webpath+psname, $ ; xsize=40.0, ysize=32.0, $
      postscript=postscript, /square, /color, /landscape

    loadct, 13, /silent
    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

;   for j = 0L, 6-1L do begin
    for j = 0L, ngalaxy-1L do begin

       specdata = read_kenn92_specfit(galaxy[j],/silent)
       wave = float(frebin(specdata[*,0],rebinpix))
       spec = float(frebin(specdata[*,1],rebinpix,/total))
       
       yrange = minmax(spec)*[0.95,1.08]
       djs_plot, [0], [0], xsty=5, ysty=5, xthick=2.0, $
         ythick=2.0, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
         charsize=2.0, charthick=2.0, position=pos[*,j], yrange=yrange, $
         xrange=[3650,6950], /noerase
;      legend, nicegalaxy[j], /left, /top, box=0, charsize=0.8, charthick=2.0

       colors = 255*findgen(rebinpix)/rebinpix
       for k = 0L, nsubpix-1L do oplot, wave[bin1[k]:bin2[k]], $
         spec[bin1[k]:bin2[k]], color=djs_mean(colors[bin1[k]:bin2[k]]), ps=10

;      for k = 0L, nsubpix-1L do oplot, wave[k*bit:(k+1)*bit-1], $
;        spec[k*bit:(k+1)*bit-1], color=djs_mean(colors[k*bit:(k+1)*bit-1]), ps=10
;      djs_oplot, wave, spec, ps=10
;      for k = 0L, rebinpix-1L do oplot, wave[k]*[1,1], spec[k]*[1,1], color=colors[k], ps=10
       
    endfor

    im_openclose, postscript=postscript, /close

    if keyword_set(postscript) then begin
       spawn, ['convert -rotate 270 '+webpath+psname+' '+webpath+pngname], /sh
       spawn, ['gzip -f '+webpath+psname], /sh
    endif

stop    
    
return
end    
