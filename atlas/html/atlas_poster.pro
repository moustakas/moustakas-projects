pro atlas_poster, postscript=postscript
; jm04may03uofa
; put the whole atlas on one page

; restore all the fitting results

    webpath = atlas_path(/web)
    atlas = read_integrated()

    srt = reverse(sort(atlas.h_alpha_ew[0])) ; reverse sort by EW(Ha)
    atlas = atlas[srt]
    
    galaxy = strcompress(atlas.galaxy,/remove)
    nicegalaxy = strtrim(atlas.nice_galaxy,2)
    ngalaxy = n_elements(galaxy)

    rebinpix = 600L
    bin1 = lindgen(rebinpix/2-1)*2
    bin2 = bin1+3L
    nsubpix = n_elements(bin1)
   
    nx = 21L
    ny = ceil(ngalaxy/float(nx))
    
    pagemaker, nx=nx, ny=ny, xspace=0.1, yspace=0.1, $;, xpage=40.0, ypage=32.0
      xmargin=0, ymargin=0, position=pos, /normal, /landscape

    psname = 'atlas_poster.ps'
    pngname = repstr(psname,'.ps','.png')

    if keyword_set(postscript) then begin
       dfpsplot, webpath+psname, /square, /color, /landscape; xsize=40.0, ysize=32.0, $
       loadct, 13, /silent
    endif else begin
       im_window, 0, /square
    endelse

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

;   for j = 0L, 6-1L do begin
    for j = 0L, ngalaxy-1L do begin

       print, format='("Galaxy ",I0,"/",I0,".",A1,$)', j+1, ngalaxy, string(13b)

       specdata = read_atlas_specfit(galaxy[j],/silent)
       wave = float(frebin(specdata[*,0],rebinpix))
       spec = float(frebin(specdata[*,1],rebinpix,/total))
       
       yrange = [min(smooth(spec,30)),max(spec)]*[0.95,1.08]
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

    if keyword_set(postscript) then begin
       dfpsclose
       spawn, ['convert -rotate 270 '+webpath+psname+' '+webpath+pngname], /sh
       spawn, ['gzip -f '+webpath+psname], /sh
    endif

stop    
    
return
end    
