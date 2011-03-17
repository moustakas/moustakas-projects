pro gto_visualize_dss, gto, ps=ps

; ###########################################################################    
; make thumbnail DSS visualizations
; ###########################################################################    

    if (n_elements(gto) eq 0L) then gto = gto_read_ancillary()
    ngalaxy = n_elements(gto)

    path = gto_path(/ancillary)
    dsspath = gto_path(/dss)
    dssfits = dsspath+strlowcase(strtrim(gto.galaxy,2))+'.fits.gz'

    if keyword_set(ps) then begin
       dfpsplot, path+'gto_dss.ps', /square, /color
       postthick = 4.0
    endif else begin
       im_window, 0, xratio=0.4, /square
       postthick = 2.0
    endelse

    for i = 0L, ngalaxy-1L do begin
       
       dssimage = readfits(dssfits[i],hdss,/silent)
       gsssextast, hdss, astr

       imsize = size(dssimage,/dimension)
       xsize = imsize[0] & xcen = xsize/2.0 & ysize = imsize[1] & ycen = ysize/2.0

       xpixscale = (astr.pltscl*1E-3*astr.xsz)/60 ; [arcmin/pixel]
       ypixscale = (astr.pltscl*1E-3*astr.ysz)/60 ; [arcmin/pixel]

       xaxis = (findgen(xsize)-xcen)*xpixscale ; centered on the image [arcsec]
       yaxis = (findgen(ysize)-ycen)*ypixscale ; centered on the image [arcsec]

       img = logscl(dssimage,exponent=1.0,negative=keyword_set(write),omin=35,omax=255)
       
       plotimage, img, /preserve_aspect, position=pos, /normal, imgxrange=minmax(xaxis), $
         imgyrange=minmax(yaxis), charsize=1.8, charthick=postthick, xthick=postthick, $
         ythick=postthick, xtitle=textoidl('\Delta\alpha [arcmin]'), $
         ytitle=textoidl('\Delta\delta [arcmin]');, color=djs_icolor('white')

       gsssadxy, astr, 15.0*im_hms2dec(gto[i].ra), im_hms2dec(gto[i].dec), xrc3, yrc3
       xrc3 = (xrc3 - xcen)*xpixscale & yrc3 = (yrc3 - ycen)*ypixscale

       if (strtrim(gto[i].galaxy,2) ne strtrim(gto[i].ned_galaxy,2)) then $
         gal = strtrim(gto[i].galaxy,2)+'='+strtrim(gto[i].ned_galaxy,2) else $
           gal = strtrim(gto[i].galaxy,2)
       legend, gal, /left, /top, box=0, charthick=postthick, charsize=1.7, $
         textcolor=djs_icolor('white')
       
       if (not keyword_set(ps)) then cc = get_kbrd(1)

    endfor
    
    if keyword_set(ps) then begin
       dfpsclose
       spawn, 'gzip -f '+path+'gto_dss.ps', /sh
    endif

return
end
