pro wavelength_redshift, postscript=postscript, blackwhite=blackwhite
; jm05jun07uofa - plot redshift versus wavelength, showing the
;                 observational accessibility of various SFR
;                 diagnostics 

    if keyword_set(postscript) then begin
       postthick = 8.0 
    endif else begin
       postthick = 2.0
       im_window, 0, xratio=0.6, /square
    endelse

    if keyword_set(blackwhite) then begin
       pspath = atlas_path(/web)+'analysis/sfrs/blackwhite/'
       psname = 'wavelength_redshift'
       color = 0L
    endif else begin
       pspath = atlas_path(/papers)+'sfrs/FIG_SFRS/'
       psname = 'wavelength_redshift'
       color = 1L
    endelse
    
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, /encapsulated, color=color

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.1, height=7.1, $
      xmargin=[1.1,0.3], ymargin=[0.3,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    zmax = 5.2
    zmin = 0.0
    dz = 0.1
    zgrid = findgen((zmax-zmin)/dz+1)*dz+zmin

    lamgrid = findgen((26E4-3500.0)/1.0+1)*1.0+3500.0
    lamgrid = lamgrid / 1E4

;   lamrest = [3727.0,4861.0,6563.0] / 1E4 ; OII, H-beta, OIII, Ha
;   linename = ['[O II]','Hb','Ha']
;   linestyle = [0,3,5]

    lamrest = [3727.0,4340.0,4861.0,6563.0] / 1E4 ; OII, H-beta, OIII, Ha
    linename = ['[O II]','Hg','Hb','Ha']
    linestyle = [0,1,3,5]

;   linestyle = [0,3,5]
;   lamrest = [3727.0,4861.0,5007.0,6563.0] / 1E4 ; OII, H-beta, OIII, Ha
;   linename = ['[O II]','Hb','[O III]','Ha']
;   linestyle = [0,2,3,5]
    
;   zgrid = alog10(1+zgrid)
;   lamgrid = alog10(lamgrid)
;   lamrest = alog10(lamrest)

;   xtitle = 'log (\lambda) ['+angstrom()+']' 
    xtitle = 'Observed Wavelength [\mu'+'m]'
;   xtitle = 'Wavelength ['+angstrom()+']'
    ytitle = ' Redshift' ; 'Redshift'
;   ytitle = 'log (1+z)' ; 'Redshift'

    xrange = [1600.0,26E3] / 1E4
;   xrange = [3000.0,30E3]
;   xrange = alog10([3000.0,30E3])
;   xrange = [min(lamrest)+min(zgrid),max(lamrest)+max(zgrid)]
    yrange = minmax(zgrid)
    
    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
      charsize=2.0, charthick=postthick, xthick=postthick, ythick=postthick, $
      xsty=1, ysty=1, xtitle=xtitle, ytitle=ytitle, position=pos[*,0];, /xlog, /ylog, $
;     ytickinterval=0.01, ytickname=['0','1','2','3','4','5'], yminor=50 ;, yticks=4
    
    for i = 0L, n_elements(lamrest)-1L do djs_oplot, lamrest[i]*(1+zgrid), $
      zgrid, line=linestyle[i], thick=postthick

; fill the optical wavelength range
    
    polyfill, [3500.0,3500.0,9500,9500] / 1E4, [!y.crange[0],!y.crange[1],$
      !y.crange[1],!y.crange[0]], color=djs_icolor('orange'), /line_fill, $
      spacing=0.1, orientation=135
    polyfill, [3500.0,3500.0,9500,9500] / 1E4, [!y.crange[0],!y.crange[1],$
      !y.crange[1],!y.crange[0]], color=djs_icolor('orange'), /line_fill, $
      spacing=0.1, orientation=45
    xyouts, 6500.0/1E4, 4.5, 'Optical', charsize=2.0, charthick=postthick, align=0.5, /data

    polyfill, [11E3,11E3,14E3,14E3] / 1E4, [!y.crange[0],!y.crange[1],$
      !y.crange[1],!y.crange[0]], color=djs_icolor('blue'), /line_fill, $
      spacing=0.1, orientation=135
    polyfill, [11E3,11E3,14E3,14E3] / 1E4, [!y.crange[0],!y.crange[1],$
      !y.crange[1],!y.crange[0]], color=djs_icolor('blue'), /line_fill, $
      spacing=0.1, orientation=45
    xyouts, 12.5E3/1E4, 4.5, 'J', charsize=2.0, charthick=postthick, align=0.5, /data

    polyfill, [15E3,15E3,18E3,18E3] / 1E4, [!y.crange[0],!y.crange[1],$
      !y.crange[1],!y.crange[0]], color=djs_icolor('dark green'), /line_fill, $
      spacing=0.1, orientation=135
    polyfill, [15E3,15E3,18E3,18E3] / 1E4, [!y.crange[0],!y.crange[1],$
      !y.crange[1],!y.crange[0]], color=djs_icolor('dark green'), /line_fill, $
      spacing=0.1, orientation=45
    xyouts, 16.5E3/1E4, 4.5, 'H', charsize=2.0, charthick=postthick, align=0.5, /data

    polyfill, [20E3,20E3,23E3,23E3] / 1E4, [!y.crange[0],!y.crange[1],$
      !y.crange[1],!y.crange[0]], color=djs_icolor('red'), /line_fill, $
      spacing=0.1, orientation=135
    polyfill, [20E3,20E3,23E3,23E3] / 1E4, [!y.crange[0],!y.crange[1],$
      !y.crange[1],!y.crange[0]], color=djs_icolor('red'), /line_fill, $
      spacing=0.1, orientation=45
    xyouts, 21.5E3/1E4, 4.5, 'K', charsize=2.0, charthick=postthick, align=0.5, /data

; legend

;   legend, textoidl(['[O II]','H\beta','H\alpha']), /right, $
    legend, textoidl(['[O II]','H\gamma','H\beta','H\alpha']), /right, $
;   legend, textoidl(['[O II]','H\beta','H\alpha']), /right, $
;   legend, textoidl(['[O II]','H\beta','[O III]','H\alpha']), /right, $
      /bottom, charsize=1.2, charthick=postthick, line=linestyle, thick=postthick, $
      clear=keyword_set(postscript), spacing=1.5, box=0

; when are each of the emission lines within the atmospheric windows?

    for iline = 0L, n_elements(linename)-1L do begin
       print, linename[iline]+':'
       print, ' Optical: <'+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),0.95),format='(F12.1)'),2)
       print, ' J-band : '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.1),format='(F12.1)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.4),format='(F12.1)'),2)
       print, ' H-band : '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.5),format='(F12.1)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),1.8),format='(F12.1)'),2)
       print, ' K-band : '+strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),2.0),format='(F12.1)'),2)+'-'+$
         strtrim(string(interpol(zgrid,lamrest[iline]*(1+zgrid),2.3),format='(F12.1)'),2)
       print
    endfor
    
    im_openclose, postscript=postscript, /close    

return
end
