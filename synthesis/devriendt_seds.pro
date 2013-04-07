pro devriendt_seds, sedcube, bandcube, ps=ps, color=color
; jm01sep6uofa
; plot the Devriendt SEDs

    path = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='plots')

    if not keyword_set(sedcube) then sedcube = read_sedcube(templates='devriendt')
    if not keyword_set(bandcube) then bandcube = read_bandcube()

    filters = ['Bw DWFS','R DWFS','K TIFKAM','IRAC 3.6','IRAC 4.5','IRAC 5.8','IRAC 8','MIPS 24','MIPS 70','MIPS 160']
    nbands = n_elements(filters)
    obands = filter_match(filters,bandcube.bandnames)

;   obands = strarr(nbands)
;   for i = 0L, nbands-1L do obands[i] = where(strmatch(bandcube.bandnames,filters[i],/fold_case) eq 1B)
    
    plotcolors = [44,30,78,83,24] ; corresponding to CZCCPLOTS
    linestyle = [0,2,3]
    
    colortable2, color
    plotfaves
;   linestyle = [0,1,2,3,0]

; ----------------------------------------------------------------------
; plot a selected set of SEDs
; ----------------------------------------------------------------------
    
    sedindx = [0,8,9]
;   sedindx = [8] ; M82
    nseds = n_elements(sedindx)

    if keyword_set(ps) then begin
       ps_open, path+'devriendt_figure', /ps_fonts, color=color
       device, /inches, /times;, xsize=7, ysize=7
    endif else window, 0, xs=450, ys=450

    xmin = 0.01
    xmax = 2E4

    ymin = 1.5
    ymax = 12.6
    
;    xmin = 0.03
;    xmax = 1000.0 ; 230000.0
;
;    ymin = 3
;    ymax = 12.5
    
    plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=[xmin,xmax], yr=[ymin,ymax], $
      /xlog, position=[0.1,0.1,0.95,0.93], $
      ytit=textoidl('log (\nu L_{\nu}) [L'+sunsymbol()+']'), $
;     xtit=textoidl('\lambda (\mu')+'m)', $
      xtit=textoidl('Wavelength (\mu')+'m)', xminor=9, $
      xthick=2.0, ythick=2.0, charsize=1.4, charthick=2.0

    for k = 0L, nseds-1L do begin

       lambda = *sedcube[sedindx[k]].lambda
       nu = 2.99793D14/lambda
       lum = nu*(*sedcube[sedindx[k]].mlum_nu)/3.826D26
       
       oplot, lambda, alog10(lum), linestyle=linestyle[k], thick=3.0, color=plotcolors[k]

    endfor

    plot, [0], [0], /noerase, /nodata, position=[0.1,0.1,0.95,0.3], $
      xsty=7, ysty=7, xtickname=replicate('',20), ytickname=replicate('',20), $
      xrange=[xmin,xmax], /xlog
    
    for i = 0L, nbands-1L do begin

       lambda0 = bandcube[obands[i]].lambda_eff*1E-4
       name = bandcube[obands[i]].mininames

       wband = *bandcube[obands[i]].wband*1E-4
       rband = *bandcube[obands[i]].rband
       rband = rband/max(rband)

       oplot, wband, rband, line=2, thick=1.5

; only overplot optical/NIR bandpass names
       
       if lambda0 lt 3.0 then $
         xyouts, [lambda0,lambda0], 1.1*[max(rband),max(rband)], name, align=0.5, $
         charsize=1.1, charthick=2.0
       
    endfor

    xyouts, [5.7,5.7], [1.1,1.1], 'IRAC', align=0.5, charsize=1.1, charthick=2.0
    xyouts, [68,68], [1.1,1.1], 'MIPS', align=0.5, charsize=1.1, charthick=2.0

    if keyword_set(ps) then ps_close

; ----------------------------------------------------------------------
; plot all the SEDs
; ----------------------------------------------------------------------

    if keyword_set(ps) then begin
       ps_open, path+'devriendt_seds', /ps_fonts, color=color
       device, /inches, /times;, xsize=7, ysize=7
    endif else window, 2, xs=450, ys=450

    plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=[xmin,xmax], yr=[ymin,ymax], $
      /xlog, position=[0.1,0.1,0.95,0.93], $
      ytit=textoidl('log (\nu L_{\nu}) [L'+sunsymbol()+']'), $
;     xtit=textoidl('\lambda (\mu')+'m)', $
    xtit=textoidl('Wavelength (\mu')+'m)', xminor=9, $
      xthick=2.0, ythick=2.0, charsize=1.4, charthick=2.0

    ntot = n_elements(sedcube)
    for k = 0L, ntot-1L do begin
       
       lambda = *sedcube[k].lambda
       nu = 2.99793D14/lambda
       lum = nu*(*sedcube[k].mlum_nu)/3.826D26
       
       cindx = [1+2*k]
       oplot, lambda, alog10(lum), thick=2.0
       
    endfor

    if keyword_set(ps) then ps_close

stop
    
    plotfaves, /restore
    

return
end
