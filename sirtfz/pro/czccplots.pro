pro czplot, colormag, nseds, zarray, ytitle, filename, plotcolors, $
            gtype, ps=ps, color=color, abmag=abmag
; create a color-redshift plot

    if not keyword_set(abmag) then ylog = 1L
    
    if keyword_set(ps) then begin
       ps_open, filename, /ps_fonts, /portrait, color=color
       device, /inches, /times, xsize=7, ysize=7
    endif else if !d.window eq -1L then window, 0, xs=550, ys=550

    if nseds gt 1L then psym = [8,lindgen(nseds-1)+4] else psym = [8]
    if nseds gt 1L then linestyle = [0,lindgen(nseds-1)+2] else linestyle = [0]
    
    for i = 0L, nseds-1L do begin ; loop on sed type

       good = where(finite(colormag[i,*]) eq 1B,ngood)
       if ngood eq 0L then message, 'Color-redshift relation not defined!'

       zcolor = colormag[i,good]
       redshift = zarray[good]

       if max(colormag,/nan) gt float(0) then yrange = [min(colormag,/nan),1.15*max(colormag,/nan)] else $
         yrange = [min(colormag,/nan),max(colormag,/nan)/1.15]
       
       if i eq 0L then begin
          
          plot, [0], [0], xr=[0,max(zarray)], yrange=yrange, /nodata, $
            charsize=1.5, charthick=2.0, xsty=3, ysty=3, ylog=ylog, $
            xtit='Redshift', ytitle=ytitle
          oplot, redshift, zcolor, thick=3.0, linestyle=linestyle[i], color=plotcolors[i];, psym=-psym[i]
          
       endif else oplot, redshift, zcolor, linestyle=linestyle[i], thick=3.0, color=plotcolors[i];, psym=-psym[i]
       
       xyouts, [0.79,0.79], [0.88,0.88]-i*0.015*nseds, gtype[i], /normal, align=0.0, $
         charsize=1.5, charthick=2.0
       plots, [0.7,0.77], [0.885,0.885]-i*0.015*nseds, linestyle=linestyle[i], $
         color=plotcolors[i], /normal, thick=3.0
;      plots, [0.885,0.885], [0.885,0.885]-i*0.015*nseds, psym=psym[i], color=plotcolors[i], /normal

    endfor 

    if keyword_set(ps) then ps_close

return
end

pro ccplot, xcolor, ycolor, nseds, xtitle, ytitle, filename, $
            plotcolors, gtype, ps=ps, color=color, abmag=abmag
; create a color-color plot

    if not keyword_set(abmag) then begin
       xlog = 1L
       ylog = 1L
    endif
    
    if keyword_set(ps) then begin
       ps_open, filename, /ps_fonts, /portrait, color=color
       device, /inches, /times, xsize=7, ysize=7
    endif else if !d.window eq -1L then window, 2, xs=550, ys=550

    good = where((xcolor gt float(0)) and (ycolor gt float(0)),ngood)
    
    xrange = fltarr(2)
    yrange = fltarr(2)

    factor = 1.2

    if min(xcolor,/nan) gt float(0) then xrange[0] = min(xcolor,/nan)/factor else xrange[0] = factor*min(xcolor,/nan)
    if max(xcolor,/nan) gt float(0) then xrange[1] = factor*max(xcolor,/nan) else xrange[1] = max(xcolor,/nan)/factor

    if min(ycolor,/nan) gt float(0) then yrange[0] = min(ycolor,/nan)/factor else yrange[0] = factor*min(ycolor,/nan)
    if max(ycolor,/nan) gt float(0) then yrange[1] = factor*max(ycolor,/nan) else yrange[1] = max(ycolor,/nan)/factor

    if nseds gt 1L then psym = [8,lindgen(nseds-1)+4] else psym = [8]
;   if nseds gt 1L then linestyle = [0,lindgen(nseds-1)+2] else linestyle = [0]
    
    for i = 0L, nseds-1L do begin

       good = where((finite(xcolor[i,*]) eq 1B) and (finite(ycolor[i,*]) eq 1B),ngood)
       if ngood eq 0L then message, 'Color-color plot not defined!'

       xxcolor = reform(xcolor[i,good])
       yycolor = reform(ycolor[i,good])

       if i eq 0L then begin
          
          plot, [xrange[0],xrange[1]], [yrange[0],yrange[1]], xrange=xrange, yrange=yrange, /nodata, $
            charsize=1.5, charthick=2.0, xsty=3, ysty=3, xtitle=xtitle, ytitle=ytitle, xlog=xlog, ylog=ylog
          oplot, xxcolor, yycolor, psym=-psym[i], thick=2.0, color=plotcolors[i], $
            linestyle=0 ; linestyle=linestyle[i]
          
       endif else oplot, xxcolor, yycolor, psym=-psym[i], thick=2.0, color=plotcolors[i], $
         linestyle=0 ; linestyle=linestyle[i]

       xyouts, 1.05*[xxcolor[0],xxcolor[0]], 1.05*[yycolor[0],yycolor[0]], gtype[i], /data

    endfor

    if keyword_set(ps) then ps_close

return
end

pro czccplots, filters, sedtype, templates=templates, zmin=zmin, $
               zmax=zmax, dz=dz, abmag=abmag, ps=ps, color=color
;+
; NAME:
;	CZCCPLOTS
;
; PURPOSE:
;	Generate color-redshift and color-color plots of model SEDs. 
;
; CALLING SEQUENCE:
;
; INPUTS:
;	sedtype   - index array of the SEDs to plot
;	filters   - filters for which to generate plots
;
; OPTIONAL INPUTS:
;	templates - SED templates
;	zmin      - minimum redshift for plotting
;	zmax      - maximum redshift
;	dz        - redshift interval
;	
; KEYWORD PARAMETERS:
;	abmag - AB magnitude color-color plots instead of flux ratios 
;	ps    - generate postscript output for each plot
;       color - generate color postscript plots (default to black
;               & white) 
;
; OUTPUTS:
;
; COMMON BLOCKS:
;	sirtf_simulations
;
; COMMENTS:
;	When specifying SED types the following is for reference:
;	Devriendt: [0,8,15] -> [E/S0,Starburst,ULIRG]
;	CWW: [0,1,3] -> [E,Sbc,Im]
;
;	This program calculates are possible color-redshift and
;	color-color combinations for the specified filters.  Careful:
;	there are [N*(N-1)/2] combinations for N filters!
;
;	Add a realistic error distribution to the color-redshift
;	diagrams to demonstrate what photometric errors do (Benitez
;	2000).
;
; EXAMPLE:
;
; PROCEDURES USED:
;	READ_MODEL_GRIDS, FILTER_MATCH(), COLORTABLE2, PS_OPEN,
;	PS_CLOSE, PLOTSYM, LEGEND, GET_ELEMENT
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 September 3, U of A
;-

    common sirtf_simulations
    common cosmology

    if n_params() ne 2L then begin
       print, 'Syntax - czccplots, filters, sedtype, [templates=], [zmin=], '+$
         '[zmax=], [dz=], abmag=abmag, ps=ps, color=color'
       return       
    endif

    if size(filters,/type) ne 7L then message, 'Filters Must be type string!' 
    sedtype = long(sedtype)

    if not keyword_set(templates) then templates = sirtf.templates
    
; restore the model grids
    
    read_model_grids, filters, mzarray, filterflux, colorflux, kcorrection, templates=templates

    obands = filter_match(filters,sirtf.bandcube.bandnames) ; observed bands
    nbands = n_elements(obands)                             ; number of filters

; plot settings and filter IDs for postscript output
    
    czpath = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='plots/czplots')
    ccpath = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='plots/ccplots')

    colortable2, c
    plotsym, 0, 1.0, /fill
    plotfaves

    plotcolors = [44,30,78,83,24]
    mini = sirtf.bandcube[obands].mininames
    ab2vega = sirtf.bandcube[obands].m_ab
    
    if max(sedtype) gt n_elements(sirtf.sedcube)-1L then message, 'Out of range SED type!'
    gtype = sirtf.sedcube[sedtype].gtype
    nseds = n_elements(sedtype)

; define the redshift array for plotting and read off the flux values
; at those redshifts through interpolation
    
    if not keyword_set(zmin) then zmin = 0.0
    if not keyword_set(zmax) then zmax = 6.0
    if not keyword_set(dz) then dz = 0.2

    zarray = (findgen((zmax-zmin+dz)/dz+1.0))*dz
    nz = n_elements(zarray)
    
; interpolate the flux at the specified redshifts (down to z=0)

    flux = fltarr(nseds,nz,nbands)
    for i = 0L, nseds-1L do for j = 0L, nbands-1L do $
      flux[i,*,j] = interpol(colorflux[sedtype[i],*,j],mzarray,zarray,/quad)

; ----------------------------------------------------------------------    
; generate a single color-redshift plot
; ----------------------------------------------------------------------    

    if nbands eq 2L then begin 

       if keyword_set(abmag) then colormag = -2.5*alog10(flux[*,*,0]/flux[*,*,1]) else $
         colormag = flux[*,*,0]/flux[*,*,1]
       
       if keyword_set(abmag) then ytitle = '-2.5 log [ F('+mini[0]+')/F('+mini[1]+') ]' else $
         ytitle = 'F('+mini[0]+')/F('+mini[1]+')'
       fname = 'cz_'+mini[0]+'_'+mini[1] ; postscript file name

;      colorvega = colormag[*,0] - (ab2vega[0]-ab2vega[1]) ; local colors
;      print, colorvega
       
       czplot, colormag, nseds, zarray, ytitle, czpath+fname, plotcolors, gtype, ps=ps, color=color
       
;      psyesno = 'y'
       psyesno = ''
       read, psyesno, prompt='Create postscript [y/n]? '
;      if (psyesno eq '') or (strlowcase(psyesno) eq 'yes') then $
       if (strlowcase(psyesno) eq 'y') then $
         czplot, colormag, nseds, zarray, ytitle, czpath+fname, plotcolors, gtype, /ps, color=color
       
; ----------------------------------------------------------------------    
; find all possible combinations of the filters
; ----------------------------------------------------------------------    

    endif else begin 

       ncombs = nbands*(nbands-1)/2     ; number of combinations
       colormag = fltarr(ncombs,nseds,nz)

; plot titles and postscript output file names
       
       ftitle = strarr(ncombs)
       fname = strarr(ncombs)
       count = 0L
       for j = 0L, nbands-2L do for k = j+1L, nbands-1L do begin
          if keyword_set(abmag) then $
            ftitle[count] = '-2.5 log [ F('+mini[j]+')/F('+mini[k]+') ]' else $
            ftitle[count] = 'F('+mini[j]+')/F('+mini[k]+')'
          fname[count] = mini[j]+'_'+mini[k]
          count = count + 1L
       endfor

; compute all the possible flux ratios
       
       for i = 0L, nseds-1L do begin 
          count = 0L
          for j = 0L, nbands-2L do for k = j+1L, nbands-1L do begin
             if keyword_set(abmag) then $
               colormag[count,i,*] = -2.5*alog10(flux[i,*,j]/flux[i,*,k]) else $
               colormag[count,i,*] = flux[i,*,j]/flux[i,*,k]
             count = count + 1L
           endfor
       endfor

; ----------------------------------------------------------------------    
; plot all the possible color-redshift combinations
; ----------------------------------------------------------------------    

       for k = 0L, ncombs-1L do begin

          czplot, reform(colormag[k,*,*],nseds,nz), nseds, zarray, ftitle[k], $
            czpath+'cz_'+fname[k], plotcolors, gtype, ps=ps, color=color, abmag=abmag

;         psyesno = 'y'
          psyesno = ''
          read, psyesno, prompt='Create postscript [y/n]? '
;         if (psyesno eq '') or (strlowcase(psyesno) eq 'yes') then $
          if (strlowcase(psyesno) eq 'y') then $
            czplot, reform(colormag[k,*,*],nseds,nz), nseds, zarray, ftitle[k], $
            czpath+'cz_'+fname[k], plotcolors, gtype, /ps, color=color, abmag=abmag

       endfor 

; ----------------------------------------------------------------------    
; plot all the possible color-color combinations
; ----------------------------------------------------------------------    

       for i = 0L, ncombs-2L do for j = i+1L, ncombs-1L do begin

          ycolor = reform(colormag[i,*,*],nseds,nz)
          xcolor = reform(colormag[j,*,*],nseds,nz)

          ccplot, xcolor, ycolor, nseds, ftitle[j], ftitle[i], ccpath+'cc_'+fname[j]+'_'+fname[i], $
            plotcolors, gtype, ps=ps, color=color, abmag=abmag

;         psyesno = 'y'
          psyesno = ''
          read, psyesno, prompt='Create postscript [y/n]? '
;         if (psyesno eq '') or (strlowcase(psyesno) eq 'yes') then $
          if (strlowcase(psyesno) eq 'y') then $
            ccplot, xcolor, ycolor, nseds, ftitle[j], ftitle[i], ccpath+'cc_'+fname[j]+'_'+fname[i], $
            plotcolors, gtype, /ps, color=color, abmag=abmag

       endfor
       
    endelse
       
stop

return
end







