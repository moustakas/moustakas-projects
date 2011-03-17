pro lyforest_plots, postscript=postscript, plotname=plotname
;+
; NAME:
;	LYFOREST_PLOTS
;
; PURPOSE:
;	Calculate the D_A and D_B constants due to ly-forest
;	absorption in the IGM at a grid of redshifts.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;	postscript - generate postscript output of the plots
;	plotname   - filename of the plot
;
; OUTPUTS:
;	If requested, generates a plot in the plots subdirectory of
;	the SIRTFz distribution.
;
; COMMENTS:
;	Takes a few minutes to generate the plots.
;
; PROCEDURES USED:
;	LYFOREST(), PS_OPEN, PS_CLOSE, TEXTOIDL(), LEGEND
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 Aug 3, U of A
;-

    wave = (findgen(1/0.01)+1.0)*0.01 ; wavelength vector (micron) [100,10000 Angstrom]
    z = (findgen(6/0.25)+1.0)*0.25        ; redshift vector [0,6]
    nz = n_elements(z)

    d_aarr = fltarr(nz)
    d_barr = fltarr(nz)

    for j = 0L, nz-1L do begin

       doit = lyforest(wave,z[j],d_a=d_a,d_b=d_b)
       d_aarr[j] = d_a
       d_barr[j] = d_b

    endfor

    if not keyword_set(plotname) then plotname = 'lyforest_plots'
    ppath = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='plots')

    if keyword_set(postscript) then begin
       ps_open, ppath+plotname, /ps_fonts, /portrait
       device, /inches, /times, xsize=7, ysize=7
    endif else window, 0, xs=450, ys=450
    
    plot, z, d_aarr, xsty=3, ysty=3, line=0, xr=[0,6], yr=[0,1], color=5, $
      xthick=2.0, ythick=2.0, thick=2.0, charsize=1.8, xtit='Redshift', $
      ytit='Fractional Absorption', charthick=2.0, tit='Lyman Forest Absorption'
    oplot, z, d_barr, line=2, color=5, thick=2.5
    legend, [textoidl('D_{A}'),textoidl('D_{B}')], line=[0,2], /left, /top, $
      box=0, charsize=1.8, charthick=2.0
    
    if keyword_set(postscript) then ps_close

return
end
