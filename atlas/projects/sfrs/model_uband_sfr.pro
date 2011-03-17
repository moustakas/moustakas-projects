;+
; NAME:
;       MODEL_UBAND_SFR
;
; PURPOSE:
;       Time-evolution of the U-band as a quantitative SFR indicator. 
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; PROCEDURES USED:
;
;
; COMMENTS:
;
;
; EXAMPLES:
;
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Mar 15, U of A
;-

pro model_uband_sfr, postscript=postscript

    datapath = atlas_path(/projects)+'sfrs/uband/'
    paperpath = atlas_path(/papers)+'sfrs/FIG_SFRS/'

    Uinfo = im_filterspecs(filterlist='bessell_U.par')

    haconst = 7.9D-42            ; K98 conversion L(Ha) --> SFR(Ha)
    loghaconst = alog10(haconst)
    hasfrconstoffset = 43.0
    Uhasfrconstoffset = 43.0

;   d4000range = [0.8,2.15]
;   sfrUharange = [41.7,45.3]

    @'xyrange_sfrs'

; read the output from MEASURE_UBAND_MODELS

    result = mrdfits(datapath+'measure_uband_models.fits',1,/silent)
    ncont = n_elements(result)
    nburst = n_elements(result[0].LU_burst)

    result_constant = mrdfits(datapath+'measure_uband_models.fits',2,/silent)
    nage = n_elements(result_constant)

; generate the plot    

    if keyword_set(postscript) then begin
       postthick = 8.0
       psname = 'd4000_vs_lu_lha_models.eps'
       dfpsplot, datapath+psname, /square, /encapsulated
;      dfpsplot, datapath+'D4000_vs_LU_LHa_models.ps', xsize=8.5, ysize=4.9
    endif else begin
       postthick = 2.0
       im_window, 0, xratio=0.6, /square
    endelse

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal
;   pagemaker, nx=1, ny=1, yspace=0, xspace=0, xmargin=[2.0,1.0], $
;     ymargin=[0.3,1.1], width=5.5, height=3.5, ypage=4.9, $
;     position=pos, /normal

    xrange = d4000range2
;   yrange = [0.1,3.9]
    yrange = sfrUharange - Uhasfrconstoffset
;   yrange = sfrUharange - Uhasfrconstoffset

    xtitle = 'D_{n}(4000)'
;   ytitle = 'log [L(U)/L(H\alpha)] '
    ytitle = 'log [10^{-'+string(hasfrconstoffset,format='(I0)')+'} L(U)/\psi(H\alpha)] ['+sfr_units()+'/erg s^{-1}]'

; plot the continuous models    
    
    plotsym, 0, 1.3, /fill
    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
      xsty=3, ysty=3, charsize=2.0, charthick=postthick, position=pos, $
      xthick=postthick, ythick=postthick, xtitle=xtitle, ytitle=ytitle

    lumratio = alog10(result.LU_cont/result.LHa_cont) - loghaconst - hasfrconstoffset
;   lumratio = alog10(result.LU_cont/result.LHa_cont)
    djs_oplot, result.D4000_cont, lumratio, ps=-8, thick=postthick

    xoffset = +[0.025,0.025,0.025,0.025]*(!x.crange[1]-!x.crange[0])
    yoffset = -[0.03,0.03,0.03,0.03]*(!y.crange[1]-!y.crange[0])
;   xfactor = [1.01,1.01,1.01,1.01]
;   yfactor = [0.92,0.90,0.88,0.88]
    
    for icont = 0L, ncont-1L do begin
       xyouts, result[icont].D4000_cont + xoffset[icont], lumratio[icont] + yoffset[icont], $
;      xyouts, result[icont].D4000_cont*xfactor[icont], yfactor[icont]*lumratio[icont], $
         strtrim(string(result[icont].tau,format='(I2)'),2), charsize=1.4, $
         charthick=postthick, align=0.5, /data
    endfor

    xyouts, 1.75, -0.3, textoidl('\tau (Gyr)'), charsize=1.4, $
      charthick=postthick, align=0.5
;   xyouts, 1.85, 1.3, textoidl('\psi\propto')+'exp(-13/'+textoidl('\tau)'), charsize=1.5, $
;   xyouts, 1.85, 1.3, textoidl('\psi\propto'+'e^{-13/\tau}'), charsize=1.6, $
    xyouts, 1.75, -0.5, textoidl('\psi\propto'+'e^{-t/\tau}'), charsize=1.6, $
      charthick=postthick, align=0.5
    
; now connect the models having fixed tau values but different "time
; since burst"

    xoffset = +[0.025,0.025,0.025,0.025]*(!x.crange[1]-!x.crange[0])
    yoffset = +[0.03,0.03,0.03,0.03]*(!y.crange[1]-!y.crange[0])
;   xfactor = [1.01,0.99,1.03]
;   yfactor = [1.05,1.07,1.02]
;   xfactor = [1.01,0.99,1.00,1.03]
;   yfactor = [1.05,1.07,1.05,1.02]
    
    plotsym, 8, 1.3, /fill
    for icont = 0L, ncont-1L do begin

       lumratio = alog10(result[icont].LU_burst/result[icont].LHa_burst) - loghaconst - hasfrconstoffset
       djs_oplot, result[icont].D4000_burst, lumratio, ps=-8, line=1, thick=postthick       

       if (icont eq 0L) then begin
          for iburst = 0L, nburst-1L do begin
             if (iburst eq 0L) then prefix = textoidl('t_{burst}=') else prefix = ''
             xyouts, result[icont].D4000_burst[iburst] + xoffset[iburst], lumratio[iburst] + yoffset[iburst], $
;            xyouts, result[icont].D4000_burst[iburst]*xfactor[iburst], yfactor[iburst]*lumratio[iburst], $
               prefix+'-'+strtrim(string(result[icont].tburst[iburst],format='(I1)'),2)+' Gyr', charsize=1.3, $
               charthick=postthick, align=0.5, /data
          endfor
       endif

    endfor

; now overplot the "continuous star formation at different ages"
; models 
    
    plotsym, 4, 1.5, /fill
    lumratio = alog10(result_constant.LU/result_constant.LHa) - loghaconst - hasfrconstoffset
    djs_oplot, result_constant.D4000, lumratio, ps=-8, thick=postthick, line=2

    xoffset = -[0.025,0.025,0.025,0.025]*(!x.crange[1]-!x.crange[0])
    yoffset = +[0.02,0.02,0.02,0.02]*(!y.crange[1]-!y.crange[0])
;   xfactor = [0.95,0.95,0.97,0.97]
;   yfactor = [1.05,1.05,1.05,1.05]
    age_string = ['0.01','0.1','0.5','10']
    
    for iage = 0L, nage-1L do begin
       xyouts, result_constant[iage].D4000 + xoffset[iage], lumratio[iage] + yoffset[iage], $
;      xyouts, result_constant[iage].D4000*xfactor[iage], yfactor[iage]*lumratio[iage], $
         age_string[iage], charsize=1.3, charthick=postthick, align=0.5, /data
    endfor

;   label = textoidl('\psi=constant')
    label = textoidl('\psi=\psi_{0}')
    xyouts, 1.05, 0.1, label, charsize=1.6, charthick=postthick, align=0.5
    xyouts, 1.05, -0.15, 'Age (Gyr)', charsize=1.4, charthick=postthick, align=0.5
;   xyouts, 1.2, 1.0, textoidl('\psi=\psi_{0}'), charsize=1.8, $
;     charthick=postthick, align=0.5

; overlay the dust extinction vector on this plot

    dust_fraction = 0.5 ; the continuum is reddened by a fractional amount of dust
    AHa = 1.0
    ebv_gas = AHa / k_lambda(6563,/odonnell)

    UHa_true =  1.5
    d4000_pos = 1.0
;   UHa_true = -1.2
;   d4000_pos = 1.82

    UHa_red = UHa_true + 0.4*ebv_gas*(dust_fraction*k_lambda(Uinfo.weff,/charlot)-k_lambda(6563,/charlot))

    arrow, d4000_pos, UHa_true, d4000_pos, UHa_red, /data, thick=postthick, hthick=postthick, hsize=-0.6
    xyouts, d4000_pos, UHa_true+0.15, textoidl('A(H\alpha) = 1'), $
      charsize=1.3, charthick=postthick, align=0.5
    
    if keyword_set(postscript) then begin
       dfpsclose
       spawn, ['/bin/cp -f '+datapath+psname+' '+paperpath], /sh
    endif
    
stop
       
return
end
