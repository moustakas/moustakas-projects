pro plot_burstmodels, postscript=postscript
; jm05oct27uofa
    
    datapath = '/home/ioannis/jobs/burstmodels/'
    paperpath = datapath

    filterlist = ['dwfs_Bw.par','bessell_R.par']
    filtinfo = im_filterspecs(filterlist=filterlist)

    d4000range = [0.8,1.9]
    bwrrange = [-2.3,0.8]

; read the output from MEASURE_BURSTMODELS

    result = mrdfits(datapath+'measure_burstmodels.fits',1,/silent)
    ncont = n_elements(result)
    nburst = n_elements(result[0].tburst)

    result_constant = mrdfits(datapath+'measure_burstmodels.fits',2,/silent)
    nage = n_elements(result_constant)

; generate the plot    

    if keyword_set(postscript) then begin
       postthick = 8.0
       psname = 'burstmodels.eps'
       dfpsplot, datapath+psname, /square, /encapsulated
    endif else begin
       postthick = 2.0
       im_window, 0, xratio=0.6, /square
    endelse

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    xrange = d4000range
    yrange = bwrrange

    xtitle = 'D_{n}(4000)'
    ytitle = 'B_{W} - R'

; plot the continuous models    
    
    plotsym, 0, 1.3, /fill
    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
      xsty=3, ysty=3, charsize=2.0, charthick=postthick, position=pos, $
      xthick=postthick, ythick=postthick, xtitle=xtitle, ytitle=ytitle

    djs_oplot, result.D4000_cont, result.bwr_rest_cont, ps=-8, thick=postthick

;   xoffset = +[0.025,0.025,0.025,0.025]*(!x.crange[1]-!x.crange[0])
;   yoffset = -[0.03,0.03,0.03,0.03]*(!y.crange[1]-!y.crange[0])
;   for icont = 0L, ncont-1L do begin
;      xyouts, result[icont].D4000_cont + xoffset[icont], lumratio[icont] + yoffset[icont], $
;        strtrim(string(result[icont].tau,format='(I2)'),2), charsize=1.4, $
;        charthick=postthick, align=0.5, /data
;   endfor

;   xyouts, 1.75, -0.3, textoidl('\tau (Gyr)'), charsize=1.4, $
;     charthick=postthick, align=0.5
;   xyouts, 1.75, -0.5, textoidl('\psi\propto'+'e^{-t/\tau}'), charsize=1.6, $
;     charthick=postthick, align=0.5
    
; now connect the models having fixed tau values but different "time
; since burst"

    xoffset = +[0.025,0.025,0.025,0.025]*(!x.crange[1]-!x.crange[0])
    yoffset = +[0.03,0.03,0.03,0.03]*(!y.crange[1]-!y.crange[0])
    
    plotsym, 8, 1.3, /fill
    for icont = 0L, ncont-1L do begin

       djs_oplot, result[icont].D4000_burst, result[icont].bwr_rest_burst, $
         ps=-8, line=1, thick=postthick       

;      if (icont eq 0L) then begin
;         for iburst = 0L, nburst-1L do begin
;            if (iburst eq 0L) then prefix = textoidl('t_{burst}=') else prefix = ''
;            xyouts, result[icont].D4000_burst[iburst] + xoffset[iburst], lumratio[iburst] + yoffset[iburst], $
;              prefix+'-'+strtrim(string(result[icont].tburst[iburst],format='(I1)'),2)+' Gyr', charsize=1.3, $
;              charthick=postthick, align=0.5, /data
;         endfor
;      endif

    endfor

; now overplot the "continuous star formation at different ages"
; models 
    
    plotsym, 4, 1.5, /fill
    djs_oplot, result_constant.D4000, result_constant.bwr_rest, ps=-8, thick=postthick, line=2

;   xoffset = -[0.025,0.025,0.025,0.025]*(!x.crange[1]-!x.crange[0])
;   yoffset = +[0.02,0.02,0.02,0.02]*(!y.crange[1]-!y.crange[0])
;   age_string = ['0.01','0.1','0.5','10']
;   
;   for iage = 0L, nage-1L do begin
;      xyouts, result_constant[iage].D4000 + xoffset[iage], lumratio[iage] + yoffset[iage], $
;        age_string[iage], charsize=1.3, charthick=postthick, align=0.5, /data
;   endfor

;   label = textoidl('\psi=\psi_{0}')
;   xyouts, 1.05, 0.1, label, charsize=1.6, charthick=postthick, align=0.5
;   xyouts, 1.05, -0.15, 'Age (Gyr)', charsize=1.4, charthick=postthick, align=0.5

stop    
    
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
