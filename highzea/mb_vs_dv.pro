pro mb_vs_dv, postscript=postscript
; jm06jul24uofa - 

    path = highzea_path(/papers)

    if keyword_set(postscript) then begin
       postthick = 5.0
       dfpsplot, path+'mb_vs_dv.ps', /color, /square
    endif else begin
       postthick = 2.0
       im_window, 0, xratio=0.5, /square
    endelse

    djs_plot, [0], [0], /nodata, xrange=[-31.0,-14.0], yrange=[5,40000], /ylog, $
      xsty=3, ysty=3, xthick=postthick, ythick=postthick, charsize=1.9, $
      charthick=postthick, xtitle='M_{B} [mag]', ytitle='- \Delta'+'v [km s^{-1}]'
    djs_oplot, !x.crange, [0,0], thick=2.0

    bals = rsex(path+'bals_trump2006_tab4_rsex.txt')

    size = n_elements(bals.m_i)
    fi = fltarr(size)
    fb = fltarr(size)
    mb = fltarr(size)
    dv = float(bals.velocity_average)

    zp = -2.5*alog10(3631.d-23) ; zeropoint for AB magnitudes
    lambda_i = 7481.            ; Angstroms
    lambda_b = 4400.            ; Angstroms
    alpha = -0.5                ; average power-law slope for quasars: f_nu ~ nu ^ alpha

    for i=0L, n_elements(bals.m_i) - 1L do begin
       fi[i] = 10.^((bals[i].m_i+zp)/(-2.5)) 
       fb[i] = fi[i] * (lambda_i / (1. + bals[i].z) / lambda_b)^alpha
       mb[i] = -2.5*alog10(fb[i]) - zp
    endfor

    indx = lindgen(n_elements(bals))
;   indx = where(dv gt -5000.0)
    x = mb[indx] - 5.0*alog10(0.75/0.70) & y = dv[indx]
    im_symbols, 105, psize=1.0, /fill, color='dark green'
    djs_oplot, x, y, ps=8

;   im_hessplot, m_b[indx], dv[indx], xrange=[-24.5,-14.0], yrange=[-5000,200], /xreverse, $
;     xsty=3, ysty=3, xthick=postthick, ythick=postthick, charsize=1.9, $
;     charthick=postthick, xtitle='M_{B} [mag]', ytitle='\Delta'+'v [km s^{-1}]'
;   djs_oplot, !x.crange, [0,0], thick=2.0

; Schwartz & Martin 2004
    
    schwartz04 = rsex(path+'04schwartz.dat')
    schwartz06 = rsex(path+'06schwartz.dat')
    schwartz = struct_append(schwartz04,schwartz06)

    indx = where((schwartz.dv lt 999.0) and (schwartz.mb gt -900.0),nindx)
    x = schwartz[indx].mb - 5.0*alog10(0.75/0.70) & y = -schwartz[indx].dv
    im_symbols, 115, psize=1.9, /fill, color='purple'
    djs_oplot, x, y, ps=8
    
;; Schwartz & Martin 2004
;    
;    schwartz04 = rsex(path+'04schwartz.dat')
;    indx = where((schwartz04.dv lt 999.0) and (schwartz04.mb gt -900.0),nindx)
;    x = schwartz04[indx].mb - 5.0*alog10(0.75/0.70) & y = schwartz04[indx].dv
;    im_symbols, 105, psize=1.5, /fill, color='dark green'
;    djs_oplot, x, y, ps=8
;    
;; Schwartz et al. 2006
;    
;    schwartz06 = rsex(path+'06schwartz.dat')
;    indx = where((schwartz06.dv lt 999.0) and (schwartz06.mb gt -900.0),nindx)
;    x = schwartz06[indx].mb - 5.0*alog10(0.75/0.70) & y = schwartz06[indx].dv
;    im_symbols, 115, psize=1.9, /fill, color='purple'
;    djs_oplot, x, y, ps=8

; Rupke et al. 2005; median (R-K) = 3.25 assume B-R = 1; convert the
; cosmology 

    rupke = rsex(path+'05rupke.dat')
    indx = where((rupke.dv lt 0.0) and (rupke.mr gt -900.0),nindx)
    x = rupke[indx].mr + 1.0 - 5.0*alog10(0.75/0.70) & y = -rupke[indx].dv
    im_symbols, 106, psize=1.2, /fill, color='blue'
    djs_oplot, x, y, ps=8
    
; Pettini et al. 2001    
    
    pettini = rsex(path+'01pettini.dat')
    indx = where((pettini.dv lt 999.0) and (pettini.mb gt -900.0),nindx)
    x = pettini[indx].mb & y = -pettini[indx].dv
    im_symbols, 108, psize=1.2, /fill, color='red'
    djs_oplot, x, y, ps=8
    
; legend

;   im_legend, ['Local starbursts','Local LIRGs/ULIRGs','High-z LBGs'], $
;     psym=[115,106,108], color=djs_icolor(['purple','blue','red']), $
;     fill=[1,1,1], charsize=1.5, charthick=postthick, /right, /bottom, thick=postthick, $
;     symsize=1.3, spacing=1.7
    
    im_legend, ['Local starbursts','LIRGs/ULIRGs <z>=0.1','LBGs <z>=3','BALs 0.5<z<4.5'], $
      psym=[115,106,108,105], color=djs_icolor(['purple','blue','red','dark green']), $
      fill=[1,1,1,1], charsize=1.4, charthick=postthick, /right, /top, thick=postthick, $
      symsize=1.3, spacing=1.7
    
    if keyword_set(postscript) then dfpsclose

stop
    
return
end
    
