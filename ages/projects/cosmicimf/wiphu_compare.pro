function wiphu_zbins, nzbins
    nzbins = 4
    zbins = replicate({zbin: 0.0, zlo: 0.0, zup: 0.0},nzbins)
    zbins.zlo = [0.05,0.20,0.35,0.50]
    zbins.zup = [0.20,0.35,0.50,0.65]
    zbins.zbin = total([[zbins.zlo],[zbins.zup]],2)/2.0
return, zbins
end

pro wiphu_compare
; jm10mar18ucsd - compare my 24-micron LF with Wiphu

    cosmicimfpath = ages_path(/projects)+'cosmicimf/'

    area = ages_survey_area()*!dtor^2
    h100 = 0.7 ; 0.45

    l24axis = im_array(6.0,15.0,0.02)
;   l24axis = im_array(7.5,12.0,0.02)
    
; read the sample, the SFRs, and the Pegase model parameters 
    allsample = read_cosmicimf_sample()
    allsfrs = read_cosmicimf_sample(/sfrs)
    
; read the redshift bins and the binsize
    binsize = 0.15
    zbins = wiphu_zbins(nzbins)

; from Wiphu's paper:     
    alpha = 0.37
    beta = 2.36
    wphistar = 1.2D-3
    wlstar = alog10([7.25,10.6,15.5,23.6]*1D9)
    wlstar_local = alog10(4.27D9)

    print, 10^(-1.6+alog10(0.66)) ; low-z rho********

; SFR conversion factor from Rieke+09
    ff = l24axis*0.0
    lo = where(10^l24axis lt 1.3d10,comp=hi)
    ff[lo] = 7.8d-10*10^l24axis[lo]
    ff[hi] = 7.8D-10*10^l24axis[hi]*(7.76D-11*10^l24axis[hi])^0.048
;   plot, l24axis, ff, /ylog

;    for ii = 0, n_elements(wlstar)-1 do begin
;       phi = lf_double_powerlaw(10^l24axis,wphistar,wlstar[ii],alpha,beta)
;       print, alog10(total(alog(10.0)*0.02*ff*phi))-alog10(0.66) ; THIS IS IT!
;    endfor
;
;stop    
;    
;    print, alog10(total(alog(10.0)*0.02*10^l24axis*phi))
;
;    print, alog10(im_integral(l24axis,ff*phi,8.0,alog10(3D10)) + $
;      im_integral(l24axis,ff*phi,alog10(3D10),15.0)) + $
;      alog10(0.66)
;
;stop    
    
; make the plot    
;   psfile = cosmicimfpath+'wiphu_compare.ps'
    im_plotconfig, 5, pos, psfile=psfile, charsize=1.6, $
      height=3.0*[1,1]
    for jj = 0, nzbins-1 do begin
       these = where((allsample.z gt zbins[jj].zlo) and $
         (allsample.z lt zbins[jj].zup) and (allsfrs.agn eq 0) and $
         (allsfrs.zmax_24 gt allsample.z) and $
         (allsample.zmax_evol gt allsample.z),ngal) ; wacky!
       sample = allsample[these]
       sfrs = allsfrs[these]

; compute Vmax       
       zmin = replicate(zbins[jj].zlo,ngal)
;      zmin = sample.zmin_evol>zbins[jj].zlo
       zmax = (sample.zmax_evol<sfrs.zmax_24)<zbins[jj].zup
       vmax = (area/3.0)*(lf_comvol(zmax)-lf_comvol(zmin))/h100^3.0 ; h=0.7
       vol = (area/3.0)*((lf_comvol(zbins[jj].zup)-$
         lf_comvol(zbins[jj].zlo)))[0]/h100^3.0 ; h=0.7

; build the Vmax-weighted 24-micron LF
       im_lf_vmax, sfrs.l24, sample.final_weight/vmax/binsize, $
         binsize=binsize, histmin=histmin, histmax=histmax, $
         binlum=binlum, phi=phi, errphi=phierr, $
         minlum=9.0, fullbin=fullbin, number=number, debug=0

; overplot Wiphu's results
       if odd(jj) then begin
          ytitle = ''
          ytickname = replicate(' ',10)
       endif else begin
          ytitle = textoidl('\Phi(L) (h_{70}^3 Mpc^{-3} dex^{-1})') 
          delvarx, ytickname 
       endelse
       if (jj lt 2) then begin
          xtitle = ''
          xtickname = replicate(' ',10)
       endif else begin
          xtitle = 'log_{10} (L_{24}/L_{\odot})'
          delvarx, xtickname 
       endelse
       djs_plot, [0], [0], /nodata, noerase=(jj gt 0), $
         position=pos[*,jj], xsty=1, ysty=1, /ylog, $
         yrange=[1E-6,2E-2], xrange=[8.5,12.1], xtickname=xtickname, $
         ytickname=ytickname, xtitle=xtitle, ytitle=ytitle
       oploterror, binlum, phi, phierr, psym=symcat(16), symsize=1.5
; overplot Wiphu's data
       if (jj eq 0) then begin
          wlum = im_array(8.1,10.7,0.2)
          wphi = [-2.44,-2.44,-2.49,-2.34,-2.47,-2.45,$
            -2.58,-2.69,-2.92,-3.08,-3.35,-3.72,-4.15,-4.33]
       endif
       if (jj eq 1) then begin
          wlum = im_array(9.1,11.3,0.2)
          wphi = [-2.88,-2.68,-2.79,-2.89,-2.97,-3.14,$
            -3.42,-3.70,-4.04,-4.40,-4.51,-5.21]
       endif
       if (jj eq 2) then begin
          wlum = im_array(9.7,11.3,0.2)
          wphi = [ -3.03,-2.88,-3.09,-3.15,-3.45,-3.81,$
            -4.15,-4.64,-4.80]
       endif
       if (jj eq 3) then begin
          wlum = im_array(10.3,11.9,0.2)
          wphi = [-3.60,-3.49,-3.59,-3.94,-4.23,-4.48,-4.97,-4.94,-5.53]
       endif
       djs_oplot, wlum, 10^wphi, psym=symcat(6,thick=5), color='orange'
       djs_oplot, l24axis, lf_double_powerlaw(l24axis,$
         wphistar,wlstar[jj],alpha,beta), color='blue'
       djs_oplot, l24axis, lf_double_powerlaw(l24axis,$
         wphistar,wlstar_local,alpha,beta), line=5, color='grey'
       legend, 'z='+strtrim(string(zbins[jj].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[jj].zup,format='(F12.2)'),2), $
         /right, /top, box=0, margin=0, charsize=1.5
;      cc = get_kbrd(1)
    endfor
    im_plotconfig, psfile=psfile, /psclose, /gzip

return
end
    
