;+
; NAME:
;       HIGHZEA_PLOT_REPEATERS
;
; PURPOSE:
;       Plot all the repeat observations to determine which objects
;       should be combined, excluded, etc.
;
; CALLING SEQUENCE:
;       highzea_plot_repeaters, /postscript
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;       After examining the postscript output from this routine,
;       generate an HIGHZEA_COMBINE_REPEATERS.TXT file and run
;       HIGHZEA_COMBINE_REPEATERS.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jul 21, U of A
;-

pro highzea_plot_repeaters, postscript=postscript
    
    analysis_path = highzea_path(/analysis)
    repeatpath = highzea_path(/spec1d)+'repeaters/'

; initialize some plotting variables

    xpage = 8.5
    ypage = 8.5
    
    pagemaker, nx=1L, ny=1L, xspace=0.0, yspace=0.0, $
      xmargin=[1.2,0.3], ymargin=[0.3,1.2], width=7.0, /normal, $
      height=7.0, position=pos, xpage=xpage, ypage=ypage

    if keyword_set(postscript) then begin
       dfpsplot, analysis_path+'highzea_plot_repeaters.ps', xsize=xpage, ysize=ypage, /color
       postthick = 5.0
       postthick2 = 2.0
    endif else begin
       im_window, 4, xratio=0.6, /square
       postthick = 2.0
       postthick2 = 1.0
    endelse

    xrange = [4000,7200]
    xtitle = 'Wavelength [\AA]'
    ytitle = 'Relative Flux [arbitrary units]'

    colorlist = ['green','red','blue','orange']
    
; read all the repeaters

    allfiles = file_search(repeatpath+'*.ms.fits',count=nall)
    allfiles = file_basename(allfiles)

    allinfo = iforage(repeatpath+allfiles)
    allgalaxy = strtrim(allinfo.galaxy,2)
;   uindx = cmset_op(uniq(allgalaxy),'OR',uniq(allinfo.aperwid)) ; 
    uindx = uniq(allgalaxy)
    galaxy = allgalaxy[uindx]
    info = allinfo[uindx]
    ngalaxy = n_elements(galaxy)

    for igalaxy = 0L, ngalaxy-1L do begin

; match on galaxy name *and* extraction aperture (this trick won't
; work for the spectral atlas)
       
;      match = where(galaxy[igalaxy] eq allgalaxy,nrepeat)
       match = where((galaxy[igalaxy] eq allgalaxy) and (info[igalaxy].aperwid eq allinfo.aperwid),nrepeat)
       if (nrepeat gt 4L) then begin
          splog, 'Number of spectra not supported!'
          return
       endif

       repeatlist = allfiles[match]
       print, repeatlist

       s1 = rd1dspec(repeatlist[0],/silent,datapath=repeatpath)
       s2 = rd1dspec(repeatlist[1],/silent,datapath=repeatpath)
       if (nrepeat ge 3L) then s3 = rd1dspec(repeatlist[2],/silent,datapath=repeatpath)
       if (nrepeat ge 4L) then s4 = rd1dspec(repeatlist[3],/silent,datapath=repeatpath)

; interpolate spectra onto a common wavelength vector to prevent
; extrapolation 

       dwave = s1.wave[1]-s1.wave[0] ; assume a common linear dispersion

       minwave = min(s1.wave)>min(s2.wave) 
       maxwave = max(s1.wave)<max(s2.wave)
       if (nrepeat ge 3L) then begin
          minwave = minwave>min(s3.wave)
          maxwave = maxwave<max(s3.wave)
       endif
       if (nrepeat ge 4L) then begin
          minwave = minwave>min(s4.wave)
          maxwave = maxwave<max(s4.wave)
       endif

       finalwave = findgen((maxwave-minwave)/dwave+1L)*dwave+minwave
       newloglam = alog10(finalwave)
       nfinalpix = n_elements(finalwave)

       dlogwave = newloglam[1]-newloglam[0]
;      binsz = min((newloglam-shift(newloglam,1))[1L:nfinalpix-1L])
       binsz = dlogwave

; interpolate the spectra       

       newflux1 = interpol(s1.spec,s1.wave,finalwave)
       newivar1 = interpol(1.0/s1.sigspec^2.0,s1.wave,finalwave)
       newsky1 = interpol(s1.sky,s1.wave,finalwave)
       
;      combine1fiber, alog10(s1.wave), s1.spec, 1.0/s1.sigspec^2, skyflux=s1.sky, $
;        newloglam=newloglam, newflux=newflux1, newivar=newivar1, newsky=newsky1, $
;        bkptbin=1.0*binsz, maxsep=2.0*binsz, binsz=binsz
;      djs_plot, s1.wave, s1.spec, ps=10
;      djs_oplot, 10^newloglam, newflux1, ps=10, color='red'
;      newsky1 = djs_maskinterp(newsky1,newivar1 eq 0,/const)
;      newivar1 = djs_maskinterp(newivar1,newivar1 eq 0,/const)
       
       newflux2 = interpol(s2.spec,s2.wave,finalwave)
       newivar2 = interpol(1.0/s2.sigspec^2.0,s2.wave,finalwave)
       newsky2 = interpol(s2.sky,s2.wave,finalwave)
       
;      combine1fiber, alog10(s2.wave), s2.spec, 1.0/s2.sigspec^2, skyflux=s2.sky, $
;        newloglam=newloglam, newflux=newflux2, newivar=newivar2, newsky=newsky2, $
;        bkptbin=1.0*binsz, maxsep=2.0*binsz, binsz=binsz
;      newsky2 = djs_maskinterp(newsky2,newivar2 eq 0,/const)
;      newivar2 = djs_maskinterp(newivar2,newivar2 eq 0,/const)

       if (nrepeat ge 3L) then begin
          newflux3 = interpol(s3.spec,s3.wave,finalwave)
          newivar3 = interpol(1.0/s3.sigspec^2.0,s3.wave,finalwave)
          newsky3 = interpol(s3.sky,s3.wave,finalwave)       

;         combine1fiber, alog10(s3.wave), s3.spec, 1.0/s3.sigspec^2, skyflux=s3.sky, $
;           newloglam=newloglam, newflux=newflux3, newivar=newivar3, newsky=newsky3, $
;           bkptbin=1.0*binsz, maxsep=2.0*binsz, binsz=binsz
;         newsky3 = djs_maskinterp(newsky3,newivar3 eq 0,/const)
;         newivar3 = djs_maskinterp(newivar3,newivar3 eq 0,/const)
       endif

       if (nrepeat ge 4L) then begin
          newflux4 = interpol(s4.spec,s4.wave,finalwave)
          newivar4 = interpol(1.0/s4.sigspec^2.0,s4.wave,finalwave)
          newsky4 = interpol(s4.sky,s4.wave,finalwave)       

;         combine1fiber, alog10(s4.wave), s4.spec, 1.0/s4.sigspec^2, skyflux=s4.sky, $
;           newloglam=newloglam, newflux=newflux4, newivar=newivar4, newsky=newsky4, $
;           bkptbin=1.0*binsz, maxsep=2.0*binsz, binsz=binsz
;         newsky4 = djs_maskinterp(newsky4,newivar4 eq 0,/const)
;         newivar4 = djs_maskinterp(newivar4,newivar4 eq 0,/const)
       endif

; compute the normalization constants       
       
       junk = im_normalize(newflux1,finalwave,normwave=5500.0,binsize=50.0,const=s1norm)
       junk = im_normalize(newflux2,finalwave,normwave=5500.0,binsize=50.0,const=s2norm)

       if (nrepeat ge 3L) then junk = im_normalize(newflux3,finalwave,normwave=5500.0,$
         binsize=50.0,const=s3norm) else s3norm = -100.0
       if (nrepeat ge 4L) then junk = im_normalize(newflux4,finalwave,normwave=5500.0,$
         binsize=50.0,const=s4norm) else s4norm = -100.0

       maxnorm = max([s1norm,s2norm,s3norm,s4norm],normindx)
;      print, repeatlist[0], s1norm, s2norm, s3norm, s4norm

       objflux = [ [newflux1*(maxnorm/s1norm)], [newflux2*(maxnorm/s2norm)] ]
       objivar = [ [newivar1*(s1norm/maxnorm)^2], [newivar2*(s2norm/maxnorm)^2] ]
       skyflux = [ [newsky1*(maxnorm/s1norm)], [newsky2*(maxnorm/s1norm)] ]

       if (nrepeat ge 3L) then begin
          objflux = [ [objflux], [newflux3*(maxnorm/s3norm)] ]
          objivar = [ [objivar], [newivar3*(s3norm/maxnorm)^2] ]
          skyflux = [ [skyflux], [newsky3*(maxnorm/s3norm)] ]
       endif
          
       if (nrepeat ge 4L) then begin
          objflux = [ [objflux], [newflux4*(maxnorm/s4norm)] ]
          objivar = [ [objivar], [newivar4*(s4norm/maxnorm)^2] ]
          skyflux = [ [skyflux], [newsky4*(maxnorm/s4norm)] ]
       endif
          
; now combine all the spectra with appropriate variance weighting

       finalflux = total(objivar*objflux,2)/total(objivar,2)
       finalivar = total(objivar,2)
       
;      inloglam = rebin(newloglam,nfinalpix,nrepeat)
;      combine1fiber, inloglam, objflux, objivar, skyflux=skyflux, $
;        newloglam=newloglam, newflux=finalflux, newivar=finalivar, $
;        newsky=newsky, bkptbin=1.2*binsz, maxsep=2.0*binsz, binsz=binsz
;      newsky = djs_maskinterp(newsky,finalivar eq 0,/const)
;      finalmask = finalivar ne 0
;      finalivar = djs_maskinterp(finalivar,finalivar eq 0,/const)

       finalsigspec = 1.0/sqrt(finalivar)
       finalmedsnr = median(finalflux/finalsigspec)
       nfinalpix = n_elements(finalflux)

; make the plot

       yrange = fltarr(2)
;      yrange[1] = max(finalflux/maxnorm)+1
       stats = im_stats(finalflux/maxnorm,sigrej=5.0)
       yrange[1] = stats.maxrej+1.2

       yrange[0] = -0.9
       if (nrepeat ge 3L) then yrange[0] = yrange[0] - (nrepeat-2.0)
       
       djs_plot, finalwave, finalflux/maxnorm+1, $
         ps=10, position=pos[*,0], thick=postthick2, $
         xthick=postthick, ythick=postthick, charthick=postthick, $
         charsize=1.5, xrange=xrange, yrange=yrange, xsty=3, ysty=3, $
         yminor=2, xtitle=xtitle, ytitle=ytitle, xminor=3, xtickinterval=1000.0
       xyouts, 5500.0, 1.6, 'S/N = '+string(finalmedsnr,format='(I0)'), $
         charsize=1.2, charthick=postthick, align=0.5, /data

       offset = 0.0
       for irepeat = 0L, nrepeat-1L do begin

          legstr = allinfo[match[irepeat]].date+', S/N = '+string(allinfo[match[irepeat]].medsnr,format='(I0)')

          if (allinfo[match[irepeat]].scanlen eq 0.0) then begin
             legstr = legstr+', AM = '+string(allinfo[match[irepeat]].airmass,format='(F5.3)')
          endif else begin
             legstr = legstr+', Aperture = '+string(allinfo[match[irepeat]].aperwid,format='(I0)')+' x '+$
               string(allinfo[match[irepeat]].scanlen,format='(I0)')+', PA = '+string(allinfo[match[irepeat]].posangle,$
               format='(I0)')
          endelse
          
          ratio = objflux[*,irepeat]/finalflux
          yesplot = where(abs(ratio) lt 1.5)

          djs_oplot, finalwave[yesplot], ratio[yesplot]-offset, $
;         djs_oplot, finalwave, objflux[*,irepeat]/maxnorm-offset, $
            color=colorlist[irepeat], ps=10, thick=postthick2
          djs_oplot, minmax(finalwave), [1.1,1.1]-offset, line=2, thick=postthick2
          djs_oplot, minmax(finalwave), [1.0,1.0]-offset, line=0, thick=postthick2
          djs_oplot, minmax(finalwave), [0.9,0.9]-offset, line=2, thick=postthick2
          xyouts, 5500.0, 1.0-offset-0.4, legstr, charsize=1.2, charthick=postthick, $
            align=0.5, /data

          offset = offset + 1.0

       endfor

       legstr = galaxy[igalaxy]
       
       legend, legstr, /left, /top, box=0, charsize=1.3, $
         charthick=postthick

       if (not keyword_set(postscript)) then cc = get_kbrd(1)

;      niceprint, wave, finalflux, finalivar, 1.0/sqrt(finalivar), long(finalmask)
;      ploterror, wave, finalflux, finalsigspec, ps=10, xsty=3, ysty=3

    endfor

    if keyword_set(postscript) then begin
       dfpsclose
       spawn, ['gzip -f '+analysis_path+'highzea_plot_repeaters.ps'], /sh
    endif
    
stop    

return
end
    
