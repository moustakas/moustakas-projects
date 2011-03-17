pro turner_combine_repeaters, wfits=wfits, postscript=postscript
; jm05jul08uofa - must be generalized to more than two repeat
;                 observations 
    
    analysis_path = atlas_path(/spec2dturner)
    datapath = atlas_path(/spec2dturner)+'spec1d/repeaters/'
    outpath = atlas_path(/spec2dturner)+'spec1d/'

    readcol, analysis_path+'turner_repeaters_list.txt', pair1, pair2, $
      outfits, galaxy, format='A,A,A,A', comment='#', delimiter=' ', /silent
    npairs = n_elements(pair1)
    
    pairs = transpose([ [pair1], [pair2] ])

; ---------------------------------------------------------------------------       
; combine the repeaters
; ---------------------------------------------------------------------------       

    xpage = 8.5
    ypage = 8.5
    
    pagemaker, nx=1L, ny=1L, xspace=0.0, yspace=0.0, $
      xmargin=[1.2,0.3], ymargin=[0.3,1.2], width=7.0, /normal, $
      height=7.0, position=pos, xpage=xpage, ypage=ypage

    if keyword_set(postscript) then begin
       dfpsplot, outpath+'turner_repeaters.ps', xsize=xpage, ysize=ypage, /color
       postthick = 8.0
       postthick2 = 5.0
    endif else begin
       im_window, 4, xratio=0.4, /square
       postthick = 2.0
       postthick2 = 2.0
    endelse

    xrange = [3500,7200]
    xtitle = 'Wavelength ['+angstrom()+']'
    ytitle = 'Relative Flux (arbitrary units)'
    
    for i = 0L, npairs-1L do begin

       s1 = rd1dspec(pairs[0,i],/silent,datapath=datapath)
       info1 = iforage(pairs[0,i],datapath=datapath)

       s2 = rd1dspec(pairs[1,i],/silent,datapath=datapath)
       info2 = iforage(pairs[1,i],datapath=datapath)
       
; normalize the spectra to the highest flux centered on 5500 Angstroms 
       
       s1spec = s1.spec & s1ivar = 1.0/s1.sigspec^2.0 & s1sky = s1.sky
       s2spec = s2.spec & s2ivar = 1.0/s2.sigspec^2.0 & s2sky = s2.sky

       junk = im_normalize(s1spec,s1.wave,normwave=5500.0,binsize=50.0,const=s1norm)
       junk = im_normalize(s2spec,s2.wave,normwave=5500.0,binsize=50.0,const=s2norm)

       maxnorm = max([s1norm,s2norm],normindx)
       if (normindx eq 0L) then begin
          s2spec = (s1norm/s2norm)*s2spec
          s2ivar = (s2norm/s1norm)^2*s2ivar
          s2sky  = (s1norm/s2norm)*s2sky
       endif else begin
          s1spec = (s2norm/s1norm)*s1spec
          s1ivar = (s1norm/s2norm)^2*s1ivar
          s1sky  = (s2norm/s1norm)*s1sky
       endelse

; interpolate both spectra onto a common wavelength vector to prevent
; extrapolation

       minwave = min(s1.wave)>min(s2.wave) 
       maxwave = max(s1.wave)<max(s2.wave)
       dwave = s1.wave[1]-s1.wave[0] ; assume a common linear dispersion

       finalwave = findgen((maxwave-minwave)/dwave+1L)*dwave+minwave
       newloglam = alog10(finalwave)

       combine1fiber, alog10(s1.wave), s1spec, s1ivar, skyflux=s1sky, $
         newloglam=newloglam, newflux=newflux1, newivar=newivar1, newsky=newsky1
       newsky1 = djs_maskinterp(newsky1,newivar1 eq 0,/const)
       newivar1 = djs_maskinterp(newivar1,newivar1 eq 0,/const)
       
       combine1fiber, alog10(s2.wave), s2spec, s2ivar, skyflux=s2sky, $
         newloglam=newloglam, newflux=newflux2, newivar=newivar2, newsky=newsky2
       newsky2 = djs_maskinterp(newsky2,newivar2 eq 0,/const)
       newivar2 = djs_maskinterp(newivar2,newivar2 eq 0,/const)

; now combine the two spectra with appropriate variance weighting

       inloglam = [ [newloglam], [newloglam] ]
       objflux = [ [newflux1], [newflux2] ]
       objivar = [ [newivar1], [newivar2] ]
       skyflux = [ [newsky1], [newsky2] ]
       
       combine1fiber, inloglam, objflux, objivar, skyflux=skyflux, $
         newloglam=newloglam, newflux=finalflux, newivar=finalivar, $
         newsky=newsky
       newsky = djs_maskinterp(newsky,finalivar eq 0,/const)
       finalmask = finalivar ne 0
       finalivar = djs_maskinterp(finalivar,finalivar eq 0,/const)

       finalsigspec = 1.0/sqrt(finalivar)
       finalmedsnr = median(finalflux/finalsigspec)
       nfinalpix = n_elements(finalflux)

; QA plot       

       yrange = [-0.8,2.8]

       djs_plot, finalwave, im_normalize(finalflux,finalwave,normwave=5500.0)+1, $
         ps=10, position=pos[*,0], thick=postthick2, $
         xthick=postthick, ythick=postthick, charthick=postthick, $
         charsize=1.5, xrange=xrange, yrange=yrange, xsty=3, ysty=3, $
         yminor=2, xtitle=xtitle, ytitle=ytitle, xminor=3, xtickinterval=1000.0
       xyouts, 5500.0, 1.6, 'S/N = '+string(finalmedsnr,format='(I0)'), $
         charsize=1.0, charthick=postthick, align=0.5, /data

       if (normindx eq 0L) then begin
          djs_oplot, s1.wave, s1spec/s1norm, ps=10, color='dark green', thick=postthick2
          xyouts, 5500.0, 0.6, info1.date+', S/N = '+string(info1.medsnr,format='(I0)'), $
            charsize=1.0, charthick=postthick, align=0.5, /data
          djs_oplot, s2.wave, s2spec/s1norm-1, ps=10, color='blue', thick=postthick2
          xyouts, 5500.0, -0.4, info2.date+', S/N = '+string(info2.medsnr,format='(I0)'), $
            charsize=1.0, charthick=postthick, align=0.5, /data
       endif else begin
          djs_oplot, s2.wave, s2spec/s2norm-1, ps=10, color='dark green', thick=postthick2
          xyouts, 5500.0, 0.6, info2.date+', S/N = '+string(info2.medsnr,format='(I0)'), $
            charsize=1.0, charthick=postthick, align=0.5, /data
          djs_oplot, s1.wave, s1spec/s2norm, ps=10, color='blue', thick=postthick2
          xyouts, 5500.0, -0.4, info1.date+', S/N = '+string(info1.medsnr,format='(I0)'), $
            charsize=1.0, charthick=postthick, align=0.5, /data
       endelse

       legend, repstr(galaxy[i],'_',' '), /left, /top, box=0, charsize=1.3, $
         charthick=postthick

       if (not keyword_set(postscript)) then cc = get_kbrd(1)

;      niceprint, wave, finalflux, finalivar, 1.0/sqrt(finalivar), long(finalmask)
;      ploterror, wave, finalflux, finalsigspec, ps=10, xsty=3, ysty=3

; write out

       if keyword_set(wfits) then begin

          finalheader = s1.header

          
          
stop          
          
       endif
          
    endfor

; ---------------------------------------------------------------------------       
; make the plot       
; ---------------------------------------------------------------------------       
       
    ncols = 2L
    nrows = ceil(npairs/float(ncols))

    xmargin = [1.2,0.3]
    ymargin = [0.3,1.2]

    width = replicate(3.5,ncols)
    height = replicate(3.5,nrows)

    xpage = total(width)+total(xmargin)
    ypage = total(height)+total(ymargin)

    pagemaker, nx=ncols, ny=nrows, xspace=0.0, yspace=0.0, $
      xmargin=xmargin, ymargin=ymargin, width=width, /normal, $
      height=height, position=pos, xpage=xpage, ypage=ypage

    if keyword_set(postscript) then begin
       dfpsplot, outpath+'turner_repeaters_paperplot.ps', xsize=xpage, ysize=ypage, /color
       postthick = 8.0 
       postthick2 = 5.0 
    endif else begin
       postthick = 2.0
       postthick2 = 2.0 
       im_window, 0, xratio=0.5, yratio=0.9
    endelse
    
    binsize = 200.0

    xrange = [3500,7200]
    yrange = [-0.6,3.0]

    for i = 0L, npairs-1L do begin

       s1 = rd1dspec(pairs[0,i],/silent,datapath=datapath)
       info1 = iforage(pairs[0,i],datapath=datapath)

       s2 = rd1dspec(pairs[1,i],/silent,datapath=datapath)
       info2 = iforage(pairs[1,i],datapath=datapath)
       
       wave = s1.wave
       flux1 = im_normalize(s1.spec,wave,normwave=5500.0)

       flux2 = interpol(s2.spec,s2.wave,wave) ; simple interpolation for the plot 
       flux2 = im_normalize(flux2,wave,normwave=5500.0)

; bin the spectra and then take the ratio

       bin1spec = im_binspec(flux1,wave,binsize=binsize,binwave=binwave,binspec_err=bin1spec_err)
       bin2spec = im_binspec(flux2,wave,binsize=binsize,binwave=binwave,binspec_err=bin2spec_err)
       nbins = n_elements(binwave)
       
       ratio = bin1spec/bin2spec-1
       ratioerr = im_compute_error(bin1spec,bin1spec_err,bin2spec,bin2spec_err,/quotient)
;      ratio = flux1/flux2 - 1.0
       junk = im_stats(100*ratio,/verbose,no_head=(i ne 0L))

;      yrange[1] = (max(flux1)>max(flux2))*1.05

       if odd(i) then begin
          ytickname = replicate(' ',10)
          ytitle = ''
       endif else begin
          delvarx, ytickname
          ytitle = 'Relative Flux (arbitrary units)'
       endelse

       if (i lt npairs-2L) then begin
          xtickname = replicate(' ',10)
          xtitle = ''
       endif else begin
          delvarx, xtickname
          xtitle = 'Wavelength ['+angstrom()+']'
       endelse
       
       djs_plot, wave, flux1, ps=10, position=pos[*,i], noerase=(i ne 0L), $
         xthick=postthick, ythick=postthick, charthick=postthick, charsize=1.5, ytickname=ytickname, $
         xtickname=xtickname, xrange=xrange, yrange=yrange, xsty=3, ysty=3, $
         yminor=2, xtitle=xtitle, ytitle=ytitle, xminor=3, xtickinterval=1000.0
       djs_oplot, wave, flux2+1, color='orange', ps=10
       djs_oplot, !x.crange, [0,0], color='grey', line=0, thick=postthick
       oploterror, binwave, ratio, replicate(binsize/2.0,nbins), ratioerr, $
         color=djs_icolor('blue'), errcolor=djs_icolor('blue'), $
         ps=3, thick=3.0, errthick=3.0
       legend, repstr(galaxy[i],'_',' '), /left, /top, box=0, charsize=1.3, charthick=postthick

       icleanup, s1 & icleanup, s2
       
    endfor

    if keyword_set(postscript) then dfpsclose

stop    

return
end
    
