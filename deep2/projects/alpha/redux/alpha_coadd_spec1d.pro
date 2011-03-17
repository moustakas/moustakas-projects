pro alpha_coadd_spec1d, mike, info, side=side, datapath=datapath, $
  qafile=qafile, fluxed=fluxed, std=std
; jm09jan06nyu - coadd the 1D spectra, both individual orders and
; repeat observations of the same object; if /FLUXED then coadded the
; fluxed 1D spectra, otherwise coadd the counts

    nobj = n_elements(info)
    if (nobj ne n_elements(mike)) then begin
       splog, 'Dimensions of INFO and MIKE do not match'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = './'

    light = 299792.458D ; [km/s]
    velpix = (side eq 1 ? 1.50d : 2.10d) * double(mike[0].rowbin) ; [km/s]
    cdelt = alog10(1.0d + velpix/light)                           ; pixel size

; loop on each object       

    if (n_elements(qafile) ne 0L) then begin
       splog, 'Building QA-plot '+qafile
       dfpsplot, qafile, /landscape, /color
       im_plotfaves, thick=5.0, /post
    endif

    if keyword_set(fluxed) then begin
       if keyword_set(std) then begin
          scale = 1D10
          ytitle = textoidl('Flux (10^{-10} '+flam_units()+')')
       endif else begin
          scale = 1D16
          ytitle = textoidl('Flux (10^{-16} '+flam_units()+')')
       endelse
       prefix = 'f.'
    endif else begin
       scale = 1.0
       ytitle = textoidl('Flux (ADU)')
       prefix = ''
    endelse

    for jj = 0L, nobj-1L do begin

; read the object structure, which contains the 1D spectra 
       objfil = mike_getfil('obj_fil',setup,$
         SUBFIL=mike[jj].img_root,/name)
       splog, 'Reading '+objfil
       objstr = xmrdfits(objfil,1,STRUCTYP='mikeobjstrct',/silent)

; read "my" object structure to identify the orders we care about 
       if keyword_set(std) then begin
          these_ordrs = objstr.order
       endif else begin
          myobjfil = repstr(objfil,'Obj_','myObj_')+'.gz'
          splog, 'Reading '+myobjfil
          myobjstr = mrdfits(myobjfil,1,/silent)
          these_ordrs = myobjstr.ordr
       endelse
       nthese_ordrs = n_elements(these_ordrs)

; store the relevant spectra of all the orders in these arrays 
       biglogwave = fltarr(5000,nthese_ordrs)
       bigflux = biglogwave*0.0
       bigivar = biglogwave*0.0
       bigskyflux = biglogwave*0.0

       spec1dfile = datapath+'spec1d/'+prefix+$
         strtrim(mike[jj].obj,2)+'.fits'
;      spec1dfile = datapath+'spec1d/'+prefix+'obj'+$
;        strtrim(mike[jj].obj,2)+'.'+mike[jj].img_root
       for oo = 0L, nthese_ordrs-1L do begin
          indx = where(objstr.order eq these_ordrs[oo])

          if keyword_set(fluxed) then begin ; spectra in flux units
             objwave = objstr[indx].wave
             objflux = objstr[indx].flux*1D-16
             objvar = (objstr[indx].sig*1D-16)^2.0
          endif else begin      ; spectra in counts
             objwave = objstr[indx].box_wv
             objflux = objstr[indx].box_fx
             objvar = objstr[indx].box_var
          endelse
;         ww = where(objwave ne 0.0) & djs_plot, objwave[ww], objflux[ww]
          
          good = where((objflux gt 0.0) and (objvar gt 0.0))
          min_logwave = alog10(min(objwave[good]))
          max_logwave = alog10(max(objwave[good]))
          logwave = dindgen((max_logwave-min_logwave)/cdelt+1.0d)*cdelt+min_logwave
          npix = n_elements(logwave)

          linterp, alog10(objwave[good]), objflux[good], logwave, flux, missing=0.0
          linterp, alog10(objwave[good]), objvar[good], logwave, var, missing=0.0
          ivar = 1.0/(var+(var eq 0.0))*(var ne 0.0)

          good = where((objwave gt 0.0) and (objvar gt 0.0) and $
            (objstr[indx].sky_wv gt 0.0) and (objstr[indx].sky gt 0.0))
          linterp, alog10(objstr[indx].sky_wv[good]), objstr[indx].sky[good], $
            logwave, skyflux, missing=0.0

          biglogwave[0:npix-1,oo] = logwave
          bigflux[0:npix-1,oo] = flux
          bigivar[0:npix-1,oo] = ivar
          bigskyflux[0:npix-1,oo] = skyflux

; write out
          spec1dfile_ordr = repstr(spec1dfile,'.fits','.'+$
            'ord'+string(objstr[indx].order,format='(I2.2)')+'.fits')
          if (not keyword_set(std)) then begin
             alpha_write_spec1d, logwave, flux, ivar, skyflux, $
               info[jj], cdelt=cdelt, outfile=spec1dfile_ordr, std=std
          endif

       endfor

; combine all the orders into a single 1D spectrum

       min_logwave = 1E6
       max_logwave = 0.0
       for oo = 0L, nthese_ordrs-1L do begin
          good = where(biglogwave[*,oo] gt 0.0)
          min_logwave = min_logwave < min(biglogwave[good,oo])
          max_logwave = max_logwave > max(biglogwave[good,oo])
       endfor

       logwave = dindgen((max_logwave-min_logwave)/cdelt+1.0d)*cdelt+min_logwave
       npix = n_elements(logwave)

       combineflux = fltarr(npix,nthese_ordrs)
       combinevar = combineflux*0.0
       combineivar = combineflux*0.0
       combineskyflux = combineflux*0.0
       for oo = 0L, nthese_ordrs-1L do begin
          good = where(biglogwave[*,oo] gt 0.0)
          thislogwave = biglogwave[good,oo]
          thisflux = bigflux[good,oo]
          thisivar = bigivar[good,oo]
          thisskyflux = bigskyflux[good,oo]
          thisvar = 1.0/(thisivar+(thisivar eq 0.0))*(thisivar ne 0.0)

          linterp, thislogwave, thisflux, logwave, flux1, missing=0.0
          linterp, thislogwave, thisvar, logwave, var1, missing=0.0
          linterp, thislogwave, thisskyflux, logwave, skyflux1, missing=0.0
          ivar1 = 1.0/(var1+(var1 eq 0.0))*(var1 ne 0.0)

          combineflux[*,oo] = flux1
          combineskyflux[*,oo] = skyflux1
          combinevar[*,oo] = var1
          combineivar[*,oo] = ivar1
       endfor

       if (nthese_ordrs eq 1L) then begin
          flux = combineflux
          var = combineivar
          ivar = combineivar
          skyflux = combineskyflux
       endif else begin
          weight = total(combineivar,2)
          notzero = where((weight ne 0L),nnotzero)
          if (nnotzero eq 0L) then message, 'This is not good'
          flux = fltarr(npix) & var = fltarr(npix)
          flux[notzero] = total(combineflux[notzero,*]*combineivar[notzero,*],2)/weight[notzero]
          var[notzero] = 1.0/weight[notzero]
          ivar = 1.0/(var+(var eq 0.0))*(var ne 0.0)
          skyflux = total(combineskyflux,2)/float(nthese_ordrs) ; simple mean
       endelse

;      if strmatch(info[jj].obj,'*044*') then stop
;      djs_plot, 10^logwave, flux, ysty=3
;      djs_oplot, 10^biglogwave[*,0], bigflux[*,0], color='red'
;      djs_oplot, 10^biglogwave[*,1], bigflux[*,1], color='blue'
;      djs_oplot, 10^biglogwave[*,2], bigflux[*,2], color='green'

; clean up cosmic rays; interpolate and set the inverse variance to
; zero; also interpolate over blank regions in the spectral orders
       if (not keyword_set(std)) then begin
          emask = emission_mask(10^logwave,z=info[jj].z,$
            /gauss,spectrum=flux,ivarspectrum=ivar,/noskymask)
          djs_iterstat, flux, mask=crmask, sigrej=5.0
          crmask[where(emask eq 0)] = 1 ; don't mask the emission lines
          domask = where(crmask eq 0,ndomask)
          flux = djs_maskinterp(flux,(crmask eq 0),logwave,/const)
;         flux = djs_maskinterp(flux,((crmask eq 0) or (flux eq 0.0)),logwave,/const)
          if (ndomask ne 0L) then begin
             ivar[domask] = 0.0
          endif
       endif       

;      if strmatch(info[jj].obj,'*044*') then stop
       
; write out
       alpha_write_spec1d, logwave, flux, ivar, skyflux, $
         info[jj], cdelt=cdelt, outfile=spec1dfile, std=std
       
; make the QA-plot          
       xrange = minmax(10^logwave)
       yrange = minmax(scale*flux)
       if (not keyword_set(std)) then begin
          oiiimask = emission_mask(10^logwave,z=info[jj].z,$
            width=15.0,linelist=[4958.91,5006.84])
          ww = where(oiiimask eq 0,nww)
          if (nww ne 0L) then yrange = [0.0,max(scale*flux[ww])] else $
            yrange = [0.0,max(scale*flux)]
       endif

;      yrange = im_max(scale*flux,sigrej=5.0,/nan)*[-0.08,1.2]
       xtitle = textoidl('Observed Wavelength (\AA)')

       if keyword_set(std) then begin
          title = strtrim(info[jj].hr,2)+'/'+strtrim(info[jj].name,2)+', '+$
            repstr(strtrim(info[jj].obj,2),'_',' ')
       endif else begin
          title = strtrim(info[jj].galaxy,2)+', '+repstr(strtrim(info[jj].obj,2),'_','')+$
            ', z='+strtrim(string(info[jj].z,format='(F12.5)'),2)
       endelse
          
       plot, 10.0^logwave, scale*flux, ps=10, /xsty, ysty=3, xrange=xrange, $
         yrange=yrange, title=title, xtitle=xtitle, ytitle=ytitle
       legend, 'Combined Spectrum', /left, /top, box=0, charsize=1.6
       if keyword_set(std) then begin ; overplot the published spectrum
          readcol, getenv('MIKE_DIR')+'/pro/Std/ESOFIL/'+$
            strtrim(info[jj].esofil,2), swave, sflux, /silent
          sflux = sflux*1D-16
          djs_oplot, swave, sflux*scale, ps=-4, color='red'
       endif else begin
          colors = ['red','blue','dark green','purple','orange']
          for oo = 0L, nthese_ordrs-1L do begin
             good = where(biglogwave[*,oo] gt 0.0)
             if (oo eq 0L) then plot, [0], [0], /nodata, ps=10, /xsty, ysty=3, $
               xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, title=title
             djs_oplot, 10.0^biglogwave[good,oo], scale*bigflux[good,oo], ps=10, $
               color=colors[oo]
             legend, 'Individual Orders', /left, /top, box=0, charsize=1.6
          endfor
       endelse
    endfor 

    if (n_elements(qafile) ne 0L) then begin
       dfpsclose
       spawn, 'gzip -f '+qafile
       im_plotfaves
    endif

return
end
    
