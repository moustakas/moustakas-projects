pro alpha_coadd_spec1d, mike, info, side=side, datapath=datapath, $
  qafile=qafile, fluxed=fluxed, std=std
; jm09jan06nyu - coadd the 1D spectra, both individual orders and
; repeat observations of the same object; if /FLUXED then coadd the
; fluxed 1D spectra, otherwise coadd the counts

    nobj = n_elements(info)
    if (nobj ne n_elements(mike)) then begin
       splog, 'Dimensions of INFO and MIKE do not match'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = './'

    pixpad = 300
    light = 299792.458D ; [km/s]
;   velpix = (side eq 1 ? 1.50d : 2.10d) * double(mike[0].rowbin) ; [km/s]
;   cdelt = alog10(1.0d + velpix/light) ; pixel size

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

    for jj = 0, nobj-1 do begin
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
;      if jj eq 1 then stop
; get the min/max wavelengths across all the orders of interest        
       minwave_all = dblarr(nthese_ordrs)
       maxwave_all = dblarr(nthese_ordrs)
       velpix_all = dblarr(nthese_ordrs)
       for oo = 0, nthese_ordrs-1 do begin
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
          good = where((objvar gt 0.0) and (objwave gt 0.0),ngood)
          minwave_all[oo] = min(objwave[good])
          maxwave_all[oo] = max(objwave[good])
          velpix_all[oo] = light*alog(maxwave_all[oo]/minwave_all[oo])/double(ngood)
       endfor
       minwave = min(minwave_all)
       maxwave = max(maxwave_all)
       velpix = djs_median(velpix_all)
       splog, minwave, maxwave, velpix_all, velpix
       
;; force the minimum wavelength for each order to be an integer
;; multiple of the minimum across all the orders, so that when we
;; combine the orders, below, we can simply shift and add, without
;; having to rebin again
;       nfact = ceil((alog(minwave_all)-alog(minwave))/(velpix/light)) ; round up
;       minwave_out = exp(alog(minwave) + nfact*(velpix/light))

       spec1dfile = datapath+'spec1d/'+prefix+$
         strtrim(mike[jj].obj,2)+'.fits'
;      spec1dfile = datapath+'spec1d/'+prefix+'obj'+$
;        strtrim(mike[jj].obj,2)+'.'+mike[jj].img_root
       for oo = 0, nthese_ordrs-1 do begin
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
          good = where((objvar gt 0.0) and (objwave gt 0.0),npix)
          srt = sort(objwave[good])
          objwave = objwave[good[srt]]
          objflux = objflux[good[srt]]
          objvar = objvar[good[srt]]

; trim PIXPAD pixels from either side
          allpix = lindgen(npix)
          good = where(allpix ge pixpad and allpix le npix-pixpad-1)
;         if strmatch(spec1dfile,'*obj_048_a*') then stop
          objwave = objwave[good]
          objflux = objflux[good]
          objvar = objvar[good]
           
          objskyflux = objstr[indx].sky
          objskywave = objstr[indx].sky_wv
          skygood = where((objskywave gt 0.0))
          skysrt = sort(objskywave[skygood])
          objskywave = objskywave[skygood[skysrt]]
          objskyflux = objskyflux[skygood[skysrt]]
          
; rebin to constant velocity pixels, conserving flux density
; (intensity)
          flux = im_log_rebin(objwave,objflux,var=objvar,/wavestrict,$
            vsc=velpix,outwave=lnwave,outvar=var,minwave=minwave,maxwave=maxwave)
          skyflux = im_log_rebin(objskywave,objskyflux,vsc=velpix,$
            /wavestrict,minwave=minwave,maxwave=maxwave)
          ivar = 1.0/(var+(var eq 0))*(var ne 0)
;         npix = n_elements(flux)

          if (oo eq 0) then begin
             bigflux = flux 
             bigivar = ivar
             bigskyflux = skyflux
          endif else begin
             bigflux = [[bigflux],[flux]]
             bigivar = [[bigivar],[ivar]]
             bigskyflux = [[bigskyflux],[skyflux]]
          endelse

; write out
          spec1dfile_ordr = repstr(spec1dfile,'.fits','.'+$
            'ord'+string(objstr[indx].order,format='(I2.2)')+'.fits')
          if (keyword_set(std) eq 0) then begin
             alpha_write_spec1d, lnwave, flux, ivar, skyflux, $
               info[jj], velpix=velpix, outfile=spec1dfile_ordr, std=std
          endif
       endfor

; combine all the orders into a single 1D spectrum
;      min_lnwave_all = dblarr(nthese_ordrs)
;      max_lnwave_all = dblarr(nthese_ordrs)
;      for oo = 0, nthese_ordrs-1 do begin
;         good = where(biglnwave[*,oo] gt 0.0)
;         min_lnwave_all[oo] = min(biglnwave[good,oo])
;         max_lnwave_all[oo] = max(biglnwave[good,oo])
;      endfor
;      min_lnwave = min(min_lnwave_all)
;      max_lnwave = max(max_lnwave_all)
;
;      npix = round((max_lnwave-min_lnwave)/(velpix/light)+1)
;      lnwave = dindgen(npix)*velpix/light+min_lnwave

;      combine1fiber, biglnwave, bigflux, bigivar, newloglam=lnwave, $
;        newflux=combineflux, newivar=combineivar

       if (nthese_ordrs eq 1) then begin
          flux = bigflux
          ivar = bigivar
          skyflux = bigskyflux
       endif else begin
          ivar = total(bigivar,2)
          flux = total(bigflux*bigivar,2)/(ivar+(ivar eq 0))*(ivar ne 0)
          skyflux = total(bigskyflux,2)/float(nthese_ordrs) ; simple mean
       endelse
;      if strmatch(spec1dfile,'*obj_048_a*') then stop


;      if strmatch(info[jj].obj,'*046*') then begin
;         dfpsclose & im_plotfaves
;         stop
;      endif
;      djs_plot, 10^lnwave, flux, ysty=3
;      djs_oplot, 10^biglnwave[*,0], bigflux[*,0], color='red'
;      djs_oplot, 10^biglnwave[*,1], bigflux[*,1], color='blue'
;      djs_oplot, 10^biglnwave[*,2], bigflux[*,2], color='green'

; clean up cosmic rays; interpolate and set the inverse variance to
; zero; also interpolate over blank regions in the spectral orders
       if (keyword_set(std) eq 0) then begin
          emask = emission_mask(exp(lnwave),z=info[jj].z,$
            /gauss,spectrum=flux,ivarspectrum=ivar,/noskymask)
          djs_iterstat, flux, mask=crmask, sigrej=5.0
          crmask[where(emask eq 0)] = 1 ; don't mask the emission lines
          domask = where(crmask eq 0,ndomask)
          flux = djs_maskinterp(flux,(crmask eq 0),lnwave,/const)
;         flux = djs_maskinterp(flux,((crmask eq 0) or (flux eq 0.0)),lnwave,/const)
          if (ndomask ne 0L) then begin
             ivar[domask] = 0.0
          endif
       endif       

;      if strmatch(info[jj].obj,'*046*') then stop
       
; write out
       alpha_write_spec1d, lnwave, flux, ivar, skyflux, $
         info[jj], velpix=velpix, outfile=spec1dfile, std=std
       
; make the QA-plot          
       xrange = minmax(exp(lnwave))
       yrange = minmax(scale*flux)

       if (keyword_set(std) eq 0) then begin
          oiiimask = emission_mask(exp(lnwave),z=info[jj].z,$
            width=2.0,linelist=[4958.91,5006.84])
          ww = where(oiiimask eq 0,nww)
          if (nww ne 0L) then $
            yrange = [im_min(scale*flux[ww],sigrej=1),im_max(scale*flux[ww],sigrej=10)] else $
              yrange = [im_min(scale*flux,sigrej=1),im_max(scale*flux,sigrej=10)]
       endif

       yrange = minmax(scale*flux)
;      yrange = im_max(scale*flux,sigrej=5.0,/nan)*[-0.08,1.2]
       xtitle = textoidl('Observed Wavelength (\AA)')

       if keyword_set(std) then begin
          title = strtrim(info[jj].hr,2)+'/'+strtrim(info[jj].name,2)+', '+$
            repstr(strtrim(info[jj].obj,2),'_',' ')
       endif else begin
          title = strtrim(info[jj].galaxy,2)+', '+repstr(strtrim(info[jj].obj,2),'_','')+$
            ', z='+strtrim(string(info[jj].z,format='(F12.5)'),2)
       endelse
          
       plot, exp(lnwave), scale*flux, ps=10, /xsty, ysty=3, xrange=xrange, $
         yrange=yrange, title=title, xtitle=xtitle, ytitle=ytitle
       legend, 'Combined Spectrum', /left, /top, box=0, charsize=1.6
       if keyword_set(std) then begin ; overplot the published spectrum
          readcol, getenv('MIKE_DIR')+'/pro/Std/ESOFIL/'+$
            strtrim(info[jj].esofil,2), swave, sflux, /silent
          sflux = sflux*1D-16
          djs_oplot, swave, sflux*scale, ps=-4, color='red'
       endif else begin
          colors = ['red','blue','dark green','purple','orange']
          for oo = 0, nthese_ordrs-1 do begin
             good = where(bigivar[*,oo] gt 0,ngood)
;            good = lindgen(n_elements(bigivar[*,oo]))
             if (oo eq 0) then plot, [0], [0], /nodata, ps=10, /xsty, ysty=3, $
               xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, title=title
             djs_oplot, exp(lnwave[good]), scale*bigflux[good,oo], ps=10, $
               color=colors[oo]
;            djs_oplot, exp(lnwave), scale*flux, psym=10, color='black'
             legend, 'Individual Orders', /left, /top, box=0, charsize=1.6
          endfor
       endelse
    endfor 

    if (n_elements(qafile) ne 0L) then begin
       dfpsclose
;      spawn, 'ps2pdf '+qafile+' '+repstr(qafile,'.ps','.pdf')+' ; \rm '+qafile, /sh
;      spawn, 'gzip -f '+qafile
       im_plotfaves
    endif

return
end
    
