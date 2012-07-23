pro alpha_final_spec1d, datapath, outpath, fluxed=fluxed
; jm09jan24nyu - coadd the final 1D spectra

crummy code    
    
    npath = n_elements(datapath)
    for ii = 0L, npath-1L do begin
       if keyword_set(fluxed) then $
         fits = file_search(datapath[ii]+'f.obj_*.fits',count=nfits) else $
           fits = file_search(datapath[ii]+'obj_*.fits',count=nfits)
       if (nfits gt 0L) then begin
          fits = fits[where(strmatch(fits,'*ord*',/fold) eq 0)]
          for jj = 0L, n_elements(fits)-1L do begin
             hdr = headfits(fits[jj])
             info1 = im_hdr2struct(hdr)
             info1 = struct_addtags({file: fits[jj]},info1)
             if (n_elements(info) eq 0L) then info = info1 else info = [info,info1]
          endfor
       endif
    endfor

;   struct_print, struct_trimtags(info,sel=['object','galaxy','z'])

    if keyword_set(fluxed) then begin
       scale = 1D16
       ytitle = textoidl('Flux (10^{-16} '+flam_units()+')')
       prefix = 'f.'
       qafile = outpath+'qa_final_spec1d_fluxed.ps'
    endif else begin
       scale = 1.0
       ytitle = textoidl('Flux (ADU)')
       prefix = ''
       qafile = outpath+'qa_final_spec1d.ps'
    endelse

    splog, 'Building QA-plot '+qafile
    dfpsplot, qafile, /landscape, /color
    im_plotfaves, thick=5.0, /post
       
    allobj = strtrim(info.object,2)
    allobj = repstr(repstr(repstr(repstr(repstr(allobj,$
      '_a',''),'_b'),'_c',''),'_d',''),'_e','')
    obj = allobj[uniq(allobj,sort(allobj))]
    nobj = n_elements(obj)

    for jj = 0, nobj-1 do begin
       these = where(obj[jj] eq allobj,nthese)
       outfile = outpath+prefix+obj[jj]+'.fits'
       outobj = allobj[these[0]]
       for kk = 0L, nthese-1L do begin
          if (kk eq 0L) then begin
             flux1 = mrdfits(info[these[kk]].file,0,hdr1,/silent)
             ivar1 = mrdfits(info[these[kk]].file,1,/silent)
             skyflux1 = mrdfits(info[these[kk]].file,2,/silent)
             logwave1 = make_wave(hdr1)
             npix = n_elements(logwave1)
; pack into a big array
             combineflux = fltarr(npix,nthese)
             combineskyflux = combineflux*0.0
             combineivar = combineflux*0.0
          endif else begin
             thisflux = mrdfits(info[these[kk]].file,0,thishdr,/silent)
             thisivar = mrdfits(info[these[kk]].file,1,/silent)
             thisskyflux = mrdfits(info[these[kk]].file,2,/silent)
             thislogwave = make_wave(thishdr)
             linterp, thislogwave, thisflux, logwave1, flux1, missing=0.0
             linterp, thislogwave, thisskyflux, logwave1, skyflu1, missing=0.0
             thisvar = 1.0/(thisivar+(thisivar eq 0.0))*(thisivar ne 0.0)
             linterp, thislogwave, thisvar, logwave1, var1, missing=0.0
             ivar1 = 1.0/(var1+(var1 eq 0.0))*(var1 ne 0.0)
          endelse
          combineflux[*,kk] = flux1
          combineskyflux[*,kk] = skyflux1
          combineivar[*,kk] = ivar1
       endfor
; combine using inverse variance weights
       if (nthese eq 1L) then begin
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
          skyflux = total(combineskyflux,2)/float(nthese) ; simple mean
       endelse

; clean up cosmic rays; interpolate and set the inverse variance to
; zero; also interpolate over blank regions in the spectral orders
       emask = emission_mask(10^logwave1,z=info[these[0]].z,$
         /gauss,spectrum=flux,ivarspectrum=ivar,/noskymask)
       djs_iterstat, flux, mask=crmask, sigrej=5.0
       crmask[where(emask eq 0)] = 1 ; don't mask the emission lines
       domask = where(crmask eq 0,ndomask)
       flux = djs_maskinterp(flux,(crmask eq 0),logwave,/const)
;      flux = djs_maskinterp(flux,((crmask eq 0) or (flux eq 0.0)),logwave,/const)
       if (ndomask ne 0L) then begin
          ivar[domask] = 0.0
       endif

;      if (outobj eq 'obj_044') then stop

; write out; possibly add the dates observed
       mkhdr, hdr, float(flux), /extend
       sxdelpar, hdr, 'COMMENT'
       sxdelpar, hdr, 'DATE'
       sxaddpar, hdr, 'OBJECT', outobj
       sxaddpar, hdr, 'GALAXY', info[these[0]].galaxy
       sxaddpar, hdr, 'EXPTIME', float(total(info[these].exptime))
       sxaddpar, hdr, 'NCOMBINE', fix(nthese)
;      sxaddpar, hdr, 'AIRMASS', float(info[these[0]].am)
       sxaddpar, hdr, 'RA', info[these[0]].ra, format='(F11.7)'
       sxaddpar, hdr, 'DEC', info[these[0]].dec, format='(F13.9)'
       sxaddpar, hdr, 'Z', info[these[0]].z
       sxaddpar, hdr, 'CTYPE1', 'LOG'
       sxaddpar, hdr, 'CRPIX1', 1.0d
       sxaddpar, hdr, 'CRVAL1', min(logwave1)
       sxaddpar, hdr, 'CDELT1', info[these[0]].cdelt1
       sxaddpar, hdr, 'DC-FLAG', 1

       splog, 'Writing '+outfile
       mwrfits, float(flux), outfile, hdr, /create
       mwrfits, float(ivar), outfile, hdr
       mwrfits, float(skyflux), outfile, hdr

; make the QA-plot          
       xrange = minmax(10^logwave1)
       oiiimask = emission_mask(10^logwave1,z=info[these[0]].z,$
         width=15.0,linelist=[4958.91,5006.84])
       yrange = [0.0,max(scale*flux[where(oiiimask eq 0)])]
;      yrange = im_max(scale*flux,sigrej=5.0)*[-0.08,1.2]

       xtitle = textoidl('Observed Wavelength (\AA)')
       title = strtrim(info[these[0]].galaxy,2)+', '+repstr(outobj,'_','')+$
         ', z='+strtrim(string(info[these[0]].z,format='(F12.5)'),2)
       colors = ['red','blue','dark green','purple','orange','cyan']
       
       plot, 10.0^logwave1, scale*flux, ps=10, /xsty, ysty=3, xrange=xrange, $
         yrange=yrange, title=title, xtitle=xtitle, ytitle=ytitle
       legend, 'Final 1D Spectrum', /left, /top, box=0, charsize=1.6
       for kk = 0L, nthese-1L do begin
          if (kk eq 0L) then plot, [0], [0], /nodata, ps=10, /xsty, ysty=3, $
            xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, title=title
          djs_oplot, 10.0^logwave1, scale*combineflux[*,kk], ps=10, $
            color=colors[kk]
       endfor
       legend, 'Individual Spectra', /left, /top, box=0, charsize=1.6
    endfor

; close the QA-plot
    
    if (n_elements(qafile) ne 0L) then begin
       dfpsclose
       spawn, 'gzip -f '+qafile
       im_plotfaves
    endif

return
end
    
