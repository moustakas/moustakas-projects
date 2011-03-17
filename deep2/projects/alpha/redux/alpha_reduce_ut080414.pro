function edit_mike_ut080414, rawmike
; jm08apr13nyu - edit the mike structure

    mike = rawmike
    mike.obj = strtrim(mike.obj,2)
    mike.obj = repstr(mike.obj,' ','_')

    mike.type = strtrim(mike.type,2)
    
; only analyze objects with 1x1 binning and the 0.7" slit; also
; only analyze the red side for now!

   good = where((mike.rowbin eq 1) and (mike.colbin eq 1) and $
     (mike.slit eq 0.7) and (mike.side ne 0),comp=bad,ncomp=nbad)
    mike[good].setup = 1
    if (nbad ne 0L) then mike[bad].flg_anly = 0

; these were water vapor observations for LCO    

    vapor = where(strmatch(mike.obj,'*water*',/fold))
    mike[vapor].flg_anly = 0

; twilight flats mis-identified as objects in the red side    
    twi = where(strmatch(mike.obj,'*twilight*',/fold))
    mike[twi].type = 'TWI' 

; rename the "internal flats" as "trace flats"; also, on the red side
; they were mis-identified as "milky flats"
    lamp = where((strmatch(mike.obj,'*qz_lamp*',/fold) eq 1) and $
      (strmatch(mike.obj,'*diffuser*',/fold) eq 0))
    mike[lamp].type = 'TFLT'

; object mis-identified as a zero
;   fixme = where(strtrim(mike.obj,2) eq '071_a')
;   mike[fixme].type = 'OBJ'

; because of the pre-processing, all the science spectra were
; mis-identified as zeros
    fixme = where(strtrim(mike.type,2) eq 'ZRO')
    mike[fixme].type = 'OBJ'
    
; only keep objects that we need    
    
    realgood = where(mike.setup eq 1 and mike.flg_anly eq 1)
    struct_print, struct_trimtags(mike[realgood],sel=['frame',$
      'img_root','obj','colbin','rowbin','type','slit','obj_id',$
      'flg_anly','setup'])

return, mike[realgood]
end

pro alpha_reduce_ut080414, init=init, edit=edit, blue=blue, flat=flat, $
  arc=arc, slitflat=slitflat, proc=proc, emlines=emlines, dotrace=dotrace, $
  skysub=skysub, extract=extract, calibrate=calibrate, combine=combine, $
  coadd=coadd, standards=standards, makesens=makesens, clobber=clobber
; jm08apr14nyu - default is to just reduce the red side
    
; alpha_reduce_ut080414, /proc, /clobber, /emlines, /dotrace, /skysub, /extract, /calibrate, /coadd
    
    datapath = getenv('DEEP2_ALPHA_DIR')+'/ut080414/'
    
    setup = 1L
    if keyword_set(blue) then side = 1L else side = 2L

; ---------------------------------------------------------------------------
; constants that will be used below

    sz_img = [2048L,4096L] ; stupidly hard-coded
    
    light = 299792.458D ; [km/s]
    deltav = 300.0      ; [km/s]
    deltav_zoom = 50.0  ; [km/s]
    snrmin = 4.0        ; minimum S/N for a line to be masked
    
    linewave = [4861.0,4958.911,5006.843] ; rest wavelength [A]
    
; ---------------------------------------------------------------------------
; pre-processing

    if keyword_set(init) then begin
       mike_strct, rawmike, /noedit, outfil='rawmike.fits'
       return
    endif

    if keyword_set(edit) then begin
       if (n_elements(rawmike) eq 0L) then $
         rawmike = mrdfits('rawmike.fits',1)

       mike = edit_mike_ut080414(rawmike)
       mike_wrstrct, mike, fits='mike.fits'

       mike_setup, mike
       mike_wrstrct, mike, fits='mike.fits'
       return
    endif
    
    if (n_elements(mike) eq 0L) then mike = mike_ar('mike.fits')

    velpix = (side eq 1 ? 1.50d : 2.10d) * double(mike[0].rowbin) ; [km/s]
    cdelt = alog10(1.0d + velpix/light) ; pixel size

; ---------------------------------------------------------------------------
; process the flats and optionally verify

    if keyword_set(flat) then begin
; (1) generate a stacked, normalized milky flat; (2) combine the trace
; flats for order and slit tracing; (3) generate a smooth model of the
; order curvature using the trace flat previously created
       if keyword_set(stepbystep) then begin
          mike_mkmflat, mike, setup, side, smooth=0, clobber=clobber
          mike_mktflat, mike, setup, side, clobber=clobber
          mike_edgeflat, mike, setup, side, inter=inter, chk=chk, clobber=clobber
       endif else begin ; do it all in one step
          mike_allflat, mike, setup, side, clobber=clobber
       endelse
       spawn, 'gv QA/Flats01/qa_trcflt_01R.ps.gz &'
;      xatv, 'Flats/Flat_R_01_M.fits' 
;      xatv, 'Flats/Flat_R_01_T.fits' 
       xatv, 'Flats/Flat_R_01_T.fits' 
       mike_chktrcflat, mike, setup, side, /nostop, /fit
    endif    
    
; ---------------------------------------------------------------------------
; process the arcs

    if keyword_set(arc) then begin
       if keyword_set(stepbystep) then begin
; pre-process the arcs and do the tracing (do not fit and do not
; generate the 2D image)
;         guess_arc = mike_getfil('guess_arc',side=side,/name,chkfil=chkfil,sz=[2048,4096])
;         mike_allarc, mike, setup, side, fits='mike.fits', $
;           clobber=clobber, chk=chk, /nowav, /noimg
       endif else begin
          mike_allarc, mike, setup, side, fits='mike.fits', $
            clobber=clobber, chk=chk, indx=21L
       endelse
    endif

; construct the slit profile
    if keyword_set(slitflat) then begin
       mike_slitflat, mike, setup, side, clobber=clobber, chk=chk
    endif

; ---------------------------------------------------------------------------
; deal with the standards    
    stdindx = where(mike.setup EQ setup AND mike.type EQ 'STD')
    std = mike[stdindx].obj_id
    nstd = n_elements(std)

    if keyword_set(standards) then begin
       if keyword_set(proc) then begin
          for kk = 0L, nstd-1L do mike_proc, mike, setup=setup, $
            obj=std[kk], side=side, clobber=clobber, /std
       endif
       if keyword_set(dotrace) then begin
          for kk = 0L, nstd-1L do mike_fntobj, mike, setup, $
            std[kk], side, chk=chk, /std
       endif
       if keyword_set(skysub) then begin
          for kk = 0L, nstd-1L do mike_skysub, mike, setup, $
            stdindx[kk], side, chk=chk, /std ; note STDINDX!
       endif
       if keyword_set(extract) then begin
          for kk = 0L, nstd-1L do mike_box, mike, setup, $
            std[kk], side, chk=chk, ochk=ochk, reschk=reschk, /std, $
            /skipskysub, nohelio=0, novac=0
       endif
; build the sensitivity function; not generalized!!
; http://www.ctio.noao.edu/spectrographs/4m_Echelle/standards.html 
       if keyword_set(makesens) then begin
          for kk = 0L, nstd-1L do begin
;            esofil = getenv('MIKE_DIR')+'/pro/Std/ESOFIL/'+'fhr3454.dat'
;            objfil = mike[stdindx[kk]].obj_fil
;            outfil = mike_getfil('sens_fil', SUBFIL=mike[stdindx[kk]].img_root,/name)
;            objstr = xmrdfits(objfil,1,STRUCTYP='mikeobjstrct',/silent)
;            x_calibstd, objstr.box_wv[0L:objstr[0].nrow-1L], objstr.box_fx[0L:objstr[0].nrow-1L], $
;              outfil, exp=mike[stdindx[kk]].exp, eso_fil=esofil, swv=swv, sfx=sfx, $
;              sens=sens, chkfit=0, /mab
             mike_calibstd, mike, stdindx[kk], esofil='fhr3454.dat' ; note STDINDX!
          endfor
       endif 
; flux-calibrate
       if keyword_set(calibrate) then begin
          fluxfil = 'Extract/Sens_mr0025.idl' ; not generalized!!
          for kk = 0L, nstd-1L do begin
             mike_flux, mike, setup, stdindx[kk], side, /std, $ ; note STDINDX!
               fluxfil=fluxfil, /boxcar
          endfor
       endif
; combine multiple exposures of the same object, which must be run on
; even a single object
;      if keyword_set(combine) then begin
;         for kk = 0L, nstd-1L do begin
;            mike_combspec, mike, setup, stdindx[kk], side, /std, /fchk ; note STDINDX!
;         endfor
;      endif 
; read the extracted 1D spectra and stitch the orders together
       if keyword_set(coadd) then begin
; read the published spectrum; NOT GENERALIZED!!
          readcol, getenv('MIKE_DIR')+'/pro/Std/ESOFIL/'+'fhr3454.dat', $
            swave, sflux, /silent
          sflux = sflux*1D-16
          
; generate a QA-plot
          qafile = datapath+'QA/qa_extract_std.ps'
          splog, 'Building QA-plot '+qafile
          dfpsplot, qafile, /landscape, /color
          im_plotfaves, thick=5.0, /post
          for kk = 0L, nstd-1L do begin
             spec1dfile = datapath+'spec1d/'+strtrim(mike[stdindx[kk]].obj,2)+$
               '.'+mike[stdindx[kk]].img_root
             objfil = mike_getfil('obj_fil',setup,$
               SUBFIL=mike[stdindx[kk]].img_root,/name)
             splog, 'Reading '+objfil
             objstr = xmrdfits(objfil,1,STRUCTYP='mikeobjstrct',/silent)

             these_ordrs = objstr.order
             nthese_ordrs = n_elements(these_ordrs)
             
             biglogwave = fltarr(5000,nthese_ordrs)
             bigflux = biglogwave*0.0
             bigivar = biglogwave*0.0
             bigskyflux = biglogwave*0.0

             for oo = 0L, nthese_ordrs-1L do begin
                indx = where(objstr.order eq these_ordrs[oo])
                spec1dfile_ordr = repstr(spec1dfile,'.fits','.'+$
                  'ord'+string(objstr[indx].order,format='(I2.2)')+'.fits')

                objwave = objstr[indx].wave
; spectra in counts
;               objflux = objstr[indx].box_fx
;               objvar = objstr[indx].box_var
; spectra in flux units
                objflux = objstr[indx].flux*1D-16
                objvar = (objstr[indx].sig*1D-16)^2.0
                
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
                flux = total(combineflux*combineivar,2)/total(combineivar,2)
                var = 1.0/total(combineivar,2)
                ivar = 1.0/var
                skyflux = total(combineskyflux,2)/float(nthese_ordrs) ; simple mean
             endelse

; write out          
             mkhdr, hdr, float(flux), /extend
             sxdelpar, hdr, 'COMMENT'
             sxdelpar, hdr, 'DATE'
;            sxaddpar, hdr, 'OBJECT', info[jj].obj
;            sxaddpar, hdr, 'GALAXY', info[jj].galaxy
;            sxaddpar, hdr, 'EXPTIME', float(info[jj].exp)
;            sxaddpar, hdr, 'AIRMASS', float(info[jj].am)
;            sxaddpar, hdr, 'RA', info[jj].ra, format='(F11.7)'
;            sxaddpar, hdr, 'DEC', info[jj].dec, format='(F13.9)'
;            sxaddpar, hdr, 'Z', info[jj].z
;            sxaddpar, hdr, 'RAWFITS', info[jj].img_root
             sxaddpar, hdr, 'CTYPE1', 'LOG'
             sxaddpar, hdr, 'CRPIX1', 1.0d
             sxaddpar, hdr, 'CRVAL1', min(logwave)
             sxaddpar, hdr, 'CDELT1', cdelt
             sxaddpar, hdr, 'DC-FLAG', 1
             
             splog, 'Writing '+spec1dfile
             mwrfits, float(flux), spec1dfile, hdr, /create
             mwrfits, float(ivar), spec1dfile, hdr
             mwrfits, float(skyflux), spec1dfile, hdr

; make the QA-plot          

             scale = 1D10
             xrange = minmax(10^logwave)
             yrange = im_max(scale*flux,sigrej=5.0)*[-0.08,1.2]
             xtitle = textoidl('Observed Wavelength (\AA)')
             ytitle = textoidl('Flux (10^{-10} '+flam_units()+')')
             title = repstr(strtrim(mike[stdindx[kk]].obj,2),'_',' ')
;            colors = ['red','blue','dark green','purple','orange']
             
             plot, 10.0^logwave, scale*flux, ps=10, /xsty, /ysty, xrange=xrange, $
               yrange=yrange, title=title, xtitle=xtitle, ytitle=ytitle
             legend, 'Combined Spectrum', /left, /top, box=0, charsize=1.6
;            for oo = 0L, nthese_ordrs-1L do begin
;               good = where(biglogwave[*,oo] gt 0.0)
;               if (oo eq 0L) then plot, [0], [0], /nodata, ps=10, /xsty, /ysty, $
;                 xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, title=title
;               djs_oplot, 10.0^biglogwave[good,oo], scale*bigflux[good,oo], ps=10;, $
;                 color=colors[oo]
;               legend, 'Individual Orders', /left, /top, box=0, charsize=1.6
;            endfor
; overplot the published spectrum
             djs_oplot, swave, sflux*scale, ps=-4, color='red'
          endfor
; close the QA-plot
          dfpsclose
          spawn, 'gzip -f '+qafile
          im_plotfaves
       endif
          
       return
    endif

; ---------------------------------------------------------------------------
; process all objects, find and trace the object, and sky-subtract 

; objects    
    objindx = where(mike.setup EQ setup AND mike.type EQ 'OBJ',nobjindx)
    objindx = objindx[uniq(mike[objindx].obj_id,sort(mike[objindx].obj_id))]
    objindx = objindx[0]
;   objindx = objindx[[0,2,3]]
    obj = mike[objindx].obj_id
    nobj = n_elements(obj)

; cross-match MIKE against the parent catalog of observed objects 
    info = alpha_match_observed_catalog(mike[objindx],$
      setup=setup,side=side)

;   niceprint, info.img_root, mike[objindx].img_root
;   niceprint, mike[objindx].frame, mike[objindx].obj, mike[objindx].obj_id

; ##################################################
; initialize the inverse variance map, overscan-subtract, and divide
; by the flat-field; the output is written in /Final
    if keyword_set(proc) then begin
       for jj = 0L, nobj-1L do mike_proc, mike, setup=setup, $
         obj=obj[jj], side=side, clobber=clobber
    endif
    
; ##################################################
; identify the orders of interest (i.e., those containing the emission
; lines)
    if keyword_set(emlines) then begin
       qafile = datapath+'QA/qa_trace_emlines.ps'
       splog, 'Building QA-plot '+qafile
       dfpsplot, qafile, /landscape, /color
       im_plotfaves, thick=5.0, /post
       for jj = 0L, nobj-1L do begin
; some input/output file names
          imgfil = mike_getfil('fin_fil',setup,/name,$
            SUBFIL=mike[objindx[jj]].img_root)
          outfil = repstr(imgfil,'.fits','_emlines.fits')
          outfil_objstr = mike_getfil('obj_fil',/name,$ 
            subfil='Extract/Obj_'+mike[objindx[jj]].img_root)
          outfil_objstr = repstr(outfil_objstr,'Obj_','myObj_')
          arcfil = strtrim(mike[objindx[jj]].arc_img,2)
; now read the processed 2D image and inverse variance map, the
; corresponding wavelength map, and the order structure
          splog, 'Reading '+imgfil
          img = xmrdfits(imgfil,0,/silent) 
          ivar = xmrdfits(imgfil,1,/silent)
          splog, 'Reading '+arcfil
          img_arc = xmrdfits(arcfil,/silent) 
          ordr_str = mike_getfil('ordr_str',setup,side=side)
; identify the orders that we care about
          sz = size(img,/dimen)
          ncol = sz[0] & nrow = sz[1]
          ordermask = x_ordermask(ncol,nrow,ordr_str,trim=0.0)
          these = where((img_arc*(ordermask gt 0.0) gt alog10(min(linewave)*(1.0+info[jj].z))) and $
            (img_arc*(ordermask gt 0.0) lt alog10(max(linewave)*(1.0+info[jj].z))),nthese)
          ordrs = ordermask[these[uniq(ordermask[these],sort(ordermask[these]))]]
          nordrs = n_elements(ordrs)
          splog, 'Identified emission lines in orders ', ordrs
; extract a 2D image centered on the order of interest and within
; +/-DELTAV km/s of each emission line, and then fit a 2D Gaussian
          line2dmodel = img*0.0
          for kk = 0L, nordrs-1L do begin
             splog, 'Working on order '+string(ordrs[kk],format='(I0)')
             indx = where(ordr_str.order eq ordrs[kk])
             for ll = 0L, n_elements(linewave)-1L do begin
                goodfit = 1
; zoom into the line to get good starting values
                these_zoom = where((ordermask eq ordr_str[indx].order) and (img_arc gt 3.0) and $
                  (img_arc gt alog10((1.0+info[jj].z)*linewave[ll]*(1.0-deltav_zoom/light))) and $
                  (img_arc lt alog10((1.0+info[jj].z)*linewave[ll]*(1.0+deltav_zoom/light))),nthese_zoom)
                if (nthese_zoom ne 0L) then begin
;                  splog, '  Fitting a 2D Gaussian to '+string(linewave[ll],format='(G0.0)')
                   xyindx_zoom = array_indices(size(img,/dim),these_zoom,/dim)
                   x1_zoom = min(xyindx_zoom[0,*]) & x2_zoom = max(xyindx_zoom[0,*])
                   y1_zoom = min(xyindx_zoom[1,*]) & y2_zoom = max(xyindx_zoom[1,*])
; setup and call MPFIT
                   parinfo = replicate({value: 1.0, fixed: 0, limited: [0,0], limits: [0.0,0.0]},7)
                   parinfo[1].limited[0] = 1      ; force the Gaussian to be positive
                   parinfo[4:5].limited[0] = 1    ; force (x0,y0) to be positive
                   ymodel_zoom = mpfit2dpeak(img[x1_zoom:x2_zoom,y1_zoom:y2_zoom],est,/tilt,$
                     /positive,/gauss,weights=ivar[x1_zoom:x2_zoom,y1_zoom:y2_zoom],$
                     perror=est_perror,parinfo=parinfo)
; now refit the line over a wider wavelength (velocity) range
                   these = where((ordermask eq ordr_str[indx].order) and (img_arc gt 3.0) and $
                     (img_arc gt alog10((1.0+info[jj].z)*linewave[ll]*(1.0-deltav/light))) and $
                     (img_arc lt alog10((1.0+info[jj].z)*linewave[ll]*(1.0+deltav/light))),nthese)
                   xyindx = array_indices(size(img,/dim),these,/dim)
                   x1 = min(xyindx[0,*]) & x2 = max(xyindx[0,*])
                   y1 = min(xyindx[1,*]) & y2 = max(xyindx[1,*])
                   parinfo.value = est
                   parinfo[4].value = est[4] + ((x1_zoom-x1)>0.0) ; account for the coordinate offset
                   parinfo[5].value = est[5] + ((y1_zoom-y1)>0.0)
                   ymodel = mpfit2dpeak(img[x1:x2,y1:y2],params,estimates=parinfo.value,$
                     parinfo=parinfo,/tilt,/positive,/gauss,weights=ivar[x1:x2,y1:y2],$
                     perror=perror)
                   xpos = params[4] + x1                            ; x-position on the full 2D image
                   ypos = params[5] + y1                            ; y-position on the full 2D image
                   if (params[1]/perror[1] gt snrmin) then begin    ; require a significant detection
                      line2dmodel[x1:x2,y1:y2] = line2dmodel[x1:x2,y1:y2] + (ymodel-params[0])
                   endif else goodfit = 0
                endif else goodfit = 0
; pack the results into a structure
                line2dfit1 = {goodfit: 0, lambda0: linewave[ll], ordr: ordr_str[indx].order, $
                  xpos: -1.0, ypos: -1.0, width: -1.0, height: -1.0, $
                  params: fltarr(7)-1.0, perror: fltarr(7)-1.0}
                if goodfit then begin
                   line2dfit1.goodfit = goodfit
                   line2dfit1.xpos = xpos
                   line2dfit1.ypos = ypos
                   line2dfit1.params = params
                   line2dfit1.perror = perror
; derive the spatial extent (width and height) of the line in the slit 
;                  atv, (ymodel-params[0]) gt 0.001, /bl
                   objpix = where((ymodel-params[0]) gt 0.001) ; 0.1% above the pedestal 
                   xyobjpix = array_indices(size(ymodel,/dim),objpix,/dim)
                   line2dfit1.width = max(xyobjpix[0,*])-min(xyobjpix[0,*])
                   line2dfit1.height = max(xyobjpix[1,*])-min(xyobjpix[1,*])
                endif
                if (n_elements(line2dfit) eq 0L) then line2dfit = line2dfit1 else $
                  line2dfit = [line2dfit,line2dfit1]
             endfor             ; close line loop
; compute the object trace and aperture width 
             ordr_center = (ordr_str[indx].rhedg-ordr_str[indx].lhedg)/2.0+ordr_str[indx].lhedg
             good = where((line2dfit.ordr eq ordr_str[indx].order) and line2dfit.goodfit,ngood)
             if (ngood gt 0L) then begin
; compute the median offset of the object from the center of the order 
                offset_from_center = djs_median(line2dfit[good].xpos-interpol(ordr_center,$
                  findgen(ordr_str[indx].ymax),line2dfit[good].ypos))
                width = djs_median(line2dfit[good].width) ; median object width [pixels]
                trace = ordr_center+offset_from_center
                trace_right = (trace+width/2.0)<ordr_str[indx].rhedg
                trace_left = (trace-width/2.0)>ordr_str[indx].lhedg
; compute the normalized aperture width, which is used by mike_fntobj
; (specifically, x_fntobj); the default values are [0.75,0.75], i.e.,
; 75% of the slit width, relative to the trace
                aper = [djs_median((trace-trace_left)/(trace-ordr_str[indx].lhedg)),$
                  djs_median((trace_right-trace)/(ordr_str[indx].rhedg-trace))]
; pack into a structure
                myobj_str1 = {ordr: ordr_str[indx].order, nline: ngood, $
                  width: width, offset_from_center: offset_from_center, $
                  trace: trace, trace_right: trace_right, $
                  trace_left: trace_left, aper: aper}
                if (n_elements(myobj_str) eq 0L) then myobj_str = myobj_str1 else $
                  myobj_str = [myobj_str,myobj_str1]
             endif
; build a QA-plot of the position of all the lines in the order
;            im_window, 0, xratio=0.7, yratio=0.6
             good = where((line2dfit.ordr eq ordr_str[indx].order) and line2dfit.goodfit,ngood)
             if (ngood ne 0L) then begin
                plot, [0], [0], /nodata, xrange=[min(ordr_str[indx].lhedg),max(ordr_str[indx].rhedg)],$
                  yrange=[ordr_str[indx].ymin,ordr_str[indx].ymax], $
                  xsty=3, ysty=1, title=file_basename(repstr(imgfil,'.fits',''))+'/Order '+$
                  string(ordr_str[indx].order,format='(I0)'), xtitle='X', ytitle='Y'
                djs_oplot, ordr_center, findgen(ordr_str[indx].ymax), line=3, thick=5.0
                djs_oplot, ordr_str[indx].rhedg, findgen(ordr_str[indx].ymax)
                djs_oplot, ordr_str[indx].lhedg, findgen(ordr_str[indx].ymax)
                for gg = 0L, ngood-1L do begin
                   plots, line2dfit[good[gg]].xpos, line2dfit[good[gg]].ypos, psym=4, $
                     sym=2.5, thick=5.0
                   xyouts, line2dfit[good[gg]].xpos*1.01, line2dfit[good[gg]].ypos*1.01, $
                     string(round(line2dfit[good[gg]].lambda0),format='(I0)'), charsize=1.3
                endfor
                djs_oplot, myobj_str1.trace, findgen(ordr_str[indx].ymax), color='red', thick=5.0
                djs_oplot, myobj_str1.trace_right, findgen(ordr_str[indx].ymax), $
                  color='blue', line=5
                djs_oplot, myobj_str1.trace_left, findgen(ordr_str[indx].ymax), $
                  color='blue', line=5
                legend, 'aper = ['+strjoin([string(myobj_str1.aper[0],format='(F5.2)'),$
                  string(myobj_str1.aper[1],format='(F5.2)')],', ')+']', /left, /bottom, $
                  box=0, charsize=1.6
             endif 
          endfor                ; close order loop
;         struct_print, line2dfit
          struct_print, struct_trimtags(myobj_str,except=['*TRACE*'])
; write out 
          splog, 'Writing '+outfil
          mwrfits, line2dmodel, outfil, /create
          mwrfits, line2dfit, outfil
          spawn, 'gzip -f '+outfil, /sh
          splog, 'Writing '+outfil_objstr
          mwrfits, myobj_str, outfil_objstr, /create
          spawn, 'gzip -f '+outfil_objstr, /sh
;         struct_print, struct_trimtags(line2dfit,except=['PERROR'])
          delvarx, line2dfit, myobj_str
       endfor ; close object loop
       dfpsclose
       spawn, 'gzip -f '+qafile
       im_plotfaves
    endif

; ##################################################
; find and trace the object in the slit; OBJAPER is the aperture to
; mask for sky subtraction (26 unbinned pixels when /STD, and 20
; pixels otherwise); FWIDTH is the fraction of the slit width to use
; when tracing (default 0.25 = 1/4 slit width)
    if keyword_set(dotrace) then begin
       for jj = 0L, nobj-1L do begin
          mike_fntobj, mike, setup, obj[jj], side, $
            objaper=objaper, fwidth=fwidth, chk=chk
          obj_strfil = mike_getfil('obj_fil',setup,$
            SUBFIL='Extract/Obj_'+mike[objindx[jj]].img_root,/name)
          obj_str = mike_getfil('obj_fil',setup,$
            SUBFIL='Extract/Obj_'+mike[objindx[jj]].img_root)
; read "my" object structure file written out by /EMLINES and
; overwrite the default trace structure; this modified trace does not
; get used by alpha_mike_skysub because we build and pass an object
; mask 
          myobj_strfil = repstr(obj_strfil,'Obj_','myObj_')+'.gz'
          myobj_str = mrdfits(myobj_strfil,1,/silent)
          for kk = 0L, n_elements(myobj_str)-1L do begin
             indx = where(myobj_str[kk].ordr eq obj_str.order)
             obj_str[indx].trace[0L:sz_img[1]-1L] = myobj_str[kk].trace
             obj_str[indx].aper = myobj_str[kk].aper
          endfor
          mwrfits, obj_str, obj_strfil, /create ; overwrite!
       endfor 
    endif

; ##################################################
; sky-subtract
    if keyword_set(skysub) then begin
       for jj = 0L, nobj-1L do begin
; read and build the emission-line mask
          imgfil = mike_getfil('fin_fil',setup,/name,$
            SUBFIL=mike[objindx[jj]].img_root)
          maskfile = repstr(imgfil,'.fits','_emlines.fits.gz')
          splog, 'Reading '+maskfile
          linemask = mrdfits(maskfile,0,/silent)
          linefit = mrdfits(maskfile,1,/silent)
          objfil = mike_getfil('obj_fil',setup,$
            SUBFIL=mike[objindx[jj]].img_root,/name)
          myobjfil = repstr(objfil,'Obj_','myObj_')+'.gz'
          splog, 'Reading '+objfil
          objstr = xmrdfits(objfil,1,STRUCTYP='mikeobjstrct',/silent)
          splog, 'Reading '+myobjfil
          myobjstr = mrdfits(myobjfil,1,/silent)
          these_ordrs = myobjstr.ordr
          ordr_str = mike_getfil('ordr_str',setup,side=side)
; in order for x_echskysub to work properly, we also have to mask out
; object pixels in every order, not just the emission lines; since we
; typically don't detect continuum, mask out some fraction of the slit
; in all the orders, but, centered on the object; this code is all in
; x_echskysub
          obj_temp = ordr_str
          aper = cmreplicate([0.2,0.2],n_elements(ordr_str)) ; +/-20% of the slit width
          slit_length = obj_temp.rhedg - obj_temp.lhedg
          obj_temp.lhedg = (objstr.trace[0:sz_img[1]-1] - $
            (aper[0,*] ## replicate(1,sz_img[1]))*slit_length/2.0) > ordr_str.lhedg
          obj_temp.rhedg = (objstr.trace[0:sz_img[1]-1] + $
            (aper[1,*] ## replicate(1,sz_img[1]))*slit_length/2.0) < ordr_str.rhedg
          mask_temp_obj = x_ordermask(sz_img[0],sz_img[1],obj_temp,trim=0.0)
; now build the emission-line mask; two possibilities: (1) mask out
; pixels that are 0.1% above the background or (2) mask out the full
; width of the slit at the position of each line to be conservative in
; the sky-subtraction 
          objmask1 = (linemask gt 0.001) ; simple method
;;        ordermask = x_ordermask(sz_img[0],sz_img[1],ordr_str,trim=0.0)
;;        objmask1 = fix(linemask*0.0)
;;        for gg = 0L, n_elements(linefit)-1L do begin
;;           if linefit[gg].goodfit then begin
;;              indx = where(ordr_str.order eq linefit[gg].ordr)
;;              y1 = floor(linefit[gg].ypos-linefit[gg].height/2.0)
;;              y2 = ceil(linefit[gg].ypos+linefit[gg].height/2.0)
;;              yy = transpose(indgen(sz_img[1]) # (intarr(sz_img[0])+1))
;;              these = where((yy ge y1) and (yy le y2) and (ordermask eq ordr_str[indx].order))
;;              objmask1[these] = objmask1[these] or 1
;;              atv, (objmask1 gt 0) or (mask_temp_obj gt 0), /bl
;;           endif
;;        endfor
          objmask = (objmask1 gt 0) or (mask_temp_obj gt 0)
;         atv, objmask, /bl
; now sky-subtract!
          alpha_mike_skysub, mike, setup, obj[jj], side, chk=chk, $
;           ordr=[these_ordrs[0],these_ordrs[n_elements(these_ordrs)-1]], $
            use_objmask=objmask
       endfor
    endif

; ##################################################
; extract 1D spectra
    if keyword_set(extract) then begin
       for jj = 0L, nobj-1L do begin
          splog, 'Performing boxcar extraction'
          objfil = mike_getfil('obj_fil',setup,$
            SUBFIL=mike[objindx[jj]].img_root,/name)
          myobjfil = repstr(objfil,'Obj_','myObj_')+'.gz'
          splog, 'Reading '+myobjfil
          myobjstr = mrdfits(myobjfil,1,/silent)
          these_ordrs = myobjstr.ordr
          nthese_ordrs = n_elements(these_ordrs)
          base_aper = myobjstr[0].aper ; aperture width figured out in /EMLINES
; do not mask cosmic rays as the algorithm in x_extechopt
          mike_box, mike, setup, obj[jj], side, chk=chk, ochk=ochk, reschk=reschk, $
            /boxonly, ordrs=these_ordrs, nohelio=0, novac=0, base_aper=base_aper, $
            /nocrmask ; do not mask cosmic rays (too simplistic!)
; this code will coadd multiple exposures of the same object, but we
; probably want to do this ourselves in post-processing          
       endfor
    endif

; ##################################################
; flux-calibrate
    if keyword_set(calibrate) then begin
       fluxfil = 'Extract/Sens_mr0025.idl' ; not generalized!!
;      fluxfil = getenv('MIKE_DIR')+'/pro/Std/Archive/sens_red'+ $
;        strtrim(mike[0].rowbin,2)+'.fits'
        for jj = 0L, nobj-1L do begin
          mike_flux, mike, setup, obj[jj], side, fluxfil=fluxfil, /boxcar
       endfor
    endif

; ##################################################
; combine multiple exposures of the same object, which must be run on
; even a single object
;   if keyword_set(combine) then begin
;      for jj = 0L, nobj-1L do begin
;         mike_combspec, mike, setup, obj[jj], side, /useboxflux
;      endfor
;   endif
    
; ##################################################
; read the extracted 1D spectra and stitch the orders together;  also
; coadd the 1D spectra of common objects  
    if keyword_set(coadd) then begin
       useflux = 1
       if useflux then begin
          scale = 1D16
          ytitle = textoidl('Flux (10^{-16} '+flam_units()+')')
          suffix = ''
       endif else begin
          scale = 1.0
          ytitle = textoidl('Flux (ADU)')
          suffix = '_unfluxed'
       endelse
; generate a QA-plot
       qafile = datapath+'QA/qa_extract'+suffix+'.ps'
       splog, 'Building QA-plot '+qafile
       dfpsplot, qafile, /landscape, /color
       im_plotfaves, thick=5.0, /post
       for jj = 0L, nobj-1L do begin
          spec1dfile = datapath+'spec1d/obj'+strtrim(info[jj].obj,2)+$
            '.'+info[jj].img_root
          objfil = mike_getfil('obj_fil',setup,$
            SUBFIL=mike[objindx[jj]].img_root,/name)
          myobjfil = repstr(objfil,'Obj_','myObj_')+'.gz'
          splog, 'Reading '+myobjfil
          myobjstr = mrdfits(myobjfil,1,/silent)
          these_ordrs = myobjstr.ordr
          nthese_ordrs = n_elements(these_ordrs)

; read the object structure and write out each order as a FITS file 
          splog, 'Reading '+objfil
          objstr = xmrdfits(objfil,1,STRUCTYP='mikeobjstrct',/silent)

          biglogwave = fltarr(5000,nthese_ordrs)
          bigflux = biglogwave*0.0
          bigivar = biglogwave*0.0
          bigskyflux = biglogwave*0.0
          
          for oo = 0L, nthese_ordrs-1L do begin
             indx = where(objstr.order eq these_ordrs[oo])
             spec1dfile_ordr = repstr(spec1dfile,'.fits','.'+$
               'ord'+string(objstr[indx].order,format='(I2.2)')+'.fits')

             objwave = objstr[indx].wave
             if useflux then begin ; spectra in flux units
                objflux = objstr[indx].flux*1D-16
                objvar = (objstr[indx].sig*1D-16)^2.0
             endif else begin ; spectra in counts
                objflux = objstr[indx].box_fx
                objvar = objstr[indx].box_var
             endelse
             
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

             mkhdr, hdr, float(flux), /extend
             sxdelpar, hdr, 'COMMENT'
             sxdelpar, hdr, 'DATE'
             sxaddpar, hdr, 'OBJECT', info[jj].obj
             sxaddpar, hdr, 'GALAXY', info[jj].galaxy
             sxaddpar, hdr, 'EXPTIME', float(info[jj].exp)
             sxaddpar, hdr, 'AIRMASS', float(info[jj].am)
             sxaddpar, hdr, 'RA', info[jj].ra, format='(F11.7)'
             sxaddpar, hdr, 'DEC', info[jj].dec, format='(F13.9)'
             sxaddpar, hdr, 'Z', info[jj].z
             sxaddpar, hdr, 'RAWFITS', info[jj].img_root
             sxaddpar, hdr, 'CTYPE1', 'LOG'
             sxaddpar, hdr, 'CRPIX1', 1.0d
             sxaddpar, hdr, 'CRVAL1', min(logwave)
             sxaddpar, hdr, 'CDELT1', cdelt
             sxaddpar, hdr, 'DC-FLAG', 1
             
             splog, 'Writing '+spec1dfile_ordr
             mwrfits, float(flux), spec1dfile_ordr, hdr, /create
             mwrfits, float(ivar), spec1dfile_ordr, hdr
             mwrfits, float(skyflux), spec1dfile_ordr, hdr
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
             flux = total(combineflux*combineivar,2)/total(combineivar,2)
             var = 1.0/total(combineivar,2)
             ivar = 1.0/var
             skyflux = total(combineskyflux,2)/float(nthese_ordrs) ; simple mean
          endelse

; write out          
          mkhdr, hdr, float(flux), /extend
          sxdelpar, hdr, 'COMMENT'
          sxdelpar, hdr, 'DATE'
          sxaddpar, hdr, 'OBJECT', info[jj].obj
          sxaddpar, hdr, 'GALAXY', info[jj].galaxy
          sxaddpar, hdr, 'EXPTIME', float(info[jj].exp)
          sxaddpar, hdr, 'AIRMASS', float(info[jj].am)
          sxaddpar, hdr, 'RA', info[jj].ra, format='(F11.7)'
          sxaddpar, hdr, 'DEC', info[jj].dec, format='(F13.9)'
          sxaddpar, hdr, 'Z', info[jj].z
          sxaddpar, hdr, 'RAWFITS', info[jj].img_root
          sxaddpar, hdr, 'CTYPE1', 'LOG'
          sxaddpar, hdr, 'CRPIX1', 1.0d
          sxaddpar, hdr, 'CRVAL1', min(logwave)
          sxaddpar, hdr, 'CDELT1', cdelt
          sxaddpar, hdr, 'DC-FLAG', 1
          
          splog, 'Writing '+spec1dfile
          mwrfits, float(flux), spec1dfile, hdr, /create
          mwrfits, float(ivar), spec1dfile, hdr
          mwrfits, float(skyflux), spec1dfile, hdr

; make the QA-plot          

          xrange = minmax(10^logwave)
          yrange = im_max(scale*flux,sigrej=5.0)*[-0.08,1.2]
          xtitle = textoidl('Observed Wavelength (\AA)')
          title = strtrim(info[jj].galaxy,2)+', obj'+repstr(strtrim(info[jj].obj,2),'_','')+$
            ', z='+strtrim(string(info[jj].z,format='(F12.5)'),2)
          colors = ['red','blue','dark green','purple','orange']
          
          plot, 10.0^logwave, scale*flux, ps=10, /xsty, /ysty, xrange=xrange, $
            yrange=yrange, title=title, xtitle=xtitle, ytitle=ytitle
          legend, 'Combined Spectrum', /left, /top, box=0, charsize=1.6
          for oo = 0L, nthese_ordrs-1L do begin
             good = where(biglogwave[*,oo] gt 0.0)
             if (oo eq 0L) then plot, [0], [0], /nodata, ps=10, /xsty, /ysty, $
               xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, title=title
             djs_oplot, 10.0^biglogwave[good,oo], scale*bigflux[good,oo], ps=10, $
               color=colors[oo]
             legend, 'Individual Orders', /left, /top, box=0, charsize=1.6
          endfor
       endfor
; close the QA-plot
       dfpsclose
       spawn, 'gzip -f '+qafile
       im_plotfaves

; now combine multiple observations of the same object

       qafile = datapath+'QA/qa_coadd'+suffix+'.ps'
       splog, 'Building QA-plot '+qafile
       dfpsplot, qafile, /landscape, /color
       im_plotfaves, thick=5.0, /post
       
       allspec1dfile = datapath+'spec1d/obj'+strtrim(info.obj,2)+$
         '.'+info.img_root
       allobj = info.parent_id
       uobj = allobj[uniq(allobj,sort(allobj))]

       for jj = 0L, n_elements(uobj)-1L do begin
          these = where((uobj[jj] eq allobj),nthese)
          outobj = string(info[these[0]].parent_id,format=('(I3.3)'))
          spec1dfile = datapath+'spec1d/obj'+outobj+'.fits'
          for kk = 0L, nthese-1L do begin
             if (kk eq 0L) then begin
                flux1 = mrdfits(allspec1dfile[these[kk]],0,hdr1,/silent)
                ivar1 = mrdfits(allspec1dfile[these[kk]],1,/silent)
                skyflux1 = mrdfits(allspec1dfile[these[kk]],2,/silent)
                logwave1 = make_wave(hdr1)
                npix = n_elements(logwave1)
; pack into a big array
                combineflux = fltarr(npix,nthese)
                combineskyflux = combineflux*0.0
                combineivar = combineflux*0.0
             endif else begin
                thisflux = mrdfits(allspec1dfile[these[kk]],0,thishdr,/silent)
                thisivar = mrdfits(allspec1dfile[these[kk]],1,/silent)
                thisskyflux = mrdfits(allspec1dfile[these[kk]],2,/silent)
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
             flux = total(combineflux*combineivar,2)/total(combineivar,2)
             var = 1.0/total(combineivar,2)
             ivar = 1.0/var
             skyflux = total(combineskyflux,2)/float(nthese) ; simple mean
          endelse

; write out
          mkhdr, hdr, float(flux), /extend
          sxdelpar, hdr, 'COMMENT'
          sxdelpar, hdr, 'DATE'
          sxaddpar, hdr, 'OBJECT', outobj
          sxaddpar, hdr, 'GALAXY', info[jj].galaxy
          sxaddpar, hdr, 'EXPTIME', float(info[jj].exp)
          sxaddpar, hdr, 'AIRMASS', float(info[jj].am)
          sxaddpar, hdr, 'RA', info[jj].ra, format='(F11.7)'
          sxaddpar, hdr, 'DEC', info[jj].dec, format='(F13.9)'
          sxaddpar, hdr, 'Z', info[jj].z
;         sxaddpar, hdr, 'RAWFITS', info[jj].img_root
          sxaddpar, hdr, 'CTYPE1', 'LOG'
          sxaddpar, hdr, 'CRPIX1', 1.0d
          sxaddpar, hdr, 'CRVAL1', min(logwave1)
          sxaddpar, hdr, 'CDELT1', cdelt
          sxaddpar, hdr, 'DC-FLAG', 1

          splog, 'Writing '+spec1dfile
          mwrfits, float(flux), spec1dfile, hdr, /create
          mwrfits, float(ivar), spec1dfile, hdr
          mwrfits, float(skyflux), spec1dfile, hdr

; make the QA-plot          

          xrange = minmax(10^logwave1)
          yrange = im_max(scale*flux,sigrej=5.0)*[-0.08,1.2]
          xtitle = textoidl('Observed Wavelength (\AA)')
          title = strtrim(info[these[0]].galaxy,2)+', obj'+outobj+$
            ', z='+strtrim(string(info[these[0]].z,format='(F12.5)'),2)
          colors = ['red','blue','dark green','purple','orange']
          
          plot, 10.0^logwave1, scale*flux, ps=10, /xsty, /ysty, xrange=xrange, $
            yrange=yrange, title=title, xtitle=xtitle, ytitle=ytitle
          legend, 'Final 1D Spectrum', /left, /top, box=0, charsize=1.6
          for kk = 0L, nthese-1L do begin
             if (kk eq 0L) then plot, [0], [0], /nodata, ps=10, /xsty, /ysty, $
               xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, title=title
             djs_oplot, 10.0^logwave1, scale*combineflux[*,kk], ps=10, $
               color=colors[kk]
          endfor
          legend, strjoin(repstr(info[these].img_root,'.fits',''),','), /left, $
            /top, box=0, charsize=1.6
       endfor
; close the QA-plot
       dfpsclose
       spawn, 'gzip -f '+qafile
       im_plotfaves
    endif
    
stop
       

return
end
    
;;; there's a bit of code here from x_extechopt to use the edges of the
;;; order as the sky pixels, but since these should all be emission-line
;;; galaxies just mask the brightest and faintest pixels
;;             if keyword_set(doit) then begin
;;                djs_iterstat, img[inorder], sigrej=2.0, invvar=ivar[inorder], mask=skymask
;;                skymask = (skymask eq 0B)                              ; 0=good, 1=bad
;;                skymask = (smooth(float(skymask),10.0,/nan,/edge) gt 0.0) ; grow the mask
;;                skymask = (smooth(float(skymask),10.0,/nan,/edge) gt 0.0)
;;                skymask = skymask eq 0B ; 1=good, 0=bad
;;                sky_repl = bspline_iterfit(ywave,img[inorder],$
;;                  invvar=ivar[inorder]*skymask,everyn=1000,/silent)
;;                skymodel = bspline_valu(ywave,sky_repl)
;;;               maxdist = max(abs(slit_frac))
;;;               edgepix = abs(slit_frac) gt maxdist-0.15
;;;               usethese = where(edgepix,nuse)
;;;               if (nuse gt 50L) then begin 
;;;                  sky_repl = bspline_iterfit(ywave[usethese],img[inorder[usethese]], $
;;;                    invvar=ivar[inorder[usethese]],everyn=50,yfit=skymodel,/silent)
;;;                  skymodel = bspline_valu(ywave,sky_repl)
;;;               endif else message, 'Should not happen'
;;                srt = sort(wave_order)
;;                djs_plot, 10^wave_order[srt], (img[inorder])[srt], ps=3, xsty=3 ;, xr=[5110,5120]
;;                djs_oplot, 10^wave_order[srt], skymodel[srt], color='red'
;;             endif
;;


; this bit of code blocks out a rectangular portion of the spectrum                   
;                     xx = findgen(x2-x1+1) # (fltarr(y2-y1+1)+1)
;                     yy = transpose(findgen(y2-y1+1) # (fltarr(x2-x1+1)+1))
;                     thismask = (xx gt (params[4]-3.0*params[2])) and (xx lt (params[4]+3.0*params[2])) and $
;                       (yy gt (params[5]-3.0*params[3])) and (yy lt (params[5]+3.0*params[3]))
;                     oiii_mask[x1:x2,y1:y2] = oiii_mask[x1:x2,y1:y2] and (thismask eq 0)

;                     oiii_mask[x1:x2,y1:y2] = oiii_mask[x1:x2,y1:y2] and ((ymodel-params[0]) lt 0.001)


; --------------------------------------------------
; test by rectifying and fitting a 2D Gaussian - this works!
;;           rimg = x_ordrectify(img,ordr_str[indx].lhedg,ordr_str[indx].rhedg)
;;           rivar = x_ordrectify(ivar,ordr_str[indx].lhedg,ordr_str[indx].rhedg)
;;           rimg_arc = x_ordrectify(img_arc,ordr_str[indx].lhedg,ordr_str[indx].rhedg)
;;           midcol = (size(rimg,/dim))[0]/2
;;           these = where((rimg_arc[midcol,*] gt alog10((1.0+info[jj].z)*5007.0*(1.0-deltav/3E5))) and $
;;             (rimg_arc[midcol,*] lt alog10((1.0+info[jj].z)*5007.0*(1.0+deltav/3E5))),nthese)
;;           yfit = mpfit2dpeak(rimg[*,these],params,/tilt,/gauss,perror=perror)
;;           mwrfits, rimg[*,these], 'junk.fits', /create
;;           mwrfits, yfit, 'junk2.fits', /create
; --------------------------------------------------
;         line2dmask = line2dmodel lt 0.001

;;    if keyword_set(extract) then begin
;;       velpix = (side eq 1 ? 1.50d : 2.10d) * double(mike[0].rowbin)
;;       cdelt = alog10(1.0d + velpix/light)
;;       for jj = 0L, nobj-1L do begin
;;          imgfil = mike_getfil('fin_fil',setup,/name,$
;;            SUBFIL=mike[objindx[jj]].img_root)
;;          arcfil = strtrim(mike[objindx[jj]].arc_img,2)
;;          myobjfil = repstr(imgfil,'.fits','_line2dmodel.fits.gz')
;;
;;          splog, 'Reading '+imgfil
;;          img = xmrdfits(imgfil,2,/silent) ; sky-subtracted image!
;;          ivar = xmrdfits(imgfil,1,/silent)
;;          splog, 'Reading '+arcfil
;;          img_arc = xmrdfits(arcfil,/silent) 
;;          ordr_str = mike_getfil('ordr_str',setup,side=side)
;;
;;          sz = size(img,/dimen)
;;          ncol = sz[0] & nrow = sz[1]
;;          ycol = dindgen(nrow)
;;          ordermask = x_ordermask(ncol,nrow,ordr_str,trim=0.0)
;;          
;;          objfil = mike_getfil('obj_fil',setup,$
;;            SUBFIL=mike[objindx[jj]].img_root,/name)
;;          obj_str = xmrdfits(objfil,1,STRUCTYP='mikeobjstrct',/silent)
;;          myobj_str = mrdfits(myobjfil,2,/silent)
;;          ordrs = myobj_str.ordr
;;          nordrs = n_elements(ordrs)
;;
;;          for kk = 0L, nordrs-1L do begin
;;             splog, 'Extracting order '+string(ordrs[kk],format='(I0)')
;;             indx = where(ordr_str.order eq ordrs[kk])
;;
;;             inorder = where((ordermask eq ordr_str[indx].order) and (img_arc gt 3.0))
;;             xstart = 1.0d*(inorder mod ncol)
;;             ystart = inorder/ncol
;;
;;; center of the order/slit
;;             slit_cen  = 0.5*(ordr_str[indx].lhedg + ordr_str[indx].rhedg)[ystart]
;;             slit_pos  = xstart - slit_cen
;;             ywave = x_qckwav(slit_pos,ystart,ordr_str[indx].arc_m,arc_slope=arc_slope)
;;; derive the mean wavelength solution along the center of the order
;;             objcen = myobj_str[kk].trace
;;;            objcen = obj_str[indx].trace[0:nrow-1] 
;;             wave_order = img_arc[inorder] ; wavelengths of every pixel in the order
;;             central_wave = extract_boxcar(img_arc,objcen,radius=0.5)
;;             wave_diff = wave_order - central_wave[ystart]
;;             zero_wave = where(wave_order LT 3.0 OR (abs(wave_diff) GT 0.0003))
;;             if (zero_wave[0] ne -1L) then begin
;;                splog, 'Rejecting ', n_elements(zero_wave), ' pixels which have '+$
;;                  'wavelengths from neighboring orders', format ='(A,I6,A)'
;;                ivar[inorder[zero_wave]] = 0
;;                wave_order[zero_wave] = 0
;;             endif
;;             raw_wave_set = bspline_iterfit(ywave,1.0d*wave_order,invvar=1.0d*$
;;               (wave_order gt 3.0),nbkpt=10L,/silent,/double)
;;             dbl_wave_set = {fullbkpt : 1.0d*raw_wave_set.fullbkpt, $
;;               bkmask: raw_wave_set.bkmask, nord: raw_wave_set.nord, $
;;               coeff: 1.0d*raw_wave_set.coeff, icoeff: 1.0d*raw_wave_set.icoeff}
;;             ybox = x_qckwav(objcen-slit_cen,ycol,ordr_str[indx].arc_m)
;;; rectify the 2D image
;;             rimg = x_ordrectify(img,ordr_str[indx].lhedg,ordr_str[indx].rhedg)
;;             rivar = x_ordrectify(ivar,ordr_str[indx].lhedg,ordr_str[indx].rhedg)
;;             rimg_arc = 10.0^x_ordrectify(img_arc,ordr_str[indx].lhedg,ordr_str[indx].rhedg)
;;; extract the 1D spectrum
;;             logwave = bspline_valu(ybox,dbl_wave_set)
;;             flux = extract_asymbox2(img,myobj_str[kk].trace_left,$
;;               myobj_str[kk].trace_right,ybox,weight_imag=ivar,f_ivar=fivar)
;;;            flux = extract_asymbox2(img,ordr_str[indx].lhedg,ordr_str[indx].rhedg,$
;;;              ybox,weight_imag=ivar,f_ivar=fivar)
;;             srt = sort(logwave) & logwave = logwave[srt] ; sort
;;             flux = flux[srt] & fivar = fivar[srt]
;;
;;;            notzero = where((fivar gt 0.0) and (flux ne 0.0))
;;;            logwave = logwave[notzero]
;;;            flux = flux[notzero] & fivar = fivar[notzero]
;;             
;;; write out, ignoring the top and bottom NCROP rows
;;             ncrop = 10L
;;
;;             spec2dfile = datapath+'spec1d/obj'+strtrim(info[jj].obj,2)+'.'+$
;;               'ord'+string(ordr_str[indx].order,format='(I2.2)')+'.2d.fits'
;;
;;             mkhdr, hdr, float(rimg[*,ncrop:nrow-ncrop-1L]), /extend
;;             sxdelpar, hdr, 'COMMENT' & sxdelpar, hdr, 'DATE'
;;             sxaddpar, hdr, 'OBJECT', info[jj].obj
;;             sxaddpar, hdr, 'GALAXY', info[jj].galaxy
;;             sxaddpar, hdr, 'EXPTIME', float(info[jj].exp)
;;             sxaddpar, hdr, 'AIRMASS', float(info[jj].AM)
;;             sxaddpar, hdr, 'RA', info[jj].ra, format='(F11.7)'
;;             sxaddpar, hdr, 'DEC', info[jj].dec, format='(F13.9)'
;;             sxaddpar, hdr, 'Z', info[jj].z
;;             sxaddpar, hdr, 'RAWFITS', info[jj].img_root
;;
;;             splog, 'Writing '+spec2dfile
;;             mwrfits, float(rimg[*,ncrop:nrow-ncrop-1L]), spec2dfile, hdr, /create
;;             mwrfits, float(rivar[*,ncrop:nrow-ncrop-1L]), spec2dfile, hdr
;;             mwrfits, float(rimg_arc[*,ncrop:nrow-ncrop-1L]), spec2dfile, hdr
;;
;;; 1D spectrum, resampled to constant log-lambda (velocity) pixels
;;
;;             spec1dfile = datapath+'spec1d/obj'+strtrim(info[jj].obj,2)+'.'+$
;;               'ord'+string(ordr_str[indx].order,format='(I2.2)')+'.fits'
;;
;;             flux_crop = flux[ncrop:nrow-ncrop-1L]
;;             fivar_crop = fivar[ncrop:nrow-ncrop-1L]
;;             logwave_crop = logwave[ncrop:nrow-ncrop-1L]
;;
;;             logwave_out = dindgen((max(logwave_crop)-min(logwave_crop))/$
;;               cdelt+1.0d)*cdelt+min(logwave_crop)
;;;            combine1fiber, logwave_crop, flux_crop, fivar_crop, $
;;;              newloglam=logwave_out, newflux=flux_out, newivar=fivar_out;, binsz=cdelt
;;             combine1fiber, logwave, flux, fivar, newloglam=logwave_out, $
;;               newflux=flux_out, newivar=fivar_out, binsz=cdelt
;;;            zero = where(fivar_out le 0.0,nzero) & if (nzero ne 0L) then message, 'Fix me'
;;;            ploterror, 10^logwave_out, flux_out, 1.0/sqrt(fivar_out), ps=10, /xsty
;;             
;;             mkhdr, hdr, float(flux_out), /extend
;;             sxdelpar, hdr, 'COMMENT' & sxdelpar, hdr, 'DATE'
;;             sxaddpar, hdr, 'OBJECT', info[jj].obj
;;             sxaddpar, hdr, 'GALAXY', info[jj].galaxy
;;             sxaddpar, hdr, 'EXPTIME', float(info[jj].exp)
;;             sxaddpar, hdr, 'AIRMASS', float(info[jj].AM)
;;             sxaddpar, hdr, 'RA', info[jj].ra, format='(F11.7)'
;;             sxaddpar, hdr, 'DEC', info[jj].dec, format='(F13.9)'
;;             sxaddpar, hdr, 'Z', info[jj].z
;;             sxaddpar, hdr, 'RAWFITS', info[jj].img_root
;;             sxaddpar, hdr, 'CTYPE1', 'LOG'
;;             sxaddpar, hdr, 'CRPIX1', 1.0d
;;             sxaddpar, hdr, 'CRVAL1', min(logwave_out)
;;             sxaddpar, hdr, 'CDELT1', cdelt
;;             sxaddpar, hdr, 'DC-FLAG', 1
;;             
;;             splog, 'Writing '+spec1dfile
;;             mwrfits, float(flux_out), spec1dfile, hdr, /create
;;             mwrfits, float(fivar_out), spec1dfile, hdr
;;;            mwrfits, float(1.0/sqrt(fivar_out)), spec1dfile, hdr
;;
;;; store the results
;;             
;;             obj_str[indx].box_wv = 10.0^logwave_out
;;             obj_str[indx].box_fx = flux_out
;;             obj_str[indx].box_var = (1.0/(fivar_out+(fivar_out eq 0.0)))*(fivar_out ne 0.0)
;;             obj_str[indx].nrow = n_elements(logwave_out)
;;             
;;;            spec1d[kk].logwave = ptr_new(logwave_out)
;;;            spec1d[kk].flux = ptr_new(flux_out)
;;;            spec1d[kk].fivar = ptr_new(fivar_out)
;;;            spec1d[kk].hdr = ptr_new(hdr)
;;;            spec1d[kk].npix = n_elements(logwave_out)
;;             
;;          endfor 
;;          
;;; now optimally combine all the orders into a single 1D spectrum
;;
;;; first we have to call/use x_echcombspec (see mike_combspec)           
;;          spec1d = {mikefspecstrct}
;;          spec1d.nexp = 1
;;          spec1d.texp = mike[objindx[jj]].exp
;;;         spec1d.nexp = nordrs
;;;         spec1d.texp[0:nordrs-1] = mike[objindx[jj]].exp
;;          copy_struct, obj_str[0], spec1d, except=['wave','fx','var','npix']
;;          these_ordrs = [ordrs[0],ordrs[nordrs-1]]
;;
;;          x_echcombspec, obj_str, spec1d, these_ordrs, 0, SIGREJ=sigrej, $
;;            ORDNM=ordnm, MCHK=1, fchk=1, NOFLUX=noflux, /useboxflux
;;          for gg = 0L, n_elements(spec1d.npix)-1L do spec1d.npix[gg] = $
;;            total(spec1d.fx[*,gg] gt 0.0)
;;
;;; define the output wavelength vector
;;
;;          min_logwave = 1E6
;;          max_logwave = 0.0
;;          for iord = 0L, nordrs-1L do begin
;;             x1 = 0L & x2 = spec1d.npix[ordrs[iord]]-1L
;;             min_logwave = min_logwave < alog10(min(spec1d.wave[x1:x2,ordrs[iord]]))
;;             max_logwave = max_logwave > alog10(max(spec1d.wave[x1:x2,ordrs[iord]]))
;;          endfor
;;
;;          final_logwave = dindgen((max_logwave-min_logwave)/cdelt+1.0d)*cdelt+min_logwave
;;          npix = n_elements(final_logwave)
;;
;;          bigflux = fltarr(npix,nordrs)
;;          bigivar = fltarr(npix,nordrs)
;;          for iord = 0L, nordrs-1L do begin
;;             x1 = 0L & x2 = spec1d.npix[ordrs[iord]]-1L
;;             logwave = alog10(spec1d.wave[x1:x2,ordrs[iord]])
;;             flux = spec1d.fx[x1:x2,ordrs[iord]]
;;             var = spec1d.var[x1:x2,ordrs[iord]]
;;             ivar = 1.0/(var+(var eq 0.0))*(var ne 0.0)
;;             linterp, logwave, flux, final_logwave, flux1, missing=0.0
;;             linterp, logwave, ivar, final_logwave, ivar1, missing=0.0
;;             bigflux[*,iord] = flux1
;;             bigivar[*,iord] = ivar1
;;          endfor
;;
;;          if (nordrs eq 1L) then begin
;;             norm = bigivar
;;             final_ferr = 1.0/(norm+(norm eq 0.0))*(norm ne 0.0)
;;             final_flux = bigflux*bigivar*final_ferr
;;          endif else begin
;;             norm = total(bigivar,2)
;;             final_ferr = 1.0/(norm+(norm eq 0.0))*(norm ne 0.0)
;;             final_flux = total(bigflux*bigivar,2)*final_ferr
;;          endelse
;;          
;;          spec1dfile_final = file_search(datapath+'spec1d/obj'+$
;;            strtrim(info[jj].obj,2)+'.ord??'+'.fits',count=nspec)
;;
;;;         ploterror, 10.0^final_logwave, final_flux, final_ferr, ps=10, xsty=3, ysty=3
;;          djs_plot, 10.0^final_logwave, final_flux, ps=10, xsty=3, ysty=3
;;          colors = ['red','blue','green','cyan','purple']
;;          for iord = 0L, nordrs-1L do $
;;            djs_oplot, spec1d.wave[0L:spec1d.npix[ordrs[iord]]-1L,ordrs[iord]], $
;;            spec1d.fx[0L:spec1d.npix[ordrs[iord]]-1L,ordrs[iord]], color=colors[iord]
;;          cc = get_kbrd(1)
;;          
;;;         alpha_mike_box, mike, setup, obj[jj], side, /chk, $
;;;           /ochk, reschk=reschk, redshift=info[jj].z, $
;;;           /novac, /nohelio, objstr=objstr1, /skipskysub
;;       endfor
;;    endif
;;       

;;; ###########################################################################
;;; now extract object spectra! - BELOW HERE IS TESTING!
;;
;;    if keyword_set(extract) then begin
;;
;;       obj_id = 1 ; obj 046_a
;;
;;       mike_allobj, mike, setup, obj_id, side, clobber=clobber, /doproc ; process
;;
;;; OBJAPER: aperture to mask for sky subtraction; FWIDTH: fraction of
;;; width used for radius in tracing (0.25 = 1/4 slit width)
;;       
;;       mike_allobj, mike, setup, obj_id, side, fwidth=0.25, objaper=[0.8,0.8], /clobber, /dofnt ; find and trace
;;
;;;      mike_allobj, mike, setup, obj_id, side, clobber=clobber, /dobox, $ ; extract
;;;        /chk, /debug, /skipskysub, highsnr=200.0, ordrs=[70L,70L], /boxonly, $
;;;        base_aper=[0.9,0.9], reschk=reschk, /nohelio, /novac
;;       mike_allobj, mike, setup, obj_id, side, /dobox, ordrs=70L, /boxonly, /chk, /ochk, /reschk, /debug
;;
;;       mike_1dspec, mike, setup, obj_id, side, /chk
;;       
;;; ---------------------------------------------------------------------------       
;;
;;       obj_id = 8L ; = 061_a
;;
;;; process, find and trace the object and sky-subtract
;;       mike_allobj, mike, setup, obj_id, side, clobber=clobber, /doproc, $
;;         /dofnt, /dosky, checkall=checkall
;;
;;       obj_nm = 'a'
;;       indx = where((mike.type EQ 'OBJ' OR mike.type EQ 'STD') $
;;         AND mike.flg_anly NE 0 AND $
;;         mike.setup EQ setup AND mike.obj_id EQ obj_id $
;;         AND mike.side EQ side,nindx)
;;       
;;;  Read in order structure
;;
;;;      mike_allobj, mike, setup, obj_id, side, clobber=clobber, /dobox, /boxonly, /novac, /nohelio
;;
;;       ordr_str = mike_getfil('ordr_str', setup, side=side)
;;       nordr = n_elements(ordr_str)
;;
;;       exp = lindgen(nindx)
;;       q = 0L ; only one exposure
;;
;;       xyoff = mike[indx[exp[q]]].arc_xyoff
;;
;;       ;; Shift ORDR_STR
;;       ordr_shift = ordr_str
;;       shft = mike_shifti(xyoff, OSTR=ordr_shift)
;;
;;; read the relevant files       
;;       objfil = mike[indx[exp[q]]].obj_fil
;;       objstr = xmrdfits(objfil, 1, STRUCTYP='mikeobjstrct', /silent)
;;       
;;       arc_img = strtrim(mike[indx[exp[q]]].arc_img,2)
;;       img_arc = xmrdfits(arc_img, /silent) 
;;
;;       imgfil = mike_getfil('fin_fil',setup,/name,$
;;         SUBFIL=mike[indx[exp[q]]].img_root)
;;       head = xheadfits(imgfil)
;;;      img = xmrdfits(imgfil, 0, /silent)    ; *not* sky-subtracted
;;       img = xmrdfits(imgfil, 2, /silent) ; sky-subtracted
;;       ivar = xmrdfits(imgfil, 1, /silent)
;;       var = 1/(ivar + (ivar EQ 0)) * (ivar GT 0)
;;; detect CRs and interpolate
;;;      ila_cosmic, img, gain=sxpar(head,'EGAIN'), masklist=crmask ; 0=good, 1=bad
;;;      cimg = djs_maskinterp(img,crmask,iaxis=iaxis)
;;;      civar = djs_maskinterp(ivar,crmask,iaxis=iaxis)
;;;      cvar = 1/(civar + (civar EQ 0)) * (civar GT 0)
;;
;;; now pull out the piece I want (this is from x_extechopt)
;;
;;       base_aper=[0.75,0.75]
;;       
;;       slit_length = ordr_shift.rhedg - ordr_shift.lhedg
;;       half_length = slit_length/2.
;;       
;;       sz = size(img,/dimen)
;;       objcen = objstr.trace[0:sz[1]-1]  
;;       
;;       left_edge  = (objcen - base_aper[0] * half_length) > ordr_shift.lhedg
;;       right_edge = (objcen + base_aper[1] * half_length) < ordr_shift.rhedg
;;       
;;       fx = extract_asymbox2(img,left_edge,right_edge)
;;       fvar  = extract_asymbox2(var, left_edge, right_edge)
;;
;;; Create the order mask (img with val = order # or -order# for gaps)
;;       ordermask = x_ordermask(sz[0], sz[1], ordr_shift, trim=msktrim)
;;       
;;       
;;stop       
;;       
;;; do everything except flux-calibrate and call MIKE_1DSPEC (note that
;;; I could use /PROCALL, but I wanted to be more explicit about each
;;; step)
;;       mike_allobj, mike, setup, obj_id, side, clobber=clobber, /doproc, $
;;         /dofnt, /dosky, /dobox, doflux=0, $
;;         /nocr, /docomb, do1d=0, checkall=checkall, /boxonly, /usebox, $
;;         noflux=noflux, /novac, /nohelio
;;
;;stop       
;;       
;;; everything except process and sky-subtraction
;;       mike_allobj, mike, setup, obj_id, side, /clobber, doproc=0, /nocr, /dofnt, dosky=0, $
;;         /dobox, /docomb, /do1d, checkall=checkall, ordrs=52, /boxonly
;;; box-car extraction
;;       mike_allobj, mike, setup, obj_id, side, /dobox, ordrs=52, /boxonly
;;; combine
;;       mike_allobj, mike, setup, obj_id, side, /docomb, ordrs=52, /usebox
;;; make the 1D spectrum; KEEPALL; do not zero out bad blaze areas  ; CRAP!
;;;      mike_allobj, mike, setup, obj_id, side, /do1d, ordrs=52, /keepall;, /chk
;;
;;; plot it
;;       obj_nm = 'a'
;;       indx = where((mike.type EQ 'OBJ' OR mike.type EQ 'STD') $
;;         AND mike.flg_anly NE 0 AND $
;;         mike.setup EQ setup AND mike.obj_id EQ obj_id $
;;         AND mike.side EQ side)
;;       subfil = strcompress(strtrim(mike[indx[0]].Obj,2),/remove_all)+obj_nm
;;
;;       if (side EQ 1) then clrc = '_b' else clrc = '_r'
;;       infil = 'FSpec/'+subfil+clrc+'.fits'
;;       print, 'mike_1dspec: Reading ', infil
;;       spec = xmrdfits(infil,1)
;;
;;       use = where(spec.phys_ordr NE 0,nuse)
;;       for ii = 0L, nuse-1L do begin & $
;;         ww = where(spec.wave[*,use[ii]] gt 0.0) & $
;;         x_splot, spec.wave[ww,use[ii]], spec.fx[ww,use[ii]], psym3=10, /block & $
;;         endfor
;;          
;;; to turn off bias subtraction call mike_proc with /skipov; in general
;;; we shouldn't use mike_allobj because not all the keywords are
;;; propagated (e.g., /skipov)
;;
;;;   for jj = 0L, nobj-1L do mike_allobj, mike, setup, obj[jj], side, clobber=clobber, $
;;;     /doproc, /dofnt, /dosky, fwidth=fwidth, objaper=objaper, checkall=checkall
;;
;;       
;;       
;;
;;stop         
;;
;;       obj_id = 1
;;       mike_allobj, mike, setup, obj_id, side, clobber=clobber, /dofnt ; find and trace
;;
;;
;;       obj_id = 5L
;;       
;;       mike_allobj, mike, setup, obj_id, side, clobber=clobber, /dofnt ; find and trace
;;       mike_allobj, mike, setup, obj_id, side, clobber=clobber, /dosky, /chk, /fchk, ord=70 ; sky-subtract
;;       mike_allobj, mike, setup, obj_id, side, clobber=clobber, /dobox, $ ; extract
;;         /chk, /debug, /skipskysub, highsnr=200.0, ordr=[70,70], /boxonly, $
;;         base_aper=[0.9,0.9], reschk=reschk, /nohelio, /novac
;;
;;       mike_box, mike, setup, obj_id, side, ordrs=66, /ochk, /chk, /reschk, $
;;         base_aper=[0.20,0.20]
;;    endif
;;
;;    if keyword_set(bobo) then begin
;;
;;       velpix = (side eq 1 ? 1.50d : 2.10d) * double(mike[0].rowbin)
;;       cdelt = alog10(1.0d + velpix/light)
;;
;;          sz = size(img,/dimen)
;;          ncol = sz[0] & nrow = sz[1]
;;          ycol = dindgen(nrow)
;;
;;; read and subtract the scattered light image
;;          scattered_light = mike_getfil('scatt_fil',setup, $
;;            subfil=mike[objindx[jj]].img_root)
;;;         img_sub = img ; no sky-subtraction
;;          img_sub = img - scattered_light
;;
;;       for jj = 0L, nobj-1L do begin
;;    
;;stop          
;;          
;;          mike_skysub, mike, setup, obj[jj], side, chk=chk
;;
;;; now sky-subtract!          
;;          for kk = 0L, nordrs-1L do begin
;;             splog, 'Sky-subtracting and extracting order '+string(ordrs[kk],format='(I0)')
;;             indx = where(ordr_str.order eq ordrs[kk])
;;; preliminaries; code largely taken from x_echskysub
;;             inorder = where((ordermask eq ordr_str[indx].order) and (img_arc gt 3.0))
;;             wave_sort = inorder[sort(img_arc[inorder])]
;;             xstart = 1.0d*(wave_sort mod ncol)
;;             ystart = wave_sort/ncol
;;             slit_length = ordr_str[indx].rhedg[ystart] - ordr_str[indx].lhedg[ystart]
;;; center of the order/slit
;;             slit_cen  = 0.5*(ordr_str[indx].lhedg + ordr_str[indx].rhedg)[ystart]
;;             slit_pos  = xstart - slit_cen
;;             ywave = x_qckwav(slit_pos,ystart,ordr_str[indx].arc_m,arc_slope=arc_slope)
;;; position along the slit [-1,1] and slit profile; do we need to
;;; normalize by the slit profile???
;;             slit_frac = frac_order(ordr_str[indx],xstart,ywave) 
;;             profile = x_slitprofile_return(slit_frac,ystart,ordr_str[indx])
;;
;;             img_norm = (img_sub)[wave_sort]/(profile+(profile eq 0))*(profile gt 0)
;;             img_norm_ivar = ivar[wave_sort]*profile^2*(profile gt 0)
;;             wave_norm =  img_arc[wave_sort]
;;             oiii_mask_norm =  oiii_mask[wave_sort]
;;
;;             sky_pix = ((ordermask[wave_sort] lt 200) and (profile gt 0.3))
;;             sky = where(sky_pix,nsky_pix)
;;
;;             evn = (1.2 * nsky_pix / (max(ystart) - min(ystart) + 1)) + 1 
;;
;;             bset_prof = bspline_iterfit(wave_norm[sky],img_norm[sky],everyn=evn,$
;;               nord=nord,upper=5.0,lower=5.0,invvar=img_norm_ivar[sky]*oiii_mask_norm[sky], $
;;               /groupbadpix,maxrej=1,bkpt=bkpt_s,maxiter=15L,outmask=pixmsk,yfit=yfit,/silent)
;;
;;             img_fit = bspline_valu(wave_norm,bset_prof)
;;             djs_plot, 10^wave_norm, img_norm, ps=10, xsty=3, ysty=3
;;             djs_oplot, 10^wave_norm, img_fit, color='red', ps=10
;;
;;             img_nosky = img_sub
;;             img_nosky[wave_norm] = img_nosky[wave_norm] - img_fit*profile
;;             
;;;            mwrfits, img_nosky, 'junk2.fits', /create
;;
;;; derive the mean wavelength solution along the center of the order
;;             objfil = mike_getfil('obj_fil',setup,$
;;               SUBFIL=mike[objindx[jj]].img_root,/name)
;;             obj_str = xmrdfits(objfil,1,STRUCTYP='mikeobjstrct',/silent)
;;             objcen = obj_str[indx].trace[0:nrow-1] 
;;             wave_order = img_arc[inorder] ; wavelengths of every pixel in the order
;;             skymask = oiii_mask[inorder]  ; mask out the [O III] lines!
;;             central_wave = extract_boxcar(img_arc,objcen,radius=0.5)
;;             wave_diff = wave_order - central_wave[ystart]
;;             zero_wave = where(wave_order LT 3.0 OR (abs(wave_diff) GT 0.0003))
;;             if (zero_wave[0] ne -1L) then begin
;;                splog, 'Rejecting ', n_elements(zero_wave), ' pixels which have '+$
;;                  'wavelengths from neighboring orders', format ='(A,I6,A)'
;;                ivar[inorder[zero_wave]] = 0
;;                wave_order[zero_wave] = 0
;;             endif
;;             raw_wave_set = bspline_iterfit(ywave,1.0d*wave_order,invvar=1.0d*$
;;               (wave_order gt 3.0),nbkpt=10L,/silent,/double)
;;             dbl_wave_set = {fullbkpt : 1.0d*raw_wave_set.fullbkpt, $
;;               bkmask: raw_wave_set.bkmask, nord: raw_wave_set.nord, $
;;               coeff: 1.0d*raw_wave_set.coeff, icoeff: 1.0d*raw_wave_set.icoeff}
;;             ybox = x_qckwav(objcen-slit_cen,ycol,ordr_str[indx].arc_m)
;;
;;             sky_repl = bspline_iterfit(ywave,img[inorder],$
;;               invvar=ivar[inorder]*skymask,everyn=1000,/silent)
;;             skymodel = bspline_valu(ywave,sky_repl)
;;
;;             srt = sort(wave_order)
;;             djs_plot, 10^wave_order[srt], (img[inorder])[srt], ps=3, xsty=3 ;, xr=[5110,5120]
;;             ww = where(skymask eq 0)
;;             djs_oplot, 10^wave_order[ww], skymodel[ww], color='cyan', ps=4, sym=2.0
;;             djs_oplot, 10^wave_order[srt], skymodel[srt], color='red'
;;
;;          endfor 
;;
;;       endfor 
;;
;;    endif
;;
