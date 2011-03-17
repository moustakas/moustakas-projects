pro alpha_trace_emlines, info, mike, qafile=qafile
; jm09jan06nyu - identify the orders of interest (i.e., those
; containing the emission lines)

    nobj = n_elements(info)
    if (nobj ne n_elements(mike)) then begin
       splog, 'Dimensions of INFO and MIKE do not match'
       return
    endif

; some parameters

    linewave = [4861.0,4958.911,5006.843] ; rest wavelength [A]
    light = 299792.458D ; [km/s]
    deltav = 300.0      ; [km/s]
    deltav_zoom = 50.0  ; [km/s]
    snrmin = 4.0        ; minimum S/N for a line to be masked
    
; loop on each object       

    if (n_elements(qafile) ne 0L) then begin
       splog, 'Building QA-plot '+qafile
       dfpsplot, qafile, /landscape, /color
       im_plotfaves, thick=5.0, /post
    endif

    for jj = 0L, nobj-1L do begin
; some input/output file names
       setup = mike[jj].setup & side = mike[jj].side
       imgfil = mike_getfil('fin_fil',setup,/name,$
         SUBFIL=mike[jj].img_root)
       outfil = repstr(imgfil,'.fits','_emlines.fits')
       outfil_objstr = mike_getfil('obj_fil',/name,$ 
         subfil='Extract/Obj_'+mike[jj].img_root)
       outfil_objstr = repstr(outfil_objstr,'Obj_','myObj_')
       arcfil = strtrim(mike[jj].arc_img,2)
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
                img_zoom = img[x1_zoom:x2_zoom,y1_zoom:y2_zoom]
                ivar_zoom = ivar[x1_zoom:x2_zoom,y1_zoom:y2_zoom]
                sz_img_zoom = size(img_zoom,/dim)
                parinfo = replicate({value: 0.0, fixed: 0, limited: [0,0], limits: [0.0,0.0]},7)
;               parinfo[0].value = 1E-3
;               parinfo[1].value = 5.0
;               parinfo[1].limited[0] = 1      ; force the Gaussian to be positive
;               parinfo[2:3].value = 0.1
;               parinfo[2:3].value = 0.1
                parinfo[4].value = sz_img_zoom[0]/2
                parinfo[5].value = sz_img_zoom[1]/2
                parinfo[4:5].limited[0] = 1    ; force (x0,y0) to be in the image
;               parinfo[4].limits[1] = sz_img_zoom[0]-1.0
;               parinfo[5].limits[1] = sz_img_zoom[1]-1.0
;               parinfo[6].value = 0.0
                ymodel_zoom = mpfit2dpeak(img_zoom,est,/tilt,/positive,parinfo=parinfo,$
                  /gauss,weights=ivar_zoom,perror=est_perror,status=status_zoom)
; now refit the line over a wider wavelength (velocity) range
                these = where((ordermask eq ordr_str[indx].order) and (img_arc gt 3.0) and $
                  (img_arc gt alog10((1.0+info[jj].z)*linewave[ll]*(1.0-deltav/light))) and $
                  (img_arc lt alog10((1.0+info[jj].z)*linewave[ll]*(1.0+deltav/light))),nthese)
                xyindx = array_indices(size(img,/dim),these,/dim)
                x1 = min(xyindx[0,*]) & x2 = max(xyindx[0,*])
                y1 = min(xyindx[1,*]) & y2 = max(xyindx[1,*])
                parinfo.value = est
;               parinfo[1:2].value = abs(parinfo[1:2].value) ; sometimes sigma is negative
                parinfo[4].value = est[4] + ((x1_zoom-x1)>0.0) ; account for the coordinate offset
                parinfo[5].value = est[5] + ((y1_zoom-y1)>0.0)
                ymodel = mpfit2dpeak(img[x1:x2,y1:y2],params,estimates=parinfo.value,$
                  parinfo=parinfo,/tilt,/positive,/gauss,weights=ivar[x1:x2,y1:y2],$
                  perror=perror,status=status)
                xpos = params[4] + x1                            ; x-position on the full 2D image
                ypos = params[5] + y1                            ; y-position on the full 2D image
; require a significant detection                
                if (params[1]/(perror[1]+(perror[1] eq 0.0))*(perror[1] ne 0.0) gt snrmin) then begin 
                   line2dmodel[x1:x2,y1:y2] = line2dmodel[x1:x2,y1:y2] + (ymodel-params[0])
                endif else goodfit = 0
             endif else goodfit = 0
; pack the results into a structure
             line2dfit1 = {goodfit: 0, lambda0: linewave[ll], ordr: ordr_str[indx].order, $
               xpos: -1.0, ypos: -1.0, width: -1.0, height: -1.0, $
               params: fltarr(7)-1.0, perror: fltarr(7)-1.0}
             if goodfit then begin
;               atv, (ymodel-params[0]) gt 0.001, /bl
                objpix = where((ymodel-params[0]) gt 0.001,nobjpix) ; 0.1% above the pedestal 
                if (nobjpix ne 0L) then begin
                   xyobjpix = array_indices(size(ymodel,/dim),objpix,/dim)
                   line2dfit1.width = max(xyobjpix[0,*])-min(xyobjpix[0,*])
                   line2dfit1.height = max(xyobjpix[1,*])-min(xyobjpix[1,*])
                   line2dfit1.goodfit = goodfit
                   line2dfit1.xpos = xpos
                   line2dfit1.ypos = ypos
                   line2dfit1.params = params
                   line2dfit1.perror = perror
                endif else stop
; derive the spatial extent (width and height) of the line in the slit 
             endif
             if (n_elements(line2dfit) eq 0L) then line2dfit = line2dfit1 else $
               line2dfit = [line2dfit,line2dfit1]
          endfor                ; close line loop
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
       endfor                   ; close order loop
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
    endfor                      ; close object loop

    if (n_elements(qafile) ne 0L) then begin
       dfpsclose
       spawn, 'gzip -f '+qafile
       im_plotfaves
    endif

return
end
    
