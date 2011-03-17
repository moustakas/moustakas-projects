pro ediscs_stitch_sensitivity, leftsens=leftsens, centersens=centersens, $
  rightsens=rightsens, datapath=datapath, outpath=outpath, suffix=suffix, $
  doplot=doplot, wfits=wfits, gzip=gzip, _extra=extra
; jm04sep12uofa
; jm04nov15uofa - only stitch together the grey-shifted sensitivity
;                 function 

    nleftsens = n_elements(leftsens)
    ncentersens = n_elements(centersens)
    nrightsens = n_elements(rightsens)
    
    if (nleftsens eq 0L) or (ncentersens eq 0L) or (nrightsens eq 0L) then begin
       print, 'Syntax - ediscs_stitch_sensitivity'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(outpath) eq 0L) then outpath = datapath
    if (n_elements(suffix) eq 0L) then suffix = '' else suffix = '_'+suffix

    psname = 'qaplot_sensitivity'+suffix+'.ps'
    sensname = 'sensitivity'+suffix+'.fits'
;   greysensname = 'sensitivity'+suffix+'_grey.fits'
    title = 'Sensitivity Function '+repstr(suffix,'_','')
    
    if keyword_set(wfits) then begin
       doplot = 0L
       postthick = 5.0 
    endif else postthick = 2.0
    
; read the three sensitivity functions

    left = irdsensfunc(leftsens,datapath=datapath)
    center = irdsensfunc(centersens,datapath=datapath)
    right = irdsensfunc(rightsens,datapath=datapath)

; create the output wavelength vector

    minwave = min(right.wave)
    maxwave = max(left.wave)
    dwave = 1.0
    senswave = findgen((maxwave-minwave)/dwave+1)*dwave+minwave
    midwave = djs_mean(senswave)

    sens = senswave*0.0
    
; figure out the overlap regions

    wright = where((right.wave gt min(center.wave)) and $
      (right.wave lt max(center.wave)),nwright)

    wleft = where((left.wave gt min(center.wave)) and $
      (left.wave lt max(center.wave)),nwleft)

; shift to the center; NOT GENERALIZED!!

    get_element, center.wave, midwave+[-25.0,+25.0], xx
    center_meansens = djs_mean(center.greysens[xx[0]:xx[1]])
    
    get_element, right.wave, midwave+[-25.0,+25.0], xx
    right_meansens = djs_mean(right.greysens[xx[0]:xx[1]])
    
    get_element, left.wave, midwave+[-25.0,+25.0], xx
    left_meansens = djs_mean(left.greysens[xx[0]:xx[1]])

    right_greyshift = center_meansens - right_meansens
    left_greyshift = center_meansens - left_meansens

; apply the grey shifts, averaging the sensitivity in the overlap
; region; finally, interpolate the final sensivity function

    linterp, right.wave, right.greysens+right_greyshift, senswave, right_temp, missing=0.0
    linterp, center.wave, center.greysens, senswave, center_temp, missing=0.0
    linterp, left.wave, left.greysens+left_greyshift, senswave, left_temp, missing=0.0

    bigsens = [ [right_temp], [center_temp], [left_temp] ]
    for ipix = 0L, n_elements(senswave)-1L do begin
       avg = where(bigsens[ipix,*] ne 0.0,navg)
       sens[ipix] = total(bigsens[ipix,*],2)/navg
    endfor

; generate the sensitivity function header

    mkhdr, senshead, sens
    sxaddpar, senshead, 'OBJECT', 'Sensitivity Function'
    sxaddpar, senshead, 'CTYPE1', 'LINEAR', ' projection type'
    sxaddpar, senshead, 'CRVAL1', minwave, ' wavelength at the reference pixel'
    sxaddpar, senshead, 'CRPIX1', 1.0, ' reference pixel number'
    sxaddpar, senshead, 'CDELT1', dwave, ' dispersion at the reference pixel'
    sxaddpar, senshead, 'CD1_1', dwave, ' dispersion in Angstroms per pixel'
    sxaddpar, senshead, 'LEFTSHFT', float(left_greyshift), ' left greyshift [mag]'
    sxaddpar, senshead, 'RGHTSHFT', float(right_greyshift), ' right greyshift [mag]'

    if keyword_set(wfits) then begin
       splog, 'Writing '+outpath+sensname+'.'
       mwrfits, sens, outpath+sensname, senshead, /create
;      mwrfits, sens, outpath+greysensname, senshead, /create
    endif
    
; generate the qaplots

    if keyword_set(doplot) or keyword_set(wfits) then begin

       if keyword_set(wfits) then dfpsplot, outpath+psname, /color, /portrait
       if keyword_set(doplot) then window, 0, xs=500, ys=500

       arm_plotconfig, ny=2, yspace=0.0, coords=pos, xmargin=[1.5,0.25], $
         ymargin=[0.5,1.1]
       
; compare the original sensitivity functions       
       
       xrange = minmax(senswave)
       yrange = [(min(left.greysens)<min(center.greysens))<min(right.greysens),$
         (max(left.greysens)>max(center.greysens))>max(right.greysens)]+[-0.1,+0.2]

       plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
         xthick=postthick, ythick=postthick, charsize=1.5, charthick=postthick, $
         xsty=3, ysty=3, xtickname=replicate(' ',10), position=pos[*,0], $
         title=title
       djs_oplot, right.wave, right.greysens, thick=postthick, color='blue', line=1
       djs_oplot, center.wave, center.greysens, thick=postthick, color='purple', line=0
       djs_oplot, left.wave, left.greysens, thick=postthick, color='red', line=2

       legend, ['Right','Center','Left'], line=[1,0,2], /right, /top, $
         box=0, color=djs_icolor(['blue','purple','red']), charsize=1.5, $
         charthick=2.0, thick=postthick
       
; compare the greyshifted sensitivity functions       
       
       yrange = minmax(sens)

       djs_plot, [0], [0], /nodata, /noerase, xrange=xrange, yrange=yrange, $
         xthick=postthick, ythick=postthick, charsize=1.5, charthick=postthick, $
         xsty=3, ysty=3, xtitle='Wavelength [\AA]', $       
         position=pos[*,1]
       djs_oplot, senswave, sens, thick=postthick, color='green'

       xyouts, pos[0,1]-0.1, pos[3,1], textoidl('2.5 log [Counts s^{-1} '+$
         '\AA^{-1}] / ['+flam_units()+']'), /normal, $
         align=0.5, orientation=90, charsize=1.5, charthick=postthick
       
       legend, 'Mean Function', line=0, /right, /top, $
         box=0, color=djs_icolor('green'), charsize=1.5, $
         charthick=2.0, thick=postthick

       if keyword_set(wfits) then dfpsclose

    endif

return
end
