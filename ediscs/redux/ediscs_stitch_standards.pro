pro ediscs_stitch_standards, leftstdlist=leftstdlist, centerstdlist=centerstdlist, $
  rightstdlist=rightstdlist, datapath=datapath, outpath=outpath, suffix=suffix, $
  extfile=extfile, aperture=aperture, doplot=doplot, wfits=wfits, gzip=gzip, $
  _extra=extra
; jm04nov19uofa

    normwave = 6780.0 ; normalization wavelength
    normwidth = 5.0   ; average within NORMWAVE+/-NORMWIDTH
    
    nleftstdlist = n_elements(leftstdlist)
    ncenterstdlist = n_elements(centerstdlist)
    nrightstdlist = n_elements(rightstdlist)
    
    if (nleftstdlist eq 0L) or (ncenterstdlist eq 0L) or (nrightstdlist eq 0L) then begin
       print, 'Syntax - ediscs_stitch_sensitivity'
       return
    endif

    if (n_elements(datapath) eq 0L) then datapath = cwd()
    if (n_elements(outpath) eq 0L) then outpath = datapath
    if (n_elements(suffix) eq 0L) then suffix = '' else suffix = '_'+suffix
    if (n_elements(extfile) eq 0L) then extfile = 'lasillaextinct.dat'
    if (n_elements(aperture) eq 0L) then aperture = 10.0

; read in the extinction curve

    extpath = getenv('ISPEC_DIR')+'/etc/'
    if (file_test(extpath+extfile,/regular) eq 0L) then begin
       splog, 'Unable to find extinction file '+extfile+'.'
       return
    endif
          
    splog, 'Reading the extinction file '+extpath+extfile+'.'
    readcol, extpath+extfile, extwave, extvals, format='D,D', /silent

; read the standard-star lists

    readcol, leftstdlist, leftlist, format='A', /silent, comment='#'
    readcol, centerstdlist, centerlist, format='A', /silent, comment='#'
    readcol, rightstdlist, rightlist, format='A', /silent, comment='#'

    nstd = n_elements(centerlist)

    if (n_elements(leftlist) ne nstd) and (n_elements(rightlist) ne nstd) then begin
       splog, 'LEFTSTDLIST, CENTERSTDLIST, and RIGHTSTDLIST have different numbers of stars.'
       return
    endif

; read the un-normalized and the normalized one-dimensional spectra    
    
    left = rd1dspec(i1dnames(leftlist,aperture=aperture),datapath=datapath)
    center = rd1dspec(i1dnames(centerlist,aperture=aperture),datapath=datapath)
    right = rd1dspec(i1dnames(rightlist,aperture=aperture),datapath=datapath)

    dleft = rd1dspec('d'+i1dnames(leftlist,aperture=aperture),datapath=datapath)
    dcenter = rd1dspec('d'+i1dnames(centerlist,aperture=aperture),datapath=datapath)
    dright = rd1dspec('d'+i1dnames(rightlist,aperture=aperture),datapath=datapath)

    leftinfo = iforage(left.specname,datapath=datapath)
    centerinfo = iforage(center.specname,datapath=datapath)
    rightinfo = iforage(right.specname,datapath=datapath)

    psname = 'qaplot_stitch_standards'+suffix+'.ps'
    
    if keyword_set(wfits) then begin
       dfpsplot, outpath+psname, /color, /landscape       
       postthick = 5.0
       doplot = 0L
    endif else begin
       window, xsize=700, ysize=600
       postthick = 2.0
    endelse
    
    pagemaker, nx=1, ny=3, yspace=0, position=pos, xmargin=[1.4,0.3], $
      ymargin=[0.9,1.4], /normal

; loop on each standard     
    
    for istd = 0L, nstd-1L do begin

       starname = strtrim(repstr(repstr(repstr(dcenter[istd].object,'center',''),'left',''),'right',''),2)
       filename = strtrim(repstr(repstr(dcenter[istd].specname,'_center',''),'.gz',''),2)

       if (istd eq 0L) then title = 'MOS Standard Stars - '+repstr(suffix,'_','') else title = ''

; create the output wavelength vector, tied to CENTER

       minwave = min(right[istd].wave)
       maxwave = max(left[istd].wave)
       dwave = arm_double(centerinfo[istd].cd1_1)
       npix = fix((maxwave-minwave)/dwave)+1L

       wave = dindgen(npix)*dwave+minwave ; output wavelength vector
       spec = wave*0.0                    ; output standard-star spectrum
       
; scale the counts in the LEFT and RIGHT spectra to the exposure time
; and airmass of the CENTER spectrum

       leftscale = (centerinfo[istd].exptime / leftinfo[istd].exptime) * $
         10D0^(0.4*interpol(extvals,extwave,left[istd].wave)*(leftinfo[istd].airmass-centerinfo[istd].airmass))
       
       left[istd].spec = left[istd].spec * leftscale
       dleft[istd].spec = dleft[istd].spec * leftscale

       rightscale = (centerinfo[istd].exptime / rightinfo[istd].exptime) * $
         10D0^(0.4*interpol(extvals,extwave,right[istd].wave)*(rightinfo[istd].airmass-centerinfo[istd].airmass))
       
       right[istd].spec = right[istd].spec * rightscale
       dright[istd].spec = dright[istd].spec * rightscale
       
; calculate the shifts required about NORMWAVE+/-NORMWIDTH
       
       get_element, dcenter[istd].wave, normwave+normwidth*[-1,1], xx
       center_meanflux = djs_mean(dcenter[istd].spec[xx[0]:xx[1]])
       
       get_element, dright[istd].wave, normwave+normwidth*[-1,1], xx
       right_meanflux = djs_mean(dright[istd].spec[xx[0]:xx[1]])
       
       get_element, dleft[istd].wave, normwave+normwidth*[-1,1], xx
       left_meanflux = djs_mean(dleft[istd].spec[xx[0]:xx[1]])

       right_shift = center_meanflux - right_meanflux
       left_shift = center_meanflux - left_meanflux

; interpolate all the spectra onto the dispersion of the CENTER
; spectrum and a common wavelength grid, simultaneously applying the
; shifts 

       linterp, dright[istd].wave, dright[istd].spec+right_shift, wave, rightspec, missing=0.0D
       linterp, dcenter[istd].wave, dcenter[istd].spec+0.0D, wave, centerspec, missing=0.0D
       linterp, dleft[istd].wave, dleft[istd].spec+left_shift, wave, leftspec, missing=0.0D

; finally, calculate the mean spectrum

       bigspec = [ [rightspec], [centerspec], [leftspec] ]
       for ipix = 0L, npix-1L do begin
          avg = where(bigspec[ipix,*] ne 0.0,navg)
          spec[ipix] = total(bigspec[ipix,*],2)/navg
       endfor

; generate the FITS header and write out the final spectrum, using
; CENTER as the base header and file name

       if (nstd gt 1L) then header = *dcenter[istd].header else header = dcenter[istd].header
       sxaddpar, header, 'NAXIS1', npix       
       sxaddpar, header, 'OBJECT', starname
       sxaddpar, header, 'CRVAL1', minwave, ' wavelength at the reference pixel'
       sxaddpar, header, 'CRPIX1', 1.0, ' reference pixel number'
       sxaddpar, header, 'CDELT1', dwave, ' dispersion at the reference pixel'
       sxaddpar, header, 'CD1_1', dwave, ' dispersion in Angstroms per pixel'
       sxaddpar, header, 'CTYPE1', "'LINEAR'", ' projection type'

;      sxaddhist, header, 'LEFTSHFT', float(left_greyshift), ' left greyshift [mag]'
;      sxaddhist, header, 'RGHTSHFT', float(right_greyshift), ' right greyshift [mag]'

       if keyword_set(wfits) then begin
          splog, 'Writing '+outpath+filename+'.'
          mwrfits, spec, outpath+filename, header, /create
          if keyword_set(gzip) then spawn, ['gzip -f '+outpath+filename], /sh
       endif
       
; generate the QA plot       

       if keyword_set(doplot) or keyword_set(wfits) then begin
          
          xrange = minmax(wave)
          yrange = [$
            (min(left[istd].spec)<min(center[istd].spec))<min(right[istd].spec),$
            (max(left[istd].spec)>max(center[istd].spec))>max(right[istd].spec)]*[0.9,1.1]
          
          djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, $
            yrange=yrange, charsize=1.5, charthick=postthick, title=title, $
            xtickname=replicate(' ',10), ytitle='Observed Counts', $
            xthick=postthick, ythick=postthick, position=pos[*,0]
          djs_oplot, center[istd].wave, center[istd].spec, ps=10, color='green', thick=2
          djs_oplot, left[istd].wave, left[istd].spec, ps=10, color='red', thick=2
          djs_oplot, right[istd].wave, right[istd].spec, ps=10, color='blue', thick=2

          legend, starname, /right, /top, charsize=1.5, charthick=postthick, box=0

          yrange = [$
            (min(dleft[istd].spec)<min(dcenter[istd].spec))<min(dright[istd].spec),$
            (max(dleft[istd].spec)>max(dcenter[istd].spec))>max(dright[istd].spec)]*[0.9,1.1]

          djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, $
            yrange=yrange, charsize=1.5, charthick=postthick, $
            xtickname=replicate(' ',10), ytitle='Normalized Counts', $
            xthick=postthick, ythick=postthick, position=pos[*,1], /noerase, yminor=3
          djs_oplot, dcenter[istd].wave, dcenter[istd].spec, ps=10, color='green', thick=2
          djs_oplot, dleft[istd].wave, dleft[istd].spec, ps=10, color='red', thick=2
          djs_oplot, dright[istd].wave, dright[istd].spec, ps=10, color='blue', thick=2

          shifts = ['Left Shift = '+strtrim(string(right_shift,format='(F12.3)'),2),$
            'Right Shift = '+strtrim(string(left_shift,format='(F12.3)'),2)]
          
          legend, shifts, /right, /top, charsize=1.5, charthick=postthick, box=0

          yrange = minmax(spec)*[0.9,1.1]
          
          djs_plot, wave, spec, xsty=3, ysty=3, xrange=xrange, ps=10, $
            yrange=yrange, charsize=1.5, charthick=postthick, thick=2, $
            xtitle='Wavelength [\AA]', ytitle='Final Counts', $
            xthick=postthick, ythick=postthick, position=pos[*,2], yminor=3, /noerase

          if keyword_set(doplot) and (nstd gt 1L) then cc = get_kbrd(1)
          
       endif
          
    endfor
       
    if keyword_set(wfits) then dfpsclose

    icleanup, left
    icleanup, center
    icleanup, right

    icleanup, dleft
    icleanup, dcenter
    icleanup, dright
    
return
end
