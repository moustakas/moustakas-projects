;+
; NAME:
;       ATLAS_COMPARE_REPEATERS
;
; PURPOSE:
;       Compare the repeat observations to make a figure for the data
;       paper.  
;
; CALLING SEQUENCE:
;       atlas_compare_repeaters, /postscript, /paper
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
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Aug 01, U of A
;-

pro atlas_compare_repeaters, result, ratiodata, postscript=postscript, paper=paper
    
    if keyword_set(paper) then postscript = 1L

    analysis_path = atlas_path(/analysis)
    repeatpath = atlas_path(/atlas1d)+'repeaters/'
    pspath = atlas_path(/papers)+'atlas/FIG_ATLAS/'

; read the repeaters text file    
    
    repeatersfile = 'atlas1d_combine_repeaters.txt'
    if (file_test(analysis_path+repeatersfile,/regular) eq 0L) then begin
       splog, 'File '+analysis_path+repeatersfile+' not found.'
       return
    endif

    repeaters = djs_readlines(analysis_path+repeatersfile)
    keep = where((strmatch(repeaters,'*#*') eq 0B) and (strcompress(repeaters,/remove) ne ''))
    repeaters = repeaters[keep]
    ngalaxy = n_elements(repeaters)

; parse the text file to determine which objects have two or more
; "good" repeat observations

    keep = lonarr(ngalaxy)+1L
    speclist = strarr(2,ngalaxy)

    for igalaxy = 0L, ngalaxy-1L do begin
       
       line = strsplit(repeaters[igalaxy],' ',/extract)
       nrepeat = long(line[0])  ; number of repeat observations

       if (nrepeat eq 1L) or strmatch(repeaters[igalaxy],'*nuclear*',/fold) then $
         keep[igalaxy] = 0L else $
         speclist[*,igalaxy] = strtrim(line[1L:2L],2)

    endfor

    speclist = speclist[*,where(keep,ngalaxy)]

; initialize some plotting variables

    pcharsize = 1.2
    lcharsize = 1.0
    
    ncols = 4L
    nrows = ceil(ngalaxy/float(ncols))

    xmargin = [1.4,0.5]
    ymargin = [0.5,1.4]

    width = replicate(3.5,ncols)
    height = replicate(3.5,nrows)

    xpage = total(width)+total(xmargin)
    ypage = total(height)+total(ymargin)

    pagemaker, nx=ncols, ny=nrows, xspace=0.0, yspace=0.0, $
      xmargin=xmargin, ymargin=ymargin, width=width, /normal, $
      height=height, position=pos, xpage=xpage, ypage=ypage

    if keyword_set(postscript) then begin
       dfpsplot, pspath+'atlas_compare_repeaters.eps', /color, xsize=8.5, ysize=11.0, /encapsulated
       postthick = 8.0 
       postthick2 = 5.0 
    endif else begin
       postthick = 2.0
       postthick2 = 2.0 
       im_window, 0, xratio=0.8, yratio=0.8
    endelse
    
    xrange = [3650,6900]
    yrange = [0.5,4.6]
    binsize = 200.0

    ytitle = 'Relative Flux (arbitrary units)'

; loop on each object    

    galaxy = strarr(ngalaxy)
    
    for igalaxy = 0L, ngalaxy-1L do begin

       s1 = rd1dspec(speclist[0,igalaxy],/silent,datapath=repeatpath)
       info1 = iforage(speclist[0,igalaxy],datapath=repeatpath)
       galaxy[igalaxy] = info1.galaxy

       s2 = rd1dspec(speclist[1,igalaxy],/silent,datapath=repeatpath)
       info2 = iforage(speclist[1,igalaxy],datapath=repeatpath)

       wave = s1.wave
       flux1 = im_normalize(s1.spec,wave,normwave=5500.0,binsize=50.0)

       flux2 = interpol(s2.spec,s2.wave,wave) ; simple interpolation for the plot 
       flux2 = im_normalize(flux2,wave,normwave=5500.0,binsize=50.0)

; bin the spectra and then take the ratio

       bin1spec = im_binspec(flux1,wave,binsize=binsize,binwave=binwave,binspec_err=bin1spec_err)
       bin2spec = im_binspec(flux2,wave,binsize=binsize,binwave=binwave,binspec_err=bin2spec_err)
       nbins = n_elements(binwave)

       ratio = bin1spec/bin2spec
       ratioerr = im_compute_error(bin1spec,bin1spec_err,bin2spec,bin2spec_err,/quotient)
;      ratio = flux1/flux2

       stats = im_stats(100*(ratio-1),no_head=(igalaxy ne 0L),verbose=0)
       keepstats = struct_trimtags(stats,select=['MIN','MAX','MEAN','MEDIAN','SIGMA','SIGMA_REJ'])
       if (igalaxy eq 0L) then result = keepstats else result = [ [result], [keepstats] ]

       if (igalaxy eq 0L) then $
         ratiodata = {binwave: float(binwave), ratio: fltarr(nbins,ngalaxy)}
       ratiodata.ratio[*,igalaxy] = 100*(ratio-1)

;      yrange[1] = (max(flux1)>max(flux2))*1.05

       if (igalaxy mod ncols) eq 0L then begin
          delvarx, ytickname
       endif else begin
          ytickname = replicate(' ',10)
          ytitle = ''
       endelse

       if (igalaxy ge ngalaxy-ncols) then begin
          delvarx, xtickname
          xtitle = 'Wavelength ['+angstrom()+']'
       endif else begin
          xtickname = replicate(' ',10)
          xtitle = ''
       endelse
       
       djs_plot, wave, flux1+1, ps=10, position=pos[*,igalaxy], noerase=(igalaxy ne 0L), $
         xthick=postthick, ythick=postthick, charthick=postthick, charsize=pcharsize, ytickname=ytickname, $
         xtickname=xtickname, xrange=xrange, yrange=yrange, xsty=3, ysty=3, $
         yminor=2, xtitle=xtitle, ytitle='', xminor=3, xtickinterval=1000.0
       djs_oplot, wave, flux2+2, color='orange', ps=10
       djs_oplot, !x.crange, [1,1], color='grey', line=0, thick=0.5
;      djs_oplot, !x.crange, [1.1,1.1], color='grey', line=2, thick=postthick
;      djs_oplot, !x.crange, [0.9,0.9], color='grey', line=2, thick=postthick
       oploterror, binwave, ratio, replicate(binsize/2.0,nbins), ratioerr, $
         color=djs_icolor('blue'), errcolor=djs_icolor('blue'), $
         ps=3, thick=3.0, errthick=3.0

       legstr = strtrim(info1.galaxy,2)
       legend, legstr, /left, /top, box=0, charsize=lcharsize, $
         charthick=postthick, clear=keyword_set(postscript)

       icleanup, s1 & icleanup, s2

       if (igalaxy eq 0L) then xyouts, 0.05, 0.535, ytitle, charsize=1.5, $
         charthick=postthick, orientation=90, align=0.5, /normal

    endfor

    result = struct_addtags(replicate({galaxy: ''},ngalaxy),result)
    result.galaxy = galaxy
    struct_print, result
    
    if keyword_set(postscript) then dfpsclose

; diagnostic plot

    im_window, 2, xr=0.4, /square
    plot, ratiodata.binwave, ratiodata.ratio[*,0], yr=[-20.0,+20.0], ps=-4, xsty=3, ysty=3
    for i = 1L, ngalaxy-1L do oplot, ratiodata.binwave, ratiodata.ratio[*,i], ps=-4

return
end
    
