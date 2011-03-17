pro nfgs_compare_spectra, atlas, nfgs, result, ratiodata, postscript=postscript
; jm04mar21uofa
; jm04dec10uofa
; jm05aug01uofa - updated

    atlaspath = atlas_path(/atlas1d)
    nfgspath = nfgs_path(/spec1d)
    pspath = atlas_path(/papers)+'atlas/FIG_ATLAS/'
    
    if (n_elements(atlas) eq 0L) then atlas = atlas_read_info()
    if (n_elements(nfgs) eq 0L) then nfgs = nfgs_read_info()

; which NFGS galaxies are in the atlas    

    raref = 15.0*im_hms2dec(atlas.ra)
    deref = im_hms2dec(atlas.dec)

    ra = 15.0*im_hms2dec(nfgs.ra)
    de = im_hms2dec(nfgs.dec)
    
    ntot = djs_angle_match(raref,deref,ra,de,dtheta=15.0/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    match = where(mindx ne -1,nmatch)
    niceprint, atlas[match].galaxy, nfgs[mindx[match]].galaxy, mdist[match]*3600.0

    splog, 'There are '+string(nmatch,format='(I0)')+' atlas galaxies in the NFGS.'

    matchatlas = atlas[match]
    matchnfgs = nfgs[mindx[match]]

; initialize some plotting variables
    
    pcharsize = 1.4
    lcharsize = 1.2

    ncols = 3L
    nrows = ceil(nmatch/float(ncols))

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
       dfpsplot, pspath+'nfgs_compare_spectra.eps', xsize=8.5, ysize=11.0, /encapsulated;, /color
       postthick = 8.0 
       postthick2 = 5.0 
    endif else begin
       postthick = 2.0
       postthick2 = 2.0 
       im_window, 0, xratio=0.5, yratio=0.9
    endelse
    
    xrange = [3650,6900]
    yrange = [0.5,4.6]
    binsize = 200.0

    ytitle = 'Relative Flux [arbitrary units]'
    
; loop on each object    

    galaxy = strarr(nmatch)
    for i = 0L, nmatch-1L do begin

       atlasfile = strtrim(matchatlas[i].drift_file,2)
       s1 = rd1dspec(atlasfile,/silent,datapath=atlaspath)
       info1 = iforage(atlasfile,datapath=atlaspath)
       galaxy[i] = repstr(repstr(strtrim(info1.galaxy,2),'NGC','NGC '),'UGC','UGC ')

       nfgsfile = strtrim(matchnfgs[i].drift_file,2)
       s2 = rd1dspec(nfgsfile,/silent,datapath=nfgspath)
       info2 = iforage(nfgsfile,datapath=nfgspath)
       
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

       stats = im_stats(100*(ratio-1),no_head=(i ne 0L),verbose=0)
       keepstats = struct_trimtags(stats,select=['MIN','MAX','MEAN','MEDIAN','SIGMA','SIGMA_REJ'])
       if (i eq 0L) then result = keepstats else result = [ [result], [keepstats] ]

       if (i eq 0L) then $
         ratiodata = {binwave: float(binwave), ratio: fltarr(nbins,nmatch)}
       ratiodata.ratio[*,i] = 100*(ratio-1)
       
;      yrange[1] = (max(flux1)>max(flux2))*1.05

       if (i mod ncols) eq 0L then begin
          delvarx, ytickname
       endif else begin
          ytickname = replicate(' ',10)
          ytitle = ''
       endelse

       if (i ge nmatch-ncols) then begin
          delvarx, xtickname
          xtitle = 'Wavelength ['+angstrom()+']'
       endif else begin
          xtickname = replicate(' ',10)
          xtitle = ''
       endelse
       
       djs_plot, wave, flux1+1, ps=10, position=pos[*,i], noerase=(i ne 0L), $
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

       legstr = [galaxy[i]+' ('+string(sxpar(s2.header,'NFGS_ID'),$
;      legstr = [strtrim(matchatlas[i].nice_galaxy,2)+' ('+string(sxpar(s2.header,'NFGS_ID'),$
         format='(I3.3)')+')']
       legend, legstr, /left, /top, box=0, charsize=lcharsize, $
         charthick=postthick, clear=keyword_set(postscript)

       icleanup, s1 & icleanup, s2

       if (i eq 0L) then xyouts, 0.06, 0.535, ytitle, charsize=1.5, $
         charthick=postthick, orientation=90, align=0.5, /normal
       
    endfor

    result = struct_addtags(replicate({galaxy: ''},nmatch),result)
    result.galaxy = galaxy
    struct_print, result

    if keyword_set(postscript) then dfpsclose

; diagnostic plot

    indx = lindgen(nmatch)
    indx = [0,1,2,4,5,6,7,8,10]
    nindx = n_elements(indx)

    im_window, 2, xr=0.4, /square
    plot, ratiodata.binwave, ratiodata.ratio[*,indx[0]], yr=[-20.0,+20.0], ps=-4, xsty=3, ysty=3
    for i = 1L, nindx-1L do oplot, ratiodata.binwave, ratiodata.ratio[*,indx[i]], ps=-4

return
end
