pro compare_kmatlas_kenn92, postscript=postscript
; jm04may06uofa
; compare overlapping spectra from Kenn92 and the K/M atlas

    pspath = atlas_path(/kenn92)+'analysis/'
    
    atlas = read_integrated(/silent)
    kenn92 = read_kenn92(/silent)

; match by coordinates

    raref = 15.0*im_hms2dec(atlas.ra)
    deref = im_hms2dec(atlas.dec)

    ra = 15.0*im_hms2dec(kenn92.ra)
    de = im_hms2dec(kenn92.dec)
    
    ntot = djs_angle_match(raref,deref,ra,de,dtheta=5.0/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    match = where(mindx ne -1,nmatch)
    niceprint, atlas[match].galaxy, kenn92[mindx[match]].galaxy, mdist[match]*3600.0

    splog, 'There are '+strn(nmatch)+' atlas galaxies in KENN92.'

    atlas = atlas[match]
    kenn92 = kenn92[mindx[match]]

    if keyword_set(postscript) then $
      dfpsplot, pspath+'kenn92_compare.ps', /color, xsize=8.5, ysize=11.0 else $
      window, 0, xs=600, ys=600*10.0/8.0

    pagemaker, nx=2, ny=nmatch/2, ypage=11.0, xpage=8.5, xspace=0.1, yspace=0.0, $
      xmargin=[0.8,0.4], ymargin=[0.1,1.2], position=pos, /normal
    xrange = [3650,7000]
    yrange = [-0.6,3.9]

    binsize = 150.0
    
    for i = 0L, nmatch-1L do begin

       s1 = rd1dspec(atlas[i].driftfile,/silent,datapath=atlas_path(/atlas1d),/normalize)
       s2 = rd1dspec(kenn92[i].specfile,/silent,datapath=atlas_path(/kenn92)+'data/',/normalize)

       flux = s1.spec
       wave = s1.wave

       linterp, s2.wave, s2.spec, wave, nflux
;      nflux = interpol(s2.spec,s2.wave,wave)

; bin the spectra and then take the ratio

       bin1spec = im_binspec(flux,wave,binsize=binsize,binwave=binwave,binspec_err=bin1spec_err)
       bin2spec = im_binspec(nflux,wave,binsize=binsize,binwave=binwave,binspec_err=bin2spec_err)
       nbins = n_elements(binwave)
       
       ratio = bin1spec/bin2spec-1
       ratioerr = im_compute_error(bin1spec,bin1spec_err,bin2spec,bin2spec_err,/quotient)
;      ratio = flux/nflux - 1.0
       junk = im_stats(100*ratio,/verbose,no_head=(i ne 0L))

;      yrange[1] = (max(flux)>max(nflux))*1.05

       if (i eq 12L) then xtitle = 'Wavelength ('+angstrom()+')' else delvarx, xtitle
       if (i lt 12L) then xtickname = replicate(' ',10) else delvarx, xtickname
       if (odd(i) eq 0L) then delvarx, ytickname else ytickname = replicate(' ',10)
       
       djs_plot, wave, flux, ps=10, position=pos[*,i], noerase=(i ne 0L), $
         xthick=5.0, ythick=5.0, charthick=5.0, charsize=1.2, ytickname=ytickname, $
         xtickname=xtickname, xrange=xrange, yrange=yrange, xsty=3, ysty=3, $
         yminor=2, xtitle=xtitle, xminor=3, xtickinterval=1000.0
       djs_oplot, wave, nflux+1, color='orange', ps=10
       djs_oplot, !x.crange, [0,0], color='grey', line=0, thick=5.0
       oploterror, binwave, ratio, replicate(binsize/2.0,nbins), ratioerr, $
         color=djs_icolor('blue'), errcolor=djs_icolor('blue'), $
         ps=3, thick=3.0, errthick=3.0
       legend, strtrim(atlas[i].nice_galaxy,2)+' - '+strtrim(kenn92[i].nice_galaxy,2), /left, $
         /top, box=0, charsize=1.3, charthick=5.0

       icleanup, s1 & icleanup, s2
       
    endfor

    xyouts, 0.05, 0.535, 'Normalized Flux Density', charsize=1.3, charthick=5.0, $
      orientation=90, align=0.5, /normal
    
    if keyword_set(postscript) then dfpsclose

return
end
    
    
