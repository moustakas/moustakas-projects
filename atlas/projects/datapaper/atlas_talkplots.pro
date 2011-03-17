pro atlas_talkplots, atlas, postscript=postscript
; jm05nov22uofa - 
    
    if (n_elements(atlas) eq 0L) then atlas = read_integrated()

    k92 = read_92kennicutt()
    nfgs = read_00jansen()
    
    talkpath = '/home/ioannis/jobs/talk/'

    if keyword_set(postscript) then begin
       postthick = 8.0
       talkcolor = 'white'
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
       talkcolor = ''
    endelse

; ---------------------------------------------------------------------------    
; B magnitude distribution    
; ---------------------------------------------------------------------------    

    psname = 'atlas_blum.ps'
    pngname = 'atlas_blum.png'
    
    if keyword_set(postscript) then dfpsplot, talkpath+psname, encapsulated=encapsulated, xsize=8.5, ysize=9.0, /color

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.0, height=7.0, $
      xmargin=[1.1,0.4], ymargin=[0.9,1.1], xpage=8.5, ypage=9.0, $
      position=pos, /normal

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

    indx = where(atlas.rc3_b_lum gt -900.0)
    x = atlas[indx].rc3_b_lum
    xabs = atlas[indx].rc3_m_b
    stats = im_stats(x,/verbose)
    
    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    xrange = [6.9,11.6]
    binsize = 0.25
;   xtitle = 'M_{B} [mag]'
;   xrange = [-12,-23]
;   binsize = 0.75

    im_plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    yrange = minmax(ybin)*[1.0,1.2]

    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=2.0, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=11, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,0], /noerase, $
      color=djs_icolor(talkcolor)
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=2.0, charthick=postthick, color=djs_icolor(talkcolor)
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=2.0, charthick=postthick, color=djs_icolor(talkcolor)

; overplot K92 and NFGS

    g = where(k92.rc3_b_lum gt -900)
    im_plothist, k92[g].rc3_b_lum, bin=binsize, /overplot, line=0, thick=postthick, $
      /halfbin, /fill, fcolor=djs_icolor('blue'), color=djs_icolor('blue')
;   im_plothist, k92[g].rc3_b_lum, bin=binsize, /overplot, line=0, thick=postthick, $
;     /halfbin, /fill, /fline, forientation=135, fcolor=djs_icolor('blue'), $
;     color=djs_icolor('blue'), fspacing=0.05

    g = where(nfgs.rc3_b_lum gt -900)
    im_plothist, nfgs[g].rc3_b_lum, bin=binsize, /overplot, line=0, thick=postthick, $
      /halfbin, /fill, /fline, forientation=45, fcolor=djs_icolor('red'), $
      color=djs_icolor('red'), fspacing=0.15

; now overplot my data    
    
    im_plothist, x, bin=binsize, /overplot, thick=postthick, /halfbin, color=djs_icolor(talkcolor)

    legend, ['Moustakas & Kennicutt (2006)','Jansen et al. (2000)','Kennicutt (1992)'], $
      /left, /top, box=0, charsize=1.7, charthick=postthick, line=[0,0,0], $
      color=djs_icolor([talkcolor,'red','blue']), thick=postthick+2, spacing=2.0, $
      textcolor=djs_icolor(talkcolor)
    
    if keyword_set(postscript) then begin
       dfpsclose
       spawn, ['convert '+talkpath+psname+' '+talkpath+repstr(psname,'.ps','.png')], /sh
    endif else cc = get_kbrd(1)

; ---------------------------------------------------------------------------    
; FIR-Optical distribution    
; ---------------------------------------------------------------------------    

    psname = 'atlas_firopt.ps'
    pngname = 'atlas_firopt.png'
    
    if keyword_set(postscript) then dfpsplot, talkpath+psname, encapsulated=encapsulated, xsize=8.5, ysize=9.0, /color

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.0, height=7.0, $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal

    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

    indx = where(atlas.l_fir_l_b gt -900)
    x = alog10(atlas[indx].l_fir_l_b)

    stats = im_stats(x,/verbose)

    xtitle = 'log L(FIR)/L(B)'
    xrange = [-2.2,2.8]
    binsize = 0.25
    
    im_plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    yrange = minmax(ybin)*[1.0,1.25]

    plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=2.0, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,0], /noerase, $
      color=djs_icolor(talkcolor)

; overplot K92 and NFGS

    g = where((k92.l_fir_l_b gt -900))
    im_plothist, alog10(k92[g].l_fir_l_b), bin=binsize, /overplot, line=0, thick=postthick, $
      /halfbin, /fill, fcolor=djs_icolor('blue'), color=djs_icolor('blue')

    g = where(nfgs.l_fir_l_b gt -900)
    im_plothist, alog10(nfgs[g].l_fir_l_b), bin=binsize, /overplot, line=0, thick=postthick, $
      /halfbin, /fill, /fline, forientation=45, fcolor=djs_icolor('red'), $
      color=djs_icolor('red'), fspacing=0.15

; now overplot my data    
    
    im_plothist, x, bin=binsize, /overplot, thick=postthick, /halfbin, color=djs_icolor(talkcolor)

    legend, ['Moustakas & Kennicutt (2006)','Jansen et al. (2000)','Kennicutt (1992)'], $
      /left, /top, box=0, charsize=1.7, charthick=postthick, line=[0,0,0], $
      color=djs_icolor([talkcolor,'red','blue']), thick=postthick+2, spacing=2.0, $
      textcolor=djs_icolor(talkcolor)

    if keyword_set(postscript) then begin
       dfpsclose
       spawn, ['convert '+talkpath+psname+' '+talkpath+repstr(psname,'.ps','.png')], /sh
    endif else cc = get_kbrd(1)

stop

; ---------------------------------------------------------------------------
; model spectra
; ---------------------------------------------------------------------------

    rootpath = atlas_path(/projects)+'simulations/nebular/'
    basepath = rootpath+'basemodels/'

    redshift = 1E-15
    mgalaxy = 1E10              ; [M_sun]
    starvdisp = 120.0           ; [km/s]
    dustratio = 0.5
    
    minwave = 1000.0            ; [rest, Angstrom]
    maxwave = 9500.0            ; [rest, Angstrom]
    fwhmres = 10.0
    ewhb = 5.0

    taumodel = 'bc03_m62_tau3.ised'

    model = im_read_bc03(isedfile=taumodel,isedpath=basepath,/silent,$
      /salpeter,minwave=minwave,maxwave=maxwave,bc03_extras=modelinfo,$
      age=10.0)

    bc03_grid2fits, model, fitsgrid, mgalaxy=mgalaxy, redshift=redshift, $
      fwhmres=fwhmres, starvdisp=starvdisp, minwave=minwave, $
      maxwave=maxwave, outfile=outfile, outpath=outpath, /espectrum, $
      hbeta_ew=ewhb, ebv_gas=0.0, dustratio=dustratio, wfits=wfits

    wave = fitsgrid.wave
    npix = n_elements(wave)

    xrange = [3500,7100]
    get_element, wave, xrange, xx

    norm = max(fitsgrid.flux[xx[0]:xx[1]])
    flux = fitsgrid.flux/norm
    eflux = fitsgrid.espectrum/norm
    cflux = fitsgrid.continuum/norm
    
    ytitle = 'Relative Amount of Light'
    xtitle = 'Wavelength'

    yrange = [0.1,1.05]

    xbalmer = [3734.0,3750.0,3771.0,3797.0,3835.0,3889.0,3970.0,4101.0,4340.0,4861.0,6563.0]
    x1balmer = xbalmer-25.0
    x2balmer = xbalmer+25.0
    xbalmer = [3850.0,4150.0,xbalmer]
    x1balmer = xbalmer-[150.0,150.0,x1balmer*0.0+25.0]
    x2balmer = xbalmer+[150.0,150.0,x2balmer*0.0+25.0]

    indexpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='etc')
    indexfile = 'indexlist.dat'
    readcol, indexpath+indexfile, licknames, w1, w2, wb1, wb2, wr1, wr2, units, $
      format='A,D,D,D,D,D,D,A', /silent, comment='#'
    x1balmer = [x1balmer,wb1]
    x2balmer = [x2balmer,wr2]
        
    xemission = [3727.0,4101.0,4340.0,4861.0,4959.0,5007.0,6548.0,6563.0,6584.0]
    x1emission = xemission-25.0
    x2emission = xemission+25.0

    bin1 = lindgen(npix/2-1)*2
    bin2 = bin1+3L
    nsubpix = n_elements(bin1)
    colors = 255*findgen(nsubpix)/float(nsubpix)

    psname = 'optical_spectrum.ps'
    if keyword_set(postscript) then dfpsplot, talkpath+psname, encapsulated=encapsulated, xsize=8.5, ysize=9.0, /color
    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

    plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=2.0, /noerase, $
      ythick=2.0, xtitle='Wavelength ['+angstrom()+']', ytitle='Relative Flux', charsize=1.7, charthick=2.0, $
      yrange=yrange, xrange=xrange, thick=2.0, color=djs_icolor(talkcolor), $
      ytickname=replicate(' ',10)
    djs_oplot, wave, cflux+eflux, ps=10, thick=0.5, color='yellow'

    if keyword_set(postscript) then begin
       dfpsclose
       spawn, ['convert '+talkpath+psname+' '+talkpath+repstr(psname,'.ps','.png')], /sh
    endif else cc = get_kbrd(1)

; all emission lines
    
    psname = 'optical_spectrum_lines.ps'
    if keyword_set(postscript) then dfpsplot, talkpath+psname, encapsulated=encapsulated, xsize=8.5, ysize=9.0, /color
    polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, color=djs_icolor('black')

    plot, [0], [0], /nodata, xsty=3, ysty=3, xthick=2.0, /noerase, $
      ythick=2.0, xtitle='Wavelength ['+angstrom()+']', ytitle='Relative Flux', charsize=1.7, charthick=2.0, $
      yrange=yrange, xrange=xrange, thick=2.0, color=djs_icolor(talkcolor), $
      ytickname=replicate(' ',10)
    djs_oplot, wave, cflux+eflux, ps=10, thick=0.5, color='cyan'
    djs_oplot, wave, cflux, ps=10, thick=0.5, color='yellow'
    for i = 0L, n_elements(x1emission)-1L do begin
       get_element, wave, [x1emission[i],x2emission[i]], xx
       djs_oplot, wave[xx[0]:xx[1]], (cflux+eflux)[xx[0]:xx[1]], ps=10, color='cyan'
    endfor
    
    if keyword_set(postscript) then begin
       dfpsclose
       spawn, ['convert '+talkpath+psname+' '+talkpath+repstr(psname,'.ps','.png')], /sh
    endif else cc = get_kbrd(1)

stop    
    
    
return
end
