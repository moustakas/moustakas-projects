pro template_tests, postscript=postscript
; jm03dec5uofa
; after running RUN_TEMPLATE_TESTS generate postscript plots for
; analysis and the data paper

; read the data

    path = atlas_path(/specfit)+'template_tests/'
    file= 'integrated_atlas_Z02'+['_04','_14','_20']

    f07 = read_integrated()
    f04 = read_integrated(linepath=path,root=file[0])
    f14 = read_integrated(linepath=path,root=file[1])
    f20 = read_integrated(linepath=path,root=file[2])

; compute some plot variables and then fill the plot structure 

    psname = 'template_tests.ps'
    
    lines = ['OII_3727','OIII_5007','H_ALPHA']
    lines = [lines,lines+'_EW']
    ntags = n_elements(lines)
    reftags = tag_names(f04)
    stags = ['x04_','x14_','x20_']
    nstruct = n_elements(stags)
    tt = {name: ''}

; emission-line fluxes and EWs    
    
    for j = 0L, nstruct-1L do begin

       case j of
          0L: f = f04
          1L: f = f14
          2L: f = f20
       endcase
       
       tags = stags[j]+lines
       for i = 0L, ntags-1L do begin

          w1 = where(lines[i] eq reftags)
          g = where(((f07.(w1))[1,*] gt 0.0) and ((f.(w1))[1,*] gt 0.0),ng)
          array = reform(((f07[g].(w1))[0,*]-(f[g].(w1))[0,*])/(f07[g].(w1))[1,*])

          ss = im_stats(array,sigrej=5.0)
          array_stats = '('+strtrim(string(ss.mean,format='(F7.3)'),2)+' \pm '+$
            strtrim(string(ss.sigma_rej,format='(F7.3)')+')',2)

          tt = create_struct(tt,tags[i],array,tags[i]+'_stat',array_stats)
          
       endfor

    endfor

; add some continuum chi2 tags

    ttchi2 = {$    
      y04_chi2:    (f04.continuum_chi2-f07.continuum_chi2), $
      y14_chi2:    (f14.continuum_chi2-f07.continuum_chi2), $
      y20_chi2:    (f20.continuum_chi2-f07.continuum_chi2)}
    tt = struct_addtags(tt,ttchi2)

; chi2 plots

    im_openclose, psname, postscript=postscript
    
    pagemaker, nx=1, ny=3, position=pos, /normal, $
      xmargin=[1.2,0.2], ymargin=[0.2,1.2], yspace=0.0
    yrange = [-1.8,1.8]
    xrange = [0.5,10.0]
    
    plotsym, 0, 1
    djs_plot, f07.continuum_chi2, tt.y04_chi2, ps=8, xrange=xrange, yrange=yrange, $
      position=pos[*,0], xthick=2.0, ythick=2.0, xtickname=replicate(' ',10), $
      charsize=1.5, charthick=2.0, ytitle='\chi^{2}(04) - \chi^{2}(07)', xsty=3, ysty=3
    oplot, !x.crange, [0,0], line=0, thick=2.0

    plotsym, 8, 1
    djs_plot, f07.continuum_chi2, tt.y14_chi2, ps=8, xrange=xrange, yrange=yrange, $
      position=pos[*,1], xthick=2.0, ythick=2.0, xtickname=replicate(' ',10), $
      charsize=1.5, charthick=2.0, /noerase, ytitle='\chi^{2}(14)-\chi^{2}(07)', xsty=3, ysty=3
    oplot, !x.crange, [0,0], line=0, thick=2.0
    
    plotsym, 3, 1, /fill
    djs_plot, f07.continuum_chi2, tt.y20_chi2, ps=8, xrange=xrange, yrange=yrange, $
      position=pos[*,2], xthick=2.0, ythick=2.0, xtitle='\chi^{2}(07)', ytitle='\chi^{2}(20)-\chi^{2}(07)', $
      charsize=1.5, charthick=2.0, /noerase, xsty=3, ysty=3
    oplot, !x.crange, [0,0], line=0, thick=2.0

    if not keyword_set(postscript) then cc = get_kbrd(1)
    
; emission-line flux plots
    
    xrange = [-2.0,2.0]
    yrange = [0,0.8]
    legstr = ['[O II] ','[O III] ','H\alpha ']      
    binsize = 0.2
    line = [0,1,2]
    
    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, xthick=2.0, $
      ythick=2.0, charsize=2.0, charthick=2.0, ysty=1, xsty=3, $
      xtitle='[F(07)-F(04)] / \delta'+'F(07)', ytitle='Fraction'
    im_plothist, tt.x04_oii_3727, bin=binsize, thick=4.0, color=djs_icolor('blue'), /overplot, /fraction, line=line[0]
    im_plothist, tt.x04_oiii_5007, bin=binsize, thick=4.0, color=djs_icolor('green'), /overplot, $
      /fraction, line=line[1], /fill, /fline, fspacing=0.1, fcolor=djs_icolor('green'), forientation=45
    im_plothist, tt.x04_h_alpha, bin=binsize, thick=4.0, color=djs_icolor('red'), /overplot, $
      /fraction, line=line[2], /fill, /fline, fspacing=0.2, fcolor=djs_icolor('red'), forientation=135
    legend, textoidl(legstr+[tt.x04_oii_3727_stat,tt.x04_oiii_5007_stat,tt.x04_h_alpha_stat]), $
      line=line, thick=4.0, spacing=2.0, color=djs_icolor(['blue','green','red']), /right, /top, box=0, $
      charsize=1.5, charthick=2.0
    
    if not keyword_set(postscript) then cc = get_kbrd(1)

    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, xthick=2.0, $
      ythick=2.0, charsize=2.0, charthick=2.0, ysty=1, xsty=3, $
      xtitle='[F(07)-F(14)] / \delta'+'F(07)', ytitle='Fraction'
    im_plothist, tt.x14_oii_3727, bin=binsize, thick=4.0, color=djs_icolor('blue'), /overplot, /fraction, line=line[0]
    im_plothist, tt.x14_oiii_5007, bin=binsize, thick=4.0, color=djs_icolor('green'), /overplot, $
      /fraction, line=line[1], /fill, /fline, fspacing=0.1, fcolor=djs_icolor('green'), forientation=45
    im_plothist, tt.x14_h_alpha, bin=binsize, thick=4.0, color=djs_icolor('red'), /overplot, $
      /fraction, line=line[2], /fill, /fline, fspacing=0.2, fcolor=djs_icolor('red'), forientation=135
    legend, textoidl(legstr+[tt.x14_oii_3727_stat,tt.x14_oiii_5007_stat,tt.x14_h_alpha_stat]), $
      line=line, thick=4.0, spacing=2.0, color=djs_icolor(['blue','green','red']), /right, /top, box=0, $
      charsize=1.5, charthick=2.0
    
    if not keyword_set(postscript) then cc = get_kbrd(1)

    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, xthick=2.0, $
      ythick=2.0, charsize=2.0, charthick=2.0, ysty=1, xsty=3, $
      xtitle='[F(07)-F(20)] / \delta'+'F(07)', ytitle='Fraction'
    im_plothist, tt.x20_oii_3727, bin=binsize, thick=4.0, color=djs_icolor('blue'), /overplot, /fraction, line=line[0]
    im_plothist, tt.x20_oiii_5007, bin=binsize, thick=4.0, color=djs_icolor('green'), /overplot, $
      /fraction, line=line[1], /fill, /fline, fspacing=0.1, fcolor=djs_icolor('green'), forientation=45
    im_plothist, tt.x20_h_alpha, bin=binsize, thick=4.0, color=djs_icolor('red'), /overplot, $
      /fraction, line=line[2], /fill, /fline, fspacing=0.2, fcolor=djs_icolor('red'), forientation=135
    legend, textoidl(legstr+[tt.x20_oii_3727_stat,tt.x20_oiii_5007_stat,tt.x20_h_alpha_stat]), $
      line=line, thick=4.0, spacing=2.0, color=djs_icolor(['blue','green','red']), /right, /top, box=0, $
      charsize=1.5, charthick=2.0

    if not keyword_set(postscript) then cc = get_kbrd(1)

; EW plots
    
    xrange = [-1.0,6.0]
    yrange = [0,0.8]
    legstr = ['EW([O II]) ','EW([O III]) ','EW(H\alpha) ']    
    binsize = 0.2
    
    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, xthick=2.0, $
      ythick=2.0, charsize=2.0, charthick=2.0, ysty=1, xsty=3, $
      xtitle='[EW(07)-EW(04)] / \delta'+'EW(07)', ytitle='Fraction'
    im_plothist, tt.x04_oii_3727_ew, bin=binsize, thick=4.0, color=djs_icolor('blue'), /overplot, /fraction, line=line[0]
    im_plothist, tt.x04_oiii_5007_ew, bin=binsize, thick=4.0, color=djs_icolor('green'), /overplot, $
      /fraction, line=line[1], /fill, /fline, fspacing=0.1, fcolor=djs_icolor('green'), forientation=45
    im_plothist, tt.x04_h_alpha_ew, bin=binsize, thick=4.0, color=djs_icolor('red'), /overplot, $
      /fraction, line=line[2], /fill, /fline, fspacing=0.2, fcolor=djs_icolor('red'), forientation=135
    legend, textoidl(legstr+[tt.x04_oii_3727_ew_stat,tt.x04_oiii_5007_ew_stat,tt.x04_h_alpha_ew_stat]), $
      line=line, thick=4.0, spacing=2.0, color=djs_icolor(['blue','green','red']), /right, /top, box=0, $
      charsize=1.5, charthick=2.0
    
    if not keyword_set(postscript) then cc = get_kbrd(1)

    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, xthick=2.0, $
      ythick=2.0, charsize=2.0, charthick=2.0, ysty=1, xsty=3, $
      xtitle='[EW(07)-EW(14)] / \delta'+'EW(07)', ytitle='Fraction'
    im_plothist, tt.x14_oii_3727_ew, bin=binsize, thick=4.0, color=djs_icolor('blue'), /overplot, /fraction, line=line[0]
    im_plothist, tt.x14_oiii_5007_ew, bin=binsize, thick=4.0, color=djs_icolor('green'), /overplot, $
      /fraction, line=line[1], /fill, /fline, fspacing=0.1, fcolor=djs_icolor('green'), forientation=45
    im_plothist, tt.x14_h_alpha_ew, bin=binsize, thick=4.0, color=djs_icolor('red'), /overplot, $
      /fraction, line=line[2], /fill, /fline, fspacing=0.2, fcolor=djs_icolor('red'), forientation=135
    legend, textoidl(legstr+[tt.x14_oii_3727_ew_stat,tt.x14_oiii_5007_ew_stat,tt.x14_h_alpha_ew_stat]), $
      line=line, thick=4.0, spacing=2.0, color=djs_icolor(['blue','green','red']), /right, /top, box=0, $
      charsize=1.5, charthick=2.0
    
    if not keyword_set(postscript) then cc = get_kbrd(1)

    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, xthick=2.0, $
      ythick=2.0, charsize=2.0, charthick=2.0, ysty=1, xsty=3, $
      xtitle='[EW(07)-EW(20)] / \delta'+'EW(07)', ytitle='Fraction'
    im_plothist, tt.x20_oii_3727_ew, bin=binsize, thick=4.0, color=djs_icolor('blue'), /overplot, /fraction, line=line[0]
    im_plothist, tt.x20_oiii_5007_ew, bin=binsize, thick=4.0, color=djs_icolor('green'), /overplot, $
      /fraction, line=line[1], /fill, /fline, fspacing=0.1, fcolor=djs_icolor('green'), forientation=45
    im_plothist, tt.x20_h_alpha_ew, bin=binsize, thick=4.0, color=djs_icolor('red'), /overplot, $
      /fraction, line=line[2], /fill, /fline, fspacing=0.2, fcolor=djs_icolor('red'), forientation=135
    legend, textoidl(legstr+[tt.x20_oii_3727_ew_stat,tt.x20_oiii_5007_ew_stat,tt.x20_h_alpha_ew_stat]), $
      line=line, thick=4.0, spacing=2.0, color=djs_icolor(['blue','green','red']), /right, /top, box=0, $
      charsize=1.5, charthick=2.0

    im_openclose, psname, postscript=postscript, /close

stop    
    
return
end    
