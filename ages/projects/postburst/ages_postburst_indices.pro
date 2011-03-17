pro ages_postburst_indices
; jm08nov24nyu - diagnostic plot showing the H-gamma and H-delta
;   Lick indices for the automatically selected sample, sorted by
;   redshift
    
    path = ages_path(/projects)+'postburst/brown/'

    indexfile = path+'indexlist_balmer.dat'
    
; read the sample

    ispec1 = read_ages(/ispec)
    readcol, path+'psb_auto_2008dec01.list', agesid, ra, dec, finalpass, $
      finalaper, format='I,D,D,I,I', comment='#', /silent
    final = speclinefit_locate(ispec1,'ages_'+string(finalpass,format='(I3.3)')+$
      '/'+string(finalaper,format='(I3.3)'))
;   final = final[29:30]
    ispec = ispec1[final]
;   ispec = ispec[sort(ispec.z)]

    ss = read_ages_specfit(ispec.galaxy)

;   psfile = path+'test.ps'
    psfile = path+'postburst_indices_08dec01.ps'
    pdffile = repstr(psfile,'.ps','.pdf')
    dfpsplot, psfile, /landscape, /color
    im_plotfaves, /postscript
;   pos = im_plotconfig(1,psfile=path+'postburst_indices.ps')
    
    for jj = 0L, n_elements(ispec)-1L do begin

       good = where((ss[*,0,jj] gt 0.0) and (ss[*,5,jj] gt 0.0),npix)
       wave = reform(ss[good,0,jj])
       ferr = reform(1.0/sqrt(ss[good,5,jj]))
       flux = reform(ss[good,1,jj]-ss[good,4,jj]) ; with emission lines
       modelflux1 = reform(ss[good,2,jj]+ss[good,3,jj])
       modelflux2 = reform(ss[good,2,jj])
       flux_nolines = reform(ss[good,1,jj]-ss[good,3,jj]) ; no emission lines

; make the plot       

       xrange = [3650,5100] ; [3650.0,6750.0] ; [4000,4500]
       get_element, wave, xrange, xx
       maxstats = im_stats(flux[xx[0]:xx[1]],sigrej=3.0)
       yrange = [-0.12*maxstats.maxrej,1.3*maxstats.maxrej]
       
       djs_plot, [0], [0], /nodata, xrange=xrange, yrange=1E17*yrange, /xsty, /ysty, $
         xtitle='Rest Wavelength (\AA)', ytitle='Flux (10^{-17} '+flam_units()+')', $
         charsize=1.8
       djs_oplot, wave, 1E17*flux, ps=10;, color='grey'
;      djs_oplot, wave, 1E17*flux_nolines, ps=10;, color='grey'
       djs_oplot, wave, 1E17*modelflux1, ps=10, color='red'
       djs_oplot, wave, 1E17*modelflux2, ps=10, color='dark green'
       djs_oplot, 4101.0*[1,1], !y.crange, line=0, thick=10.0, color='blue'
       djs_oplot, 4340.0*[1,1], !y.crange, line=0, thick=10.0, color='orange'
       djs_oplot, 5577.0*[1,1]/(1.0+ispec[jj].z), !y.crange, line=0, thick=10.0, color='dark green'
       legend, [repstr(ispec[jj].galaxy,'_',' '),$
         'z = '+strtrim(string(ispec[jj].z,format='(F12.3)'),2)], $
         /left, /top, box=0, charsize=1.8, /clear
       legend, textoidl(['H\delta \lambda4101','H\gamma \lambda4340','OI \lambda5577']), $
         /right, /top, box=0, /clear, charsize=1.8, $
         textcolor=djs_icolor(['blue','orange','dark green'])

; make an inset plot for H-alpha
       
       xrange = 6563+[-75.0,+75.0]
       get_element, wave, xrange, xx
       maxstats = im_stats(flux[xx[0]:xx[1]],sigrej=5.0)
       yrange = [-0.05*maxstats.maxrej,1.1*maxstats.max]
       djs_plot, [0], [0], /nodata, xrange=xrange, yrange=1E17*yrange, /xsty, /ysty, $
         xtitle='\lambda (\AA)', ytitle='F_{\lambda}', xtickinterval=50.0, $
         charsize=1.2, position=[0.6,0.25,0.85,0.45], /noerase
       djs_oplot, wave, 1E17*flux, ps=10;, color='red'
;      djs_oplot, wave, 1E17*flux_nolines, ps=10;, color='grey'
       djs_oplot, wave, 1E17*modelflux1, ps=10, color='red'
       djs_oplot, wave, 1E17*modelflux2, ps=10, color='dark green'
       
;      doit = spectral_indices(wave,thisflux,indexfile=indexfile,$
;        /nobreaks,position=pos,/postscript)
; index plot with the lines subtracted
       doit = spectral_indices(wave,flux_nolines,indexfile=path+'indexlist_hdelta.dat',$
         /nobreaks,/silent,/postscript);,position=pos[*,0])
       djs_oplot, 5577.0*[1,1]/(1.0+ispec[jj].z), !y.crange, line=0, thick=10.0, color='dark green'
       doit = spectral_indices(wave,flux_nolines,indexfile=path+'indexlist_hgamma.dat',$
         /nobreaks,/silent,/postscript);,position=pos[*,0])
       djs_oplot, 5577.0*[1,1]/(1.0+ispec[jj].z), !y.crange, line=0, thick=10.0, color='dark green'
; index plot with the lines subtracted
;      doit = spectral_indices(wave,flux,indexfile=path+'indexlist_hdelta.dat',$
;        /nobreaks,/silent,/postscript);,position=pos[*,0])
;      djs_oplot, 5577.0*[1,1]/(1.0+ispec[jj].z), !y.crange, line=0, thick=10.0, color='red'
;      doit = spectral_indices(wave,flux,indexfile=path+'indexlist_hgamma.dat',$
;        /nobreaks,/silent,/postscript);,position=pos[*,0])
;      djs_oplot, 5577.0*[1,1]/(1.0+ispec[jj].z), !y.crange, line=0, thick=10.0, color='red'
       
    endfor

    dfpsclose
    spawn, 'ps2pdf '+psfile+' '+pdffile
    spawn, 'pdftk '+pdffile+' cat 1-endW output /tmp/junk.pdf'
    spawn, '/bin/mv -f /tmp/junk.pdf '+pdffile
;   spawn, '/bin/rm -f '+psfile
;   cleanplot, /silent
    im_plotfaves
    
stop    
    
return
end
    
