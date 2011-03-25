pro qaplot_ages_zptoffsets
; jm10aug31ucsd - compare the observed and synthetic magnitudes
; jm11feb09ucsd - updated to just use the K-correct results
    
    zptpath = ages_path(/mycatalogs)+'zptoffsets/'
    qapath = zptpath

; --------------------------------------------------
; show the dispersion in the measurements and the residuals vs
; magnitude and redshift  
    zmin = 0.05
    zmax = 0.70
    levels = [0.1,0.25,0.5,0.75,0.9,0.95]
    npix = 20
    rrange1 = 2.9*[-1,1]

; read the zeropoint tweaks and the K-corrections
    zptfile = zptpath+'ages_zptoffsets.fits.gz'
    datafile = repstr(zptfile,'.fits.gz','_data.fits.gz')
    zpt = mrdfits(zptfile,1)
    zptoffset = total(zpt.zptoffset,2)
    
    kcorr_before = mrdfits(datafile,1)
    kcorr_after = mrdfits(datafile,2)
    ngal = n_elements(kcorr_before)

    filt = (bootes_filterlist())[0:7] ; no irac
    band = ['U','B_{W}','R','I','z','J','H','K_{s}']
    nfilt = n_elements(filt)

; make the plot
    psfile = qapath+'qa_zptoffsets_ages.ps'
    im_plotconfig, 4, pos, psfile=psfile, xmargin=[1.4,0.3], $
      height=[3.0,2.5,2.5], yspace=[0.0,1.0], width=6.8, $
      charsize=1.6, ymargin=[0.5,1.1]
    
    for jj = 0, nfilt-1 do begin
       for pp = 0, 1 do begin
          if (pp eq 0) then begin
             kcorr = kcorr_before
             zptlabel = [band[jj],'Observed Photometry']
          endif else begin
             kcorr = kcorr_after
             zptlabel = [band[jj],'\Delta'+'Zpt = '+strtrim(string(zptoffset[jj],format='(F12.3)'),2)]
          endelse
          
          good = where((kcorr.k_maggies[jj] gt 0),ngood)
          zobj = kcorr[good].k_zobj
          xx = -2.5*alog10(kcorr[good].k_maggies[jj])
          yy = -2.5*alog10(kcorr[good].k_bestmaggies[jj]) ; k-correct
          zptlabel[0] = zptlabel[0]+': '+im_string_stats(yy-xx,ndecimal=3)
          
          if (pp eq 0) then begin
             mrange1 = [im_min(xx,sigrej=10)<im_min(yy,sigrej=5),$
               (im_max(xx,sigrej=5)>im_max(yy,sigrej=5))<30]
             rrange1 = im_max(yy-xx,sigrej=5)*[-1,1]
          endif
          
          hogg_scatterplot, xx, yy, position=pos[*,0], xsty=1, ysty=1, $
            xrange=mrange1, yrange=mrange1, /outliers, outcolor=djs_icolor('blue'), $
            levels=levels, xtitle='', xtickname=replicate(' ',10), $
            ytitle=band[jj]+' (synth, AB mag)', xnpix=npix, ynpix=npix, yminor=2
          djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
          im_legend, ['N='+string(ngood,format='(I0)')+'/'+$
            string(ngal,format='(I0)')], /right, /bottom, $
            box=0, charsize=1.3, margin=0

          im_legend, [zptlabel], /left, /top, box=0, charsize=1.3, margin=0

; residual vs mag
          hogg_scatterplot, xx, yy-xx, position=pos[*,1], /noerase, $
            xsty=1, ysty=1, xrange=mrange1, yrange=rrange1, /outliers, $
            levels=levels, outcolor=djs_icolor('blue'), $
            xtitle=band[jj]+' (observed, AB mag)', ytitle='Residuals (mag)', $
            xnpix=npix, ynpix=npix, /internal, yminor=2
          djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'

; residuals vs redshift
          hogg_scatterplot, zobj, yy-xx, position=pos[*,2], /noerase, $
            xsty=3, ysty=1, xrange=[zmin,zmax], yrange=rrange1, $
            levels=levels, xtitle='Redshift', ytitle='Residuals (mag)', $
            /outliers, outcolor=djs_icolor('blue'), xnpix=npix, ynpix=npix, $
            /internal, yminor=2
          djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
       endfor  
    endfor 
    im_plotconfig, /psclose, psfile=psfile, /gzip

return
end
