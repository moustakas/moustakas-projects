; jm00nov13uofa
pro plot_smflf, name, lfunction, phi, log=log

	if keyword_set(log) then begin
            lf = alog10(float(phi)>1.)
            ytitle = 'log N(m)'
        endif else begin
            lf = phi
            ytitle = 'N(m)'
        endelse

        nstars = n_elements(lfunction.mags)

        plot, lfunction.mag, lf, line=0, yminor=3, color=7, $
          xsty=3, ysty=3, xminor=3, tit=name+' Luminosity Function', $
          ytitle=ytitle, xtit='I', yrange = [0,max(lf)*1.1]
        legend, ['Stars '+strn(nstars), 'Binsize '+strn(lfunction.binsize,format='(F5.2)')], $
          charsize = 1.5, textcolor=[1,1], box=0

return
end
