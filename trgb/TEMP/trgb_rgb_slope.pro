pro make_lfunction, imags, binsize, minmag, maxmag, hist, xmag

        hist = histogram(imags,binsize=binsize,min=minmag)
        nhist = n_elements(hist)
        xnorm = max(findgen(nhist))/(maxmag-minmag) ; magnitude array
        xmag = findgen(nhist)/xnorm+18.

return
end

; jm00june20ucb

; compare the slope of the RGB population in Leo I, Sextans B and
; UGC07577


pro trgb_rgb_slope

	ps_open, '/deepscr1/ioannis/slope'

        !x.ticklen=0.03
        p_thick = !p.thick & chars = !p.charsize
        charth = !p.charthick & xth = !x.thick & yth = !y.thick
        !p.thick = 2 & !p.charsize=2 & !p.charthick=2 & !x.thick=2 & !y.thick=2

	colortable1

        binsize = 0.1
        minmag = 18.
        maxmag = 26.

	cd, '/deepscr1/ioannis/trgb/ugc07577'
        info = sread('/deep1/ioannis/trgb/hst_object.dat')
        indx = where(strupcase(info.object) eq 'UGC07577')

        readfast, 'ugc07577_starlist.dat', data, ncols=5
        imags=data[3,*] - info[indx].a_i        

;       window, 0, xs=800, ys=800
        make_lfunction, imags, binsize, minmag, maxmag, lfunction, xmag

        !p.multi = [0,1,3]
        plot, xmag, lfunction, ps=10, line=0, color=7, $
          xsty=3, ysty=3, /ylog, xr=[20,27], xminor=2, $
          tit='UGC 07577 Luminosity Function', $
          ytit = 'Number Density', xtit='Magnitude'

; --------------------------------------------------------------------------------
        
        info = sread('/deep1/ioannis/trgb/keck_object.dat')

        cd, '/deep3/marc/trgb/data/22dec97/sextansb'
        indx = where(strupcase(info.object) eq 'SEXTANSB')

        readfast, 'SextansB_IVmags.dat', data, skip=2, ncols=8
        imags=data[3,*] - info[indx].a_i        

        make_lfunction, imags, binsize, minmag, maxmag, lfunction, xmag

        plot, xmag, lfunction, ps=10, line=0, color=5, $
          xsty=3, ysty=3, /ylog, xr=[20,27], xminor=2, $
          tit='Sextans B Luminosity Function', $
          ytit = 'Number Density', xtit='Magnitude'
        
; --------------------------------------------------------------------------------
        
        cd, '/deep3/marc/trgb/data/23dec97/leoi'
        indx = where(strupcase(info.object) eq 'LEOI')

        readfast, 'leoi_IVmags.dat', data, skip=2, ncols=8
        imags=data[3,*] - info[indx].a_i        

        make_lfunction, imags, binsize, minmag, maxmag, lfunction, xmag

        plot, xmag, lfunction, ps=10, line=0, color=5, $
          xsty=3, ysty=3, /ylog, xr=[19,26], xminor=2, $
          tit='Leo I Luminosity Function', $
          ytit = 'Number Density', xtit='Magnitude'
        
        !p.multi=0

        !x.ticklen=0
        !p.thick = p_thick & !p.charsize=chars & !p.charthick=charth
        !x.thick = xth & !y.thick=yth

        ps_close

stop
return
end
