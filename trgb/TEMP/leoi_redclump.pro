function rcfunc, mag, p
;	rgbfit = p[0] + p[1]*(mag-p[5]) + p[2]*(mag-p[5])^2
	rcfit = (p[0]/p[1]*sqrt(2.*!pi)) * exp(-0.5*((mag-p[2])/p[1])^2)
        fit = rcfit ; + rgbfit
        return, fit
end

; jm00june23ucb

; program to find the distance to leoi I from the red clump

pro leoi_redclump

        !x.ticklen=0.03
        p_thick = !p.thick & chars = !p.charsize
        charth = !p.charthick & xth = !x.thick & yth = !y.thick
        !p.thick = 2 & !p.charsize=1.8 & !p.charthick=2 & !x.thick=2 & !y.thick=2

	colortable1

        cd, '/deep3/marc/trgb/data/23dec97/leoi'
        datapath = '/deep3/marc/trgb/data/23dec97/leoi'

        info = sread('/deep1/ioannis/trgb/keck_object.dat')
        indx = where(strupcase(info.object) eq 'LEOI')
            
        filename = 'leoi_IVmags.dat'
        readfast, filename, data, skip=2, ncols=8

        nstars = n_elements(data[0,*])

        imags = data[3,*] - info[indx].a_i   ; apply extinction correction
        ierr = data[4,*]
        color = data[7,*] - info[indx].e_v_i ; color excess correction

        colorstart = -1.
        colorend = 2.
        magstart = 20.
        magend = 22.
        
        window, 0, xs=550, ys=550
        plot, color, imags, ps=3, color=7, xsty=3, ysty=3, syms=1.2, $
          xtit = 'V-I', ytit='I', xr=[-1,4], yr=[max(imags),min(imags)], $
          title = 'Leo I CMD'
        xyouts, [0.25,0.25], [0.83,0.83], strn(nstars)+' Stars', /normal, $
          charthick=2, charsize=1.5, color=1
        oplot, [colorstart,colorstart], [!y.crange[0],!y.crange[1]], line=2, color=2, thick=2.5
        oplot, [colorend,colorend], [!y.crange[0],!y.crange[1]], line=2, color=2, thick=2.5
        oplot, [!x.crange[0],!x.crange[1]], [magstart,magstart], line=2, color=2, thick=2.5
        oplot, [!x.crange[0],!x.crange[1]], [magend,magend], line=2, color=2, thick=2.5

        rc = where((color gt colorstart) and (color lt colorend) and $
                    (imags gt magstart) and (imags lt magend))

; create the luminosity function

        binsize = 0.05

        lfunction = histogram(imags[rc],binsize=binsize,min=magstart,max=magend)
;       lfuncerr = histogram(ierr[rc],binsize=binsize,min=magstart,max=magend)
        nhist = n_elements(lfunction)
        xnorm = max(findgen(nhist))/(magend-magstart) ; magnitude array
        xmag = findgen(nhist)/xnorm+magstart

        window, 2, xs=450, ys=450
        plot, xmag, lfunction, ps=10, line=0, color=7, $
          xsty=3, ysty=3, $; xr=[magstart,magend], $
          yminor=1, tit='Leo I Red Clump LF', $
          ytit = 'Number Density', xtit='Magnitude'

;       rcguess = fltarr(6)
        rcguess = fltarr(3)
;        rcguess[0] = 1.
;        rcguess[1] = 1.
;        rcguess[2] = 1.
        rcguess[0] = 80.	; amplitude
        rcguess[1] = 0.5	; dispersion
        rcguess[2] = 21.	; centroid

        err = float(lfunction-lfunction)+1.
        rclump = mpfitfun('rcfunc',xmag,lfunction,err,rcguess,/quiet)

;	rgb = rclump[0] + rclump[1]*(xmag-rclump[5]) + rclump[2]*(xmag-rclump[5])^2
	gauss = (rclump[0]/rclump[1]*sqrt(2.*!pi)) * exp(-0.5*((xmag-rclump[2])/rclump[1])^2)
        model = gauss ; + rgb
        oplot, xmag, model, line=0, color=5, thick=2

        mabs_rc = -0.23	; red clump absolute magnitude

        modulus = rclump[2]-mabs_rc
        distance = 10^((modulus-25.)/5.)
        
        print, 'Distance = '+strn(distance)
        
        !x.ticklen=0
        !p.thick = p_thick & !p.charsize=chars & !p.charthick=charth
        !x.thick = xth & !y.thick=yth

        stop

return
end

