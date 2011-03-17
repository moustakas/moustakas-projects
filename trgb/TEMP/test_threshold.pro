; jm00june14ucb

; test what changing the threshold has on the luminosity function of
; the stars found.

pro test_threshold

	path = '/deepscr1/ioannis/trgb/ugc07577/CHIP3/TEST/'
        obj = 'ugc07577_I1_3'

        openr, lun1, path+obj+'.ap_2', /get_lun
        openw, lun2, path+'lfunction_2.dat', /get_lun
        
        temp = strarr(4)
        data = fltarr(4,2)
        
        readf, lun1, temp	; read the header

        while not eof(lun1) do begin
            readf, lun1, data
            printf, lun2, data[0,1] ; write the magnitude of the star
        endwhile

        free_lun, lun1, lun2

; --------

        openr, lun3, path+obj+'.ap_20', /get_lun
        openw, lun4, path+'lfunction_20.dat', /get_lun

        readf, lun3, temp

        while not eof(lun3) do begin
            readf, lun3, data
            printf, lun4, data[0,1]
        endwhile

        free_lun, lun3, lun4

; --------

        openr, lun5, path+obj+'.ap_5', /get_lun
        openw, lun6, path+'lfunction_5.dat', /get_lun

        readf, lun5, temp

        while not eof(lun5) do begin
            readf, lun5, data
            printf, lun6, data[0,1]
        endwhile

        free_lun, lun5, lun6

; read the luminosity functions and plot them

        read1col, path+'lfunction_2.dat', lfunc2
        read1col, path+'lfunction_5.dat', lfunc5
        read1col, path+'lfunction_20.dat', lfunc20

        colortable1
        window, 0, xs=450, ys=450
        plothist, lfunc2, xsty=3, ysty=3, line=0, thick=3, color=3, xtit='Magnitude', $
          ytit='Number of Stars', tit='Stellar Luminosity Function'
        plothist, lfunc5, line=1, thick=2, color=7, /overplot
        plothist, lfunc20, line=2, thick=2, color=5, /overplot
        legend, ['Threshold = 2', 'Threshold = 5', 'Threshold = 20'], line=[0,1,2], color=[3,7,5], thick=2, $
          /right

stop
return
end
