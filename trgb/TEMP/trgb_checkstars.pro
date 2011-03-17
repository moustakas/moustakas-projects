; jm00mar22ucb
; modified: jm00may24ucb

pro trgb_checkstars, image, masterdir, apdata, newdata, chimatch
; this program updates the .exc file created by matchstars and creates
; a new exclude file called .excnew.  interactive.

	window, 0, xs=450, ys=450
	window, 1, xs=450, ys=450
        window, 2, xs=450, ys=450
        window, 3, xs=450, ys=450
        plotsym, 0, 1, /fill
;       loadct, 3

; read in the image with the neighbors subtracted

        rawim = readfits(masterdir+'/'+image+'.fits',head)
	im = readfits(image+'.sub.fits')

        sz = size(rawim)
        nx = sz[1]
        ny = sz[2]

        aprad = [4, 6, 8, 10, 12, 14, 16, 18, 20]	; aperture radii

        bad = ' '
        keep = lonarr(n_elements(apdata[*,0]))

        print & print, 'To remove a star from the list type x, or'
        print, 'press RETURN to go on to the next star.' & print

        for k = 0L, n_elements(apdata[*,0])-1L do begin

            print
            print, 'ALLFRAME chi: ', strn(chimatch[k])
            print

            xc = apdata[k,1]
            yc = apdata[k,2]

; cut out a section of the image by 50 pixels

            imsub = imcut(im,xc,yc,50)
;           imsub = hist_equal(imsub)
;           imsub = alog10(imsub)
;           im = alog10(im)
;           rawim = alog10(rawim)
            
            wset, 0	; subimage
            plotimage, bytscl(imsub)
            tvcircle, 4, 25, 25, /data	;  4-pixel aperture
            tvcircle, 20, 25, 25, /data ; 20-pixel aperture

            wset, 1	; star-subtracted image
            display, im, min=0, max=max(im)
            tvcircle, 20, xc, yc, /data ; 20-pixel aperture

            wset, 2	; curve of growth
            plot, aprad, apdata[k,3:11], ps=8, xsty=3, ysty=3, $
              xtit='Aperture Radii', ytit='Relative Magnitude', $
              tit='Curve of Growth', thick=2, xthick=2, ythick=2, $
              charsize=2, charthick=2

            wset, 3	; raw image
            display, rawim, min=0, max=max(rawim)
            tvcircle, 20, xc, yc, /data ; 20-pixel aperture

            read, bad
            if strupcase(bad) eq 'X' then keep[k] = -1L else keep[k] = k
        endfor

        keep = keep[where(keep ne -1L)]

        newdata = apdata[keep,*] ; good photometry

        wdelete, 0
        wdelete, 1
        wdelete, 2
        wdelete, 3

return
end
