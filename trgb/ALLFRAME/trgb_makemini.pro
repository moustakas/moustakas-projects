pro trgb_makemini, objname
;+
; NAME:
;	TRGB_MAKEMINI
;
; PURPOSE:
;	Interactively select a miniature photometric starlist.
;
; INPUTS:
;	objname : galaxy name
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; PROCEDURES USED:
;	RDTXT()
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 May 25, UCB, based on B. Mendez's code
;-

	spawn, ['pwd'], datapath	; current directory
        datapath = datapath[0]

        daohead = rdtxt(datapath+'/'+objname+'.mag',nlines=2) ; header
	readfast, datapath+'/'+objname+'.mag', data, skip=3, ncols=9
        mag = data[3,*]

        nstars = n_elements(mag)

        xc = data[1,*]
        yc = data[2,*]

	okay = 'Q'
        while strupcase(okay) eq 'Q' do begin

            window, 0, xs=450, ys=450
            plot, xc, yc, ps=3, xsty=1, ysty=1, $
              xtit='X-coordinate', ytit='Y-coordinate', $
              thick=2, ythick=2, xthick=2, charsize=1.5, $
              charthick=2, tit=strupcase(objname)

            print & print, 'Please select a region to create a mini starlist.' & print

            button = 0
            repeat begin
                print & print, 'Click on the lower left corner of your box. '
                cursor, x1, y1, /data, /down
                button = !mouse.button
            endrep until (button eq 1 or button eq 2)
            
            repeat begin
                print & print, 'And now click on the upper right corner to close the box.'
                cursor, x2, y2, /data, /down
                button = !mouse.button
            endrep until (button eq 1 or button eq 2)
            
            flag = bytarr(nstars)
            for i = 0, nstars-1 do $
              if (xc[i] gt x1 and xc[i] lt x2 and yc[i] gt y1 and yc[i] lt y2) then $
              flag[i] = 1B

            ndata = data[*,where(flag EQ 1B)]

            oplot, ndata[1,*], ndata[2,*], ps=2

            hist = histogram(mag,bin=0.5,min=5.)	; mag histogram
            nrhist = n_elements(hist) 			; number of elements
            marray = findgen(nrhist)/2.+5. 		; magnitude array

            window, 2, xs=450, ys=450
            plot, marray, hist, xtit='Relative Magnitude', $
              ytit='Number of Stars', thick=2, ythick=2, $
              xthick=2, charsize=1.5, charthick=2, $
              tit=strupcase(objname)+' Magnitude Histogram', $
              xsty=3, ysty=3
            xyouts, [0.45, 0.45], [0.85, 0.85], $
              '# Stars: '+strn(n_elements(ndata[1,*])), $
              /normal, align=0.5, charsize=2, charthick=2

            nmag = ndata[3,*]

            nhist = histogram(nmag,bin=0.5,min=5.)
            nnrhist = n_elements(nhist)
            nmarray = findgen(nnrhist)/2.+5.

            oplot, nmarray, nhist

            read, okay, prompt='Press ENTER to continue or Q to re-select a mini starlist: '

        endwhile

        outfile = objname+'_mini.mag'
        print
        print, 'Writing ', outfile+'.'

        frmt =  '(1x,I5,1x,F8.3,1x,F8.3,1x,F8.3,1x,F8.4,1x,F8.2,1x,F8.0,1x,F8.2,1x,F8.3)'
        openw, lun1, outfile, /get_lun
        printf, lun1, daohead[0]
        printf, lun1, daohead[1]
        printf, lun1, ' '
        for k = 0, n_elements(ndata[0,*])-1L do $
          printf, lun1, ndata[*,k], format=frmt
        free_lun, lun1

return
end
