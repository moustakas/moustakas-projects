; jm00june15ucb
; routine to read in column-formated data very quickly.  at this point
; the routine is not very intelligent.  everything is assumed floating point.

pro readfast, filename, data, header, skipline=skipline, ncols=ncols

	on_error, 2 ; return to user

;	t1 = systime(1)

        nrows = 0L
        spawn, ['wc -l '+filename], wcout	; file size
        reads, wcout, nrows

        openr, lun1, filename, /get_lun

        if keyword_set(skipline) then begin

            nrows = nrows-skipline 
            header = strarr(skipline)
            readf, lun1, header	; read the header

        endif

        data = fltarr(ncols,nrows)

        readf, lun1, data       ; read the data
        free_lun, lun1

;	print, systime(1)-t1

return
end


; routine to combine the final calibrated starlists from all four
; wfpc2 chips

; should be run from the master object directory

; the master lists for each chip should be called objname_IVmags.dat
; (see trgb_hst_calibrate)

; if vdata is *not* set then the routine assumes there is no v-band data

pro trgb_hst_combine_starlists2, objname

;	spawn, ['clear']

	cd, '/deepscr1/ioannis/trgb/'+objname
        datapath = '/deepscr1/ioannis/trgb/'+objname

        obj_vdata = ['UGC03476','UGC03755','IC342'] ; HST objects with color
        check = where(strupcase(obj_vdata) eq strupcase(objname),count)
        if count gt 0L then ncol = 8 else ncol = 5

  ;      pcscale = 0.045	/60. ; PC plate scale (arcmin/pixel)
  ;      wfcscale = 0.1 /60. ; WFC plate scale (arcmin/pixel)
  ;      rotchips = [0,1,2,3] ; ROTATE values for each chip
	pcscale=1.
	wfcscale=1.
	rotchips=[0,0,0,0]
; do not create a coordinate system for WFPC2 (1600x1600)

; chip 1

        readfast, datapath+'/CHIP1/'+objname+'_IVmags.dat', chip1, skip=2, ncols=ncol
        chip1[1:2,*] = chip1[1:2,*] * pcscale-1.
	nchip1=chip1[1,*]*0
        print,min(chip1[1:2,*]),max(chip1[1:2,*])
; chip 2

        readfast, datapath+'/CHIP2/'+objname+'_IVmags.dat', chip2, skip=2, ncols=ncol
        chip2[1:2,*] = chip2[1:2,*] * wfcscale-1.
        ;chip2[1,*] = -temp2[1,*]
        ;chip2[2,*] = temp2[0,*]
        nchip2=chip2[1,*]*0+1
; chip 3

        readfast, datapath+'/CHIP3/'+objname+'_IVmags.dat', chip3, skip=2, ncols=ncol
        chip3[1:2,*] = chip3[1:2,*] * wfcscale-1.
	nchip3=chip3[1,*]*0+2
; chip 4

        readfast, datapath+'/CHIP4/'+objname+'_IVmags.dat', chip4, skip=2, ncols=ncol
        chip4[1:2,*] = chip4[1:2,*] * wfcscale-1.
        ;chip4[1,*] = temp4[1,*]
        ;chip4[2,*] = -temp4[0,*]
	nchip4=chip4[1,*]*0+3
        starlist = [[chip1],[chip2],[chip3],[chip4]]
	chipnums=[[nchip1],[nchip2],[nchip3],[nchip4]]
	help,chipnums
        nout = n_elements(starlist[0,*])
        print, 'Total number of stars: ', nout

;       plot, starlist[1,*], starlist[2,*], ps=3, xsty=1, ysty=1

        openw, lun1, datapath+'/'+objname+'_starlist2.dat', /get_lun

        if count gt 0L then $
          printf, lun1, '#  ID   Xcenter   Ycenter    I       Ierr      V       Verr     V-I   Chip no.' else $
          printf, lun1, '#  ID   Xcenter   Ycenter    I       Ierr   Chip no.' 
        printf, lun1, ' '
        for k = 0L, nout-1L do printf, lun1, starlist[*,k],chipnums[k], format = '(1x,I5,7F9.3,I5)'
        free_lun, lun1

return
end





