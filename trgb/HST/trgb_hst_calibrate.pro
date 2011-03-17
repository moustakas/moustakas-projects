function newttrgb, x
	common trgb, imag, vmag, iapcor, vapcor
        return, [ imag - 0.063*(x[1]-x[0]) + 0.025*(x[1]-x[0])^2 + iapcor, $
                  vmag - 0.052*(x[1]-x[0]) + 0.027*(x[1]-x[0])^2 + vapcor]
;        return, [ imag - 0.063*(x[1]-x[0]) + iapcor, $
;                  vmag - 0.052*(x[1]-x[0]) + vapcor]
end

; jm00june10ucb

; apply calibrations to the data

; keywords:
;        vdata: indicates that v-band data are present
;        chip:  pc, wfc1, wfc2, wfc3
 

pro calibrate, chip, data, itrue, vtrue, vdata=vdata

	common trgb, imag, vmag, iapcor, vapcor

        case strupcase(chip) of 

            'PC': begin

                iapcor = -0.17
                vapcor = -0.17
                zmagi = 21.585
                zmagv = 22.469
                
            endcase
            
            'WFC2': begin

                iapcor = -0.00
                vapcor = -0.04
                zmagi = 21.594
                zmagv = 22.478

            endcase

            'WFC3': begin

                iapcor = 0.1
                vapcor = -0.1
                zmagi = 21.596
                zmagv = 22.480
                
            endcase

            'WFC4': begin
                
                iapcor = 0.01
                vapcor = 0.00
                zmagi = 21.568
                zmagv = 22.452
                        
            endcase

            else: print, 'You need to indicate a WFPC2 chip!'

        endcase
            
        imags = data[3,*]
        
        if keyword_set(vdata) then begin
            
            vmags = data[5,*]
            
            x = dblarr(2)
            itrue = imags-imags
            vtrue = vmags-vmags
            
            for k = 0L, n_elements(imags)-1L do begin

                imag = imags[k]	- (25.-zmagi) ; instrumental magnitude
                vmag = vmags[k] - (25.-zmagv)

                x[0] = imags[k] ; initial guess is the instrumental magnitude
                x[1] = vmags[k]

                for h = 0L, 9L do x = newttrgb(x) ; iterate to solve

;                true = newton(x, 'newttrgb', itmax=50, $
;                               stepmax = 100.0, tolf=1.E-10, $
;                               tolmin=1.E-6, tolx=1.E-7, /double)
                itrue[k] = x[0]
                vtrue[k] = x[1]

            endfor

        endif else begin

; no color information: use 1.42 as the TRGB color (Key Project value)

            itrue = imags - (25. - zmagi) - 0.063*1.3 + 0.025*1.3^2 + iapcor
            
        endelse

;        window, 0, xs=450, ys=450
;        plothist, itrue, bin=0.1, col=7, xr=[22,32], thick=2, line=0
;        plothist, imags, /overplot, bin=0.1, thick=2, color=5, line=2
;        legend, ['Calibrated Mags', 'Instrumental Mags'], line=[0,2], color=[7,5], thick=2
;        stop

return
end

function wfpc2_cte, alfdata
; jm00june15ucb
; apply the charge-transfer efficiency (CTE) correction

	xcen = alfdata[1,*]	; x-centroid
	ycen = alfdata[2,*]	; y-centroid
	mags = alfdata[3,*]	; magnitude
        sky  = alfdata[5,*]	; modal sky value

; convert DAOPHOT magnitude to ADU then electrons

        gain = 7D

        flux = gain*10D^(0.4D*(double(25)-mags))

        skyflux = gain*sky
        for j = 0L, n_elements(skyflux[0,*])-1L do $
          if skyflux[0,j] lt 0. then skyflux[0,j] = 0.

        xnorm = xcen/400D	; stellar position in the 800x800 image
        ynorm = ycen/400D

        starval = alog10(flux)-4D
        skyval = alog10(1.+skyflux)-1D

; CTE coefficients

        A = 0.67 & B = -0.88 & C = -0.75
        D = -1.63 & E = -0.25 & F = -0.09

        xcte = exp(A + C*starval + E*skyval)
        ycte = exp(B + D*starval + F*skyval)
        
        trueflux = flux * (1.+0.01*ycte*ynorm+0.01*xcte*xnorm)

        truemags = mags - 2.5*alog10(1.+0.01*ycte*ynorm+0.01*xcte*xnorm)

;        colortable1
;        plothist, truemags, line=0, bin=0.1, color=3, thick=2, $
;          xtit='Instrumental Magnitude', ytit='Number of Stars', $
;          tit='CTE Effect'
;        plothist, mags, /overplot, line=2, color=7, bin=0.1, thick=2
;        legend, ['CTE Corrected Mags', 'Instrumental Mags'], line=[0,2], color=[3,7], thick=2
;
;; residual plot
;        plot, truemags, mags-truemags, ps=3, xsty=1, ysty=1, yr=[0,0.2], color=5, $
;          thick=2
;        diff  = mags-truemags
;        indx = where((truemags gt 19.) and (truemags lt 20))
;        print, median(diff[indx])
;        stop

return, truemags
end


pro trgb_hst_calibrate, objname, chip
;+
; NAME:
;	TRGB_HST_CALIBRATE
;
; PURPOSE:
;	Calibrate HST/WFPC2 ALLFRAME photometry onto a standard
;	system.  Applies the Holtzmann 1995 zero-points, corrects for
;	the CTE effect, applies exposure-time corrections, and
;	averages all ALLFRAME magnitudes.
;
; INPUTS:
;	objname : galaxy name
;	chip	; pc, wfc1, wfc2, wfc3
;
; PROCEDURES USED:
;	RDTXT(), READFAST, HEADFITS(), SXPAR(), TRGB_ROBUST()

; MODIFICATION HISTORY:
;	Written, J. Moustakas, 2000 June 9, UC Berkeley
;-

        npar = n_params()
        if npar ne 1 then begin
            print, 'Syntax: trgb_hst_calibrate, objname, chip'
            return
        endif
        
        spawn, ['pwd'], datapath
        datapath = datapath[0]

        imnames = rdtxt(datapath+'/image.names')        
        nr = n_elements(imnames)

        i_indx = where(strpos(imnames,'_I') gt 0,ni)
        if ni gt 0 then iband = imnames[i_indx]	else begin
            print & print, 'No I-band frames found!'
            return
        endelse

        v_indx = where(strpos(imnames,'_V') gt 0,nv)
        if nv gt 0 then vband = imnames[v_indx] else begin
            print & print, 'No V-band frames.' & print
        endelse
        
; apply the CTE correction to all the stars and create new alf's

        for j = 0L, nr-1L do begin
            readfast, datapath+'/'+imnames[j]+'.alf', alfdata, header, skip=3, ncols=9

            truemags = wfpc2_cte(alfdata)
            newdata = [alfdata[0:2,*],truemags[0,*],alfdata[4:8,*]]

            openw, lun1, datapath+'/'+imnames[j]+'.alf_cte', /get_lun
            for i = 0L, 2L do printf, lun1, header[i]
            for k = 0L, n_elements(newdata[0,*])-1L do $
              printf, lun1, newdata[*,k], format = '(1x,I5,8F9.3)'
            free_lun, lun1
	endfor

        trgb_makemch, objname, /alf, /hst ; update the transformation file
        trgb_daomaster, objname, nr, /verbose ; create a new raw file

; read the new alfs and robustly average all the frames
            
        readfast, objname+'.raw', rawdata, ncol=2*ni+(2L*nv+4L)+1L, skip=3
        nstars = n_elements(rawdata[0,*]) ; total number of stars

; I-band
        
        if ni eq 1L then begin	; only one I-band frame

            imags = rawdata[5,*]   ; raw I-band data
            ierr = rawdata[6,*] ; raw I-band errors

            ihdr = headfits(iband[0]+'.fits')
            itime = sxpar(ihdr,'EXPTIME')	; exposure time

            imags = imags + 2.5*alog10(itime)

        endif else begin

            icols = 2*indgen(ni)+3L 	  ; I-band magnitude columns
            iraw = rawdata[icols,*] 	  ; raw I-band data
            ierrraw = rawdata[icols+1L,*] ; raw I-band errors
        
; calibrate

            for q = 0L, ni-1L do begin
                
                ihdr = headfits(iband[q]+'.fits')
                itime = sxpar(ihdr,'EXPTIME') ; exposure time
                
                iraw[q,*] = iraw[q,*] + 2.5*alog10(itime)
                
            endfor
            
            imags = fltarr(nstars)	; average the data
            ierr = fltarr(nstars)
            for i = 0L, nstars-1L do imags[i] = trgb_robust(iraw[*,i],ierrraw[*,i])
            for i = 0L, nstars-1L do ierr[i] = 1./(sqrt(total(1./(ierrraw[*,i])^2)))

;; rms dispersion
;
;            window, 0, xs=450, ys=450
;            window, 10, xs=450, ys=450
;            for h = 1L, ni-1L do begin
;                wset, 0
;                plothist, imags-7.7, bin=0.1, col=7, xr=[22,32]-7.7, thick=2, line=0
;                plothist, iraw[0,*]-7.7, /overplot, bin=0.1, thick=2, color=5, line=2
;                legend, ['Calibrated Mags', 'Instrumental Mags'], line=[0,2], color=[7,5], thick=2
;
;                wset, 10
;                plot, iraw[0,*], imags-iraw[0,*], ps=2, xsty=1, ysty=1, syms=0.3, yr=[-1,1], $
;                  xr=[22,30]
;                oplot, iraw[h,*], imags-iraw[h,*], ps=2, syms=0.2
;            endfor
; stop
        endelse

; V-band

        if nv gt 0L then begin

            if nv eq 1L then begin
                
                vcols = 2*ni+5L 	   ; V-band columns
                vmags =  rawdata[vcols,*]  ; raw V-band data
                verr = rawdata[vcols+1L,*] ; raw V-band magnitude errors

                vhdr = headfits(vband[0]+'.fits')
                vtime = sxpar(vhdr,'EXPTIME') ; exposure time

                vmags = vmags + 2.5*alog10(vtime)

            endif else begin
                
                vcols = 2*ni+2*indgen(nv)+3L
                vraw =  rawdata[vcols,*]
                verrraw = rawdata[vcols+1L,*]

; add the exposure time
                
                for q = 0L, nv-1L do begin
                    
                    vhdr = headfits(vband[q]+'.fits')
                    vtime = sxpar(vhdr,'EXPTIME') ; exposure time
                    
                    vraw[q,*] = vraw[q,*] + 2.5*alog10(vtime)
                
                endfor

                vmags = fltarr(nstars)
                verr = fltarr(nstars)
                for j = 0L, nstars-1L do vmags[j] = trgb_robust(vraw[*,j],verrraw[*,j])
                for j = 0L, nstars-1L do verr[j] = 1./(sqrt(total(1./(verrraw[*,j])^2)))

            endelse

; flag and filter bad stars

            flag = bytarr(nstars)+1B

            for k = 0L, nstars-1L do begin
                if (imags[k] gt 30.) then flag[k] = 0B ; I mag > 30
                if (vmags[k] gt 30.) then flag[k] = 0B ; V mag > 30
                if (ierr[k] gt 1.) then flag[k] = 0B   ; I mag error > 1
                if (verr[k] gt 1.) then flag[k] = 0B   ; V mag error > 1
                if (rawdata[1,k] lt 1. or rawdata[2,k] lt 1.) then flag[k] = 0B ; too close to the chip edge
                if (rawdata[1,k] gt 800-1. or rawdata[2,k] gt 800-1.) then flag[k] = 0B
            endfor

            data = [[rawdata[0:2,*]],transpose(imags),transpose(ierr),$
                    transpose(vmags),transpose(verr),transpose(vmags-imags)]
            data = data[*,where(flag eq 1B)]

; apply HST/WFPC2 calibrations

            calibrate, chip, data, itrue, vtrue, /vdata
            dataout = [data[0:2,*],itrue,data[4,*],vtrue,data[6,*],vtrue-itrue]
        
        endif else begin

            flag = bytarr(nstars)+1B

            for k = 0L, nstars-1L do begin
                if (imags[k] gt 30.) then flag[k] = 0B ; I mag > 30
                if (ierr[k] gt 1.) then flag[k] = 0B   ; error > 1
                if (rawdata[1,k] lt 1. or rawdata[2,k] lt 1.) then flag[k] = 0B ; too close to the chip edge
                if (rawdata[1,k] gt 800-1. or rawdata[2,k] gt 800-1.) then flag[k] = 0B
            endfor

            data = [[rawdata[0:2,*]],transpose(imags),transpose(ierr)]
            data = data[*,where(flag eq 1B)]

; apply HST/WFPC2 calibrations

            calibrate, chip, data, itrue
            dataout = [data[0:2,*],itrue,data[4,*]]

        endelse

        nout = n_elements(dataout[0,*])

; write the final I- and V-band magnitudes, errors, and V-I colors to a file:

        print & print, 'Writing '+objname+'_IVmags.dat . . .' & print

        openw, lun1, datapath[0]+'/'+objname+'_IVmags.dat', /get_lun

        if nv gt 0L then $
          printf, lun1, '#  ID   Xcenter   Ycenter    I       Ierr      V       Verr     V-I' else $
          printf, lun1, '#  ID   Xcenter   Ycenter    I       Ierr'
        printf, lun1, ' '
        for k = 0L, nout-1L do printf, lun1, dataout[*,k], format = '(1x,I5,7F9.3)'

        free_lun, lun1

;        window, 2, xs=450, ys=450
;        plot, dataout[3,*], dataout[4,*], ps=3, yr=[0,1], $
;          xtit='I', ytit='Magnitude Error', tit=strupcase(objname)+' Diagnostic', $
;          thick=2, charthick=1.5, charsize=1.5, xthick=1.5, ythick=1.5, $
;          xsty=3, ysty=3

;       if nv gt 0L then trgb_cmd, objname, dataout ; color-magnitude diagram
;       trgb_lfunction, objname, transpose(dataout), /plotit ; luminosity function

;; some testing
;
;        window, 0, xs=450, ys=450
;        colortable1
;        plothist, rawdata[3,*], bin=0.1, col=7, xr=[15,25], thick=2, line=0
;        plothist, rawdata[3,where(flag eq 1B)], bin=0.1, col=5, thick=2, line=2, $
;          /overplot
;        legend, ['Instrumental Mags', 'Filtered Stars'], line=[0,2], color=[7,5], thick=2

;        window, 2, xs=450, ys=450
;        plot, ierr[where(flag eq 1B)], abs(iraw[0,where(flag eq 1B)]-imags[where(flag eq 1B)]), $
;          ps=3, color=3, yr=[0,0.5], xsty=3, ysty=3
;        window, 0, xs=450, ys=450
;        newvar = abs(iraw[0,where(flag eq 1B)]-imags[where(flag eq 1B)])/ierr[where(flag eq 1B)]
;        plothist, newvar, color=5, bin=0.1, xr=[-0.2,5]

return
end
