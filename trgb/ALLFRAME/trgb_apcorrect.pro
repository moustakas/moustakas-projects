pro find_apstars, alfdata
; selection criteria (dolphin 2000):  -0.5<sharp<0.5, chi<2.5, 4.3" =
; 20 pixels away from any other centroid (star)

	nstars = n_elements(alfdata[0,*])

        id = alfdata[0,*]
        xc = alfdata[1,*]
	yc = alfdata[2,*]

        sharp = where((alfdata[7,*] gt -0.5) and (alfdata[7,*] lt 0.5))
        chi = where(alfdata[6,sharp] lt 1.5,ngood)

        print, 'Stars that satisfy sharpness and chi criterion: '+strn(ngood)+'.' & print

        for j = 1L, 10L do begin	; increase the radius criterion
       
            print & print, 'Current aperture radius = '+strn(2.*j)+'.'
            
            for k = 0L, ngood-1L do begin

                flag = where((abs(xc[k]-xc) lt 2.*j) and (abs(yc[k]-yc) lt 2.*j),count)
                match = where(flag eq k,mcount)
                if (mcount eq -1L) and (count le 10) then begin
                    print, 'Star = '+strn(k)+'.'
                    print, 'Radius = '+strn(2.*j)+'.'
                    print, 'Near stars = '+strn(count)+'.'
                endif

            endfor

        endfor

        plot, xc, yc, ps=3, xsty=3, ysty=3
        oplot, xc[flag], yc[flag], ps=2, syms=2, color=2
        tvcircle, 10, xc[k], yc[k], col=5, /data

stop

;; use the nstars from the middle of the list as aperture correction stars
;
;        nstars = 200L
;        ntot = n_elements(alfdata[0,*])
;        nhalf = ntot/2L
;
;        print & print, 'Writing '+objname+'.exc.'
;
;        openw, lun1, objname+'.exc', /get_lun
;        printf, lun1, daohead[0] ; daophot header
;        printf, lun1, daohead[1]
;        printf, lun1, ' '
;        for j = nhalf, nhalf+nstars-1L do printf, lun1, alfdata[0:2,j], format='(1x,I5,2F9.3)'
;        free_lun, lun1
;
;

return
end

pro trgb_apcorrect, objname, verbose=verbose, nosub=nosub
;+
; NAME:
;	TRGB_APCORRECT
;
; PURPOSE:
;	Determines the value of the aperture correction (flux
;	difference between a 20 pixel circular aperture and the flux
;	from psf photometry) for a single image.
;
; INPUTS:
;	The user should be in an otherwise empty subdirectory with
;	the following file:
;		master.dir  : Directory where the fits image and its
;			      associated .psf and .alf files reside.
;
; OUTPUTS:
;	Creates the DAOPHOT .opt files for this image; a .lst file,
;	a coordinate file of the aperture correction stars; a .exc
;	file, a coordinate file of the .lst stars with the ALLFRAME
;	position and magnitudes; a file called image.apcor containing
;	aperture photometry on the aperture correction stars;
;	image.sub.fits, a image with all but the aperture correction
;	stars subtracted; and a file called aperture.correction which
;	contains the results.
;
; PROCEDURE:
;	The program uses IRAF's image.pst file as a list of bright,
;	isolated, unsaturated stars with which to determine the
;	aperture correction (these are the aperture correction stars).
;	The file is converted to a .lst file for use by DAOPHOT.
;	These stars are then matched by ID numbers with the ALLFRAME
;	photometry file which contains superior centroid positions to
;	create an exclude (.exc) file.  Using the ALLFRAME (.alf) file,
;	SUBSTAR is called to subtract all the neighbor stars in the
;	image.  Next, 20 pixel aperture photometry is done on the
;	aperture correction stars, and the difference of the resulting
;	magnitude with the ALLFRAME magnitude is computed.
;
; PROCEDURES USED:
;	TRGB_MAKEOPT, IRAF2DAO, DAOHEADER, READFAST(), READ_DAOPHOT, 
;
; RESTRICTIONS:
;	Any known old .coo.1 or .pst files should be placed in an
;	OLD/ subdirectory or deleted.  Deletes all files created by a
;	previous run of this code (to prevent DAOPHOT from crashing). 
;
; MODIFICATION HISTORY:
;	Written, J. Moustakas, 2000 March 19-22
; 	modified jm00may24ucb
;-

	npar = n_params()
	if npar ne 1 then begin
            print, 'Syntax: apcorrect, objname, verbose=verbose, nosub=nosub'
            return
        endif

	spawn, ['pwd'], datapath	; current directory
        datapath = datapath[0]
        masterdir = rdtxt(datapath+'/master.dir') ; master object directory
        masterdir = masterdir[0]

; create the environment setting and temporary alias to the image

        setenv, 'dir='+masterdir
        spawn, ['alias dir cd $dir']

; delete any old .apcor and .sub.fits files from a previous run of
; this program to prevent DAOPHOT from crashing

        spawn, ['\rm '+objname+'.apcor']
        if not keyword_set(nosub) then spawn, ['\rm '+objname+'.sub.fits']

; create the .opt files for DAOPHOT appropriate to this image

        spawn, ['find '+masterdir+' -name "'+objname+'.pars" -print'], pname
        nold = where(strpos(pname,'OLD') eq -1,count)
        if count ne 0 then pname = pname[nold]
        pname = pname[0]

; find the IRAF psf

        spawn, ['find '+masterdir+' -name "'+objname+'.psf.fits" -print'], irafpsf
        nold = where(strpos(irafpsf,'OLD') eq -1,count)
        if count ne 0 then irafpsf = irafpsf[nold]
        irafpsf = irafpsf[0]

; pop out the fits header for this image

        fhead = headfits(masterdir+'/'+objname+'.fits')
        psfhead = headfits(irafpsf)

; create the .opt files for DAOPHOT appropriate to this image

        print, 'Creating the .opt files.'
        trgb_makeopt, pname, psfhead

        daoheader, pname, fhead, daohead, /lst	; create a daophot header

; sometimes the magerr and the skyval overlap because of the way
; DAOPHOT writes files, so use the filter to figure out how many
; columns to read.

        filter = sxpar(fhead,'FILTER')
        if strupcase(strn(filter)) eq 'I' then cols = 8 else cols = 9

        readfast, masterdir+'/'+objname+'.alf', alfdata, skip=3, ncols=cols

; select aperture correction stars and write them to an exclude file

        find_apstars, alfdata

        if not keyword_set(nosub) then begin

; SUBSTAR script:
            openw, lun1, 'substar_script', /get_lun
            printf, lun1, 'ATTACH dir:'+objname+'.fits'	; image
            printf, lun1, 'SUBSTAR'
            printf, lun1, 'dir:'+objname+'.psf' ; psf
            printf, lun1, 'dir:'+objname+'.alf' ; allframe photometry file
            printf, lun1, 'Y'   		; answer yes to exclude a starlist
            printf, lun1, objname+'.exc' 	; exclude the psf stars
            printf, lun1, objname+'.sub.fits' 	; output subtracted image
            printf, lun1, 'EXIT' 		; exit
            free_lun, lun1

            print & print, 'Subtracting away the neighbor stars . . .' & print
            if keyword_set(verbose) then begin
                spawn, ['daophot < substar_script']
            endif else begin
                spawn, ['daophot < substar_script > substar_log']
                print, '. . . finished.' & beep
            endelse

        endif

; APERTURE PHOTOMETRY (20 pixels) on the aperture correction stars 
        openw, lun2, 'aphot_script', /get_lun
        printf, lun2, 'ATTACH '+objname+'.sub.fits' ; neighbor-subtracted image
;       printf, lun2, 'ATTACH '+objname+'.fits'     ; raw image
        printf, lun2, 'PHOT'
        printf, lun2, 'photo.opt'		; photometry file
        printf, lun2, 'A2 = 6'			; multiple aperture photometry
        printf, lun2, 'A3 = 8'			; generate the curve of growth
        printf, lun2, 'A4 = 10'
        printf, lun2, 'A5 = 12'
        printf, lun2, 'A6 = 14'
        printf, lun2, 'A7 = 16'
        printf, lun2, 'A8 = 18'
        printf, lun2, 'A9 = 20'
        printf, lun2, 'AA = 25'
        printf, lun2, 'AB = 30'
        printf, lun2, 'AC = 35'
        printf, lun2, ' '
        printf, lun2, objname+'.exc' 		; coordinate file
        printf, lun2, objname+'.apcor'		; output aperture photometry file
        printf, lun2, 'EXIT'
        free_lun, lun2

        print & print, 'Doing aperture photometry on the starlist.'
        if keyword_set(verbose) then begin
            spawn, ['daophot < aphot_script']
        endif else spawn, ['daophot < aphot_script > aphot_log']

; read in the aperture photometry file and reorganize it

	naps = 11L  ; naps is the number of apertures
	read_daophot, objname+'.apcor', apdata, naps=naps

; plot (psf_mag - aperture_mag) for the aperture correction stars

        psfmag = alfdata[3,nhalf:nhalf+nstars-1L]
        apmag = apdata[3:3+naps-1L,*]

        delta = apmag-apmag
        meandel = fltarr(naps)
        meddel = meandel-meandel

        for k = 0L, naps-1L do begin

            delta[k,*] = apmag[k,*] - psfmag[0,*]
            goodap = where(delta[k,*] lt 50.)

            meandel[k] = mean(delta[k,goodap])
            meddel[k] = median(delta[k,goodap])

        endfor

        medap = fltarr(naps)
        for j = 0L, naps-1L do medap[j] = median(apmag[j,*])
        meanap = fltarr(naps)
        for j = 0L, naps-1L do meanap[j] = mean(apmag[j,*])

        apertures = [2,6,8,10,12,14,16,18,20,25,30]
            
        colortable1
        plotsym, 0, 1, /fill

        window, 0, xs=450, ys=450
        plot, apertures, meandel, ps=8, xsty=3, ysty=3, color=3, $
          yr=[min(meandel)<min(meddel),max(meandel)>max(meddel)]
        oplot, apertures, meddel, ps=8, color=5

        window, 2, xs=450, ys=450
        plot, psfmag, ps=8, col=1, yr=[12,18]
        for k = 0L, naps-1L do oplot, apmag[k,*], ps=8, color=k
;       oplot, apmag[0,*], ps=8, col=7       

        window, 3, xs=450, ys=450
        plot, apertures, medap, ps=8, color=5, xsty=3, ysty=3


;        plot, apertures, delta[*,0], ps=8, xsty=3, ysty=3, color=3, $
;          yrange=[-5,5]
;        for j = 1L, nstars-1L do oplot, apertures, delta[*,j], ps=8, color=3
;        oplot, [!x.crange[0],!x.crange[1]], [0,0], line=2, color=2 

        stop       

;        corr = newdata[*,11]-magalf[match] 	; aperture correction
;        openw, lun4, masterdir+'/'+objname+'.apc', /get_lun
;        printf, lun4, '# Aperture correction for '+objname
;        printf, lun4, ' '
;        printf, lun4, ' Mean    '+string(mean(corr),format='(F7.4)')
;        printf, lun4, ' Median  '+string(median(corr),format='(F7.4)') 
;        printf, lun4, ' StDev   '+string(stdev(corr),format='(F7.4)')
;        free_lun, lun4            
;
;        print, 'Aperture Correction:'
;        print, '--------------------'
;        print, 'Mean   : ', mean(corr)
;        print, 'Median : ', median(corr)
;        print, 'StDev  : ', stdev(corr)
;        print
;
;        window, 0, xs=450, ys=450
;        plothist, corr, bin=stdev(corr)<0.2, xsty=3, ysty=3, xtit='Aperture Correction', ytit='Number', $
;          thick=2, charsize=1.3, charthick=1.3, xthick=2, ythick=2
;        xyouts, [0.75,0.75], [0.85,0.85], '# Stars: '+strn(n_elements(corr)), $
;          /normal, align=0.5, charsiz=1.8, charthick=1.8
;
;; create a postscript file of the histogram of aperture corrections
;        ps_open, 'apcorhist', /ps_fonts, /portrait
;        device, xs=5, ys=5, /inches, /times
;        plothist, corr, bin=stdev(corr)<0.2, xsty=3, ysty=3, xtit='Aperture Correction', ytit='Number', $
;          thick=2, charsize=1.3, charthick=1.3, xthick=2, ythick=2
;        xyouts, [0.75,0.75], [0.85,0.85], '# Stars: '+strn(n_elements(corr)), $
;          /normal, align=0.5, charsiz=1.8, charthick=1.8
;        ps_close

stop        
        
; destroy the temporary alias and environment variable

        spawn, ['unalias dir']
        spawn, ['unsetenv dir']

return
end