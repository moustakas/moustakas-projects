pro trgb_calibrate, objname
;+
; NAME:
;	TRGB_CALIBRATE
;
; PURPOSE:
;	Calibrate the raw magnitudes output from ALLFRAME onto a
;	standard system for the KECK data.
;
; INPUTS:
;	The user must specify objname, the name of the object being
;	analyzed.  A file called image.names must exist in the
;	directory from which this program is run, and must contain a
;	list, one per line, of all the frames for the current object.
;	The user must be in the directory where all the FITS images
;	reside.  
;
; OUTPUTS:
;	Creates a file containing the standardized I and V magnitudes
;	and errors and V-I colors for all the stars in the starlist.
;
; PROCEDURE:
;	Applies aperture corrections (TRGB_APCORRECT), exposure-time 
;	corrections, and color-corrections.
;
; PROCEDURES USED:
;	TRGB_DATAPATH(), RDTXT(), READFAST(), ROBUST(), READFITS(),
;	SKY
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 April 14, UCB, based on Bryan Mendez's
;	CALIBRATE. 
;-

	paths = trgb_datapath()

	objdata = ['holmbergii','holmbergix','i342','ngc1560','ngc2366',$
                   'ngc2903','ngc2976','ngc3109','sextansb']
        check = where(strupcase(objdata) eq strupcase(objname),count)
        if count gt 0L then date = '22dec97/' else date = '23dec97/'
        
        pushd, paths[2]+date+strlowcase(objname)
        datapath = paths[2]+date+strlowcase(objname)
        
        imnames = rdtxt(datapath+'/image.names')        
        iband = imnames[where(strpos(imnames,'_I') gt 0)]  ; I-band frames
        vband = imnames[where(strpos(imnames,'_V') gt 0)]  ; V-band frames

        ni = n_elements(iband)
        nv = n_elements(vband)

; read in the aperture corrections

       iapcor = fltarr(ni)-fltarr(ni)
;      iapcor = fltarr(ni)
;      for j = 0L, ni-1L do begin
;          readcol, iband[j]+'.apc', cor, skipline=2, /silent, format='X,F7.4'
;          iapcor[j] = cor[1]
;      endfor
       vapcor = fltarr(nv)-fltarr(nv)
;      vapcor = fltarr(nv)
;      for j = 0L, nv-1L do begin
;          readcol, vband[j]+'.apc', cor, skipline=2, /silent, format='X,F7.4'
;          vapcor[j] = cor[1]
;      endfor

; transformation equation coefficients (from the standards):

; v_instr = V_abs + v1 + v2*XV + v3*(V-I) + v4*(V-I)*XV
; i_instr = I_abs + i1 + i2*XI + i3*(V-I) + i4*(V-I)*XI

        v1 = -2.598
        v2 = 0.12
        v3 = 0.0
        v4 = 0.0

        i1 = -2.438
        i2 = 0.07
        i3 = 0.0
        i4 = 0.0
        
        iter =  5  ; number of iterations for robust averaging

; read the ALLSTAR/ALLFRAME photometry

        print & print, 'Reading in the data . . . ' & print

	alfile = objname+'.raw'	; allframe/daomaster photometry file

        readfast, alfile, rawdata, ncol=2*ni+(2*nv+4)+1, skip=3
        rawdata = transpose(rawdata)
stop        
        nstars = n_elements(rawdata[*,0])	; total number of stars

        icols = 2*indgen(ni)+3L	 	; I-band magnitude columns
        iraw = rawdata[*,icols]		; raw I-band data
        ierrraw = rawdata[*,icols+1L]	; raw I-band errors
        
        vcols = 2*ni+2*INDGEN(nv)+3L	; V-band columns
        vraw =  rawdata[*,vcols]	; raw V-band data
        verrraw = rawdata[*,vcols+1L]	; raw V-band magnitude errors
        
; add exposure time normalization and aperture correction, then apply
; the inverse transformation to the standard system
        
        print & print, 'Applying the corrections . . . ' & print

        for j = 0L, ni-1L do begin	; I-band
            ihdr = headfits(iband[j]+'.fits')
            itime = sxpar(ihdr,'EXPOSURE')	; exposure time
            airmass = sxpar(ihdr,'AIRMASS')	; airmass

            iraw[*,j] = iraw[*,j] + 2.5*alog10(itime) + iapcor[j] - i1 - i2*airmass
        endfor

        for k = 0L, nv-1L do begin	; V-band
            vhdr = headfits(vband[k]+'.fits')
            vtime = sxpar(vhdr,'EXPOSURE')
            airmass = sxpar(vhdr,'AIRMASS')

            vraw[*,k] = vraw[*,k] + 2.5*alog10(vtime) + vapcor[k] - v1 - v2*airmass
        endfor

; robustly average multiple frames in the same filter

        print & print, 'Averaging the data . . . ' & print

        if (ni ne 1L) then begin
            imags = fltarr(nstars) & ierr = fltarr(nstars)
            for i = 0L, nstars-1L do imags[i] = robust(iraw[i,*],ierrraw[i,*],iter)
            for i = 0L, nstars-1L do ierr[i] = median(ierrraw[i,*])
        endif else begin
            imags = iraw
            ierr = ierrraw
        endelse

        if (nv ne 1L) then begin
            vmags = fltarr(nstars) & verr = fltarr(nstars)
            for j = 0L, nstars-1L do vmags[j] = robust(vraw[j,*],verrraw[v,*],iter)
            for j = 0L, nstars-1L do verr[j] = median(verrraw[j,*])
        endif else begin
            vmags = vraw
            verr = verrraw
        endelse

        data = [[rawdata[*,0:2]],[imags],[ierr],[vmags],[verr],[vmags-imags]]

; read in the ALLFRAME-subtracted image of the frame used to create
; the master starlist (assumed the first I-band frame)

        image = readfits(iband[0]+'j.fits',header)

        imsize = size(image)
        xsize = imsize[1]
        ysize = imsize[2]

        flag = bytarr(nstars)+1B

; flag and filter bad stars

        for l = 0L, nstars-1L do begin
            if (data[l,3] gt 30.) then flag[l] = 0B	; I mag > 30
            if (data[l,5] gt 30.) then flag[l] = 0B	; V mag > 30
            if (data[l,4] gt 1.) then flag[l] = 0B	; I mag error > 1
            if (data[l,6] gt 1.) then flag[l] = 0B	; V mag error > 1
            if (data[l,1] lt 1. or data[l,2] lt 1.) then flag[l] = 0B	; too close to the chip edge
            if (data[l,1] gt xsize-1. or data[l,2] gt ysize-1.) then flag[l] = 0B
        endfor

        goodata = data[where(flag eq 1B),*]
        ngood = n_elements(goodata[*,0])

; create a mask from the subtracted image

        print & print, 'Creating an image mask . . . '

        med = median(image,18)	; median filter
        diff = image-med
        
        sky, diff, dum, skysig, /silent
        masq = abs(diff) gt abs(4.*skysig) ; 4-sigma clip

        mask = median(masq,8)	; 1 is bad & 0 is good
        mask = mask eq 0B	; 0 is bad & 1 is good

; flag as bad any stars whose coordinates are near the masked points

        flagmask = bytarr(ngood)+1B
        for j = 0L, ngood-1L do $
          if (mask[goodata[j,1]-1L,goodata[j,2]-1L] eq 0B) then flagmask[j]=0B

        dataout = goodata[where(flagmask eq 1B),*]
        nout = n_elements(dataout[*,0])

; write the final I- and V-band magnitudes, errors, and V-I colors to a file:

        print & print, 'And writing '+objname+'_IVmags.dat . . .' & print

        openw, lun1, datapath+'/'+objname+'_IVmags.dat', /get_lun
        printf, lun1, '#  ID   Xcenter   Ycenter    I       Ierr      V       Verr     V-I'
        printf, lun1, ' '
        for k = 0L, nout-1L do $
          printf, lun1, dataout[k,*], format = '(1x,I6,7F9.3)'
        free_lun, lun1

        print & print, '. . . done.'

        popd
        
return
end
