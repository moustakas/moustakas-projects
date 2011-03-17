;+
; NAME: 
;   image_align
; PURPOSE: 
;   To register and combine images through a median filter
; DESCRIPTION:
; CATEGORY:
;   CCD data processing
; CALLING SEQUENCE:
;   IMAGE_ALIGN,datacube,datafile,mask,out,outmask,[/FITS,seq_num0,header]
; INPUTS:
;   datacube  - name of the array stack of images
;   datafile  - name of the file produced by DIMAGE that contains 
;               stellar positions
;   mask      - bad pixel mask
; OPTIONAL INPUT PARAMETERS:
; KEYWORD INPUT PARAMETERS:
;   FITS      - If this is set each file will be written out to a fits file
;        ***If Keyword FITS is set these must also be set.*** 
;   seq_num0  - Sequence # of 1st image.
;   header    - name of header array
; OUTPUTS:
;   out     - name of the output idl array
;   outmask - name of output masks
; KEYWORD OUTPUT PARAMETERS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; EXAMPLE(s):
; MODIFICATION HISTORY:
;   06/18/98 - Bryan Mendez, UC Berkeley, Department of Astronomy
;-
PRO IMAGE_ALIGN,datacube,datafile,mask,out,outmask,$
		FITS=fits,seqnum0,header
         
;------------Read coordinate file---------------------------------------

temp=SIZE(datacube)
num_frames=temp(3)

junk=''
OPENR,lun,datafile,/GET_LUN
READF,lun,junk,FORMAT='(A80)'
n=0
WHILE (NOT EOF(lun)) DO BEGIN
	READF,lun,x,y,FORMAT='(6x,f0,6x,f0)'  
	n=n+1
ENDWHILE 
CLOSE,lun

num_stars=n/num_frames

xi=FLTARR(num_stars*num_frames)
yi=FLTARR(num_stars*num_frames)
OPENR,lun,datafile,/GET_LUN
READF,lun,junk,FORMAT='(A80)'
i=0
WHILE (NOT EOF(lun)) DO BEGIN

	READF,lun,x,y,FORMAT='(6x,f0,6x,f0)'  
	xi(i)=x
	yi(i)=y

	i=i+1

ENDWHILE 
CLOSE,lun

;----------compute shifts-------------------------------------------
xn=FLTARR(num_stars,num_frames)
yn=FLTARR(num_stars,num_frames)

FOR j=0,num_frames-1 DO BEGIN

	FOR k=0,num_stars-1 DO BEGIN

	xn(k,j)=xi(k+j*num_stars)
	yn(k,j)=yi(k+j*num_stars)

	ENDFOR
ENDFOR

shiftx=FLTARR(num_stars,num_frames-1)
shifty=FLTARR(num_stars,num_frames-1)
FOR l=0,num_stars-1 DO BEGIN

	FOR m=1,num_frames-1 DO BEGIN

	shiftx(l,m-1)=xn(l,0)-xn(l,m)
	shifty(l,m-1)=yn(l,0)-yn(l,m)

	ENDFOR
ENDFOR

sx=FLTARR(num_frames-1)
sy=FLTARR(num_frames-1)
FOR n=0,num_frames-2 DO BEGIN

	sx(n)=ROUND(MEAN(shiftx(*,n)))
	PRINT,'Rounded Average Shift in X for image',n+2,'=',sx(n)
	sy(n)=ROUND(MEAN(shifty(*,n)))
	PRINT,'Rounded Average Shift in Y for image',n+2,'=',sy(n)

ENDFOR

;----check size of shifts and pad by appropriate amount-------------------

mmsx=MINMAX(sx)
maxsx=MAX(ABS(mmsx))
mmsy=MINMAX(sy)
maxsy=MAX(ABS(mmsy))
padx=2*maxsx
pady=2*maxsy

bigarray=FLTARR(1600+padx,2048+pady,num_frames)
bigmask=BYTARR(1600+padx,2048+pady,num_frames)
left=FIX(maxsx)
right=FIX(maxsx+1600-1)
bottom=FIX(maxsy)
top=FIX(maxsy+2048-1)
bigarray[left:right,bottom:top,0]=datacube[*,*,0]
bigmask[left:right,bottom:top,0]=mask

	FOR m=0,num_frames-2 DO BEGIN

		bigarray[left+sx(m):right+sx(m),bottom+sy(m):top+sy(m),m+1]= $
                                                 datacube[*,*,m+1]
		bigmask[left+sx(m):right+sx(m),bottom+sy(m):top+sy(m),m+1]= $
						mask

	ENDFOR

IF (MAX(sx) GT 0) THEN lborder=left+MAX(sx) ELSE lborder=left

IF (MAX(sy) GT 0) THEN bborder=bottom+MAX(sy) ELSE bborder=bottom

IF (MIN(sx) LT 0) THEN rborder=right+MIN(sx) ELSE rborder=right

IF (MIN(sy) LT 0) THEN tborder=top+MIN(sy) ELSE tborder=top


xsize=rborder-lborder+1
ysize=tborder-bborder+1

out=FLTARR(xsize,ysize,num_frames)
outmask=BYTARR(xsize,ysize,num_frames)	
	FOR n=0,num_frames-1 DO BEGIN

	out[*,*,n]=bigarray[lborder:rborder,bborder:tborder,n]
	outmask[*,*,n]=bigmask[lborder:rborder,bborder:tborder,n]

	ENDFOR


;----Optional write out of array to separate fits files-----------------

IF KEYWORD_SET(fits) THEN BEGIN
	data_size=SIZE(out)
	num_frame=data_size(3)
	seq=STRARR(num_frame)
	seq_num=INDGEN(num_frame)+seqnum0
	seq_num1=STRING(FIX(seq_num))
        v=0
	ni=0

        FOR i=0,num_frame-1 DO BEGIN

;		fi=SXPAR(header[*,i],'REDFILT')
		fi=SXPAR(header[*,i],'FILTER')
		fl=STRPOS(fi,'V')
		IF (fl EQ 0) THEN v=v+1 ELSE ni=ni+1

	ENDFOR

	FOR i=0,num_frame-1 DO BEGIN
	
		image=out[*,*,i]
		h=header[*,i]
;		SXADDPAR,h,'GAIN',1.97,/PDU
	;	SXADDPAR,h,'RDNOISE',6.3,/PDU
		;SXADDPAR,h,'GDATAMAX',64000,/PDU
		filter=SXPAR(h,'FILTER')
;		SXADDPAR,h,'FILTER',filter,/PDU
;		SXDELPAR,h,'REDFILT'
		date=SXPAR(h,'DATE')
;		SXADDPAR,h,'DATE-OBS',date,/PDU
;		SXDELPAR,h,'DATE'
		epoch=SXPAR(h,'EQUINOX')
;		SXADDPAR,h,'EPOCH',epoch,/PDU
;		SXDELPAR,h,'EQUINOX'
;		SXADDPAR,h,'OBSERVAT','keck'
;		saturated=WHERE(image GT 65535.0,c5)
;		IF (c5 NE 0) THEN image(saturated)=65535.0
;		image=image-32768
;		image=FIX(image)
;		d_max=MAX(image)+32768
;		d_min=MIN(image)+32768
;		SXADDPAR,h,'DATAMIN',d_min,/PDU
;		SXADDPAR,h,'DATAMAX',d_max,/PDU
;		SXADDPAR,h,'BSCALE',1.0,/PDU
;		SXADDPAR,h,'BZERO',32768.0,/PDU
		history='Overscan Corrected, Cropped, Dark subtracted, Flat Fielded, and Registered'
                sxdelpar, h, 'HISTORY'
		SXADDHIST,history,h,/PDU
		root=SXPAR(h,'TARGNAME')
		name=''
		fflag=STRPOS(filter,'I')
		IF (fflag EQ 0) THEN name=root+'_'+filter+seq_num1(i)+'.fits'
		IF (fflag EQ 0) AND (ni EQ 1) THEN name=root+'_'+filter+'.fits'
		IF (fflag EQ -1) AND (v GT 1) THEN name=root+'_'+filter+seq_num1(i)+'.fits' 
		IF (fflag EQ -1) AND (v EQ 1) THEN name=root+'_'+filter+'.fits'
		name=STRCOMPRESS(name,/REMOVE_ALL)
		PRINT,'Writing fits file...',name
		WRITEFITS,name,image,h

		h_size=SIZE(h)
		header[0:h_size(1)-1,i]=h
		

	ENDFOR


ENDIF


END



