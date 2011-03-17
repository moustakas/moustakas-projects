;+
; NAME: 
;  lris_proc
; PURPOSE: 
;  Process Keck LRIS image data, in 2 amplifier mode.
; DESCRIPTION:
; CATEGORY:
;  CCD data processing
; CALLING SEQUENCE:
;  LRIS_PROC,num_frame,seq_num0,darkarr,flatarr,imgodf,header, $
;            [/FITS,RVAL=rval]
; INPUTS:
;   num_frame - number of frames to be reduced
;   seq_num0  - sequence number of the first lris****.fits file (input should 
;                be 4 characters long and files must be in a sequence)
;   darkarr   - name of the idl dark frame to be used
;   flatarr   - name of the idl flat frame to be used
; OPTIONAL INPUT PARAMETERS:
; KEYWORD INPUT PARAMETERS:
;   FITS      - If this keyword is set then the reduced images will be 
;                 written out to a fits file as well as an idl array. 
;   RVAL    - Set this keyword to use predetermined rscale value in OVERSCAN
; OUTPUTS:
;   imgodf    - name of the output image array
;   header    - name of output header array
; KEYWORD OUTPUT PARAMETERS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; MODIFICATION HISTORY:
;   07/24/98 - Bryan Mendez, UC Berkeley, Department of Astronomy
;   
;-
PRO LRIS_PROC,num_frame,seq_num0,darkarr,flatarr,imgodf,header, $
              FITS=fits, RVAL=rval


seq=STRARR(num_frame)
seq_num=INDGEN(num_frame)+seq_num0
seq_num1=STRING(FIX(seq_num))

infile=''

header=STRARR(100,num_frame)
img=FLTARR(2250,2048,num_frame)
imgo=FLTARR(1600,2048,num_frame)
imgodf=FLTARR(1600,2048,num_frame)
FOR i=0,num_frame-1 DO BEGIN

	IF (seq_num1(i) lt 10) THEN BEGIN
	rootname='lris000'
	ENDIF

	IF (seq_num1(i) ge 10 and seq_num1(i) lt 100) THEN BEGIN
	rootname='lris00'
	ENDIF

	IF (seq_num1(i) ge 100 and seq_num1(i) lt 1000) THEN BEGIN
	rootname='lris0'
	ENDIF

	infile=rootname+seq_num1(i)+'.fits'
	infile=STRCOMPRESS(infile,/REMOVE_ALL)
	PRINT,'Reading...',infile
	img[*,*,i]=READFITS(infile,hdr,/SILENT)
	IF KEYWORD_SET(rval) THEN BEGIN
        imgo[*,*,i]=OVERSCAN(img[*,*,i],2094,2166,2174,2244,250,1849,0,2047,RSCALE=rval)
	ENDIF ELSE BEGIN 
	imgo[*,*,i]=OVERSCAN(img[*,*,i],2094,2166,2174,2244,250,1849,0,2047)
	ENDELSE 

	exptime=SXPAR(hdr,'EXPOSURE')
	exptime=FLOAT(exptime)
	imgodf[*,*,i]=(imgo[*,*,i]-(darkarr*exptime))/flatarr 
	hdr_size=SIZE(hdr)
	header(0:hdr_size(1)-1,i)=hdr

ENDFOR



IF KEYWORD_SET(fits) THEN BEGIN

	FOR i=0,num_frame-1 DO BEGIN
	
	image=imgodf[*,*,i]
	h=header[*,i]
	SXADDPAR,h,'GAIN',1.97,/PDU
	SXADDPAR,h,'RDNOISE',6.3,/PDU
	SXADDPAR,h,'GDATAMAX',64000,/PDU
	filter=SXPAR(h,'REDFILT')
	SXADDPAR,h,'FILTER',filter,/PDU
	SXDELPAR,h,'REDFILT'
	date=SXPAR(h,'DATE')
	SXADDPAR,h,'DATE-OBS',date,/PDU
	SXDELPAR,h,'DATE'
	epoch=SXPAR(h,'EQUINOX')
	SXADDPAR,h,'EPOCH',epoch,/PDU
	SXDELPAR,h,'EQUINOX'
	SXADDPAR,h,'OBSERVAT','keck'
	saturated=WHERE(image GT 65535.0,c1)
	IF (c1 NE 0) THEN image(saturated)=65535.0
	image=image-32768
	image=FIX(image)
	d_max=MAX(image)+32768
	d_min=MIN(image)+32768
	SXADDPAR,h,'DATAMIN',d_min
	SXADDPAR,h,'DATAMAX',d_max
	SXADDPAR,h,'BSCALE',1.00000000,/PDU
	SXADDPAR,h,'BZERO',32768.00000000,/PDU
	history='Overscan Corrected, Cropped, Dark subtracted, Flat Fielded'
	SXADDHIST,history,h,/PDU
	root=SXPAR(h,'TARGNAME')
	name=''
	name=root+'_'+filter+seq_num1(i)+'.fits'
	name=STRCOMPRESS(name,/REMOVE_ALL)
	PRINT,'Writing fits file...',name
	WRITEFITS,name,image,h

	h_size=SIZE(h)
	header[0:h_size(1)-1,i]=h


	ENDFOR

ENDIF



END

