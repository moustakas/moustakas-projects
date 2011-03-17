; jm00june19ucb

; program to create symbolic links to the data

pro create_links, objname, datafiles, maskfiles, datanames, masknames, idata=idata, vdata=vdata

	ndata = n_elements(datafiles)

        datanames = strarr(ndata)
        masknames = strarr(ndata)

        if keyword_set(idata) then suffix = 'I'
        if keyword_set(vdata) then suffix = 'V'
        
        print & print, 'Creating the symbolic links for the '+suffix+'-band data . . . '

	if ndata eq 1L then begin	; only one image

            spawn, ['ln -s '+datafiles[j]+' '+objname+'_'+suffix+'.fits'] 
            spawn, ['ln -s '+maskfiles[j]+' '+objname+'_'+suffix+'_dq.fits']

            datanames[0] = objname+'_'+suffix+'.fits'
            masknames[0] = objname+'_'+suffix+'_dq.fits'

        endif else begin		; multiple images

            for j = 0L, ndata-1L do begin

                spawn, ['ln -s '+datafiles[j]+' '+objname+'_'+suffix+strn(1L+j)+'.fits']
                spawn, ['ln -s '+maskfiles[j]+' '+objname+'_'+suffix+strn(1L+j)+'_dq.fits']

                datanames[j] = objname+'_'+suffix+strn(1L+j)+'.fits'
                masknames[j] = objname+'_'+suffix+strn(1L+j)+'_dq.fits'

                rootname = strmid(datanames[j],0,strpos(datanames[j],'.fits'))

            endfor

        endelse

return
end

; jm00june16ucb

; structure for the data, the bad pixel mask and the data header

pro wfpc2_structure, datanames, masknames, wfpc2_struc, idata=idata, vdata=vdata

	if keyword_set(idata) then str_name = 'WFPC2 I-band Data Structure'
	if keyword_set(vdata) then str_name = 'WFPC2 V-band Data Structure'

	ndata = n_elements(datanames)

        template = {name: str_name, $
                    image: fltarr(800,800,ndata), $
                    mask: bytarr(800,800,ndata), $
                    header: strarr(500,ndata)}
        wfpc2_struc = replicate(template,4) 	; array of structures for all four chips

; fill the structure

        for i = 0L, ndata-1L do begin ; loop on the images

            imdata = readfits(datanames[i],dhead) ; data
            immask = readfits(masknames[i]) 	  ; data-quality image (initial mask)

            for indx = 0L, 3L do begin ; loop on the chip number

                wfpc2_struc[indx].image[*,*,i] = imdata[*,*,indx]
                wfpc2_struc[indx].mask[*,*,i] = byte(immask[*,*,indx])
                wfpc2_struc[indx].header[0:n_elements(dhead)-1L,i] = dhead
                    
            endfor

        endfor

return
end

pro wfpc2_redux, struc, vigfits, pixfits

	nsize = size(struc.image)
        ndata = nsize[3]

        for k = 0L, ndata-1L do begin

            struc.image[*,*,k] = struc.image[*,*,k] * pixfits
            struc.mask[*,*,k] = struc.mask[*,*,k] + vigfits
            struc.mask[*,*,k] = struc.mask[*,*,k] eq 0B

            (struc.image[*,*,k])[where(struc.mask[*,*,k] eq 0B)] = -10.E10 ; mark the bad data            

        endfor
            
return
end

pro wfpc2_out, datanames, struc, idata=idata, vdata=vdata

        spawn, ['pwd'], dir
        dir = dir[0]

        ndata = n_elements(datanames)

        if keyword_set(idata) then psf = 'newi'
        if keyword_set(vdata) then psf = 'newv'
        
        for j = 0L, ndata-1L do begin

            rootname = strmid(datanames[j],0,strpos(datanames[j],'.fits'))
            
            for indx = 0L, 3L do begin

                path = dir+'/CHIP'+strn(1L+indx)+'/'
            
                writefits, path+rootname+'_'+strn(1L+indx)+'.fits', struc[indx].image[*,*,j], struc[indx].header[*,j]
            
                spawn, ['ln -s /deepscr1/ioannis/trgb/PSFS/'+psf+strn(1L+indx)+'.psf '$
                        +path+rootname+'_'+strn(1L+indx)+'.psf']

            endfor

        endfor

return
end

pro trgb_hst_redux, objname, datapath
;+
; NAME:
;	TRGB_HST_REDUX
;
; PURPOSE:
;	This program is the master HST/WFPC2 reduction program.  It:
;	(1) creates symbolic links to the data files; (2) applies the
;	data-quality mask; (3) flags cosmic rays; (4) creates a mean,
;	cosmic-ray rejected image from which a starlist can be made;
;	(5) creates a subdirectory structure; and (6) writes out the
;	fits data files by chip number.   
;
; INPUTS:
;	objname  : galaxy name
;	datapath : full path name to the raw data (the program creates
;		   symbolic links to the raw data)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; COMMON BLOCKS:
;
; RESTRICTIONS:
;	Only understands V- and I-band WFPC2 data.
;
; PROCEDURES USED:
;	HEADFITS(), SXPAR(), TRGB_HST_CRREJ, WRITEFITS
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 7, UCB
;-

	npar = n_params()
        if npar ne 2 then begin
            print, 'Syntax: trgb_hst_redux, objname, datapath'
            return
        endif

        spawn, ['pwd'], dir
        dir = dir[0]

; create the directory structure

        print & print, 'Creating the directory structure . . . '

        spawn, ['mkdir -m 0755 -p CHIP1/STARLIST']
        spawn, ['mkdir -m 0755 -p CHIP2/STARLIST']
        spawn, ['mkdir -m 0755 -p CHIP3/STARLIST']
        spawn, ['mkdir -m 0755 -p CHIP4/STARLIST']

; create master.dir in the appropriate directories

        for indx = 0L, 3L do begin
            openw, lun1, dir+'/CHIP'+strn(1L+indx)+'/master.dir', /get_lun
            printf, lun1, dir+'/CHIP'+strn(1L+indx)
            openw, lun2, dir+'/CHIP'+strn(1L+indx)+'/STARLIST/master.dir', /get_lun
            printf, lun2, dir+'/CHIP'+strn(1L+indx)
            free_lun, lun1, lun2
        endfor

	data = findfile(datapath+'/*_c0f.fits')
        mask = findfile(datapath+'/*_c1f.fits')	; data-quality images
        nr = n_elements(data)
        
        if n_elements(mask) ne nr then begin
            print & print, 'The data files do not match the data-quality files!' & return
        endif

; read in the vignetting mask and the pixel-area mask

        vigfits = readfits('/deepscr1/ioannis/trgb/vignetting.fits')
        pixfits = readfits('/deepscr1/ioannis/trgb/pixarea.fits')

; distinguish the data by filter name

        filter = strarr(nr)
        for k = 0L, nr-1L do begin
            header = headfits(data[k])
            filter[k] = sxpar(header,'FILTNAM1')
        endfor

; --------------------------------------------------------------------------------
; begin reductions:  I-band
; --------------------------------------------------------------------------------

        for i = 0L, n_elements(filter)-1L do filter[i] = strn(filter[i])

        i_indx = where(filter eq 'F814W', ni)
        if ni eq 0L then begin

            print & print, 'There are no I-band images!'
            return

        endif else begin

            print & print, 'Beginning I-band reductions . . . '

; create symbolic links to the data

            create_links, objname, data[i_indx], mask[i_indx], idatanames, masknames, /idata ; I-band

            wfpc2_structure, idatanames, masknames, iband, /idata ; fill the data structure

; apply the wfpc2 reductions

            print & print, 'Applying the WFPC2 reductions to the I-band data . . . '
            wfpc2_redux, iband, vigfits, pixfits

; create a cosmic ray mask, flag the CRs in the data image and write
; out a mean, CR-rejected image for each chip 

            for indx = 0L, 3L do begin

                trgb_hst_crrej, iband[indx], cr_mask, meanim 

                iband[indx].image[where(cr_mask eq 0B)] = -10.E10
                    
                print & print, 'Creating the mean I-band image for CHIP'+strn(1L+indx)+' . . .'
                mean_name = dir+'/CHIP'+strn(1L+indx)+'/'+objname+'_Im_'+strn(1L+indx)+'.fits'
                writefits, mean_name, meanim, iband[indx].header[*,0]
                
                spawn, ['ln -s /deepscr1/ioannis/trgb/PSFS/newi'+strn(1L+indx)+'.psf '$
                        +dir+'/CHIP'+strn(1L+indx)+'/'+objname+'_Im_'+strn(1L+indx)+'.psf']

            endfor

; write each image to separate directories and create symbolic links
; to the psf's 
                    
            wfpc2_out, idatanames, iband, /idata

        endelse

; --------------------------------------------------------------------------------
; V-band reductions
; --------------------------------------------------------------------------------

        v_indx = where(filter eq 'F555W', nv)
        if nv eq 0L then begin
            
            print & print, 'There are no V-band images.'
            
        endif else begin

            print & print, 'Beginning V-band reductions . . . '

            create_links, objname, data[v_indx], mask[v_indx], vdatanames, masknames, /vdata

            wfpc2_structure, vdatanames, masknames, vband, /vdata

            print & print, 'Applying the WFPC2 reductions to the V-band data . . . '
            wfpc2_redux, vband, vigfits, pixfits

            for indx = 0L, 3L do begin

                trgb_hst_crrej, vband[indx], cr_mask, meanim 

                vband[indx].image[where(cr_mask eq 0B)] = -10.E10
                    
                print & print, 'Creating the mean V-band image for CHIP'+strn(1L+indx)+' . . .'
                mean_name = dir+'/CHIP'+strn(1L+indx)+'/'+objname+'_Vm_'+strn(1L+indx)+'.fits'
                writefits, mean_name, meanim, vband[indx].header[*,0]
                
                spawn, ['ln -s /deepscr1/ioannis/trgb/PSFS/newv'+strn(1L+indx)+'.psf '$
                        +dir+'/CHIP'+strn(1L+indx)+'/'+objname+'_Vm_'+strn(1L+indx)+'.psf']

            endfor

            wfpc2_out, vdatanames, vband, /vdata

        endelse

; create image.names

        if nv gt 0L then names = [idatanames,vdatanames] else names = idatanames
        ncount = n_elements(names)
        for indx = 0L, 3L do begin
            openw, lun10, dir+'/CHIP'+strn(1L+indx)+'/image.names', /get_lun
            for h = 0L, ncount-1L do printf, lun1, $
              strmid(names[h],0,strpos(names[h],'.fits'))+'_'+strn(1L+indx)
            free_lun, lun10
        endfor

        print & print, 'All finished!' & print

return
end


