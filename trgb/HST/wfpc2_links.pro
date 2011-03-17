pro create_links, objname, datafiles, maskfiles, datanames, masknames, idata=idata, vdata=vdata

	ndata = n_elements(datafiles)

        datanames = strarr(ndata)
        masknames = strarr(ndata)

        if keyword_set(idata) then suffix = 'I'
        if keyword_set(vdata) then suffix = 'V'
        
        print
        print, 'Creating the symbolic links for the '+suffix+'-band data . . . '
        print

	if ndata eq 1L then begin	; only one image

            datanames[0] = objname+'_'+suffix+'.fits'
            masknames[0] = objname+'_'+suffix+'_dq.fits'

            print, '   '+datanames[0]+'.'
            print, '   '+masknames[0]+'.'
            print

            spawn, ['ln -s '+datafiles[0]+' '+objname+'_'+suffix+'.fits'] 
            spawn, ['ln -s '+maskfiles[0]+' '+objname+'_'+suffix+'_dq.fits']

        endif else begin		; multiple images

            for j = 0L, ndata-1L do begin

                datanames[j] = objname+'_'+suffix+strn(1L+j)+'.fits'
                masknames[j] = objname+'_'+suffix+strn(1L+j)+'_dq.fits'

                print, '   '+datanames[j]+'.'
                print, '   '+masknames[j]+'.'
                print

                spawn, ['ln -s '+datafiles[j]+' '+objname+'_'+suffix+strn(1L+j)+'.fits']
                spawn, ['ln -s '+maskfiles[j]+' '+objname+'_'+suffix+strn(1L+j)+'_dq.fits']

                rootname = strmid(datanames[j],0,strpos(datanames[j],'.fits'))

            endfor

        endelse

return
end

pro wfpc2_links, objname, datapath
;+
; NAME:
;	WFPC2_LINKS
;
; PURPOSE:
;	Create symbolic links from any subdirectory to raw WFPC2
;	data.  Essentially, renames the raw WFPC2 fits files to a file
;	name which contains the object name and the filter, and
;	numbers the data files sequentially.
;
; CALLING SEQUENCE:
;	wfpc2_links, objname, datapath
;
; INPUTS:
;	objname  : string name of the galaxy
;	datapath : full string datapath to the raw data
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	Writes symbolic links to the WFPC2 data in the current working
;	directory.  
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;	This routine cannot overwrite existing symbolic links.
;	Currently, the routine only understands I- and V-band data.  
;
; PROCEDURE:
;	Delete any pre-existing symbolic links.  Determine the name of
;	the object and the full datapath to the raw WFPC2 data, and
;	call this routine.
;
; EXAMPLE:
;	wfpc2_links, 'ngc2903', '/deepscr1/marc/trgb/hst/n2903'
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 30, UCB
;-

	npar = n_params()
        if npar ne 2 then begin
            print, 'Syntax: wfpc2_links, objname, datapath'
            return
        endif

	data = findfile(datapath+'/*_c0f.fits')
        mask = findfile(datapath+'/*_c1f.fits')	; data-quality images
        nr = n_elements(data)
        
        if n_elements(mask) ne nr then begin
            print & print, 'The data files do not match the data-quality files!' & return
        endif

; distinguish the data by filter name

        filter = strarr(nr)
        for k = 0L, nr-1L do begin
            header = headfits(data[k])
            filter[k] = sxpar(header,'FILTNAM1')
        endfor

        for i = 0L, n_elements(filter)-1L do filter[i] = strn(filter[i])

        i_indx = where(filter eq 'F814W', ni)
        if ni eq 0L then begin

            print & print, 'There are no I-band images.'
            return

        endif else begin

            create_links, objname, data[i_indx], mask[i_indx], idatanames, masknames, /idata ; I-band

        endelse
            
        v_indx = where(filter eq 'F555W', nv)
        if nv eq 0L then begin
            
            print & print, 'There are no V-band images.' & print
            
        endif else begin

            create_links, objname, data[v_indx], mask[v_indx], vdatanames, masknames, /vdata

        endelse
            
return
end


