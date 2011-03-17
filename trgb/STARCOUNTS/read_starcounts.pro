pro read_starcounts;, keck=keck, hst=hst
;+
; NAME:
;	READ_STARCOUNTS
;
; PURPOSE:
;	Routine to read in the differential star counts files for a
;	given (l,b) and to organize them.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; COMMON BLOCKS:
;
; RESTRICTIONS:
;	Assumes a datapath for now.
;
; PROCEDURE:
; 	The files are organized as:
;
;		MAG.BIN.  logNTOTAL  DISK    BULGE     ARMS    RING	HALO
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 28, UCB
;	jm00nov12uofa, added hst and keck keywords
;-

;        on_error, 2
;        if (not keyword_set(hst)) and (not keyword_set(keck)) then $
;          message, 'Please specify a KECK or HST keyword.'

	path = '/deepscr1/ioannis/STARCOUNTS/'

        flist = findfile(path)
        good = where(strpos(flist,'all') ge 0,fcount)
        flist = flist[good]

        ifiles = flist[where(strpos(flist,'id') ge 0,icount)]
        vfiles = flist[where(strpos(flist,'vd') ge 0,vcount)]

;       hstfov = 5.66^2	 ; HST FOV (field of view) in arcminutes
        hstfov = 3*(1.33^2) + 0.6133^2
        keckfov = [42., 39., 36., 36., 42., 36., 35., 34., 42., 33., 42., 36.]

        template = {name:	'', $
                    filename:	strarr(1), $
                    magbin:	ptr_new(1D), $
                    logcounts:	ptr_new(1D), $
                    counts:	ptr_new(1D)}
        starcounts = replicate(template,icount)
        
; read in the I-band data

;       if keyword_set(hst) then begin
        
            for k = 0L, icount-1L do begin

;               readfast, path+ifiles[k], istars, ncols=7
                read_data, path+ifiles[k], istars, ncol=7, /silent
                istars = transpose(istars)

                starcounts[k].filename = ifiles[k]
                starcounts[k].magbin = ptr_new(istars[0,*])
                
                logcounts = istars[1,*] ;/(4. * 60.^2) * hstfov^2
                starcounts[k].logcounts = ptr_new(logcounts)
                
                temp = (10^logcounts)/(4. * 60.^2) * hstfov
                starcounts[k].counts = ptr_new(temp)
                
            endfor

;       endif

;       if keyword_set(keck) then begin
        
            for k = 12, 23 do begin

                read_data, path+ifiles[k], istars, ncol=7, /silent
                istars = transpose(istars)
                
                starcounts[k].filename = ifiles[k]
                starcounts[k].magbin = ptr_new(istars[0,*])
                
                logcounts = istars[1,*]
                starcounts[k].logcounts = ptr_new(logcounts)  
                
                temp = (10^logcounts)/(4. * 60.^2) * keckfov[k-12]
                starcounts[k].counts = ptr_new(temp)
                
            endfor

;       endif
        
        window, 0, xs=450, ys=450
        colortable1
        plot, *starcounts[0].magbin, *starcounts[0].counts, $
          xsty=3, ysty=3, xr=[18,27], yr=[0,55], color=1, line=1, $
          xtitle='I', ytitle='Counts/quarter-mag/HST FOV', title='Star Counts'
;       for j = 1L, icount-1L do oplot, *starcounts[j].magbin, *starcounts[j].logcounts, ps=8
        for j = 1, 11 do oplot, *starcounts[j].magbin, *starcounts[j].counts, line = 1
        for j = 12, 23 do oplot, *starcounts[j].magbin, *starcounts[j].counts, line = 2
        
	ptr_free, starcounts[*].magbin
	ptr_free, starcounts[*].logcounts
	ptr_free, starcounts[*].counts
        
return
end

