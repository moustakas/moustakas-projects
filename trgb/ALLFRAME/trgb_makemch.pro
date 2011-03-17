pro trgb_makemch, objname, ap=ap, als=als, alf=alf, mini=mini, hst=hst
;+
; NAME:
;	TRGB_MAKEMCH
;
; PURPOSE:
;	Create a DAOPHOT .mch (transformation) file.  Will create a
;	transformation file for .ap, .als and .alf photometry files.
;
; INPUTS:
;	objname	: galaxy name
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; RESTRICTIONS:
;	The HST data needs the CTE-corrected alf photometry files. 
;
; PROCEDURES USED:
;	RDTXT()
;	
; MODIFICATION HISTORY:
;	John Moustakas, 2000 May 25, UCB
;-

	npar = n_params()
        if npar ne 1 then begin
            print, 'Syntax: trgb_makemch, objname, ap=ap, als=als, alf=alf'
            return
        endif

	spawn, ['pwd'], datapath	; current directory
        datapath = datapath[0]

        imnames = rdtxt(datapath+'/image.names') ; read in the image names
        nr = n_elements(imnames)		 ; number of images

 ; for the mini starlist we want to isolate the I-band images

        if keyword_set(mini) then begin
            iband = imnames[where(strpos(imnames,'_I') gt 0)] ; I-band frames
            ni = n_elements(iband)
            indx = ni
            file = objname+'_mini.mch'
        endif else begin
            indx = nr
            file = objname+'.mch'
        endelse
            
        mchstr = strarr(indx)

        if keyword_set(ap) then begin ; no mch file exists: create one
            
            for i = 0L, indx-1L do $
              mchstr[i] = "'"+imnames[i]+'.ap'+"  '"+$
              '   0.000     0.000   1.00000   0.00000   0.00000   1.00000     0.000'
            print
            print, 'Next exit IDL and run DAOMASTER, '
            print, 'giving it this newly created file.'
            print

        endif else begin 
            if keyword_set(als) then begin ; assume a mch file exists

                spawn, ['find '+datapath+' -name "'+objname+'.mch" -print'], mchfile
                mchtext = rdtxt(mchfile) ; read the file
                for i = 0L, indx-1L do mchstr[i] = "'"+imnames[i]+'.als'+strmid(mchtext[i],rstrpos(mchtext[i],"'"))
                
            endif else begin 
                if keyword_set(alf) then begin ; same as above

                    spawn, ['find '+datapath+' -name "'+objname+'.mch" -print'], mchfile
                    mchtext = rdtxt(mchfile) ; read the file
                    if keyword_set(hst) then $
                      for j = 0L, indx-1L do mchstr[j] = "'"+imnames[j]+'.alf_cte'+$
                      strmid(mchtext[j],rstrpos(mchtext[j],"'")) else $
                      for j = 0L, indx-1L do mchstr[j] = "'"+imnames[j]+'.alf'+$
                      strmid(mchtext[j],rstrpos(mchtext[j],"'"))

                endif else begin
                    
                    print
                    print, 'Please specify which transformation file you want '
                    print, 'to create (ap, als, alf) using the keywords.'
                    return
                    
                endelse
            endelse
        endelse
        
        print & print, 'Writing '+file & print
        openw, lun1, file, /get_lun
        for k = 0L, indx-1L do printf, lun1, mchstr[k]
        free_lun, lun1       

return
end
