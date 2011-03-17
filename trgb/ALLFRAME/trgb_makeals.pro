pro trgb_makeals, hst=hst
;+
; NAME:
;	TRGB_MAKEALS
;
; PURPOSE:
;	Create the two-line header .als files needed by ALLFRAME. 
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;	hst : specify HST data
;
; OUTPUTS:
;
; PROCEDURES USED:
;	RDTXT(), HEADFITS(), DAOHEADER
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 May 1, UCB
;	modified, jm00may25ucb
;-

	spawn, ['pwd'], datapath	; current directory
        datapath = datapath[0]

        imnames = rdtxt(datapath+'/image.names') ; read in the image names

        for i = 0L, n_elements(imnames)-1L do begin  ; loop on the images

            image = imnames[i]

            if keyword_set(hst) then begin

                daohead = strarr(2)
                frmt = '(2x,I1,1x,I4,1x,I4,1x,F7.1,1x,F7.1,2x,F6.1,3x,f5.2,2x,f5.2,3x,f5.2,3x,f5.2)'

                daohead[0] = rdtxt('/deep1/ioannis/trgb/daohead.sample_als')
                daohead[1] = string([1],[800],[800],[-109.9],[8000.0],$
                                    [6.16],[3.00],[7.00],[0.75],[1.75],format=frmt)

            endif else begin

                fhead = headfits(datapath+'/'+image+'.fits') ; fits header
                spawn, ['find '+datapath+' -name "'+image+'.pars" -print'], parfile
                parfile = parfile[0]
                
                daoheader, parfile, fhead, daohead, /als ; create the als header

            endelse
            
            daofile = datapath+'/'+image+'.als'

            print               ; write the als file
            print, 'Writing '+daofile+'.'
            openw, lun1, daofile, /get_lun
            printf, lun1, daohead[0] ; insert the daophot header
            printf, lun1, daohead[1] ; insert the header values
            printf, lun1, ' '        ; insert a blank line
            print
            free_lun, lun1

        endfor

return
end

