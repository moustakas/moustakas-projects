pro trgb_hst_master, galaxy_list
;+
; NAME:
;	TRGB_HST_MASTER
;
; PURPOSE:
;	Master HST/WFPC2 ALLFRAME program.
;
; INPUTS:
;	galaxy_list : name of a text file containing the galaxy names to
;		      be analyzed, one per line
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; PROCEDURES USED:
;	RDTXT(), TRGB_DATAPATH(), TRGB_HST_PRE_STARLIST, TRGB_MAKEMCH,
;	TRGB_HST_DAOMASTER, TRGB_STARLIST, TRGB_ALLFRAME,
;	TRGB_HST_CALIBRATE, TRGB_HST_COMBINE_STARLISTS 
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 19, UCB
;-


; read in the galaxy list

	glist = rdtxt(galaxy_list)
        ngals = n_elements(glist)

        paths = trgb_datapath()
        datapath = paths[0]

        chip_names = ['pc','wfc2','wfc3','wfc4']

        for k = 0L, ngals-1L do begin	; loop on the galaxy list

            objname = glist[k]

            print & print, 'Object: '+objname+'.' & print

            pushd, datapath+objname

            for j = 1L, 4L do begin     ; loop on the wfpc2 chips

                cd, datapath+objname+'/CHIP'+strn(j)

                if j eq 1L then pcamera = 1 else pcamera = 0

                trgb_hst_pre_starlist, objname, /mean, pc=pcamera;, /verbose
                trgb_hst_pre_starlist, objname, pc=pcamera;, /verbose
                    
                trgb_makemch, objname, /ap
                trgb_hst_daomaster, objname
                    
                cd, datapath+objname+'/CHIP'+strn(j)+'/STARLIST'

                meanim = objname+'_Im_'+strn(j)
                trgb_starlist, meanim, /hst, pc=pcamera, /verbose

                trgb_allframe, objname, /verbose, /hst, /nocheck

                chip = chip_names[j-1L]
                trgb_hst_calibrate, objname, chip

            endfor

            trgb_hst_combine_starlists, objname

            popd
	
        endfor

return
end
