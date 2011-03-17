pro trgb_hst_redux_wrapper, galaxy_list, galaxy_dirs
;+
; NAME:
;	TRGB_HST_REDUX_WRAPPER
;
; PURPOSE:
;	Run TRGB_HST_REDUX on a galaxy list.
;
; INPUTS:
;	galaxy_list : name of a text file containing the galaxy names to
;		      be analyzed, one per line
;	galaxy_dirs : name of a text file containing the full path to
;		      the raw data of each corresponding galaxy in
;		      galaxy_list 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	See TRGB_HST_REDUX
;
; COMMON BLOCKS:
;
; PROCEDURE:
;	Loops on each galaxy and calls TRGB_HST_REDUX.
;
; PROCEDURES USED:
;	RDTXT(), TRGB_DATAPATH(), TRGB_HST_REDUX
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 18, UCB
;
;-

	npar = n_params()
        if npar ne 2 then begin
            print
            print, 'Please specify the full path to the list'
            print, 'of galaxy names and galaxy directories.'
            return
        endif

; read in the galaxy list

	glist = rdtxt(galaxy_list)
        gdirs = rdtxt(galaxy_dirs)
        ngals = n_elements(glist)

        paths = trgb_datapath()

        for k = 0L, ngals-1L do begin

            objname = glist[k]
            gdir = gdirs[k]

            cd, paths[2]+objname

            print & print, 'Reducing '+objname+'.'

            trgb_hst_redux, objname, gdir

        endfor

        print & print, 'Done.'

return
end
