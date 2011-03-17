pro multiphot_script_wrapper, galaxy_list, galaxy_path
;+
; NAME:
;	MULTIPHOT_SCRIPT_WRAPPER
;
; PURPOSE:
;	Run MULTIPHOT_SCRIPT on a galaxy list.
;
; CALLING SEQUENCE:
;	multiphot_script_wrapper, galaxy_list.
;
; INPUTS:
;	galaxy_list : name of a text file containing the galaxy names to
;                     be analyzed, one per line
;	galaxy_path : name of a text file containing the full path to
;		      the raw data of the galaxy names given in
;		      galaxy_list 
;
; OUTPUTS:
;	See MULTIPHOT_SCRIPT
;
; PROCEDURE:
;	Loops on each galaxy and calls WFPC2_LINKS and
;	MULTIPHOT_SCRIPT. 
;
; EXAMPLE:
;	multiphot_script_wrapper,
;	'/deepscr1/ioannis/trgb/galaxy_list', '/deepscr1/ioannis/trgb/galaxy_dir'
;
; PROCEDURES USED:
;	TRGB_DATAPATH, RDTXT(), MULTIPHOT_SCRIPT, WFPC2_LINKS
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 August 3, UCB
;
;-

	if ((not keyword_set(galaxy_list)) and (not keyword_set(galaxy_dir))) then begin
            print, 'Syntax: multiphot_script_wrapper, galaxy_list, galaxy_dir'
            return
        endif

        paths = trgb_datapath()

        glist = rdtxt(galaxy_list)
        gdir = rdtxt(galaxy_path)

        ngals = n_elements(glist)
        ndirs = n_elements(gdir)
        if ngals ne ndirs then begin
            print, 'The text files need to match up!'
            return
        endif
        
        for k = 0L, ngals-1L do begin ; make all the links first
            
            objname = glist[k]
            dir = gdir[k]

            pushd, paths[0]+objname
            wfpc2_links, objname, dir
            popd

        endfor

        for k = 0L, ngals-1L do begin ; run multiphot

            objname = glist[k]

            print & print, 'Running MULTIPHOT on '+objname+' . . .'
            multiphot_script, objname, /log, /nocheck

        endfor

        print & print, 'All done!'

return
end
