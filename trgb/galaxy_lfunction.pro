pro galaxy_lfunction, objname, binsize=binsize, minmag=minmag, maxmag=maxmag, $
                      hst=hst, halo=halo, core=core, ccut=ccut, log=log, smf=smf, $
                      help=help

;+
; NAME:
;	GALAXY_LFUNCTION
;
; PURPOSE:
;	Generate a galaxy luminosity function, ready to write to postscript.
;
; INPUTS:
;	objname : string galaxy name
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;	hst : keyword specifying HST data
;
; OUTPUTS:
;
; PROCEDURES USED:
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 September 5, UofA
;-


	!x.ticklen=0.05
        p_thick = !p.thick & chars = !p.charsize
        charth = !p.charthick & xth = !x.thick & yth = !y.thick
        !p.thick = 3 & !p.charsize=2. & !p.charthick=3 & !x.thick=3 & !y.thick=3

        if keyword_set(help) then begin
            print
            print, 'Syntax : galaxy_lfunction, objname, binsize=binsize, '
            print, 'minmag=minmag, maxmag=maxmag, hst=hst, halo=halo, '
            print, 'core=core, ccut=ccut, log=log, smf=smf, help=help'
            print
            return
        endif

        colortable1

        trgb_readata, objname, datapath, data, infobase, hst=hst, halo=halo, core=core, ccut=ccut

        if not keyword_set(binsize) then binsize = 0.05
        if not keyword_set(minmag) then minmag = min(data[3,*])
        if not keyword_set(maxmag) then maxmag = max(data[3,*])

        trgb_lfunction, data[3,*], data[4,*], lf, binsize=binsize, minmag=minmag, maxmag=maxmag

        window, 0, xs=500, ys=500 ; linear scaling
        if keyword_set(smf) then begin
            error = trgb_error_func(objname,datapath,lf.mag,minmag=minmag,$ ; smoothed error function
                                    maxmag=maxmag,hst=hst,halo=halo,core=core)
            error = 3.*error
            smf_lfunction, lf, error, phi
            plot_smflf, infobase.truename, lf, phi, log=log
        endif else begin
            plot_lfunction, infobase.truename, lf, log=log
        endelse

        okay = 'N'
        read, okay, prompt='Generate a postscript file (Y/[N])? '

        if strupcase(okay) eq 'Y' then begin
            path = trgb_datapath()
            path = path[3]      ; plot subdirectory
            
            if keyword_set(smf) then begin ; Sakai method
                
                if keyword_set(log) then begin ; log plot
                    
                    if (keyword_set(ccut) and keyword_set(halo)) then $
                      fname = path+objname+'_halo_ccut_log_smf' else $
                      if (keyword_set(ccut) and keyword_set(core)) then $
                      fname = path+objname+'_core_ccut_log_smf' else $
                      if keyword_set(ccut) then fname = path+objname+'_ccut_log_smf' else $
                      if keyword_set(halo) then fname = path+objname+'_halo_log_smf' else $
                      if keyword_set(core) then fname = path+objname+'_core_log_smf' else $
                      fname = path+objname+'_log_smf'
                    
                endif else begin ; linear plot
                    
                    if (keyword_set(ccut) and keyword_set(halo)) then $
                      fname = path+objname+'_halo_ccut_lin_smf' else $
                      if (keyword_set(ccut) and keyword_set(core)) then $
                      fname = path+objname+'_core_ccut_lin_smf' else $
                      if keyword_set(ccut) then fname = path+objname+'_ccut_lin_smf' else $
                      if keyword_set(halo) then fname = path+objname+'_halo_lin_smf' else $
                      if keyword_set(core) then fname = path+objname+'_core_lin_smf' else $
                      fname = path+objname+'_lin_smf'
                    
                endelse
                
            endif else begin
                
                if keyword_set(log) then begin ; log plot
                    
                    if (keyword_set(ccut) and keyword_set(halo)) then $
                      fname = path+objname+'_halo_ccut_log' else $
                      if (keyword_set(ccut) and keyword_set(core)) then $
                      fname = path+objname+'_core_ccut_log' else $
                      if keyword_set(ccut) then fname = path+objname+'_ccut_log' else $
                      if keyword_set(halo) then fname = path+objname+'_halo_log' else $
                      if keyword_set(core) then fname = path+objname+'_core_log' else $
                      fname = path+objname+'_log'
                    
                endif else begin ; linear plot
                    
                    if (keyword_set(ccut) and keyword_set(halo)) then $
                      fname = path+objname+'_halo_ccut_lin' else $
                      if (keyword_set(ccut) and keyword_set(core)) then $
                      fname = path+objname+'_core_ccut_lin' else $
                      if keyword_set(ccut) then fname = path+objname+'_ccut_lin' else $
                      if keyword_set(halo) then fname = path+objname+'_halo_lin' else $
                      if keyword_set(core) then fname = path+objname+'_core_lin' else $
                      fname = path+objname+'_lin'
                    
                endelse
                
            endelse
            
            ps_open, fname, /ps_fonts
            device, /times, /inches
            if keyword_set(smf) then plot_smflf, infobase.truename, lf, phi, log=log else $
              plot_lfunction, infobase.truename, lf, log=log
            ps_close

        endif
            
        !x.ticklen=0
        !p.thick = p_thick & !p.charsize=chars & !p.charthick=charth
        !x.thick = xth & !y.thick=yth

return
end
