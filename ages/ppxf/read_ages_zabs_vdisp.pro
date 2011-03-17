function read_ages_zabs_vdisp, pass, zabsvdisppath=zabsvdisppath
; jm09dec13ucsd - read the output from AGES_GET_ZABS_VDISP and perform
; some tweaks to the measured quantities

    if (n_elements(pass) eq 0) then begin
       doc_library, 'read_ages_zabs_vdisp'
       return, -1
    endif
    
    if (n_elements(zabsvdisppath) eq 0) then begin
       version = ages_version(/ppxf_specfit)
       zabsvdisppath = ages_path(/ppxf)+'fluxed/zabs_vdisp/'+version+'/'
    endif

    zabsvdispfile = zabsvdisppath+'zabs_vdisp_'+$
      string(pass,format='(I0)')+'.fits.gz'
    splog, 'Reading '+zabsvdispfile
    zabsvdisp = mrdfits(zabsvdispfile,1)
    nobj = n_elements(zabsvdisp)
       
; clean up the ZABS and VDISP values 
    zabsvdisp = struct_addtags(zabsvdisp,$
      replicate({zabs_raw: 0.0, vdisp_raw: 0.0},nobj))

    zabsvdisp.zabs_raw = zabsvdisp.zabs
    errzero = where(zabsvdisp.zabs_err eq 0.0,nerrzero)
    if (nerrzero ne 0) then zabsvdisp[errzero].zabs = zabsvdisp[errzero].z ; AGES
    zabsvdisp.vdisp_raw = zabsvdisp.vdisp
    zabsvdisp.vdisp = ages_ppxf_vdisp_good(zabsvdisp) ; cull out crap

; some redshifts require additional tweaks; these estimates were
; derived "by-hand" by calling AGES_GET_ZABS_VDISP with parameters
; that were different from the default values
    case strtrim(pass,2) of
       '114': begin
; this object is a broad-line AGN with very different continuum and
; emission-line redshifts
          fix = where(zabsvdisp.aper eq 4,nfix)
          zabsvdisp[fix].zabs = 0.35119579
          zabsvdisp[fix].zabs_raw = 0.35119579
          zabsvdisp[fix].zabs_err = 1.6632965e-05
       end
       '314': begin
; the original AGES redshift is more correct
          fix = where(zabsvdisp.aper eq 154,nfix)
          zabsvdisp[fix].zabs_raw = zabsvdisp[fix].zabs
          zabsvdisp[fix].zabs = zabsvdisp[fix].z
          zabsvdisp[fix].zabs_err = 0.0
       end
       else: 
    endcase

return, zabsvdisp
end
