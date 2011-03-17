function read_ages_spec1d, pass1, unfluxed=unfluxed, $
  silent=silent, fix_redshifts=fix_redshifts
; jm09nov13ucsd - read all the observed spectra (after PCA
;   sky-subtraction) associated with a given PASS number (e.g., 101) 
; jm09dec09ucsd - added FIX_REDSHIFTS keyword to tweak or repair wrong
;   redshifts in AGES

    if (n_elements(pass1) eq 0) then begin
       doc_library, 'read_ages_spec1d'
       return, -1
    endif
    if (n_elements(pass1) ne 1) then begin
       splog, 'PASS must be a scalar'
       return, -1
    endif

    base_spec1dpath = ages_path(/spec1d)
    if keyword_set(unfluxed) then dir = 'unfluxed' else dir = 'fluxed'
    spec1dpath = base_spec1dpath+dir+'/after_skysubpca/'

    pass = string(pass1,format='(I3.3)')
    file = spec1dpath+'ages_'+pass+'.fits.gz'
    if (file_test(file,/reg) eq 0) then begin
       splog, 'Pass file '+file+' not found'
       return, -1
    endif
    splog, 'Reading '+file
    spec1d = mrdfits(file,1,silent=silent)

; fix/tweak some redshifts
    if keyword_set(fix_redshifts) then begin
       case strtrim(pass1,2) of
          '602': begin
             fix = where(spec1d.aper eq 238,nfix)
             if (nfix ne 0) then spec1d.z[fix] = 0.693 ; z=0.69-->0.693
          end 
          else: 
       endcase
    endif
    
return, spec1d
end
