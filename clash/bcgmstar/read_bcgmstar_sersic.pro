function read_bcgmstar_sersic, cluster, results=results, radius=radius, $
  model=model, band=band
; jm14aug12siena - read the appropriate Sersic-fitting results, which
; varies from cluster to cluster
    
; note that the individual structures cannot be 'stacked' because of a
; varying number of bandpasses (the covariance matrices are all
; different sizes)    
    
    if n_elements(cluster) eq 0 then begin
       splog, 'CLUSTER input required'
       return, -1
    endif

    sersicpath = bcgmstar_path(/sersic)
    
    case cluster of
; double-Sersic 
       'a209': suffix = '-alphabetazero'
       'a2261': suffix = '-alphabetazero'
; single-Sersic 
       'macs2129': suffix = '-alphabetazero'
       'macs1149': suffix = '-alphabetazero'
       'rxj2129': suffix = '-alphabetazero'
       'macs1720': suffix = '-alphabetazero'
       else: suffix = ''
    endcase

    sersic = mrdfits(sersicpath+cluster+'-allsersic'+suffix+'.fits.gz',1,/silent)
    if arg_present(results) then results = mrdfits(sersicpath+cluster+$
      '-allsersic'+suffix+'-results.fits.gz',1,/silent)

; chose one or more bandpasses
    if n_elements(band) ne 0 then begin
       match, strtrim(sersic.band,2), band, m1, m2
       if n_elements(m2) ne n_elements(band) then message, 'Unknown bandpass(es)!'
       srt = sort(m2) & m1 = m1[srt] & m2 = m2[srt]
       sersic = sersic[m1]
       if arg_present(results) then results = results[m1]
    endif
    
; evaluate the model
    if n_elements(radius) ne 0 then begin
       nfilt = n_elements(sersic)
       model = fltarr(n_elements(radius),nfilt)
       for ib = 0, nfilt-1 do begin
; double-Sersic
          if cluster eq 'a209' or cluster eq 'a2261' then begin
             model[*,ib] = bcgmstar_sersic2_func(radius,params=sersic[ib],/allbands)
          endif else begin
; single-Sersic
             model[*,ib] = bcgmstar_sersic_func(radius,params=sersic[ib],/allbands)
          endelse 
       endfor 
    endif 
    
return, sersic
end
    
