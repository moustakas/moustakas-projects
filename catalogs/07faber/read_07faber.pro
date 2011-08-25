function read_07faber, zbin=zbin, blue=blue, red=red, combo17=combo17
; jm11aug25ucsd - read the Faber+07 B-band luminosity functions; both
; have h=0.7 and Vega mags

    if (n_elements(zbin) eq 0) then zbin = 1
    path = getenv('CATALOGS_DIR')+'/07faber/'

    if keyword_set(combo17) then begin
       sample = 'a'
       if keyword_set(blue) then sample = 'b'
       if keyword_set(red) then sample = 'r'
       case zbin of
          1: zsuffix = '0.2_0.4'
          2: zsuffix = '0.4_0.6'
          3: zsuffix = '0.6_0.8'
          4: zsuffix = '0.8_1.0'
          5: zsuffix = '1.0_1.2'
       endcase
       file = path+'combo_'+sample+'_'+zsuffix+'.dat'
       readcol, file, absmag, phi, phierr, format='F,F,X,X,X,F', /silent
       absmag = absmag - 0.10 ; Vega --> AB from Table 1 in Willmer+06
    endif else begin
       sample = 'all'
       if keyword_set(blue) then sample = 'blue'
       if keyword_set(red) then sample = 'red'
       case zbin of
          1: zsuffix = '0.20_0.40'
          2: zsuffix = '0.40_0.60'
          3: zsuffix = '0.60_0.80'
          4: zsuffix = '0.80_1.00'
          5: zsuffix = '1.00_1.20'
       endcase
       file = path+'deep2_optimal_'+sample+'_vmax_'+zsuffix+'.da'
       readcol, file, absmag, phi, phierr, ngal, format='F,X,X,F,F,L', /silent
       absmag = absmag - 0.11 ; Vega --> AB from Table 1 in Willmer+06
    endelse
    nmag = n_elements(absmag)
    
    lf = replicate({absmag: 0.0, phi: 0.0, phierr: 0.0},nmag)
    lf.absmag = absmag
    lf.phi = phi
    lf.phierr = phierr
    
return, lf
end
    
    
