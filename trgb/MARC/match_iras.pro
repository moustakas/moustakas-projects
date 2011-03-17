
pro match_iras

	readfast, '/deep0/marc/trgb/vpylm_neighbor_05.11', iras, ncol=11
        readfast, '/deep1/ioannis/trgb/data_object.txt', gal, ncol=2

        niras = n_elements(iras[0,*])
        ngal = n_elements(gal[0,*])
        match = lonarr(niras)

;        for k = 0L, niras-1L do begin
;            temp = where((iras[0,k] gt round(gal[0,*]-1.)) and (iras[0,k] lt round(gal[0,*]+1.)) $
;                             and (iras[1,k] gt round(gal[1,*]-1.)) and (iras[1,k] lt round(gal[1,*]+1.)),mcount)
;           if mcount gt 1L then match[k] = -1L else match[k] = temp
;    endfor

        for j = 0L, ngal-1L do lmatch = where(round(gal[0,j]) eq round(iras[0,*]))
        for k = 0L, ngal-1L do bmatch = where(round(gal[1,k]) eq round(iras[1,*]))

        prin

        good = match[where(match ne -1L)]
        irasgood = iras[0:4,good]
        galgood = gal[*,good]

        srt = sort(gal[0,good])

        print, [irasgood[0:1,srt],galgood[0:1,srt]]
        
        
        
stop
        
return
end
