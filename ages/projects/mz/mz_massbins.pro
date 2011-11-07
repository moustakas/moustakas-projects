function mz_massbins, nmassbins, rev=rev, masscut=masscut, masskeep=masskeep
; jm11may08ucsd - note the REVERSE!

    mm = 'log(M/M_{'+sunsymbol()+'})'

    minmass = 9.1D
    maxmass = 11.1D
    binsize = 0.2D
    nmassbins = long((maxmass-minmass)/binsize)

    massbins = replicate({massbin: 0.0, lomass: 0.0, himass: 0.0, label: '', $
      color: 'black', psym: 6, symsize: 1.0, line: 0},nmassbins)
    
    massbins.lomass = findgen(nmassbins)*binsize+minmass
    massbins.himass = massbins.lomass+binsize
    massbins.massbin = massbins.lomass+binsize/2.0
;   niceprint, lomass, himass    

    massbins.label = string(massbins.lomass,format='(F4.1)')+'<'+mm+'<'+$
      string(massbins.himass,format='(F4.1)')

    if n_elements(masscut) eq 0 then masscut = min(massbins.lomass)
    masskeep = where(massbins.lomass ge masscut)
    
    if keyword_set(rev) then begin
       massbins = massbins[masskeep]
       massbins.massbin = reverse(massbins.massbin)
       massbins.lomass = reverse(massbins.lomass)
       massbins.himass = reverse(massbins.himass)
       massbins.label = reverse(massbins.label)
       masskeep = reverse(masskeep)
    endif
    
;   nmassbins = 5
;   massbins = replicate({massbin: 0.0, lomass: 0.0, himass: 0.0, label: '', $
;     color: '', psym: -1, symsize: -1, line: 0},nmassbins)
;   massbins.lomass = reverse([9.0,9.5,10.0,10.5,11.0])
;   massbins.himass = reverse([9.2,10.0,10.5,11.0,12.0])
;   massbins.massbin = reverse([9.25,9.75,10.25,10.75,11.25]) ; note the last bin
;   massbins.label = reverse(['9.0<'+mm+'<9.5','9.5<'+mm+'<10',$
;     '10<'+mm+'<10.5','10.5<'+mm+'<11',mm+'>11'])

;   massbins.color = ['black','forest green','firebrick','navy','orange']
;   massbins.psym = [6,5,9,16,15]
;   massbins.symsize = [2.5,2.9,2.5,2.5,2.5]
;   massbins.line = [0,3,5,4,1]
    
return, massbins
end
    
