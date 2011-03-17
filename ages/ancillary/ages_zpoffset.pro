function ages_zpoffset, nozpoffset=nozpoffset
; jm10jan07ucsd - zeropoint offsets derived in BOOTES_ZEROPOINTS 
    zpoffset = [0.0,0.0,bootes_zpoffset()]
    if keyword_set(nozpoffset) then zpoffset = zpoffset*0.0
return, zpoffset
end
