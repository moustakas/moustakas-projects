function atlas_read_leda
; jm05jul21uofa - read the LEDA output for the spectral atlas

    ledapath = atlas_path(/ned)
    ledaname = 'atlas_leda.dat'

    leda = parse_leda(ledaname,ledapath=ledapath)
    
return, leda
end
    
