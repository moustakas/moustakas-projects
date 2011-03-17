function nfgs_read_leda
; jm05jul21uofa - read the LEDA output for the NFGS galaxies

    ledapath = nfgs_path(/ned)
    ledaname = 'nfgs_leda.dat'

    leda = parse_leda(ledaname,ledapath=ledapath)
    
return, leda
end
    
