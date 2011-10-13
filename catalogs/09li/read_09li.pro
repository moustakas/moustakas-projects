function read_09li, nolog=nolog, h100=h100, salpeter=salpeter
; jm11oct10ucsd - read the Li & White (2009) stellar mass functions;
; h=1, Omega0=0.25, OmegaL=0.75, Chabrier IMF

    if n_elements(h100) eq 0 then h100 = 0.7
    
    path = getenv('CATALOGS_DIR')+'/09li/'
    file = path+'massfun.dr72bbright0'

    data = rsex(file)
    ndata = n_elements(data)
    
    mf = replicate({mass: 0.0, phi: 0.0, phierr: 0.0},ndata)
    mf.mass = data.meanmass+2*alog10(1.0/h100)
    mf.phi = data.phi*(h100/1.0)^3
    mf.phierr = data.phierr*(h100/1.0)^3

    if keyword_set(salpeter) then mf.mass += 0.26 ; Chabrier-->Salpeter
    
    if (keyword_set(nolog) eq 0) then begin
       mf.phierr = mf.phierr/mf.phi/alog(10)
       mf.phi = alog10(mf.phi)
    endif

return, mf
end
