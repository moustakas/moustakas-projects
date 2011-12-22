function read_08baldry, blue=blue, red=red, nolog=nolog, $
  h100=h100, salpeter=salpeter
; jm11dec20ucsd - read the Baldry et al. 2008 stellar mass functions;
; h=0.7, Omega0=0.3, OmegaL=0.7, Kroupa IMF

    if n_elements(h100) eq 0 then h100 = 0.7
    
    file = getenv('CATALOGS_DIR')+'/08baldry/gsmf-BGD08.txt'
    data = rsex(file)

    good = where(data.ngal gt 1,ndata)
    mf = replicate({mass: 0.0, phi: 0.0, phierr: 0.0},ndata)
    mf.mass = data[good].mass
    mf.phi = data[good].phi
    mf.phierr = data[good].phierr

    if keyword_set(salpeter) then mf.mass += 0.26 ; Chabrier-->Salpeter
    
    if (keyword_set(nolog) eq 0) then begin
       mf.phierr = mf.phierr/mf.phi/alog(10)
       mf.phi = alog10(mf.phi)
    endif

return, mf
end
