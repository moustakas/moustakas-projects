function read_12baldry, red=red, blue=blue, nolog=nolog, h100=h100, salpeter=salpeter
; jm12jan01ucsd - read the Baldry+12 stellar mass functions;
; h=0.7, Omega0=0.3, OmegaL=0.7, Chabrier IMF

    oldh100 = 0.7
    if n_elements(h100) eq 0 then h100 = 0.7
    
    path = getenv('CATALOGS_DIR')+'/12baldry/'

    file = path+'mfall.sex'
    if keyword_set(red) then file = path+'mfred.sex'
    if keyword_set(blue) then file = path+'mfblue.sex'

    mf = rsex(file)
    if keyword_set(red) eq 0 and keyword_set(blue) eq 0 then factor = 1D-3 else factor = 1D
    mf.phi *= factor*(h100/oldh100)^3
    mf.phierr *= factor*(h100/oldh100)^3

    if keyword_set(salpeter) then mf.mass += 0.26 ; Chabrier-->Salpeter
    
    if (keyword_set(nolog) eq 0) then begin
       mf.phierr = mf.phierr/mf.phi/alog(10)
       mf.phi = alog10(mf.phi)
    endif

return, mf
end
