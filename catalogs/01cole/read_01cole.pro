function read_01cole, nolog=nolog, h100=h100, salpeter=salpeter
; jm11dec20ucsd - read the Cole+2001 stellar mass functions;
; h=1, Salpeter IMF

    oldh100 = 1.0
    if n_elements(h100) eq 0 then h100 = 0.7
    
    file = getenv('CATALOGS_DIR')+'/01cole/cole01_table4.sex'
    data = rsex(file)
    ndata = n_elements(data)
    
    mf = replicate({mass: 0.0, phi: 0.0, phierr: 0.0},ndata)
    mf.mass = data.mass+2*alog10(oldh100/h100)
    mf.phi = data.phi*(h100/oldh100)^3
    mf.phierr = data.phierr*(h100/oldh100)^3

    if (keyword_set(salpeter) eq 0) then mf.mass -= 0.26 ; Salpeter-->Chabrier
    
    if (keyword_set(nolog) eq 0) then begin
       mf.phierr = mf.phierr/mf.phi/alog(10)
       mf.phi = alog10(mf.phi)
    endif

return, mf
end
