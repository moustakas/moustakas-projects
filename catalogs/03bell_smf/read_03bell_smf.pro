function read_03bell_smf, blue=blue, red=red, nolog=nolog, $
  h100=h100, salpeter=salpeter
; jm11dec20ucsd - read the Bell et al. 2003 stellar mass functions;
; h=1, Omega0=0.3, OmegaL=0.7, diet Salpeter IMF

; just does the g-band MFs for now (K-band also available)
    
    oldh100 = 1.0
    if n_elements(h100) eq 0 then h100 = 0.7
    
    path = getenv('CATALOGS_DIR')+'/03bell_smf/'
    file = path+'gmf.txt'
    if keyword_set(blue) then file = path+'gmflatecol.txt'
    if keyword_set(red) then file = path+'gmfearlycol.txt'
    
    data = rsex(file)
    ndata = n_elements(data)

    mf = replicate({mass: 0.0, phi: 0.0, phierr: 0.0},ndata)
    mf.mass = data.mass-0.1+2*alog10(oldh100/h100) ; diet Salpeter-->Chabrier

    mf.phi = data.phi*(h100/oldh100)^3
    mf.phierr = (data.phiup-data.philo)/2.0*(h100/oldh100)^3

    if keyword_set(salpeter) then mf.mass += 0.26 ; Chabrier-->Salpeter
    
    if (keyword_set(nolog) eq 0) then begin
       mf.phierr = mf.phierr/mf.phi/alog(10)
       mf.phi = alog10(mf.phi)
    endif

return, mf
end
