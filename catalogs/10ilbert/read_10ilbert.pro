function read_10ilbert, zbin, sf=sf, active=active, interm=interm, $
  quiescent=quiescent, sty=sty, nolog=nolog, h100=h100, salpeter=salpeter
; jm11dec20ucsd - read the Ilbert et al. 2010 stellar mass functions;
; h=0.7, Omega0=0.3, OmegaL=0.7, Chabrier IMF

; for SF combine active and intermediate
    if keyword_set(sf) then begin
       act = read_10ilbert(zbin,/active,/nolog)
       int = read_10ilbert(zbin,/interm,/nolog)
       sf = act
       sf.mass = (act.mass+int.mass)/2.0
       sf.phi = act.phi+int.phi
       sf.phierr = sqrt(act.phierr^2+int.phierr^2)
       sf.phierr = sf.phierr/sf.phi/alog(10)
       sf.phi = alog10(sf.phi)

       good = where(sf.phierr lt 10)
       sf = sf[good]
       return, sf
    endif
    
    if n_elements(h100) eq 0 then h100 = 0.7

    if keyword_set(sty) then method = 'STY' else method = 'Vmax'
    type = 'All'
    if keyword_set(active) then type = 'SB' ; actively star-forming 
    if keyword_set(interm) then type = 'Spi2' ; intermediate
    if keyword_set(quiescent) then type = 'Ell' ; quiescent

; redshift bins:
;   zbin0 = 0.05-0.2
;   zbin1 =  0.2-0.4
;   zbin2 =  0.4-0.6
;   zbin3 =  0.6-0.8
;   zbin4 =  0.8-1.0
;   zbin5 =  1.0-1.2
;   zbin6 =  1.2-1.5
;   zbin7 =  1.5-2.0
    
    if n_elements(zbin) eq 0 then message, 'ZBIN input required!'
    if (zbin lt 1) or (zbin gt 8) then message, '0<=ZBIN<=7'

    file = getenv('CATALOGS_DIR')+'/10ilbert/MF_'+method+'_K_'+type+strtrim(zbin,2)+'.dat'
    splog, file
    readcol, file, mass, phi, philo, phiup, format='F,F,F,F', /silent, comment='#'
    ndata = n_elements(mass)

    mf = replicate({mass: 0.0, phi: 0.0, phierr: 0.0},ndata)
    mf.mass = mass
    mf.phi = phi
    mf.phierr = (philo+phiup)/2.0

    if keyword_set(salpeter) then mf.mass += 0.26 ; Chabrier-->Salpeter
    
    if keyword_set(nolog) then begin
       mf.phi = 10D^mf.phi
       mf.phierr = mf.phierr*mf.phi*alog(10)
    endif

return, mf
end
