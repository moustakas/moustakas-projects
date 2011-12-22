function read_08perez, zbin, nolog=nolog, h100=h100, salpeter=salpeter
; jm11dec20ucsd - read the Perez-Gonzalez et al. 2008 stellar mass
; functions; h=0.7, Omega0=0.3, OmegaL=0.7, Salpeter IMF

; use the I-band selected SMF    
    
    if n_elements(h100) eq 0 then h100 = 0.7
    
    file = getenv('CATALOGS_DIR')+'/08perez-gonzalez/08perez-gonzalez.fits.gz'
    data = mrdfits(file,1,/silent)
    ndata = n_elements(data)

    if n_elements(zbin) ne 0 then begin
       case zbin of
          1: these = where(data.z0 eq 0.2,ndata)
          2: these = where(data.z0 eq 0.4,ndata)
          3: these = where(data.z0 eq 0.6,ndata)
          4: these = where(data.z0 eq 0.8,ndata)
          else: message, 'ZBIN must be 1-4!'
       endcase
    endif else these = lindgen(ndata)

    mf = replicate({mass: 0.0, phi: 0.0, phierr: 0.0},ndata)
    mf.mass = data[these].logm

    mf.phi = data[these].logphi2
    mf.phierr = data[these].e_logphi2
    good = where(finite(mf.phi),ngood)
    mf = mf[good]
    
    if (keyword_set(salpeter) eq 0) then mf.mass -= 0.26 ; Salpeter-->Chabrier
    
    if keyword_set(nolog) then begin
       mf.phi = 10D^mf.phi
       mf.phierr = mf.phierr*mf.phi*alog(10)
    endif

return, mf
end
