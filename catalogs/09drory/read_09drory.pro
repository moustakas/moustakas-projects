function read_09drory, zbin, sf=sf, quiescent=quiescent, nolog=nolog, $
  h100=h100, salpeter=salpeter
; jm11dec20ucsd - read the Drory et al. 2009 stellar mass functions;
; h=0.7, Omega0=0.3, OmegaL=0.7, Chabrier IMF

    if n_elements(h100) eq 0 then h100 = 0.7
    
    file = getenv('CATALOGS_DIR')+'/09drory/cosmos_mf_2009.dat'

    data = rsex(file)
    ndata = n_elements(data)

    if n_elements(zbin) ne 0 then begin
       case zbin of
          1: these = where(data.zmin eq 0.2D,ndata)
          2: these = where(data.zmin eq 0.4D,ndata)
          3: these = where(data.zmin eq 0.6D,ndata)
          4: these = where(data.zmin eq 0.8D,ndata)
          else: message, 'ZBIN must be 1-4!'
       endcase
    endif else these = lindgen(ndata)

    mf = replicate({mass: 0.0, phi: 0.0, phierr: 0.0},ndata)
    mf.mass = data[these].logm

    mf.phi = data[these].mftot
    mf.phierr = data[these].mftot_err
    if keyword_set(sf) then begin
       mf.phi = data[these].mfblue
       mf.phierr = data[these].mfblue_err
    endif
    if keyword_set(quiescent) then begin
       mf.phi = data[these].mfred
       mf.phierr = data[these].mfred_err
    endif
    good = where(mf.phi lt 0D,ngood)
    mf = mf[good]
    
    if keyword_set(salpeter) then mf.mass += 0.26 ; Chabrier-->Salpeter
    
    if keyword_set(nolog) then begin
       mf.phi = 10D^mf.phi
       mf.phierr = mf.phierr*mf.phi*alog(10)
    endif

return, mf
end
