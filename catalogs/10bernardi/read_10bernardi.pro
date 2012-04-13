function read_10bernardi, nolog=nolog, h100=h100, salpeter=salpeter
; jm12apr11ucsd - read the Bernardi+2010 stellar mass functions (Fig
; 22 in their paper) h=0.7, Chabrier IMF

    oldh100 = 0.7
    if n_elements(h100) eq 0 then h100 = 0.7
    
    file = getenv('CATALOGS_DIR')+'/10bernardi/10bernardi_fig22_smf.dat'
    data = rsex(file)
    ndata = n_elements(data)
    
    mf = replicate({mass: 0.0, phi: 0.0, phierr: 0.0},ndata)
    mf.mass = data.mass+2*alog10(oldh100/h100)
    mf.phi = data.phi+3.0*alog10(h100/oldh100)
    mf.phierr = total([[data.phierr_lo],[data.phierr_up]],2)/2.0 ; average
    mf.phierr = mf.phierr+3*alog10(h100/oldh100)

    if (keyword_set(salpeter) eq 0) then mf.mass -= 0.26 ; Salpeter-->Chabrier
    
return, mf
end
