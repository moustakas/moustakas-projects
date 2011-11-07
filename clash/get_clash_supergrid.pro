function get_clash_supergrid, supergrid, nsuper=nsuper, superstring=superstring, $
  bcg=bcg, arc=arc
; jm11oct14ucsd - read the supergrid parameter file    
    supergrid_paramfile = getenv('CLASH_DIR')+'/clash_supergrid.par'
    super = yanny_readone(supergrid_paramfile)
    
; get the requested supergrids; default is to read the arc supergrids,
; unless /BCG
    keep = lindgen(n_elements(super))
    if keyword_set(bcg) then keep = where(super.supergrid ge 10)
    if keyword_set(arc) then keep = where(super.supergrid lt 10)
    super = super[keep]
    
    if (n_elements(supergrid) ne 0) then begin
       match2, super.supergrid, supergrid, m1, m2
       if (total(m2 eq -1) ne 0) then message, 'Unknown supergrid!'
       match, super.supergrid, supergrid, m1, m2
       srt = sort(m2) & m1 = m1[srt] & m2 = m2[srt]
       super = super[m1]
    endif
    
    superstring = string(super.supergrid,format='(I2.2)')
    nsuper = n_elements(super)
return, super
end
