function get_maskpops_supergrid, supergrid, nsuper=nsuper, superstring=superstring
; jm12may08ucsd - read the supergrid parameter file    
    supergrid_paramfile = getenv('CLASH_DIR')+'/maskpops/maskpops_supergrid.par'
    super = yanny_readone(supergrid_paramfile)
    
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
