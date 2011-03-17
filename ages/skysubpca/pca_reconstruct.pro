;*** PCA_RECONSTRUCT
; AIM: perform reconstruction of a spectrum 
; INPUTS: pcs and espectra
; OUTPUT: reconstructed spectrum
;******************************************************************


FUNCTION PCA_RECONSTRUCT, pcs, espec

nbin = (size(espec))[1]
nrecon = n_elements(pcs)

if (size(espec))[2] ne nrecon then message,'pcs or espec array wrong size'

return, total(rebin(Transpose(pcs),nbin,nrecon)*espec,2)

END
