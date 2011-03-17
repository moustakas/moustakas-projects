pro plot_ages_template_basis
; jm08sep18nyu 

    aa = read_ages(/ispec)

    path = im_primus_path()
    templatefile = 'ages_basis_080918.fits.gz'
    basis = mrdfits(path+templatefile,1)
    basis = basis[1:n_elements(basis)-1] ; exclude the red template
    nbasis = n_elements(basis)
    
    ww = lindgen(nbasis)
    for jj = 0, nbasis-1 do ww[jj] = where(basis[jj].ages_id eq aa.ages_id)
    niceprint, aa[ww].ages_id, basis.ages_id
;   ww = ww[0:2]

    ages_display_spectrum, aa[ww], plottype=2, $
      nsmooth=0, labeltype=3, /postscript, /pdf, $
      psname=path+'ages_template_basis_080918.ps'

stop    
    
return
end
