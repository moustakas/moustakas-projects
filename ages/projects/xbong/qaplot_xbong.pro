pro qaplot_xbong
; jm09dec16ucsd

    path = ages_path(/projects)+'xbong/'
    bong = rsex(path+'xbong_xray.txt')
    bonggal = string(bong.pass,format='(I3.3)')+'/'+$
      string(bong.aper,format='(I3.3)')
    all = read_ages_gandalf(/ppxf)
    allgal = string(all.pass,format='(I3.3)')+'/'+$
      string(all.aper,format='(I3.3)')

    indx = where_array(bonggal,allgal)
    data = all[indx]
    specfit = read_ages_gandalf_specfit(data)
    
    im_mwrfits, data, path+'xbong_gandalf.fits', /clobber
    qaplot_ages_gandalf_specfit, data, specfit, $
      psfile=path+'qaplot_xbong.ps'
    
stop    
    
return
end
    
