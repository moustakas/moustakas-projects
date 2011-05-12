function mz_massbins, nmassbins
; jm11may08ucsd - note the REVERSE!
    nmassbins = 5
    massbins = replicate({massbin: 0.0, lomass: 0.0, himass: 0.0, label: '', $
      color: '', psym: -1, symsize: -1, line: 0},nmassbins)
    massbins.lomass = reverse([9.0,9.5,10.0,10.5,11.0])
    massbins.himass = reverse([9.5,10.0,10.5,11.0,12.0])
    massbins.massbin = reverse([9.25,9.75,10.25,10.75,11.25]) ; note the last bin

    mm = 'log(M/M_{'+sunsymbol()+'})'
    massbins.label = reverse(['9.0<'+mm+'<9.5','9.5<'+mm+'<10',$
      '10<'+mm+'<10.5','10.5<'+mm+'<11',mm+'>11'])

    massbins.color = ['black','forest green','firebrick','navy','orange']
    massbins.psym = [6,5,9,16,15]
    massbins.symsize = [2.5,2.9,2.5,2.5,2.5]
    massbins.line = [0,3,5,4,1]
    
return, massbins
end
    
