pro plotbcgmstar_montage
; jm14may15siena - build a 5x3 panel color montage of the BCGs (see
; PLOTBCGSFHS_MONTAGE) 

    coloroutpath = bcgsfhs_path(/colormosaics)
    paperpath = bcgmstar_path(/paper)

    sample = read_bcgsfhs_sample()
    sample = sample[sort(sample.mvir)]
    cluster = strtrim(sample.shortname,2)
    ncl = n_elements(sample)

    outfile = bcgmstar_path()+'bcg-mstar-finalmontage.png'
    infile = strjoin(file_search(coloroutpath+cluster+'/'+cluster+'-image-cutout.png'),' ')
    cmd = 'montage -bordercolor white -borderwidth 1 '+ $
      '-tile 5x3 -geometry +0+0 -quality 100 '+$ ; -resize 1024x1024 '+$
      infile+' '+outfile
;   cmd = 'montage '+infile+' '+outfile
    splog, cmd
;   spawn, cmd, /sh ; this doesn't work!

    cmd = 'convert -scale 50% '+outfile+' '+paperpath+$
      'bcg-mstar-finalmontage-small.png'
    splog, cmd
;   spawn, cmd ; this doesn't work!

stop    
    
return
end
    

