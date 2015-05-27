pro A2109_gallery

    ff = ['g','r','z']
    for ii = 0, 2 do begin
       im = mrdfits('A2109-'+ff[ii]+'.fits',0,hdr)
       sxaddpar, hdr, 'FILTER', ff[ii]
       mwrfits, im, 'A2109-image-'+ff[ii]+'.fits', hdr, /create
    endfor

    name = 'A2109'
    infile = name+'-image.in'
    openw, lun, infile, /get_lun
    printf, lun, 'B'
    printf, lun, name+'-image-g.fits'
    printf, lun, ''
    printf, lun, 'G'
    printf, lun, name+'-image-r.fits'
    printf, lun, ''
    printf, lun, 'R'
    printf, lun, name+'-image-z.fits'
    printf, lun, ''
    printf, lun, 'indir /Users/ioannis/tmp/'
    printf, lun, 'outname '+name+'-SDSS'
    printf, lun, 'outdir .'
    printf, lun, 'noiselum 0.15'
    printf, lun, 'show 0'
    printf, lun, 'legend 0'
    printf, lun, 'testfirst 0'
    free_lun, lun

    cmd = 'python '+getenv('IMPY_DIR')+'/image/trilogy.py '+infile
    splog, cmd
    spawn, cmd, /sh

return
end

