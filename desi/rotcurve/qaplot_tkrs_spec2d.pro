pro qaplot_tkrs_spec2d

    outpath = getenv('IM_PROJECTS_DIR')+'/desi/'
    cat = mrdfits(outpath+'sim_zcat.fits.gz',1)
    ngal = n_elements(cat)

    nplot = 10
    these = shuffle_indx(ngal,num=nplot)
    srt = sort(cat[these].z)

    cat = cat[these[srt]]
    prefix = string(cat.id,format='(I7.7)')

    psfile = 'qaplot_tkrs_spec2d.ps'
    ps_start, psfile
    pos = cglayout([1,nplot],ygap=1)
    for ii = 0, nplot-1 do begin
       spec2d = mrdfits(outpath+'sim_'+prefix[ii]+'.fits.gz',0,hdr,/silent)
;      spec2d = spec2d/max(spec2d)
       cgimage, spec2d, noerase=ii gt 0, position=pos[*,ii];, $
;        stretch='Log', /negative, minvalue=median(spec2d)
    endfor
    ps_end, /png, gs_path='/usr/local/bin/', /delete_ps

stop    
    
return
end
    
