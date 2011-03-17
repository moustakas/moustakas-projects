pro render_hbhg_plot, dust, nodust, pos=pos, label=label

    lineratio, dust, 'h_beta', 'h_gamma', '', '', hbhg, hbhg_err, $
      snrcut=0.0, /nolog, index=indx
    ebv = get_ebv(hbhg,decrement_err=hbhg_err,/hbhg,ebv_err=ebv_err)

    ploterror, nodust[indx].ebv_hahb, ebv, $
      nodust[indx].ebv_hahb_err, ebv_err, position=pos[*,0], $
      psym=6, xsty=1, ysty=1, xrange=[0,3], yrange=[0,3], $
      xtitle=textoidl('E(B-V) [H\alpha/H\beta]'), $
      ytitle=textoidl('E(B-V) [H\beta/H\gamma]')
    djs_oplot, !x.crange, !y.crange, line=0, thick=4
    legend, label, /right, /bottom, box=0, margin=0

    ww = where(dust[indx].h_beta_ew[0] gt 15)
    splog, median(ebv[ww]-nodust[indx[ww]].ebv_hahb)
    ploterror, dust[indx].h_beta_ew[0], ebv-nodust[indx].ebv_hahb, $
      dust[indx].h_beta_ew[1], sqrt(nodust[indx].ebv_hahb_err^2+ebv_err^2), $
      position=pos[*,1], /noerase, psym=6, xsty=1, ysty=1, /xlog, yrange=[-3,3], $
      xtitle=textoidl('EW(H\beta) (\AA)'), $
      ytitle='Residuals (mag)', xrange=[1,300]
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=4

return
end    

pro plotsings_check_hbhg
; jm10jul02ucsd - plot some example PPXF fitting results for the paper 

    nuclear = read_sings_log12oh_samples(/nuclear)
    drift20 = read_sings_log12oh_samples(/drift20)
    drift56 = read_sings_log12oh_samples(/drift56)

    nuclear_nodust = read_sings_log12oh_samples(/nodust_nuclear)
    drift20_nodust = read_sings_log12oh_samples(/nodust_drift20)
    drift56_nodust = read_sings_log12oh_samples(/nodust_drift56)

    psfile = sings_path(/projects)+'log12oh/check_hbhg.ps'
    im_plotconfig, 6, pos, psfile=psfile, yspace=1.0, $
      height=[3.5,2.5]
    render_hbhg_plot, nuclear, nuclear_nodust, pos=pos, label='Nuclear'
    render_hbhg_plot, drift20, drift20_nodust, pos=pos, label='Circumnuclear'
    render_hbhg_plot, drift56, drift56_nodust, pos=pos, label='Radial Strip'
    im_plotconfig, psfile=psfile, /psclose, /gzip

return
end
    
