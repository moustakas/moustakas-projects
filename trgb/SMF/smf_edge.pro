pro smf_edge, lf, phi, resplin, resplog, noise, trgblin, trgblog

	resplin = resplin/noise
        resplog = resplog*noise

;       maxlin = max(smooth(resplin,4),linindx)
;       maxlog = max(smooth(resplog,4),logindx)

        maxlin = max(resplin,linindx)
        maxlog = max(resplog,logindx)

        trgblin = (lf.mag)[linindx]
        trgblog = (lf.mag)[logindx]

return
end
