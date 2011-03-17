function trgb_lfmodel, e, x, cut, beta, alpha=alpha, nosmooth=nosmooth

	model = 10.^(alpha*x)*(x ge 0.) + 10.^(-cut+beta*x)*(x lt 0)

        if keyword_set(nosmooth) then return, model else begin ; boxcar smooth

            gsmooth = model-model
            for i = 0L, n_elements(gsmooth[*,0])-1L do $
              gsmooth[i,*] = total(model[e.eminus[i]:e.eplus[i],*],1)/e.nsmooth[i] 

            gsmooth[e.a1[0]:e.a1[e.x2-1]] = model[e.a1[0]:e.a1[e.x2-1]] ; fix the endpoint kink

            return, gsmooth

        endelse

end
