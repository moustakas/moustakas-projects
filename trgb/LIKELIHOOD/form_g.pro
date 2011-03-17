; mdjune2000
function form_g, x, cut, alpha=alpha, beta=beta, gamma=gamma
; function to evaluate a distribution function consisting of three linear pieces.  

	model = 10.^(alpha*x)*(x ge 0.) + $
          10.^(gamma*x)*(x lt 0 and x ge -cut/gamma) + $
          10.^(-cut)*10.^(beta*(x+cut/gamma))*(x lt -cut/gamma)

return, model
end
