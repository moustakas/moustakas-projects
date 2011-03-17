
function doublet, x, p, dp, inst_res = inst_res, line_ctr = line_ctr, dr=dr

;------------------------------
; Following Rupke et al 2005

n_doublets = n_elements(dr)

yval = fltarr(n_elements(x)) + 1.0
cspeed = 2.99792e5

for ii = 0, n_doublets - 1 do begin

  x1   = p[0 + ii*4]  ; Angstroms           
  x2   = p[0 + ii*4] + line_ctr[1 + ii*2] - line_ctr[0 + ii*2]
  b    = p[1 + ii*4]  ; km/s
  tau0 = p[2 + ii*4]  ; 
  cf   = p[3 + ii*4]

  tau1 = tau0 * exp(-1*(x - x1)^2 / (x1 * b/cspeed)^2)
  tau2 = tau0 * dr[ii] * exp(-1*(x - x2)^2 / (x2 * b/cspeed)^2)
  yval = yval * (1 - cf + cf * exp(-tau1 - tau2))

endfor

return, yval

end
