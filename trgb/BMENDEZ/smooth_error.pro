FUNCTION SMOOTH_ERROR, err

GETELEMENT_VECTOR, err, MAX(err), e2
truncerr = err[0:e2]
smooth = SMOOTH2(truncerr, 4)
err[0:e2] = smooth


  return, err
end

