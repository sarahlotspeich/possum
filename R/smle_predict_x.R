smle_predict_x = function(row_data, bspline_coeff) {
  row_bspline = row_data[rep(1, nrow(bspline_coeff)), 
                         colnames(bspline_coeff)]
  xhat = sum(row_bspline * bspline_coeff)
  return(xhat)
}