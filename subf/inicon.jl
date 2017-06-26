function inicon(dx::Float64)
  sqr2 = dx * sqrt(2.0);
  sqr5 = dx * 0.5 * sqrt(5.0);
  sqr10 = dx * (1.0/3.0) * sqrt(10.0);
  sqr13 = dx * (1.0/3.0) * sqrt(13.0);
  sqr17 = dx * (1.0/4.0) * sqrt(17.0);
  a4 = dx * (5.0/12.0);

  c = zeros(Float64,10);

  c[1] = dx;                             # 1st order
  c[2] = sqrt(2.0) * dx;                 # 2nd order
  c[3] = 0.5 * sqrt(5.0) * dx;           # 3rd order
  c[4] = (1.0/3.0) * sqrt(10.0) * dx;    # 4th order
  c[5] = (1.0/3.0) * sqrt(13.0) * dx;    # 4th order
  c[6] = (1.0/6.0) * sqrt(13.0) * dx;    # 4th order
  c[7] = (1.0/4.0) * sqrt(17.0) * dx;    # 5th order
  c[8] = (5.0/12.0) * dx * 3.0;          # 5th order
  c[9] = (5.0/12.0) * dx;                # 5th order
  c[10] = (5.0/12.0) * dx *2.0;          # 5th order

  return sqr2,sqr5,sqr10,sqr13,sqr17,a4,c
end
