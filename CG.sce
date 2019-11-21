clear
function [y] = mat_vec_mul(matvecmul_vec)

//matric-vector multiplication for A_D
for i = 1:n
    y(i) = A_D(i) * matvecmul_vec(i);
end

//matric-vector multiplication for A_L
for i = 1:n
    for j = row_ptr(i):row_ptr(i+1)-1
        y(i) = y(i) + val(j) * matvecmul_vec(col_ind(j));
    end
end

//matric-vector multiplication for A_LT
for i = 1:n
    for j = row_ptr(i):row_ptr(i+1)-1
        y(col_ind(j)) = y(col_ind(j)) + val(j) * matvecmul_vec(i);
    end
end
endfunction

if ( exists('func1') == 0 )
  exec('func1.sci');
end
if ( exists('func2') == 0 )
  exec('func2.sci');
end
if ( exists('GenLS') == 0 )
  exec('GenLS.sci');
end
itr = 0;
N = 11;
[A_D, val, col_ind, row_ptr, b] = GenLS(N);
n = length(A_D);

//initilaze vector x
for i = 1:n
    x(i) = 0.0;
end

//CG method
[y] = mat_vec_mul(x);
r0 = b - y;
p = r0;
//main loop of CG method
while 1
    itr = itr + 1;
    [q] = mat_vec_mul(p);
    alpha = (r0' * r0) / (p' * q);
    x = x + alpha * p;
    rk = r0 - alpha * q;
    if (norm(rk, 2) / norm(b, 2) <= 1.0E-12) then
        break;
    end
    beta = (rk' * rk) / (r0' * r0);
    p = rk + beta * p;
    r0 = rk;
end

disp(itr);