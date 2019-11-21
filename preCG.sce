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

function [L,D] = IC0(AD, AL, col_ind, row_ptr)
    n = length(AD);
    nz = length(AL);
    D = AD;
    L = zeros(nz,1);
    for i=1:n
        w = zeros(i-1,1);
        for j=row_ptr(i):row_ptr(i+1)-1
            w(col_ind(j)) = AL(j);
            for k=row_ptr(col_ind(j)):row_ptr(col_ind(j)+1)-1
                w(col_ind(j)) = w(col_ind(j)) - L(k) * w(col_ind(k));
            end
            L(j) = w(col_ind(j)) / D(col_ind(j));
        end
        for j=row_ptr(i):row_ptr(i+1)-1
            D(i) = D(i) - L(j) * w(col_ind(j));
        end
    end
endfunction

function z = LDLTsolve(L, D, r, col_ind, row_ptr)
    n = length(D);
    z = r;
    for i=1:n
        for j=row_ptr(i):row_ptr(i+1)-1
            z(i) = z(i) - L(j) * z(col_ind(j));
        end
    end
    for i=1:n
        z(i) = z(i) / D(i);
    end
    for i=n:-1:1
        for j=row_ptr(i+1)-1:-1:row_ptr(i)
            z(col_ind(j)) = z(col_ind(j)) - L(j) * z(i);
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

//initilaze vector
for i = 1:n
    x(i) = 1.0;
    zk(i) = 0.0;
    rk(i) = 0.0;
    zk(i) = 0.0;
    r0(i) = 0.0;
    z0(i) = 0.0;
end

//Preconditioned CG method with IC(0)
[y] = mat_vec_mul(x);
r0 = b - y;
[L,D] = IC0(A_D, val, col_ind, row_ptr);
z0 = LDLTsolve(L, D, r0, col_ind, row_ptr);
p = z0;
//main loop of CG method
while 1
    itr = itr + 1;
    [q] = mat_vec_mul(p);
    alpha = (r0' * z0) / (p' * q);
    x = x + alpha * p;
    rk = r0 - alpha * q;
    if (norm(rk, 2) / norm(b, 2) <= 1.0E-12) then
        break;
    end
    zk = LDLTsolve(L, D, rk, col_ind, row_ptr);
    beta = (rk' * zk) / (r0' * z0);
    p = zk + beta * p;
    r0 = rk;
    z0 = zk;
end

disp(itr);