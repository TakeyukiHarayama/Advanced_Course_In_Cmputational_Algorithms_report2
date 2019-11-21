function [AD, AL, col_ind, row_ptr,b] = GenLS(N)
  AD      =  ones((N-1)^2,1) * 4;
  AL      = -ones(2*N^2-6*N+4,1);
  col_ind = zeros(2*N^2-6*N+4,1);
  row_ptr = zeros((N-1)^2+1,1);
  b       = zeros((N-1)^2,1);
  
  if ( exists('func1') == 0 )
    exec('func1.sci');
  end
  if ( exists('func2') == 0 )
    exec('func2.sci');
  end
  
  row_ptr(1) = 1;
  k=1;
  l=1;

  //----
  j = 1;
  i = 1;
  b(k) = func2(i,j-1,N) + func2(i-1,j,N);
  row_ptr(k+1) = row_ptr(k);
  k = k + 1;
  
  //----
  for i=2:N-2
    col_ind(l) = func1(i-1,j,N);
    l = l + 1;
    
    b(k) = func2(i,j-1,N);
    row_ptr(k+1) = row_ptr(k) + 1;
    k = k + 1;
  end

  //----
  i = N - 1;
  col_ind(l) = func1(i-1,j,N);
  l = l + 1;
  
  b(k) = func2(i,j-1,N) + func2(i+1,j,N);
  row_ptr(k+1) = row_ptr(k) + 1;
  k = k + 1;

  for j=2:N-2
    //----
    i = 1;

    col_ind(l) = func1(i,j-1,N);
    l = l + 1;
    
    b(k) = func2(i-1,j,N);    
    row_ptr(k+1) = row_ptr(k) + 1;
    k = k + 1;
    
    //----
    for i=2:N-2
      col_ind(l  ) = func1(i,j-1,N);
      col_ind(l+1) = func1(i-1,j,N);
      l = l + 2;
      
      b(k) = 0;
      row_ptr(k+1) = row_ptr(k) + 2;
      k = k + 1;
    end
    
    //----
    i = N - 1;
    col_ind(l  ) = func1(i,j-1,N);
    col_ind(l+1) = func1(i-1,j,N);
    l = l + 2;
    
    b(k) = func2(i+1,j,N);
    row_ptr(k+1) = row_ptr(k) + 2;
    k = k + 1;
  end

  //----
  j = N - 1;
  i = 1;
  col_ind(l) = func1(i,j-1,N);
  l = l + 1;
  
  b(k) = func2(i-1,j,N) + func2(i,j+1,N);
  row_ptr(k+1) = row_ptr(k) + 1;
  k = k + 1;
  
  //----
  for i=2:N-2
    col_ind(l) = func1(i,j-1,N);
    col_ind(l+1) = func1(i-1,j,N);
    l = l + 2;
    
    b(k) = func2(i,j+1,N);
    row_ptr(k+1) = row_ptr(k) + 2;
    k = k + 1;
  end
  
  //----
  i = N - 1;
  col_ind(l) = func1(i,j-1,N);
  col_ind(l+1) = func1(i-1,j,N);
  
  b(k) = func2(i+1,j,N) + func2(i,j+1,N);
  row_ptr(k+1) = row_ptr(k) +2;
endfunction
