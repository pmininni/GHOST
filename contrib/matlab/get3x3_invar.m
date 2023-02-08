% Get tensor invariants for 3x3 tensors, given 1d array with
% tensor components in the order:
%
%    [bI bII bIII] = get3x3_invar(b,1)
%
%    a11, a12, a13, a22, a23, a33
%
%    a11, a12, a13, a21, a22, a23, a31, a32, A33
%
function [bI bII bIII] = get3x3_invar(b,itrans)

  if length(b) ~= 6 & length(b) ~= 9
    error('Invalid 1d array length');
  end

  if length(b) == 6
    A(1,1) = b  (1); A(1,2) = b  (2); A(1,3) = b(3);
    A(2,1) = A(1,2); A(2,2) = b  (4); A(2,3) = b(5);
    A(3,1) = A(1,3); A(3,2) = A(2,3); A(3,3) = b(6);
  else
    A(1,1) = b  (1); A(1,2) = b  (2); A(1,3) = b(3);
    A(2,1) = b  (4); A(2,2) = b  (5); A(2,3) = b(6);
    A(3,1) = b  (7); A(3,2) = b  (8); A(3,3) = b(9);
  end

  if 0
    eigA = eig(A);
    bI = sum(eigA);
    bII = eigA(1)*eigA(2) + eigA(1)*eigA(3) + eigA(2)*eigA(3);
    bIII = eigA(1)*eigA(2)*eigA(3);
  end


  if nargin > 1 & itrans > 0
    A = A';
  end

  bI = A(1,1)+A(2,2)+A(3,3);

  bII = 0;
  for i = 1:3
    for k = 1:3
      bII = bII + A(i,k)*A(k,i);
    end
  end
  bII = bII/2;

  bIII = 0;
  for i = 1:3
    for j = 1:3
      for k = 1:3
        bIII = bIII + A(i,k)*A(k,j)*A(j,i);
      end
    end
  end
  bIII = bIII/3;
