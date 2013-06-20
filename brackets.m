% The brackets <N|psi> are sum(nq*evan(f,a),2)
%    <N|Dpsi> is [nq*evan(f,a), ndqr*evan(f,a)]

global N aop nhn nq ndqr

aop = diag(sqrt(1:N),1);
nhn = ((0:N).*(-1:N-1)).';
nq = diag(1./sqrt(factorial(0:N)));
ndqr = diag(sqrt((1:N)./factorial(0:N-1)),-1);
