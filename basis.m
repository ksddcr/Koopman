function  psi = basis(x,p)

    psi = zeros(p+1,1);
    % Polynomial basis
    for i = 1:p+1
        psi(i) = x.^(i-1);
    end
%     %Hermite polynomial
%     psi(1) = x.^0;
%     psi(2) = x;
%     for i = 3:p+1
%         psi(i) = x.*psi(i-1) - (i-1)*psi(i-2);
%     end
end