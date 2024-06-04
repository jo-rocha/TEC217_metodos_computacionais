% Definicao da matrix
A = [5 1 1; 3 4 1; 3 3 6];
b = [5; 6; 0];
x0 = [0; 0; 0];

% tolerancia
tol = 0.05;
max_iter = 1000;

function x = jacobi(A, b, x0, tol, max_iter)
    fprintf('ITERACAO\t X1\t\t X2\t\t X3\t\t ERRO\n');
    n = length(b);
    x = x0;
    x_new = zeros(n, 1);
    iter = 0;
    while iter < max_iter
        for i = 1:n
            sigma = 0;
            for j = 1:n
                if j != i
                    sigma = sigma + A(i, j) * x(j);
                end
            end
            x_new(i) = (b(i) - sigma) / A(i, i);
        end
        err = (norm(x_new - x) / norm(x_new));
        fprintf('%d\t\t %.6f\t %.6f\t %.6f\t %.6f\n', iter, x_new(1), x_new(2), x_new(3), err);
        if err < tol
            x = x_new;
            break;
        end
        x = x_new;
        iter = iter + 1;
    end
end

x = jacobi(A, b, x0, tol, max_iter);

disp('Solução:');
disp(x);
