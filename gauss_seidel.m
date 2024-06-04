% Definicao da matrix
A = [1 0.5 -0.1; 0.2 1 -0.2; -0.1 -0.2 1];
b = [0.2; -2; 1];
x0 = [0; 0; 0];

% tolerancia
tol = 0.05;
max_iter = 1000;

function x = gauss_seidel(A, b, x0, tol, max_iter)
    n = length(b);
    x = zeros(n, 1);
    fprintf('ITERACAO\t X1\t\t X2\t\t X3\t\t ERRO\n');

    for k = 1:max_iter
      x_old = x;

      for i = 1:n
        sum1 = 0;
        sum2 = 0;

        for j = 1:i-1
          sum1 += A(i, j) * x(j);
        end

        for j = i+1:n
          sum2 += A(i, j) * x_old(j);
        end

        x(i) = (b(i) - sum1 - sum2) / A(i, i);
      end

      err = norm(x - x_old) / norm(x);
      fprintf('%d\t\t %.6f\t %.6f\t %.6f\t %.6f\n', k, x(1), x(2), x(3), err);

      if err < tol
        disp('---------------------------------------');
        fprintf('%d iterações para convergir\n', k);
        return;
      end
    end

    fprintf('Número máximo de iterações %d\n', max_iter);
end

function sassenfeld(A)
    % criterio de Sassenfeld
    n = size(A, 1);
    B = zeros(n,1);
    for i = 1:n
      sum_a = 0;
      sum_b = 0;
      for j = 1:i-1
        sum_a = sum_a + abs(A(i,j)) * B(j);
      end
      for j = i+1:n
        sum_b = sum_b + abs(A(i,j));
      end
      sassenfeld(i) = (sum_a + sum_b) / abs(A(i,i));
    end
    maxim = max(sassenfeld);
    if maxim >= 1
      error('O critério de Sassenfeld não foi satisfeito');
    else
      disp('---------------------------------------');

      fprintf('B = %d\t O critério de Sassenfeld foi satisfeito\n', maxim);

      disp('---------------------------------------');
    end
end


x = gauss_seidel(A, b, x0, tol, max_iter);
sassenfeld(A);

disp('Solução:');
disp(x);
