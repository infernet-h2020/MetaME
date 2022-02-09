function [orth_matrix, idx] = orth_norm(idxBM, basis)


    fprintf('Performing Gram Schmidt method...');
    aux = basis(:,idxBM);
    basis(:,idxBM) = basis(:,1);
    basis(:,1) = aux;
    orth_matrix = mgrscho(basis);
    idx = [];
    N = size(orth_matrix,2);
    for i = 1:N-1
        for j = i+1:N
            if(abs(dot(orth_matrix(:,i),orth_matrix(:,j))) > 1e-3)
                'not orthogonal'
                dot(orth_matrix(:,i),orth_matrix(:,j))
                [i j]
                pause
            else
                if i == 1
                    idx = cat(1, j, idx);
                end
            end
        end
    end
    aux = orth_matrix(:,idxBM);
    orth_matrix(:,idxBM) = orth_matrix(:,1);
    orth_matrix(:,1) = aux;
   fprintf( 'done\n');
end
