function markov = MarkovParam(Ff,Gg,Cc,Dd,n_markov,modelname)
    markov_params = zeros(n_markov, size(Cc, 1), size(Gg, 2));

    markov_params(1, :, :) = Dd;  
    for k = 1:n_markov
        markov_params(k, :, :) = Cc * (Ff^(k-1)) * Gg;  
    end
    
    if ~isempty(modelname)
        figure;
        for i = 1:size(Cc, 1)
            for j = 1:size(Gg, 2)
                subplot(size(Cc, 1), size(Gg, 2), (i-1)*size(Gg, 2) + j);
                stem(0:n_markov-1, squeeze(markov_params(:, i, j)));
                title(sprintf('Markov Parameter G_{%d%d}', i, j));
                xlabel('Time step (k)');
                ylabel('Markov Parameter');
            end
        end
        sgtitle(modelname)
    end
    
    markov = cell(n_markov, 1);

    for i = 1:n_markov
        markov{i} = squeeze(markov_params(i, :, :)); 
    end
end