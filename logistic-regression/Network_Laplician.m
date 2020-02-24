function Network_Laplician(I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate graph using Erdos-Renyi model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while true
    pErdosRenyi = 0.1;
    p_data_gen = 1 - sqrt(1 - pErdosRenyi);    
    Adj = rand(I);
    idx1 = (Adj >= p_data_gen);
    idx2 = (Adj < p_data_gen);
    Adj(idx1) = 0;
    Adj(idx2) = 1;
    
    %Adj = randraw('bernoulli',pErdosRenyi,I,I);
    NotI = ~eye(I);
    Adj = Adj.*NotI; % eliminate zeros on diag
    Adj = or(Adj,Adj'); % symmetrize, undirected
    degree = diag(sum(Adj));  %Degree matrix
    StandardLap = degree - Adj;
    lambda = sort(eig(StandardLap));
    if lambda(2) > 0.1 %&& lambda(2) < 0.5
        fprintf(['The Erdos-Renyi graph is generated. Algebraic Connectivity: ', num2str(lambda(2)),'\n']);
        break;
    end
    
%     testAdj = (eye(I)+Adj)^I; % check if G is connected
%     
%     if ~any(any(~testAdj))
%         fprintf('The ER graph is connected\n');
%         break;
%     end
end

%%%% MH weight %%%%
A = zeros(I);
for i=1:I
    i_link=find(Adj(i,:)>0);
    for j=1:I
        if i~=j && sum(find(j==i_link))>0
            A(i,j)=1/(max(degree(i,i),degree(j,j))+1);
        end
    end
end
Lap = diag(sum(A)) - A; 
save('laplacian.mat', 'Lap')
end