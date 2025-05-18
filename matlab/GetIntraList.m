function pairs = GetIntraList(params)

Ns = params.Ns;
Nc = params.Nc;

%a = repmat(1:Ns,Ns,1); b = a'; I = find(a > b); a = a(I); b = b(I);
%cell_pair = [b,a];
% Initialize pairs list
pairs = zeros(Nc*(Ns^2-Ns)/2, 2);

%pairs = 

% Counter for pairs
count = 0;

% Iterate over each column
for c = 1:Nc
    % Iterate over each combination of pairs within the column
    for i = 1:Ns-1
        for j = i+1:Ns
            % Increment counter
            count = count + 1;
            % Convert i,j to linear indices
            idx_i = sub2ind([Ns Nc], i, c);
            idx_j = sub2ind([Ns Nc], j, c);
            % Assign linear indices to the pairs list
            pairs(count, :) = [idx_i, idx_j];
        end
    end
end