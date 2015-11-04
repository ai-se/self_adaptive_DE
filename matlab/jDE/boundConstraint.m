function vi=boundConstraint(vi, lu);

[NP, D] = size(vi);  % the population size and the problem's dimension

xl = repmat(lu(1, :), NP, 1);
xu = repmat(lu(2, :), NP, 1);

%% check the lower bound
pos = vi < xl;
vi(pos) = 2 .* xl(pos) - vi(pos);
pos_=vi(pos) > xu(pos);
vi(pos(pos_)) = xu(pos(pos_));

%% check the upper bound
pos = vi > xu;
vi(pos) = 2 .* xu(pos) - vi(pos);
pos_=vi(pos) < xl(pos);
vi(pos(pos_)) = xl(pos(pos_));