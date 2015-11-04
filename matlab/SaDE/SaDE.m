%**************************************************************************************************
%Reference:  A. K. Qin, V. L. Huang, and P. N. Suganthan,¡°Differential evolution
%                     algorithm with strategy adaptation for global numerical optimization,¡±
%                     IEEE Trans. Evolut. Comput., vol. 13, no. 2, pp. 398¨C417, Apr. 2009.
%
% Note: We obtained the MATLAB source code from the authors, and did some
%           minor revisions in order to solve the 25 benchmark test functions,
%           however, the main body was not changed.
%**************************************************************************************************

clc;
clear all;
tic;

format long;
format compact;

'SaDE'

% Choose the problems to be tested. Please note that for test functions F7
% and F25, the global optima are out of the initialization range. For these
% two test functions, we do not need to judge whether the variable violates
% the boundaries during the evolution after the initialization.
problemSet = [1 : 6 8 : 10];
for prolemIndex = 1 : 9

    problem = problemSet(prolemIndex)
    
    % Define the dimension of the problem
    D = 30;

    NP = 50;

    switch problem

        case 1

            % lu: define the upper and lower bounds of the variables
            lu = [-100 * ones(1, D); 100 * ones(1, D)];
            % Load the data for this problem
            load sphere_func_data
            A = []; M = []; a = []; alpha = []; b = [];

        case 2

            lu = [-100 * ones(1, D); 100 * ones(1, D)];
            load schwefel_102_data
            A = []; M = []; a = []; alpha = []; b = [];

        case 3

            lu = [-100 * ones(1, D); 100 * ones(1, D)];
            load high_cond_elliptic_rot_data
            A = []; a = []; alpha = []; b = [];

            if D == 2, load elliptic_M_D2,
            elseif D == 10, load elliptic_M_D10,
            elseif D == 30, load elliptic_M_D30,
            elseif D == 50, load elliptic_M_D50,
            end

        case 4

            lu = [-100 * ones(1, D); 100 * ones(1, D)];
            load schwefel_102_data
            A = []; M = []; a = []; alpha = []; b = [];

        case 5

            lu = [-100 * ones(1, D); 100 * ones(1, D)];
            load schwefel_206_data
            M = []; a = []; alpha = []; b = [];

        case 6

            lu = [-100 * ones(1, D); 100 * ones(1, D)];
            load rosenbrock_func_data
            A = []; M = []; a = []; alpha = []; b = [];

        case 7

            lu = [0 * ones(1, D); 600 * ones(1, D)];
            load griewank_func_data
            A = []; a = []; alpha = []; b = [];

            if D == 2, load griewank_M_D2,
            elseif D == 10, load griewank_M_D10,
            elseif D == 30, load griewank_M_D30,
            elseif D == 50, load griewank_M_D50,
            end

        case 8

            lu = [-32 * ones(1, D); 32 * ones(1, D)];
            load ackley_func_data
            A = []; a = []; alpha = []; b = [];

            if D == 2, load ackley_M_D2,
            elseif D == 10, load ackley_M_D10,
            elseif D == 30, load ackley_M_D30,
            elseif D == 50, load ackley_M_D50,
            end

        case 9

            lu = [-5 * ones(1, D); 5 * ones(1, D)];
            load rastrigin_func_data
            A = []; M = []; a = []; alpha = []; b = [];

        case 10

            lu = [-5 * ones(1, D); 5 * ones(1, D)];
            load rastrigin_func_data
            A = []; a = []; alpha = []; b = [];
            if D == 2, load rastrigin_M_D2,
            elseif D == 10, load rastrigin_M_D10,
            elseif D == 30, load rastrigin_M_D30,
            elseif D == 50, load rastrigin_M_D50,
            end

        case 11

            lu = [-0.5 * ones(1, D); 0.5 * ones(1, D)];
            load weierstrass_data
            A = []; a = []; alpha = []; b = [];
            if D == 2, load weierstrass_M_D2, ,
            elseif D == 10, load weierstrass_M_D10,
            elseif D == 30, load weierstrass_M_D30,
            elseif D == 50, load weierstrass_M_D50,
            end

        case 12

            lu = [-pi * ones(1, D); pi * ones(1, D)];
            load schwefel_213_data
            A = []; M = []; o = [];

        case 13

            lu = [-3 * ones(1, D); 1 * ones(1, D)];
            load EF8F2_func_data
            A = []; M = []; a = []; alpha = []; b = [];

        case 14

            lu = [-100 * ones(1, D); 100 * ones(1, D)];
            load E_ScafferF6_func_data
            if D == 2, load E_ScafferF6_M_D2, ,
            elseif D == 10, load E_ScafferF6_M_D10,
            elseif D == 30, load E_ScafferF6_M_D30,
            elseif D == 50, load E_ScafferF6_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 15

            lu = [-5 * ones(1, D); 5 * ones(1, D)];
            load hybrid_func1_data
            A = []; M = []; a = []; alpha = []; b = [];

        case 16

            lu = [-5 * ones(1, D); 5 * ones(1, D)];
            load hybrid_func1_data
            if D == 2, load hybrid_func1_M_D2,
            elseif D == 10, load hybrid_func1_M_D10,
            elseif D == 30, load hybrid_func1_M_D30,
            elseif D == 50, load hybrid_func1_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 17

            lu = [-5 * ones(1, D); 5 * ones(1, D)];
            load hybrid_func1_data
            if D == 2, load hybrid_func1_M_D2,
            elseif D == 10, load hybrid_func1_M_D10,
            elseif D == 30, load hybrid_func1_M_D30,
            elseif D == 50, load hybrid_func1_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 18

            lu = [-5 * ones(1, D); 5 * ones(1, D)];
            load hybrid_func2_data
            if D == 2, load hybrid_func2_M_D2,
            elseif D == 10, load hybrid_func2_M_D10,
            elseif D == 30, load hybrid_func2_M_D30,
            elseif D == 50, load hybrid_func2_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 19

            lu = [-5 * ones(1, D); 5 * ones(1, D)];
            load hybrid_func2_data
            if D == 2, load hybrid_func2_M_D2,
            elseif D == 10, load hybrid_func2_M_D10,
            elseif D == 30, load hybrid_func2_M_D30,
            elseif D == 50, load hybrid_func2_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 20

            lu = [-5 * ones(1, D); 5 * ones(1, D)];
            load hybrid_func2_data
            if D == 2, load hybrid_func2_M_D2,
            elseif D == 10, load hybrid_func2_M_D10,
            elseif D == 30, load hybrid_func2_M_D30,
            elseif D == 50, load hybrid_func2_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 21

            lu = [-5 * ones(1, D); 5 * ones(1, D)];
            load hybrid_func3_data
            if D == 2, load hybrid_func3_M_D2,
            elseif D == 10, load hybrid_func3_M_D10,
            elseif D == 30, load hybrid_func3_M_D30,
            elseif D == 50, load hybrid_func3_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 22

            lu = [-5 * ones(1, D); 5 * ones(1, D)];
            load hybrid_func3_data
            if D == 2, load hybrid_func3_HM_D2,
            elseif D == 10, load hybrid_func3_HM_D10,
            elseif D == 30, load hybrid_func3_HM_D30,
            elseif D == 50, load hybrid_func3_HM_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 23

            lu = [-5 * ones(1, D); 5 * ones(1, D)];
            load hybrid_func3_data
            if D == 2, load hybrid_func3_M_D2,
            elseif D == 10, load hybrid_func3_M_D10,
            elseif D == 30, load hybrid_func3_M_D30,
            elseif D == 50, load hybrid_func3_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 24

            lu = [-5 * ones(1, D); 5 * ones(1, D)];
            load hybrid_func4_data
            if D == 2, load hybrid_func4_M_D2,
            elseif D == 10, load hybrid_func4_M_D10,
            elseif D == 30, load hybrid_func4_M_D30,
            elseif D == 50, load hybrid_func4_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 25

            lu = [2 * ones(1, D); 5 * ones(1, D)];
            load hybrid_func4_data
            if D == 2, load hybrid_func4_M_D2,
            elseif D == 10, load hybrid_func4_M_D10,
            elseif D == 30, load hybrid_func4_M_D30,
            elseif D == 50, load hybrid_func4_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

    end

    outcome = [];  % record the best results

    %Main body which was provided by the authors

    time = 1;

    numst = 4;
    
    % The total number of runs
    totalTime = 1;

    while time <= totalTime

        aaaa = cell(1, numst);

        learngen = 50;

        lpcount = [];
        npcount = [];

        %Record the number of success or failure
        ns = [];
        nf = [];

        %Record the success rate
        pfit = ones(1, numst);

        %Record the median of CR
        ccm = 0.5 * ones(1, numst);

        %-----Initialize population and some arrays-------------------------------
        pop = zeros(NP, D); %initialize pop to gain speed
        XRRmin = repmat(lu(1, :), NP, 1);
        XRRmax = repmat(lu(2, :), NP, 1);
        rand('seed', sum(100 * clock));
        pop = XRRmin + (XRRmax - XRRmin) .* rand(NP, D);

        popold   = zeros(size(pop));   % toggle population
        val    = zeros(1, NP);      % create and reset the "cost array"
        DE_gbest  = zeros(1, D);      % best population member ever
        nfeval = 0;           % number of function evaluations

        %------Evaluate the best member after initialization----------------------

        val  = benchmark_func(pop, problem, o, A, M, a, alpha, b);
        [DE_gbestval, ibest] = min(val);
        DE_gbest = pop(ibest, :);

        %     ibest  = 1;            % start with first population member
        %     val(1)  = benchmark_func(pop(ibest, :), problem, o, A, M, a, alpha, b);
        %     DE_gbestval = val(1);         % best objective function value so far
        %     nfeval  = nfeval + 1;
        %     for i = 2:NP            % check the remaining members
        %       val(i) = benchmark_func(pop(i, :), problem, o, A, M, a, alpha, b);
        %       nfeval  = nfeval + 1;
        %       if (val(i) < DE_gbestval)      % if member is better
        %         ibest  = i;         % save its location
        %         DE_gbestval = val(i);
        %       end
        %     end
        %     DE_gbest = pop(ibest, :);     % best member of current iteration

        %------DE-Minimization---------------------------------------------
        %------popold is the population which has to compete. It is--------
        %------static through one iteration. pop is the newly--------------
        %------emerging population.----------------------------------------

        pm1 = zeros(NP, D);        % initialize population matrix 1
        pm2 = zeros(NP, D);        % initialize population matrix 2
        pm3 = zeros(NP, D);        % initialize population matrix 3
        pm4 = zeros(NP, D);        % initialize population matrix 4
        pm5 = zeros(NP, D);        % initialize population matrix 5
        bm  = zeros(NP, D);        % initialize DE_gbestber matrix
        ui  = zeros(NP, D);        % intermediate population of perturbed vectors
        mui = zeros(NP, D);        % mask for intermediate population
        mpo = zeros(NP, D);        % mask for old population
        rot = (0 : 1 : NP-1);        % rotating index array (size NP)
        rt  = zeros(NP);         % another rotating index array
        a1  = zeros(NP);         % index array
        a2  = zeros(NP);         % index array
        a3  = zeros(NP);         % index array
        a4  = zeros(NP);         % index array
        a5  = zeros(NP);         % index array
        ind = zeros(4);

        iter = 1;
        nfeval = NP;
        while nfeval < 10000 * D

            popold = pop;          % save the old population
            ind = randperm(4);        % index pointer array
            a1  = randperm(NP);       % shuffle locations of vectors
            rt = rem(rot + ind(1), NP);     % rotate indices by ind(1) positions
            a2  = a1(rt + 1);         % rotate vector locations
            rt = rem(rot + ind(2), NP);
            a3  = a2(rt + 1);
            rt = rem(rot + ind(3), NP);
            a4  = a3(rt + 1);
            rt = rem(rot + ind(4), NP);
            a5  = a4(rt + 1);

            pm1 = popold(a1, :);       % shuffled population 1
            pm2 = popold(a2, :);       % shuffled population 2
            pm3 = popold(a3, :);       % shuffled population 3
            pm4 = popold(a4, :);       % shuffled population 4
            pm5 = popold(a5, :);       % shuffled population 5

            bm = repmat(DE_gbest, NP, 1);

            if (iter >= learngen)
                for i = 1:numst
                    if  ~isempty(aaaa{i})
                        ccm(i) = median(aaaa{i}(:, 1));
                        d_index = find(aaaa{i}(:, 2) == aaaa{i}(1, 2));
                        aaaa{i}(d_index, :) = [];
                    else
                        ccm(i) = rand;
                    end
                end
            end

            for i = 1 : numst
                cc_tmp = [];
                for k = 1 : NP
                    tt = normrnd(ccm(i), 0.1);
                    while tt > 1 | tt < 0
                        tt = normrnd(ccm(i), 0.1);
                    end
                    cc_tmp = [cc_tmp; tt];
                end
                cc(:, i) = cc_tmp;
            end

            % Stochastic universal sampling
            rr = rand;
            spacing = 1/NP;
            randnums = sort(mod(rr : spacing : 1 + rr - 0.5 * spacing, 1));

            normfit = pfit / sum(pfit);
            partsum = 0;
            count(1) = 0;
            stpool = [];

            for i = 1 : length(pfit)
                partsum = partsum + normfit(i);
                count(i + 1) = length(find(randnums < partsum));
                select(i, 1) = count(i + 1) - count(i);
                stpool = [stpool; ones(select(i, 1), 1) * i];
            end
            stpool = stpool(randperm(NP));

            for i = 1 : numst
                atemp = zeros(1, NP);
                aaa{i} = atemp;
                index{i} = [];
                if ~isempty(find(stpool == i))
                    index{i} = find(stpool == i);
                    atemp(index{i}) = 1;
                    aaa{i} = atemp;
                end
            end

            aa = zeros(NP, D);
            for i = 1 : numst
                aa(index{i}, :) = rand(length(index{i}), D) < repmat(cc(index{i}, i), 1, D);      % all random numbers < CR are 1, 0 otherwise
            end
            mui = aa;

            % jrand
            dd = ceil(D * rand(NP, 1));
            for kk = 1 : NP
                mui(kk, dd(kk)) = 1;
            end
            mpo = mui < 0.5;         % inverse mask to mui

            for i = 1 : numst
                %-----------jitter---------
                F = [];
                m = length(index{i});
                F = normrnd(0.5, 0.3, m, 1);
                F = repmat(F, 1, D);
                if i == 1
                    ui(index{i}, :) = pm3(index{i}, :) + F .* (pm1(index{i}, :) - pm2(index{i}, :));     % differential variation
                    ui(index{i}, :) = popold(index{i}, :) .* mpo(index{i}, :) + ui(index{i}, :) .* mui(index{i}, :);   % crossover
                end
                if i == 2
                    ui(index{i}, :) = popold(index{i}, :) + F .* (bm(index{i}, :)-popold(index{i}, :)) + F .* (pm1(index{i}, :) - pm2(index{i}, :) + pm3(index{i}, :) - pm4(index{i}, :));    % differential variation
                    ui(index{i}, :) = popold(index{i}, :) .* mpo(index{i}, :) + ui(index{i}, :) .* mui(index{i}, :);   % crossover
                end
                if i == 3
                    ui(index{i}, :) = pm5(index{i}, :) + F .* (pm1(index{i}, :) - pm2(index{i}, :) + pm3(index{i}, :) - pm4(index{i}, :));    % differential variation
                    ui(index{i}, :) = popold(index{i}, :) .* mpo(index{i}, :) + ui(index{i}, :) .* mui(index{i}, :);   % crossover
                end
                if i == 4
                    ui(index{i}, :) = popold(index{i}, :) + rand .* (pm5(index{i}, :)-popold(index{i}, :)) + F .* (pm1(index{i}, :) - pm2(index{i}, :));
                end
            end

            for i = 1 : NP
                outbind = find(ui(i, :) < lu(1, :));
                XRmin = lu(1, :);
                XRmax = lu(2, :);
                if size(outbind, 2) ~= 0
                    ui(i, outbind) = XRmin(outbind) + (XRmax(outbind) - XRmin(outbind)) .* rand(1, size(outbind, 2));
                end
                outbind = find(ui(i, :) > lu(2, :));
                if size(outbind, 2) ~= 0
                    ui(i, outbind) = XRmin(outbind) + (XRmax(outbind) - XRmin(outbind)) .* rand(1, size(outbind, 2));
                end
            end

            lpcount = zeros(1, numst);
            npcount = zeros(1, numst);

            tempval = benchmark_func(ui, problem, o, A, M, a, alpha, b);  % check cost of competitor
            nfeval  = nfeval + NP;

            for i = 1 : NP

                if (tempval(i) <= val(i)) % if competitor is better than value in "cost array"

                    pop(i, :) = ui(i, :);  % replace old vector with new one (for new iteration)
                    val(i)  = tempval(i);  % save value in "cost array"

                    tlpcount = zeros(1, numst);
                    for j = 1 : numst
                        temp = aaa{j};
                        tlpcount(j) = temp(i);
                        if tlpcount(j) == 1
                            aaaa{j} = [aaaa{j}; cc(i, j) iter];
                        end
                    end
                    lpcount = [lpcount; tlpcount];

                else

                    tnpcount = zeros(1, numst);
                    for j = 1:numst
                        temp = aaa{j};
                        tnpcount(j) = temp(i);
                    end
                    npcount = [npcount; tnpcount];

                end

            end %---end for imember = 1:NP

            ns = [ns; sum(lpcount, 1)];
            nf = [nf; sum(npcount, 1)];

            if iter >= learngen,
                for i = 1 : numst
                    if (sum(ns(:, i)) + sum(nf(:, i))) == 0
                        pfit(i) = 0.01;
                    else
                        pfit(i) = sum(ns(:, i)) / (sum(ns(:, i)) + sum(nf(:, i))) + 0.01;
                    end
                end
                if ~isempty(ns), ns(1, :) = [];  end
                if ~isempty(nf), nf(1, :) = [];  end
            end
            iter = iter + 1;

        end

        outcome = [outcome min(val)];

        time = time + 1;

    end

    sort(outcome)
    mean(outcome)
    std(outcome)

end
toc;
