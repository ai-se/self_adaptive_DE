%**************************************************************************************************
%Reference:  C. Garcia-Martinez, M. Lozano, F. Herrera, D. Molina, and A.M. Sanchez,
%                   "Global and local real-coded genetic algorithms based on parent-centric
%                      crossover operators," European Journal of Operatoral Research, vol. 185,
%                      pp. 1088¨C1113, 2008.
%
% Note: We obtained the C source code from the authors, and developed the
%           MATLAB source code based on the C source code.
%**************************************************************************************************

clc;
clear all;
tic;

'GL-25'

format long;
format compact;

% Choose the problems to be tested. Please note that for test functions F7
% and F25, the global optima are out of the initialization range. For these
% two test functions, we do not need to judge whether the variable violates
% the boundaries during the evolution after the initialization.
problemSet = [1 : 6 8 : 24];
for problemIndex = 1 : 23

    problem = problemSet(problemIndex)

    % Define the dimension of the problem
    n = 30;

    switch problem

        case 1

            % lu: define the upper and lower bounds of the variables
            lu = [-100 * ones(1, n); 100 * ones(1, n)];
            % Load the data for this problem
            load sphere_func_data
            A = []; M = []; a = []; alpha = []; b = [];

        case 2

            lu = [-100 * ones(1, n); 100 * ones(1, n)];
            load schwefel_102_data
            A = []; M = []; a = []; alpha = []; b = [];

        case 3

            lu = [-100 * ones(1, n); 100 * ones(1, n)];
            load high_cond_elliptic_rot_data
            A = []; a = []; alpha = []; b = [];

            if n == 2, load elliptic_M_D2,
            elseif n == 10, load elliptic_M_D10,
            elseif n == 30, load elliptic_M_D30,
            elseif n == 50, load elliptic_M_D50,
            end

        case 4

            lu = [-100 * ones(1, n); 100 * ones(1, n)];
            load schwefel_102_data
            A = []; M = []; a = []; alpha = []; b = [];

        case 5

            lu = [-100 * ones(1, n); 100 * ones(1, n)];
            load schwefel_206_data
            M = []; a = []; alpha = []; b = [];

        case 6

            lu = [-100 * ones(1, n); 100 * ones(1, n)];
            load rosenbrock_func_data
            A = []; M = []; a = []; alpha = []; b = [];

        case 7

            lu = [0 * ones(1, n); 600 * ones(1, n)];
            load griewank_func_data
            A = []; a = []; alpha = []; b = [];

            c = 3;
            if n == 2, load griewank_M_D2,
            elseif n == 10, load griewank_M_D10,
            elseif n == 30, load griewank_M_D30,
            elseif n == 50, load griewank_M_D50,
            end

        case 8

            lu = [-32 * ones(1, n); 32 * ones(1, n)];
            load ackley_func_data
            A = []; a = []; alpha = []; b = [];

            if n == 2, load ackley_M_D2,
            elseif n == 10, load ackley_M_D10,
            elseif n == 30, load ackley_M_D30,
            elseif n == 50, load ackley_M_D50,
            end

        case 9

            lu = [-5 * ones(1, n); 5 * ones(1, n)];
            load rastrigin_func_data
            A = []; M = []; a = []; alpha = []; b = [];

        case 10

            lu = [-5 * ones(1, n); 5 * ones(1, n)];
            load rastrigin_func_data
            A = []; a = []; alpha = []; b = [];
            if n == 2, load rastrigin_M_D2,
            elseif n == 10, load rastrigin_M_D10,
            elseif n == 30, load rastrigin_M_D30,
            elseif n == 50, load rastrigin_M_D50,
            end

        case 11

            lu = [-0.5 * ones(1, n); 0.5 * ones(1, n)];
            load weierstrass_data
            A = []; a = []; alpha = []; b = [];
            if n == 2, load weierstrass_M_D2, ,
            elseif n == 10, load weierstrass_M_D10,
            elseif n == 30, load weierstrass_M_D30,
            elseif n == 50, load weierstrass_M_D50,
            end

        case 12

            lu = [-pi * ones(1, n); pi * ones(1, n)];
            load schwefel_213_data
            A = []; M = []; o = [];

        case 13

            lu = [-3 * ones(1, n); 1 * ones(1, n)];
            load EF8F2_func_data
            A = []; M = []; a = []; alpha = []; b = [];

        case 14

            lu = [-100 * ones(1, n); 100 * ones(1, n)];
            load E_ScafferF6_func_data
            if n == 2, load E_ScafferF6_M_D2, ,
            elseif n == 10, load E_ScafferF6_M_D10,
            elseif n == 30, load E_ScafferF6_M_D30,
            elseif n == 50, load E_ScafferF6_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 15

            lu = [-5 * ones(1, n); 5 * ones(1, n)];
            load hybrid_func1_data
            A = []; M = []; a = []; alpha = []; b = [];

        case 16

            lu = [-5 * ones(1, n); 5 * ones(1, n)];
            load hybrid_func1_data
            if n == 2, load hybrid_func1_M_D2,
            elseif n == 10, load hybrid_func1_M_D10,
            elseif n == 30, load hybrid_func1_M_D30,
            elseif n == 50, load hybrid_func1_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 17

            lu = [-5 * ones(1, n); 5 * ones(1, n)];
            load hybrid_func1_data
            if n == 2, load hybrid_func1_M_D2,
            elseif n == 10, load hybrid_func1_M_D10,
            elseif n == 30, load hybrid_func1_M_D30,
            elseif n == 50, load hybrid_func1_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 18

            lu = [-5 * ones(1, n); 5 * ones(1, n)];
            load hybrid_func2_data
            if n == 2, load hybrid_func2_M_D2,
            elseif n == 10, load hybrid_func2_M_D10,
            elseif n == 30, load hybrid_func2_M_D30,
            elseif n == 50, load hybrid_func2_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 19

            lu = [-5 * ones(1, n); 5 * ones(1, n)];
            load hybrid_func2_data
            if n == 2, load hybrid_func2_M_D2,
            elseif n == 10, load hybrid_func2_M_D10,
            elseif n == 30, load hybrid_func2_M_D30,
            elseif n == 50, load hybrid_func2_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 20

            lu = [-5 * ones(1, n); 5 * ones(1, n)];
            load hybrid_func2_data
            if n == 2, load hybrid_func2_M_D2,
            elseif n == 10, load hybrid_func2_M_D10,
            elseif n == 30, load hybrid_func2_M_D30,
            elseif n == 50, load hybrid_func2_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 21

            lu = [-5 * ones(1, n); 5 * ones(1, n)];
            load hybrid_func3_data
            if n == 2, load hybrid_func3_M_D2,
            elseif n == 10, load hybrid_func3_M_D10,
            elseif n == 30, load hybrid_func3_M_D30,
            elseif n == 50, load hybrid_func3_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 22

            lu = [-5 * ones(1, n); 5 * ones(1, n)];
            load hybrid_func3_data
            if n == 2, load hybrid_func3_HM_D2,
            elseif n == 10, load hybrid_func3_HM_D10,
            elseif n == 30, load hybrid_func3_HM_D30,
            elseif n == 50, load hybrid_func3_HM_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 23

            lu = [-5 * ones(1, n); 5 * ones(1, n)];
            load hybrid_func3_data
            if n == 2, load hybrid_func3_M_D2,
            elseif n == 10, load hybrid_func3_M_D10,
            elseif n == 30, load hybrid_func3_M_D30,
            elseif n == 50, load hybrid_func3_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 24

            lu = [-5 * ones(1, n); 5 * ones(1, n)];
            load hybrid_func4_data
            if n == 2, load hybrid_func4_M_D2,
            elseif n == 10, load hybrid_func4_M_D10,
            elseif n == 30, load hybrid_func4_M_D30,
            elseif n == 50, load hybrid_func4_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

        case 25

            lu = [2 * ones(1, n); 5 * ones(1, n)];
            load hybrid_func4_data
            if n == 2, load hybrid_func4_M_D2,
            elseif n == 10, load hybrid_func4_M_D10,
            elseif n == 30, load hybrid_func4_M_D30,
            elseif n == 50, load hybrid_func4_M_D50,
            end
            A = []; a = []; alpha = []; b = [];

    end

     % Record the best results
     outcome = []; 

    % Main body which is developed according to the C source code provided
    % by the authors
    time = 1;

    % The total number of runs
    totalTime = 1;

    while time <= totalTime

        rand('seed', sum(100 * clock));

        popsize = 400;

        % Initialize the male group of the global search
        GM = repmat(lu(1, :), popsize, 1) + rand(popsize, n) .* (repmat(lu(2, :) - lu(1, :), popsize, 1));
        GM_fit = benchmark_func(GM, problem, o, A, M, a, alpha, b);

        % Initialize the female group of the global search
        [fitval, fitindex] = sort(GM_fit);

        GF = GM(fitindex(1 : 200), :);
        GF_fit = GM_fit(fitindex(1 : 200), :);

        % Record the number of offspring generated by the female parent
        numTimesMother = zeros(200, 1);

        % Judge the global search or the local search: j = 0 mean the global
        % search and j = 1 means the local search
        j = 0;

        FES = popsize;
        while FES < n * 10000

            if FES > n * 10000 * 0.25 & j == 0

                popsize = 100;

                [fitval, fitindex] = sort(GM_fit);
                % Initialize the male group of the local search
                GM = GM(fitindex(1 : 100), :);
                GM_fit = GM_fit(fitindex(1 : 100), :);

                % Initialize the female group of the local search
                GF = GM(1 : 5, :);
                GF_fit = GM_fit(1 : 5, :);

                % Updat the number of offspring generated by the female
                % parent
                numTimesMother = zeros(5, 1);

                j = 1;

            end

            % Uniform fertility selection
            [Timesval, Timesindex] = min(numTimesMother);

            mother = GF(Timesindex, :);
            mother_fit = GF_fit(Timesindex, :);

            numTimesMother(Timesindex) = numTimesMother(Timesindex) + 1;

            % Negative assortative mating
            temp = floor(rand(1, 5) * popsize) + 1;
            temp_father = GM(temp, :);
            temp_father_fit = GM_fit(temp, :);

            Edistance = sum((ones(5, 1) * mother-temp_father) .^ 2, 2);

            [Edistanceval, Edistanceindex] = max(Edistance);

            father = temp_father(Edistanceindex, :);
            father_fit = temp_father_fit(Edistanceindex, :);

            %PBX-alfa
            temp = abs(mother - father);
            L = max(lu(1, :), mother - temp .* 0.8);
            U = min(lu(2, :), mother + temp .* 0.8);
            offspring = L + rand(1, n) .* (U-L);

            % Evaluation
            offspring_fit = benchmark_func(offspring, problem, o, A, M, a, alpha, b);

            % Replace the worst individual in the population
            if offspring_fit < mother_fit

                [fitval, fitindex] = max(GM_fit);
                GM(fitindex, :) = mother;
                GM_fit(fitindex, :) = mother_fit;

                GF(Timesindex, :) = offspring;
                GF_fit(Timesindex, :) = offspring_fit;

                numTimesMother(Timesindex, :) = 0;

            elseif offspring_fit < max(GM_fit)

                [fitval, fitindex] = max(GM_fit);
                GM(fitindex, :) = offspring;
                GM_fit(fitindex, :) = offspring_fit;

            end

            FES = FES + 1;

        end

        outcome = [outcome min(GF_fit)];

        time = time + 1;

    end

    sort(outcome)
    mean(outcome)
    std(outcome)

end

toc;