%**************************************************************************************************
%Author: Yong Wang
%Last Edited: July 1, 2010
%Email: ywang@csu.edu.cn
%Reference: Differential Evolution with Composite Trial Vector Generation Strategies
%                                                       and Control Parameters
%                            IEEE Transactions on Evolutionary Computation, Accepted
%**************************************************************************************************

clc;
clear all;
tic;

format long;
format compact;

'Adaptive CoDE'

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
            % Load the data for this test function
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

    % Main body

    popsize = 30;

    % Learning period
    LEP = 50;

    time = 1;

    % The total number of runs
    totalTime = 1;

    while time <= totalTime

        rand('seed', sum(100 * clock));

        % Initialize the main population
        p = repmat(lu(1, :), popsize, 1) + rand(popsize, n) .* (repmat(lu(2, :) - lu(1, :), popsize, 1));

        % Evaluate the objective function values
        fit = benchmark_func(p, problem, o, A, M, a, alpha, b);

        for j = 1 : 3

            % Success memory for each trial vector generation strategy
            sucMemo{ j } = zeros(LEP, 3);
            % Failure memory for each trial vector generation strategy
            failMemo{ j } = zeros(LEP, 3);

        end

        % The number of function evaluations (FES)
        FES = popsize;

        gen = 1;

        while FES < n * 10000

            if gen < LEP

                % The control parameter settings have the same probability
                % to be selected for each trial vector generation strategy.
                paraProb = ones(3, 3) * 1/3;

            else

                for k = 1 : 3

                    % Compute the success rate of each control parameter
                    % setting for each trial vector generation strategy
                    paraSR = sum(sucMemo{ k }) ./ (sum(sucMemo{ k }) + sum(failMemo{ k })) + 0.01;
                    % Normalized the success rates
                    paraProb(k, :) = paraSR ./ sum(paraSR);

                end

            end

            pTemp = p;
            fitTemp = fit;
            
            % uSet: the set of trial vectors
            uSet = zeros(3 * popsize, n);

            % numP{ i }: at each generation, record which control parameter
            % setting is used to produce the trial vector for each target
            % vector, please note that three trial vectors are produced for
            % each target vector
            numP = cell(popsize, 1);
            for i = 1 : popsize
                numP{ i } = zeros(3);
            end
            % numS: at each generation, record which control parameter
            % setting is used to produce the trial vector entering the next
            % population successfully
            numS = zeros(3);

            for i = 1 : popsize

                % Generate the trail vectors
                [u, tempPara] = reproduce(p, lu, i, popsize, n, paraProb);

                uSet(i * 3-2 : 3 * i, :) = u;

                FES = FES + 3;

                for k = 1 : 3

                    % Judge which control parameter setting is used to
                    % produce the trial vector for each target vector,
                    % please note that three trial vectors are generated
                    numP{ i }(k, 1) = numP{ i }(k, 1) + length(find(tempPara(k, 1) == 1));
                    numP{ i }(k, 2) = numP{ i }(k, 2) + length(find(tempPara(k, 1) == 2));
                    numP{ i }(k, 3) = numP{ i }(k, 3) + length(find(tempPara(k, 1) == 3));

                end

            end

            % Evaluate the trial vectors
            fitSet = benchmark_func(uSet, problem, o, A, M, a, alpha, b);

            for i = 1 : popsize

                % Find the best trial vector of the three trial vectors
                % generated for each target vector
                % Minindex: denotes which trial vector generation strategy
                % is used to generate the best trial vector
                [minVal, minIndex] = min(fitSet(3 * i - 2 : 3 * i, :));
                bestInd = uSet(3 * (i - 1) + minIndex, :);
                bestIndFit = fitSet(3 * (i - 1) + minIndex, :);

                % Selection
                if fit(i) >= bestIndFit

                    pTemp(i, :) = bestInd;
                    fitTemp(i, :) = bestIndFit;

                    % Judge which control parameter setting is used to
                    % generate the best trial entering the next population
                    temp = find(numP{ i }(minIndex, :) == 1);
                    numS(minIndex, temp) = numS(minIndex, temp) + 1;

                end

            end

            p = pTemp;
            fit = fitTemp;

            % Update the success memory
            for k = 1 : 3

                sucMemo{ k }(1, :) = [];
                sucMemo{ k }(LEP, :) = numS(k, :);

            end

            % Record the total number of each control parameter setting
            % used at each generation for each trial vector generation
            % strategy
            totalNum = zeros(3);
            for i = 1 : popsize

                totalNum = totalNum + numP{ i };

            end

            % Update the failure memory
            for k = 1 : 3

                failMemo{ k }(1, :) = [];
                failMemo{ k }(LEP, :) = totalNum(k, :) - numS(k, :);

            end

            gen = gen + 1;

        end

        outcome = [outcome min(fit)];

        time = time + 1;

    end

    sort(outcome)
    mean(outcome)
    std(outcome)

end

toc;
