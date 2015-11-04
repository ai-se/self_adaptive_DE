%**************************************************************************************************
%Reference:  J. J. Liang, A. K. Qin, P. N. Suganthan, and S. Baskar, "Comprehensive
%                     learning particle swarm optimizer for global optimization of multimodal
%                     functions," IEEE Trans. Evolut. Comput., vol. 10, no. 3, pp. 281¨C295, 2006.
%
% Note: We obtained the MATLAB source code from the authors, and did some
%           minor revisions in order to solve the 25 benchmark test functions,
%           however, the main body was not changed.
%**************************************************************************************************

clc;
clear all
tic

'CLPSO'

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

    % Main body which was provided by the authors

    me = 300000/40;
    ps = 40;
    maxFES = 300000;

    outcome = [];
    
    % The total number of runs
    totalTime = 1;

    for time = 1 : totalTime

        rand('state', sum(100 * clock));

        t = 0 : 1 / (ps - 1) : 1;
        t = 5 .* t;
        Pc = 0.0 + (0.5 - 0.0) .* (exp(t) - exp(t(1))) ./ (exp(t(ps)) - exp(t(1)));
        m = 0 .* ones(ps, 1);
        iwt = 0.9 - (1 : me) * (0.7 / me);
        cc = [1.49445 1.49445]; %acceleration constants
        mv = 0.2 * (lu(2, :) - lu(1, :));
        VRmin = repmat(lu(1, :), ps, 1);
        VRmax = repmat(lu(2, :), ps, 1);
        Vmin = repmat(-mv, ps, 1);
        Vmax = -Vmin;

        pos = VRmin + (VRmax - VRmin) .* rand(ps, n);

        e = benchmark_func(pos, problem, o, A, M, a, alpha, b);

        fitcount = ps;
        vel = Vmin + 2 .* Vmax .* rand(ps, n); %initialize the velocity of the particles

        pbest = pos;
        pbestval = e;  %initialize the pbest and the pbest's fitness value

        [gbestval, gbestid] = min(pbestval);
        gbest = pbest(gbestid, :); %initialize the gbest and the gbest's fitness value
        gbestrep = repmat(gbest, ps, 1);

        stay_num = zeros(ps, 1);

        ai = zeros(ps, n);
        f_pbest = 1 : ps;
        f_pbest = repmat(f_pbest', 1, n);
        for k = 1 : ps

            ar = randperm(n);
            ai(k, ar(1 : m(k))) = 1;
            fi1 = ceil(ps * rand(1, n));
            fi2 = ceil(ps * rand(1, n));
            fi = (pbestval(fi1) < pbestval(fi2))' .* fi1 + (pbestval(fi1) >= pbestval(fi2))' .* fi2;
            bi = ceil(rand(1, n) - 1 + Pc(k));

            if bi == zeros(1, n),
                rc = randperm(n);
                bi(rc(1)) = 1;
            end

            f_pbest(k, :) = bi .* fi + (1 - bi) .* f_pbest(k, :);

        end

        stop_num = 0;

        for i = 2 : me

            valid = [];
            for k = 1 : ps

                if stay_num(k) >= 5

                    stay_num(k) = 0;
                    ai(k, :) = zeros(1, n);
                    f_pbest(k, :) = k .* ones(1, n);
                    ar = randperm(n);
                    ai(k, ar(1 : m(k))) = 1;
                    fi1 = ceil(ps * rand(1, n));
                    fi2 = ceil(ps * rand(1, n));
                    fi = (pbestval(fi1) < pbestval(fi2))' .* fi1 + (pbestval(fi1) >= pbestval(fi2))' .* fi2;
                    bi = ceil(rand(1, n) - 1 + Pc(k));

                    if bi == zeros(1, n),
                        rc = randperm(n);
                        bi(rc(1)) = 1;
                    end

                    f_pbest(k, :) = bi .* fi + (1 - bi) .* f_pbest(k, :);

                end

                for dimcnt = 1 : n

                    pbest_f(k, dimcnt) = pbest(f_pbest(k, dimcnt), dimcnt);

                end

                temp = pos(k, :);
                temp_ = temp + vel(k, :);

                aa(k, :) = cc(1) .* (1 - ai(k, :)) .* rand(1, n) .* (pbest_f(k, :) - pos(k, :)) + cc(2) .* ai(k, :) .* rand(1, n) .* (gbestrep(k, :) - pos(k, :));
                vel(k, :) = iwt(i) .* vel(k, :) + aa(k, :);
                vel(k, :) = (vel(k, :) > mv) .* mv + (vel(k, :) <= mv) .* vel(k, :);
                vel(k, :) = (vel(k, :) < (-mv)) .* (-mv) + (vel(k, :) >= (-mv)) .* vel(k, :);
                pos(k, :) = pos(k, :) + vel(k, :);

                if (sum(pos(k, :) > VRmax(k, :)) + sum(pos(k, :) < VRmin(k, :))) == 0;
                    valid = [valid k];
                    fitcount = fitcount + 1;
                end

            end

            if ~isempty(valid)

                e(valid, 1) = benchmark_func(pos(valid, :), problem, o, A, M, a, alpha, b);

                for k = 1 : length(valid)

                    tmp = (pbestval(valid(k)) <= e(valid(k)));
                    if tmp == 1
                        stay_num(valid(k)) = stay_num(valid(k)) + 1;
                    end
                    temp = repmat(tmp, 1, n);
                    pbest(valid(k), :) = temp .* pbest(valid(k), :) + (1 - temp) .* pos(valid(k), :);
                    pbestval(valid(k)) = tmp .* pbestval(valid(k)) + (1 - tmp) .* e(valid(k)); %update the pbest
                    if pbestval(valid(k)) < gbestval
                        gbest = pbest(valid(k), :);
                        gbestval = pbestval(valid(k));
                        gbestrep = repmat(gbest, ps, 1); %update the gbest
                    end

                end

            end

            if fitcount >= maxFES
                break;
            end
            if (i == me) & (fitcount < maxFES)
                i = i-1;
            end

        end

        outcome = [outcome min(pbestval)];

    end

    sort(outcome)
    mean(outcome)
    std(outcome)

end
toc