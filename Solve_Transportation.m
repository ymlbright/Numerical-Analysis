function [baseMatrix, rmin] = Solve_Transportation( startPoint, endPoint, priceMatrix)
% [plain, rmin] = Solve_Transportation Problem( startPoint, endPoint, priceMatrix)
% 使用表上作业法计算运输问题
%参数：
% startPoint    运输起点供应量
% endPoint      运输终点需求量
% priceMatraix  运输至各个目的地的费用，一行表示一个运输起点，每一列代表运输到对应终点的费用
%返回值：
% baseMatrix    运输方案，一行表示一个运输起点，每一列代表运输到对应终点的运输量
% rmin          最小运输花费
%测试:
% Solve_Transportation([14 27 19],[22 13 12 13],[6 7 5 3;8 4 2 7;5 9 10 6])
% Solve_Transportation([14 27 19],[22 13 12 18],[6 7 5 3;8 4 2 7;5 9 10 6])
% Solve_Transportation([35,25,15,40],[25,20,10,25,10,15,10],[6 5 2 6 3 6 3;3 7 5 8 6 9 2;4 8 6 5 5 8 5; 7 4 4 7 4 7 4;])

    % 初始化变量
    startPointLength = length(startPoint);
    endPointLength = length(endPoint);
    startSum = sum(startPoint);
    endSum = sum(endPoint);
    if startSum < endSum
        startPointLength = startPointLength + 1;
        startPoint(startPointLength) = endSum - startSum;
        priceMatrix(startPointLength,:) = 0;
    elseif startSum > endSum
        endPointLength = endPointLength + 1;
        startPoint(endPointLength) = startSum - endSum;
        priceMatrix(:,endPointLength) = 0;
    end
    init_stack(startPointLength+endPointLength);
    baseMatrix = -1*ones(startPointLength, endPointLength);

    % 西北角法选基
    i = 1; j = 1;
    while i+j<=startPointLength+endPointLength
        if startPoint(i) < endPoint(j)
            baseMatrix(i,j) = startPoint(i);
            endPoint(j) = endPoint(j) - startPoint(i);
            i = i + 1;
        else
            while j<=endPointLength && endPoint(j) <= startPoint(i)
                baseMatrix(i,j) = endPoint(j);
                startPoint(i) = startPoint(i) - endPoint(j);
                j = j + 1;
            end
        end
    end

    count = 0;
    % 运输表迭代
    checkMatrix = calc_potential(baseMatrix, priceMatrix); % 计算势
    [C,I] = max(checkMatrix);
    [C,J] = max(C);
    while C > 0
        count = count + 1;
        if count > 10000; disp('可能发生了退化！'); break; end
        [path, vmin] = find_cycle(I(J), J, baseMatrix); % 寻找闭合回路
        firstZeroFlag = true;
        i = 1;
        while path(i,1)~= 0
            if mod(i,2)
                baseMatrix(path(i,1),path(i,2)) = baseMatrix(path(i,1),path(i,2)) - vmin;
                if baseMatrix(path(i,1),path(i,2)) == 0 && firstZeroFlag; firstZeroFlag = false; baseMatrix(path(i,1),path(i,2)) = -1; end
            else
                if baseMatrix(path(i,1),path(i,2)) == -1; baseMatrix(path(i,1),path(i,2)) = baseMatrix(path(i,1),path(i,2)) + vmin + 1; else baseMatrix(path(i,1),path(i,2)) = baseMatrix(path(i,1),path(i,2)) + vmin; end
            end
            i = i + 1;
        end
        checkMatrix = calc_potential(baseMatrix, priceMatrix);
        [C,I] = max(checkMatrix);
        [C,J] = max(C);
    end
    
    % 计算最小值
    baseMatrix(baseMatrix==-1) = 0;
    rmin = sum(sum(baseMatrix.*priceMatrix));
end

function [checkMatrix] = calc_potential(baseMatrix, priceMatrix)
    clear()
    [m,n] = size(baseMatrix);
    u = zeros(m,1);
    v = zeros(n,1);
    checkMatrix = zeros(m, n);
    for i = 1:n
        if baseMatrix(1,i) ~= -1; push(1,i,0); checkMatrix(1,i)=1; end
    end
    % 计算位势
    while ~empty()
        [x,y,f] = pop();
        if f    % scan line
            u(x) = priceMatrix(x,y) - v(y);
            for i = 1:n
                if baseMatrix(x,i) ~= -1 && ~checkMatrix(x,i); push(x,i,0); checkMatrix(x,i)=1; end
            end
        else    % scan column
            v(y) = priceMatrix(x,y) - u(x);
            for i = 1:m
                if baseMatrix(i,y) ~= -1 && ~checkMatrix(i,y); push(i,y,1); checkMatrix(i,y)=1; end
            end
        end
    end
    % 计算检验数
    for i = 1:m
        for j = 1:n
            if checkMatrix(i,j); checkMatrix(i,j) = 0; else checkMatrix(i,j) = u(i) + v(j) - priceMatrix(i,j); end
        end
    end
end

function [path, vmin] = find_cycle(I, J, baseMatrix)
    clear()
    [m,n] = size(baseMatrix);
    path = zeros(m+n,2);
    vmin = Inf;
    flagMatrix = zeros(m, n);
    flagMatrix(I,J) = 1;
    continueIndex = 1;
    finishFlag = false;
    push(I,J,0); % index x, index y, type
    while ~empty()
        [x, y, t] = top();
        foundFlag = false;
        if t % scan line
            for i = continueIndex:n
                if baseMatrix(x,i) ~= -1 && ~flagMatrix(x,i)
                    foundFlag = true;
                    break;
                elseif stack_size() ~= 1 && x == I && i == J && finish_check()
                    finishFlag = true;
                    break;
                end
            end
            if foundFlag
                continueIndex = 1;
                flagMatrix(x,i) = 1;
                push(x,i,0);
                continue;
            elseif ~finishFlag % not found in line, try to traceback and check finishFlag
                [x, y, ~] = pop();
                flagMatrix(x,y) = 0;
                continueIndex = x+1;
            end
        else % scan column
            for i = continueIndex:m
                if baseMatrix(i,y) ~= -1 && ~flagMatrix(i,y)
                    foundFlag = true;
                    break;
                elseif stack_size() ~= 1 && i == I && y == J && finish_check()
                    finishFlag = true;
                    break;
                end
            end
            if foundFlag
                continueIndex = 1;
                flagMatrix(i,y) = 1;
                push(i,y,1);
                continue;
            elseif stack_size() == 1 && x == I && y == J && finish_check()
                continueIndex = 1;
                [x, y, ~] = pop();
                push(x,y,1);
            elseif ~finishFlag % not found in line, try to traceback and check finishFlag
                [x, y, ~] = pop();
                flagMatrix(x,y) = 0;
                continueIndex = y+1;
            end
        end
        if finishFlag
            i = 1;
            while ~empty()
                [path(i,1),path(i,2),~] = pop();
                if mod(i,2) && vmin > baseMatrix(path(i,1),path(i,2)) && baseMatrix(path(i,1),path(i,2)) ~= -1; vmin = baseMatrix(path(i,1),path(i,2)); end
                i = i + 1;
            end
            return
        end
    end
    disp('ERROR: Can not find a cycle in baseMatrix!')
    disp('Press Ctrl+C to terminate programe.')
    pause
end

function [r] = finish_check()
    global STACK;
    global STACK_END;
    if STACK(1,1) == STACK(2,1) && STACK(1,1) == STACK(STACK_END,1)
        r = ~((STACK(2,2)>STACK(1,2)&&STACK(STACK_END,2)>STACK(1,2))||(STACK(2,2)<STACK(1,2)&&STACK(STACK_END,2)<STACK(1,2)));
    elseif STACK(1,2) == STACK(2,2) && STACK(1,2) == STACK(STACK_END,2)
        r = ~((STACK(2,1)>STACK(1,1)&&STACK(STACK_END,1)>STACK(1,1))||(STACK(2,1)<STACK(1,1)&&STACK(STACK_END,1)<STACK(1,1)));
    else
        r = 1;
    end
end

function init_stack(M)
    global STACK;
    global STACK_LENGTH;
    global STACK_END;
    STACK = zeros(M,4);
    STACK_LENGTH = M;
    STACK_END = 0;
end

function push(x, y, f)
    global STACK;
    global STACK_LENGTH;
    global STACK_END;
    if STACK_END+1 <= STACK_LENGTH
        STACK_END = STACK_END + 1;
        STACK(STACK_END,1) = x;
        STACK(STACK_END,2) = y;
        STACK(STACK_END,3) = f;
    else
        disp('Warning: StackOverflow!');
    end
end

function [x, y, f] = pop()
    global STACK;
    global STACK_END;
    if STACK_END > 0
        x = STACK(STACK_END,1);
        y = STACK(STACK_END,2);
        f = STACK(STACK_END,3);
        STACK_END = STACK_END - 1;
    else
        disp('Warning: Stack is empty!');
    end
end

function clear()
    global STACK_END;
    STACK_END = 0;
end

function [r] = empty()
    global STACK_END;
    r = STACK_END == 0;
end

function [x, y, f] = top()
    global STACK;
    global STACK_END;
    if STACK_END > 0
        x = STACK(STACK_END,1);
        y = STACK(STACK_END,2);
        f = STACK(STACK_END,3);
    else
        disp('Warning: Stack is empty!');
    end
end

function [r] = stack_size()
    global STACK_END;
    r = STACK_END;
end