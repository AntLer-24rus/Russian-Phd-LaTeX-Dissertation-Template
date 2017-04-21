function Main_functions()
% close all;
clear;clc;cla
t = [];
y = [];

load experiment_data;
y = ex.getRadius(0.5).getStep(10).getCut(1).rawVertical(2,1:end-1);
t = ex.getRadius(0.5).getStep(10).getCut(1).rawVertical(1,1:end-1);

% y = y - offset;
res = DropGrossError(y, 1);% 0 - îòáğàñûâàòü òîëüêî ãğóáûå
                           % 1 - ñïğàøèâàòü îá îòáğîñå ñğåäíèõ
% plot(t,res)
                           
                           
offset = findOffSet(y);
offset_drop = findOffSet(y(~(res==1)));
% y = y - offset;

plot([0 max(t)],[offset offset], '--k','lineWidth', 1)  %Ëèíèÿ íóëÿ
hold on
plot(t, y)                              %Ñûğûå äàííûå
plot(t(~(res==1)), y(~(res==1)), 'r')       %Ñûğûå äàííûå ñ îòáğîøåííûìè ãğóáûìè


plot(t(res==1), y(res==1), 'sr','MarkerFaceColor','r')        %Îøèáî÷íûå èçìåğåíèÿ (ãğóáûå)
plot(t(res==3), y(res==3), 'sg','MarkerFaceColor','g')        %Îøèáî÷íûå èçìåğåíèÿ (ïîäîçğåíèå)
% plot(t(res==3), y(res==3), 'sg','MarkerFaceColor','g')        %Íåëüçÿ îòêèäûâàòü

% [xFit_a, yFit_a] = createFit(t(~((res==1) | (res==3))), y(~((res==1) | (res==3))), 0.95, 1000);
% plot(xFit_a, yFit_a, 'y', 'lineWidth', 3)                   %Èíòåğïîëÿöèÿ ñî âñåìåè îòáğîøåííûìè òî÷êàìè

[xFit_d, yFit_d] = createFit(t(~(res==1)), y(~(res==1)), 0.95, 1000);
plot(xFit_d, yFit_d, 'r', 'lineWidth', 3)                   %Èíòåğïîëÿöèÿ ñ îòáğîøåííûìè ãğóáûìè
axis([0 max(t) min(y) max(y)])
grid on

set(gca, 'XTick', linspace(0, max(t), 10), 'XMinorTick', 'on')
xlabel('Âğåìÿ (ñåê)')
set(gca, 'YTick', linspace(min(y), max(y), 10), 'YMinorTick', 'on')
ylabel('Èçìåğåííîå çíà÷åíèå (Â)')
% difN = sort([max(yFit_d), max(yFit_a)]);


% fprintf('Ğàçíèöà ìàêñèìóìîâ â %%                 %7.3f\n', 100 - ((difN(1)*100)/difN(2)));
fprintf('Âûáîğî÷íîå ñğåäíåå                     %7.3f\n', mean(y));
sum_y_m = sum((y - mean(y)).^2);
D = sum_y_m/length(y);
fprintf('Äèñïåğñèÿ (âòîğîé öåíòğàëüíûé ìîìåíò)  %7.3f ÑÊÎ %7.4f\n', D, sqrt(D));
N = sum_y_m/(length(y) - 1);
fprintf('Íåñìåùåííàÿ îöåíêà                     %7.3f ÑÊÎ %7.4f\n', N, sqrt(N));
Kv = sqrt(N)/mean(y);
fprintf('Êîıôôèöåíò âàğèàöèè                     %7.3f â %% %7.3f\n', Kv, Kv * 100);

function out = DropGrossError(y_in, f)
p_l = 5;    l=1;
p_r = 0.1;  r=2;

n = length(y_in);
out = zeros(1, n);
y_error = false(1, n);
y_warn = false(1, n);
y_good = false(1, n);
y_idx = 1:n;
e = true;
while e
    y_idx_step = y_idx(~y_error);
    y_step = y_in(~y_error);
    
    n = length(y_step);
    t_tab = abs((tinv([p_l p_r] ./ 100, n - 2) .* sqrt(n-1)) ./ sqrt(n - 2 + tinv([p_l p_r] ./ 100, n - 2) .^ 2));
    
    y_t = zeros(1, 2);
    ind_y_t = zeros(1, 2);
    
    [y_t(1), ind_y_t(1)] = max(y_step);
    [y_t(2), ind_y_t(2)] = min(y_step);
    
    [d_max, d_max_ind] = max(abs(y_t-mean(y_step)));    
    d_max_ind = ind_y_t(d_max_ind);
    
    t = d_max/std(y_step);
    
    e = t > t_tab(r);
    w = t_tab(l) < t & t < t_tab(r);
    g = t < t_tab(l);
    
    if w && f
        hFq = figure;
        plot(y_step);hold on
        plot(d_max_ind, y_step(d_max_ind), 'sr','MarkerFaceColor','r')
        axis([0 n min(y_step)-2 max(y_step)+2])
        choice = ...
            questdlg('Óäàëèòü îòìå÷åííó íà ãğàôèêå òî÷êó?', ...
                          'Óäàëåíèå ãğóáûõ òî÷åê', ...
	                      'Äà','Íåò','Íåò');
        close(hFq)   
        e = strcmp('Äà', choice);
    end
    
    y_error(y_idx_step(d_max_ind)) = e;% | (w & f);
    y_warn(y_idx_step(d_max_ind)) = w;
    y_good(y_idx_step(d_max_ind)) = g;    
end
% out(y_warn) = 3;
out(y_error) = 1;
% out(y_error & y_warn) = 2;

% out(y_good) = 4;


function offset = clearOffSet(y)
    %% 
    
    %% Time specifications:
    Fs = 100;                      % samples per second
    N = size(y,2);
    %% Fourier Transform:
    X = fftshift(fft(x));
    %% Frequency specifications:
    dF = Fs/N;                      % hertz
    f = -Fs/2:dF:Fs/2-dF;           % hertz
    %% Plot the spectrum:
    figure;
    stem(f(end/2+1:5:end),abs(X(end/2+1:5:end))/N,'k');
    axis square
    axis([-1 50 0 max(abs(X)/N)])
    xlabel('Frequency (in hertz)');
    ylabel('Ñèëà (êÍ)')
    xlabel('×àñòîòà (Ãö)')
    grid on
    box off
    set(gca, 'TickDir', 'both')
    set(gca, 'XMinorTick', 'on')
    set(gca, 'YMinorTick', 'on')
    set(gca, 'XTick', 0:5:50)
    set(gca, 'YTick', linspace(0,max(abs(X)/N), 10))

    ac = abs(X(end/2+1))/N;
    X(end/2+1) = 0;
    ix = ifft(ifftshift(X));
    %%
    figure
    hRawS = plot(t,x,'k');
    axis square
    axis([0 max(t) min(x) max(x)])
    hold on
    grid on
    box off
    set(gca, 'TickDir', 'both')
    set(gca, 'XMinorTick', 'on')
    set(gca, 'YMinorTick', 'on')
    ylabel('Ñèëà (êÍ)')
    xlabel('Âğåìÿ (ñåê.)')
    hFS = plot(t,ix, '--k');
    plot(t, zeros(1,length(t)),'--k')
    set(gca,'XTick',[0:0.5:4, max(t)])
    set(gca,'YTick',linspace(min(x), max(x), 10))
    legend([hRawS,hFS],'C ïîñòîÿííîé ñîñòàâëÿşùåé','Áåç ïîñòîÿííîé ñîñòàâëÿşùåé')
    fprintf('Ïîñòîÿííàÿ ñîñòàâëÿşùàÿ %.4f\n',(ac) * sign(mean(x)))
    x_wc = ix;
    fprintf('RMS_raw %.3f\n', sqrt(sum(x.^2)/length(x)))
    fprintf('RMS %.3f\n', sqrt(sum(x_wc.^2)/length(x_wc)))
    fprintf('Mean %.3f\n', mean(x_wc))
    fprintf('Max_raw %.3f\n', max(x))
    fprintf('Max %.3f\n', max(x_wc))

function [xFit, yFit] = createFit(t, y, Smoothing, points)
x = t;
[xData, yData] = prepareCurveData(x, y);
ft = fittype('smoothingspline');
opts = fitoptions('Method', 'SmoothingSpline');
opts.SmoothingParam = Smoothing;

fitresult = fit( xData, yData, ft, opts );
xFit = linspace(0, max(x), points);
yFit = feval(fitresult, xFit)';